#include "arg_parse.hpp"
#include "models.hpp"

#include "estimate_hmm_parameter.hpp"
#include "create_tracks.hpp"
#include "ecm.hpp"
#include "newick.hpp"

#include "run.hpp"

#include <omp.h>

#include <sstream>

struct TrackCLIParams
{
    unsigned threads = omp_get_max_threads();
    float phylo_threshold = 0.1;
    bool phylo_smooth = false;
    bool phylo_raw = true;
    bool phylo_power = true;
    std::string output_path = "";
};

void merge_job_output_files(const std::string & base_name, const unsigned jobs)
{
    std::string job_output_file = base_name + ".0";
    if (rename(job_output_file.c_str(), base_name.c_str()) == -1)
    {
        printf("Could not rename file %s to %s (error code: %d, %s)\n", job_output_file.c_str(), base_name.c_str(), errno, strerror(errno));
        exit(-1);
    }

    FILE *merged_file = fopen(base_name.c_str(), "ab");
    for (unsigned job_id = 1; job_id < jobs; ++job_id)
    {
        job_output_file = base_name + "." + std::to_string(job_id);
        append(merged_file, job_output_file.c_str());
        unlink(job_output_file.c_str());
    }
    fclose(merged_file);
}

void run_tracks(const std::string & alignment_path, const Model & model, const TrackCLIParams & params)
{
    unsigned jobs = (params.threads > 1) ? params.threads * 10 : 1;

    std::vector<Data> data(params.threads);

    // store output next to alignment file if no output_path was specified
    std::string output_folder = params.output_path;
    if (output_folder == "")
    {
        size_t pos = alignment_path.find_last_of('/');
        if (pos == std::string::npos)
            output_folder = "./";
        else
            output_folder = alignment_path.substr(0, pos);
    }

    // prepare alignment
    std::vector<alignment_t> alignments;
    const uint16_t leaves = (model.phylo_array.size() + 1)/2;
    for (unsigned i = 0; i < params.threads; ++i)
    {
        alignments.emplace_back(leaves);
    }

    std::vector<std::vector<double> > lpr_per_codon(params.threads), bls_per_bp(params.threads);

    parallel_maf_reader maf_rd(alignment_path.c_str(), jobs, &model.seqid_to_phyloid, true);
    jobs = maf_rd.get_jobs(); // maybe file is too small and a smaller number of jobs is used

    #pragma omp parallel for num_threads(params.threads) schedule(dynamic, 1) default(none) shared(jobs, alignments, maf_rd, data, model, params, output_folder, lpr_per_codon, bls_per_bp)
    for (unsigned job_id = 0; job_id < jobs; ++job_id)
    {
        unsigned thread_id = omp_get_thread_num();
        auto & aln = alignments[thread_id];
        // TODO: merge arrays file_range_pos and _end are thread-safe for cache locality. should be thread-safe

        std::vector<double> scores;
        std::vector<scored_region> region;

        const std::string filename_power = output_folder + "/PhyloCSFpower.wig." + std::to_string(job_id);
        FILE *file_power = fopen(filename_power.c_str(), "w");

        FILE *file_score_raw[6];
        FILE *file_score[6];

        for (uint8_t i = 0; i < 6; ++i)
        {
            char strand = (i < 3) ? '+' : '-';
            unsigned frame = (i % 3) + 1;

            if (params.phylo_raw)
            {
                const std::string filename_score_raw = output_folder + "/PhyloCSFRaw" + std::string(1, strand) + std::to_string(frame) + ".wig." + std::to_string(job_id);
                file_score_raw[i] = fopen(filename_score_raw.c_str(), "w");
                if (file_score_raw[i] == NULL)
                {
                    printf("Error creating file!");
                    exit(1);
                }
            }

            if (params.phylo_smooth)
            {
                const std::string filename_score = output_folder + "/PhyloCSF" + std::string(1, strand) + std::to_string(frame) + ".wig." + std::to_string(job_id);
                file_score[i] = fopen(filename_score.c_str(), "w");
                if (file_score[i] == NULL)
                {
                    printf("Error creating file!");
                    exit(1);
                }
            }
        }

        maf_rd.skip_partial_alignment(aln, job_id);

        while (maf_rd.get_next_alignment(aln, job_id))
        {
//            printf("%10ld\t10%ld\n", aln.start_pos, aln.seqs[0].size());

            // on first iteration, compute bls scores (used by all 6 frames then!)
            bls_per_bp[thread_id].clear();
            compute_bls_score<true>(model.phylo_tree, aln, model, bls_per_bp[thread_id]);

            const uint64_t orig_start_pos = aln.start_pos;

            for (char strand = '+'; strand <= '-'; strand += 2)
            {
                for (unsigned frame = 1; frame <= 3; ++frame)
                {
                    const uint8_t file_index = frame - 1 + (strand == '+' ? 0 : 3);
                    aln.update_seqs(orig_start_pos, strand, frame);

                    lpr_per_codon[thread_id].clear();
                    run_tracks(data[thread_id], model, aln, lpr_per_codon[thread_id]/*, bls_per_bp[thread_id]*/);
                    data[thread_id].clear();

                    if (strand == '-')
                    {
                        std::reverse(lpr_per_codon[thread_id].begin(), lpr_per_codon[thread_id].end());

                        const size_t excess_basepairs = aln.length() % 3;
                        aln.start_pos += excess_basepairs;
                    }

                    if (params.phylo_power && frame == 3 && strand == '+')
                    {
                        // since we iterate over a codon array, there must be 3 bp for each codon
                        // the last 0-2 remaining basepairs in the bls array do not have a codon entry
                        assert(lpr_per_codon[thread_id].size() * 3 <= bls_per_bp[thread_id].size());

                        fprintf(file_power, "fixedStep chrom=%s start=%ld step=3 span=3\n", aln.chrom.c_str(), aln.start_pos);
                        for (uint32_t pos = frame - 1; pos + 2 < bls_per_bp[thread_id].size(); pos += 3)
                        {
                            const float bls_codon_avg = (bls_per_bp[thread_id][pos]
                                                         +  bls_per_bp[thread_id][pos + 1]
                                                         +  bls_per_bp[thread_id][pos + 2]) / 3.0;

                            my_fprintf(file_power, "%.4f", bls_codon_avg);
                        }
                    }

                    size_t bls_pos = 0;
                    // compute offset
                    if (strand == '+')
                        bls_pos += (frame - 1);
                    else
                        bls_pos += aln.length() % 3;

                    int64_t prevPos = -4;
                    int64_t startBlockPos = aln.start_pos;
                    for (uint32_t xx = 0; xx < lpr_per_codon[thread_id].size(); ++xx, bls_pos += 3)
                    {
                        const float bls_codon_sum = bls_per_bp[thread_id][bls_pos]
                                                    + bls_per_bp[thread_id][bls_pos + 1]
                                                    + bls_per_bp[thread_id][bls_pos + 2];
                        if (bls_codon_sum < 0.1f * 3)
                        {
                            if (params.phylo_smooth && scores.empty())
                                startBlockPos = aln.start_pos + ((xx + 1) * 3);
                            continue;
                        }

                        int64_t newPos = aln.start_pos + (xx * 3);

                        if (prevPos + 3 != newPos)
                        {
                            if (params.phylo_raw)
                                fprintf(file_score_raw[file_index], "fixedStep chrom=%s start=%ld step=3 span=3\n", aln.chrom.c_str(), newPos);

                            if (params.phylo_smooth && !scores.empty())
                            {
                                fprintf(file_score[file_index], "fixedStep chrom=%s start=%ld step=3 span=3\n", aln.chrom.c_str(), startBlockPos);

                                process_scores(model._hmm, scores, startBlockPos, region, SCORE_CODON);

                                for(size_t i = 0; i < region.size(); i++)
                                {
                                    my_fprintf(file_score[file_index], "%.3f", region[i].log_odds_prob);
                                }

                                scores.clear();
                                region.clear();
                                startBlockPos = aln.start_pos + (xx * 3);
                            }
                        }

                        prevPos = newPos;

                        if (params.phylo_raw)
                            my_fprintf(file_score_raw[file_index], "%.3f", lpr_per_codon[thread_id][xx]);

                        if (params.phylo_smooth)
                            scores.push_back(lpr_per_codon[thread_id][xx]);
                    }

                    if (params.phylo_smooth && !scores.empty())
                    {
                        fprintf(file_score[file_index], "fixedStep chrom=%s start=%ld step=3 span=3\n", aln.chrom.c_str(), startBlockPos);
                        process_scores(model._hmm, scores, startBlockPos, region, SCORE_CODON);

                        for(size_t i = 0; i < region.size(); i++)
                        {
                            my_fprintf(file_score[file_index], "%.3f", region[i].log_odds_prob);
                        }

                        scores.clear();
                        region.clear();
                    }
                }

                // compute reverse complement for neg. strand
                for (auto & seq : aln.seqs)
                {
                    std::reverse(seq.begin(), seq.end());
                    for (uint64_t j = 0; j < seq.size(); ++j)
                    {
                        seq[j] = complement(seq[j]);
                    }
                }
            }

            for (auto & seq : aln.seqs)
                seq = "";
        }

        if (params.phylo_power)
            fclose(file_power);

        for (uint8_t i = 0; i < 6; ++i)
        {
            if (params.phylo_raw)
                fclose(file_score_raw[i]);

            if (params.phylo_smooth)
                fclose(file_score[i]);
        }
    }

    // merge power file
    if (params.phylo_power)
        merge_job_output_files(output_folder + "/PhyloCSFpower.wig", jobs);

    // merge 6 frame files for scores and raw scores
    for (char strand = '+'; strand <= '-'; strand += 2)
    {
        for (unsigned frame = 1; frame <= 3; ++frame)
        {
            if (params.phylo_raw)
                merge_job_output_files(output_folder + "/PhyloCSFRaw" + std::string(1, strand) + std::to_string(frame) + ".wig", jobs);

            if (params.phylo_smooth)
                merge_job_output_files(output_folder + "/PhyloCSF" + std::string(1, strand) + std::to_string(frame) + ".wig", jobs);
        }
    }
}

int main_tracks(int argc, char **argv)
{
    ArgParse args("phylocsf++ tracks",
                  "Computes PhyloCSF score tracks for each codon and all 6 frames from alignments \n"
                  "in FASTA or MAF files. Output is written to wig files. Optionally a PhyloCSF \n"
                  "power track containing the branch length scores can written to a wig file \n"
                  "(confidence of the PhyloCSF scores).");

    TrackCLIParams params;

    std::string model_list = get_list_of_models();

    args.add_option("output-phylo", ArgParse::Type::BOOL, "Compute all 6 smoothened PhyloCSF tracks. Requires coding exon coordinates and genome length. Default: " + std::to_string(params.phylo_smooth), ArgParse::Level::GENERAL, false);
    args.add_option("output-raw-phylo", ArgParse::Type::BOOL, "Compute all 6 unsmoothened PhyloCSF tracks. Default: " + std::to_string(params.phylo_raw), ArgParse::Level::GENERAL, false);
    args.add_option("output-power", ArgParse::Type::BOOL, "Output confidence track (branch length score). Default: " + std::to_string(params.phylo_power), ArgParse::Level::GENERAL, false);
    args.add_option("power-threshold", ArgParse::Type::FLOAT, "Minimum confidence score to output PhyloCSF score. Default: " + std::to_string(params.phylo_threshold), ArgParse::Level::GENERAL, false);
    args.add_option("genome-length", ArgParse::Type::INT, "Total genome length (needed for --output-phylo).", ArgParse::Level::GENERAL, false);
    args.add_option("coding-exons", ArgParse::Type::STRING, "BED-like file (chrom name, strand, phase, start, end) with coordinates of coding exons (needed for --output-phylo).", ArgParse::Level::GENERAL, false);
    args.add_option("threads", ArgParse::Type::INT, "Parallelize scoring of multiple MSAs in a file using multiple threads. Default: " + std::to_string(params.threads), ArgParse::Level::GENERAL, false);
    args.add_option("output", ArgParse::Type::STRING, "Path where tracks in wig format will be written. If it does not exist, it will be created. Default: output files are stored next to the input files.", ArgParse::Level::GENERAL, false);
    args.add_option("mapping", ArgParse::Type::STRING, "If the MSAs don't use common species names (like Human, Chimp, etc.) you can pass a two-column tsv file with a name mapping.", ArgParse::Level::GENERAL, false);
    args.add_option("species", ArgParse::Type::STRING, "Comma separated list of species to reduce <model> to a subset of species to improve running time, e.g., \"Human,Chimp,Seaturtle\"", ArgParse::Level::GENERAL, false);
    args.add_option("model-info", ArgParse::Type::STRING, "Output the organisms included in a specific model. Included models are: " + model_list + ".", ArgParse::Level::GENERAL, false);
    // TODO: add dry-run option (to check whether all species names can be mapped and it doesn't fail after 2 days of runtime)
    // or even better: continue where the tool left off!

    // TODO: ignore sequences not occuring in model (instead of failing)

    args.add_positional_argument("model", ArgParse::Type::STRING, "Path to parameter files, or one of the following predefined models: " + model_list + ".", true);
    args.add_positional_argument("alignments", ArgParse::Type::STRING, "One or more files containing multiple sequence alignments. Formats: MAF and multi FASTA. Multiple MSAs can be stored in a single file separated by empty lines.", true);
    args.parse_args(argc, argv);

    // additional help strings
    if (args.is_set("model-info"))
    {
        const std::string & model_name = args.get_string("model-info");
        return print_model_info(model_name);
    }

    args.check_args(); // check here because if "--model-info" is set, we don't want to require mandatory arguments

    // retrieve flags/options that were set
    if (args.is_set("threads"))
        params.threads = args.get_int("threads");
    if (args.is_set("output-phylo"))
        params.phylo_smooth = args.get_bool("output-phylo");
    if (args.is_set("output-raw-phylo"))
        params.phylo_raw = args.get_bool("output-raw-phylo");
    if (args.is_set("power-threshold"))
        params.phylo_power = args.get_bool("power-threshold");

    // this has to hold true: phylo_smooth => (args.is_set("genome-length") && args.is_set("coding-exons"))
    if (params.phylo_smooth && (!args.is_set("genome-length") || !args.is_set("coding-exons")))
    {
        print_error_msg("For smoothened tracks (--output-phylo) you need to provide --genome-length and --coding-exons.");
        return -1;
    }

    if (args.is_set("output"))
    {
        params.output_path = args.get_string("output");
        // create a directory if it doesn't exist yet
        if (create_directory(params.output_path))
            print_info_msg("Created the output directory.");
    }

    if (args.is_set("mapping"))
    {
        std::string mapping_file = args.get_string("mapping");
        update_sequence_name_mapping(mapping_file);
    }

    if (!args.is_set_positional("model") || !args.is_set_positional("alignments"))
    {
        print_error_msg("No model or alignments provided.");
        return -1;
    }

    const std::string alignment_path = args.get_positional_argument("alignments");

    // load and prepare model
    uint64_t genome_length = 0;
    if (args.is_set("genome-length"))
        genome_length = args.get_int("genome-length");

    std::string coding_exons_path = "";
    if (args.is_set("coding-exons"))
        coding_exons_path = args.get_string("coding-exons");

    Model model;
    load_model(model, args.get_positional_argument("model"), args.get_string("species"),
               params.phylo_smooth, genome_length, coding_exons_path);

    // run for every alignment file
    for (uint16_t i = 1; i < args.positional_argument_size(); ++i)
    {
        // printf("%d, %s\n", i, args.get_positional_argument(i).c_str());
        run_tracks(args.get_positional_argument(i).c_str(), model, params);
    }

    return 0;
}

// genome_length        = 2870184193
// model_str            = 100vertebrates
// selected_species_str = Rat,Frog_X._tropicalis,Zebrafish,Platypus,Chicken,Turkey,Panda,Mouse,Prairie_vole,Rhesus,Rabbit,Cat,Opossum,Cow,Human,Chimp,Dog,Guinea_pig
// aln_path             = /ccb/salz4-2/pocki/PhyloCSF++data/rn6/rn6.20way.chr20.maf
// hmm_data_path        = /ccb/salz4-2/pocki/PhyloCSF++data/rn6/RatCodingExons.txt
// output_folder        = /ccb/salz4-2/pocki/PhyloCSF++data/rn6/out_chr20_30_times_100jobs

// genome_length        = 1065365434
// model_str            = 53birds
// selected_species_str =
// aln_path             = /home/chris/dev-uni/PhyloCSF++/phylo_galgal/galGal6_chrM.maf
// hmm_data_path        = /home/chris/dev-uni/PhyloCSF++/phylo_galgal/ChickenCodingExonsV2.txt
// output_folder        = /home/chris/dev-uni/PhyloCSF++/galgal_out

// genome_length        = 12157105
// model_str            = 7yeast
// selected_species_str =
// aln_path             = /home/chris/dev-uni/PhyloCSF++/phylo_yeast/chrI.maf
// hmm_data_path        = /home/chris/dev-uni/PhyloCSF++/phylo_yeast/YeastCodingExons.txt
// output_folder        = /home/chris/dev-uni/PhyloCSF++/yeast_out_test