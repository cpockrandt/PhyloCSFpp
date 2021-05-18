#include "common.hpp"

#include "arg_parse.hpp"
#include "models.hpp"

#include "estimate_hmm_parameter.hpp"
#include "create_tracks.hpp"
#include "ecm.hpp"
#include "newick.hpp"

#include "run.hpp"

#include <sstream>

struct TrackCLIParams
{
    unsigned threads = omp_get_max_threads();
    float phylo_threshold = 0.1;
    bool phylo_smooth = false;
    bool phylo_raw = true;
    bool phylo_power = true;
    bool phylo_bed_file = false;
    bool smooth_wo_wig = false;
    std::string output_path;
};

void merge_job_output_files(const std::string & base_name, const unsigned jobs, const bool final_file_exists)
{
    std::string job_output_file;

    // the first file processed: rename .wig.0 to .wig and append .wig.1, .wig.2, etc. to it
    // for every subsequent file: append .wig.0, .wig.1, .wig.2, etc. to (already) existing .wig

    if (!final_file_exists)
    {
        job_output_file = base_name + ".0";
        if (rename(job_output_file.c_str(), base_name.c_str()) == -1)
        {
            printf("Could not rename file %s to %s (error code: %d, %s)\n", job_output_file.c_str(), base_name.c_str(), errno, strerror(errno));
            exit(-1);
        }
    }

    FILE *merged_file = fopen(base_name.c_str(), "ab");
    unsigned job_id = (final_file_exists) ? 0 : 1;
    for (; job_id < jobs; ++job_id)
    {
        job_output_file = base_name + "." + std::to_string(job_id);
        append(merged_file, job_output_file.c_str());
        unlink(job_output_file.c_str());
    }
    fclose(merged_file);
}

void run_tracks(const std::string & alignment_path, const Model & model, const TrackCLIParams & params,
                const uint32_t file_id, const uint32_t files, std::vector<std::vector<bool> > & species_seen_in_alignment)
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
    maf_rd.setup_progressbar(file_id, files);

    #pragma omp parallel for num_threads(params.threads) schedule(dynamic, 1)
    for (unsigned job_id = 0; job_id < jobs; ++job_id)
    {
        unsigned thread_id = omp_get_thread_num();
        auto & aln = alignments[thread_id];

        std::vector<double> scores;
        std::vector<scored_region> region;
        std::vector<scored_bed_region> bedregions;

        FILE *file_power = NULL;

        if (params.phylo_power)
        {
            const std::string filename_power = output_folder + "/PhyloCSFpower.wig." + std::to_string(job_id);
            file_power = fopen(filename_power.c_str(), "w");
            if (file_power == NULL)
            {
                printf(OUT_ERROR "Error creating file!\n" OUT_RESET);
                exit(1);
            }
        }

        FILE *file_score_raw[6];
        FILE *file_score[6];
        FILE *file_score_bed[6];

        for (uint8_t i = 0; i < 6; ++i)
        {
            const char strand = (i < 3) ? '+' : '-';
            const unsigned frame = (i % 3) + 1;

            if (params.phylo_raw)
            {
                const std::string filename_score_raw = output_folder + "/PhyloCSFRaw" + std::string(1, strand) + std::to_string(frame) + ".wig." + std::to_string(job_id);
                file_score_raw[i] = fopen(filename_score_raw.c_str(), "w");
                if (file_score_raw[i] == NULL)
                {
                    printf(OUT_ERROR "Error creating file!\n" OUT_RESET);
                    exit(1);
                }
            }

            if (params.phylo_smooth)
            {
                const std::string filename_score = output_folder + "/PhyloCSF" + std::string(1, strand) + std::to_string(frame) + ".wig." + std::to_string(job_id);
                file_score[i] = fopen(filename_score.c_str(), "w");

                if (file_score[i] == NULL)
                {
                    printf(OUT_ERROR "Error creating wig file!\n" OUT_RESET);
                    exit(1);
                }
                const std::string filename_bed = output_folder + "/PhyloCSF" + std::string(1, strand) + std::to_string(frame) + "Regions.bed." + std::to_string(job_id);
                file_score_bed[i] = fopen(filename_bed.c_str(), "w");
                if (file_score_bed[i] == NULL) {
                    printf(OUT_ERROR "Error creating bed file!\n" OUT_RESET);
                    exit(1);
                }
            }
        }

        maf_rd.skip_partial_alignment(aln, job_id);

        while (maf_rd.get_next_alignment(aln, job_id, &species_seen_in_alignment[thread_id]))
        {
            maf_rd.print_progress();

            // on first iteration, compute bls scores (used by all 6 frames then!)
            bls_per_bp[thread_id].clear();
            compute_bls_score<true>(model.phylo_tree, aln, model, bls_per_bp[thread_id]);

            if (params.phylo_power)
            {
                // skip the first 0-2 bases as in aln.update_seqs(orig_start_pos, strand = '+', frame = 3);
                int64_t skip_bases = static_cast<int64_t>(3 - aln.start_pos) % 3;
                if (skip_bases < 0)
                    skip_bases += 3;

                // do not print header if no values are following
                if (2 + 2 < bls_per_bp[thread_id].size())
                    fprintf(file_power, "fixedStep chrom=%s start=%" PRId64 " step=3 span=3\n", aln.chrom.c_str(), aln.start_pos + skip_bases);

                for (uint32_t pos = 2; pos + 2 < bls_per_bp[thread_id].size(); pos += 3)
                {
                    const float bls_codon_avg = (bls_per_bp[thread_id][pos]
                                                 +  bls_per_bp[thread_id][pos + 1]
                                                 +  bls_per_bp[thread_id][pos + 2]) / 3.0;

                    my_fprintf(file_power, "%.4f", bls_codon_avg);
                }
            }

            const uint64_t orig_start_pos = aln.start_pos;

            for (char strand = '+'; strand <= '-'; strand += 2)
            {
                for (unsigned frame = 1; frame <= 3; ++frame)
                {
                    const uint8_t file_index = (frame - 1) + (strand == '+' ? 0 : 3);
                    aln.update_seqs(orig_start_pos, strand, frame);

                    lpr_per_codon[thread_id].clear();
                    run_tracks(data[thread_id], model, aln, lpr_per_codon[thread_id]);
                    data[thread_id].clear();

                    if (strand == '-')
                    {
                        std::reverse(lpr_per_codon[thread_id].begin(), lpr_per_codon[thread_id].end());

                        // aminoacids are translated from right to left, i.e., remove the remaining 0-2 dna bases from
                        // the left, i.e., increase the start_pos
                        aln.start_pos += aln.length() % 3;
                    }

                    int64_t prevPos = -4;
                    int64_t startBlockPos = aln.start_pos;

                    // since we iterate over a codon array, there must be 3 bp for each codon
                    // the last 0-2 remaining basepairs in the bls array do not have a codon entry
                    assert(lpr_per_codon[thread_id].size() * 3 <= bls_per_bp[thread_id].size());

                    if (strand == '+') {
                        for (size_t xx = 0, bls_pos = aln.skip_bases; xx < lpr_per_codon[thread_id].size(); ++xx, bls_pos += 3)
                        {
                            const float bls_codon_sum = bls_per_bp[thread_id][bls_pos]
                                                        + bls_per_bp[thread_id][bls_pos + 1]
                                                        + bls_per_bp[thread_id][bls_pos + 2];
                            if (bls_codon_sum < params.phylo_threshold * 3)
                            {
                                if (params.phylo_smooth && scores.empty())
                                    startBlockPos = aln.start_pos + ((xx + 1) * 3);
                                continue;
                            }

                            int64_t newPos = aln.start_pos + (xx * 3);

                            if (prevPos + 3 != newPos)
                            {
                                if (params.phylo_raw)
                                    fprintf(file_score_raw[file_index], "fixedStep chrom=%s start=%" PRId64 " step=3 span=3\n", aln.chrom.c_str(), newPos);

                                if (params.phylo_smooth && !scores.empty())
                                {
                                    if (!params.smooth_wo_wig) {
                                        fprintf(file_score[file_index],
                                                "fixedStep chrom=%s start=%" PRId64 " step=3 span=3\n", aln.chrom.c_str(),
                                                startBlockPos);
                                    }

                                    process_scores(model._hmm, scores, startBlockPos, region, bedregions);

                                    if (!params.smooth_wo_wig) {
                                        for (size_t i = 0; i < region.size(); i++) {
                                            my_fprintf(file_score[file_index], "%.3f", region[i].log_odds_prob);
                                        }
                                    }
                                    if (params.phylo_bed_file) {
                                        for (size_t i = 0; i < bedregions.size(); i++) {
                                            fprintf(file_score_bed[file_index],
                                                    "%s\t%" PRIu32 "\t%" PRIu32 "\t%s:%" PRIu32 "-%" PRIu32 "\t0\t+\t%" PRIu32 "\t%" PRIu32 "\t%" PRIu32 ",%" PRIu32 ",%" PRIu32 "\n",
                                                    aln.chrom.c_str(),
                                                    bedregions[i].region_start, bedregions[i].region_end, aln.chrom.c_str(),
                                                    bedregions[i].region_start + 1,
                                                    bedregions[i].region_end, bedregions[i].region_start,
                                                    bedregions[i].region_end,
                                                    bedregions[i].color, bedregions[i].color, bedregions[i].color);
                                        }
                                    }

                                    scores.clear();
                                    region.clear();
                                    bedregions.clear();
                                    startBlockPos = aln.start_pos + (xx * 3);
                                }
                            }

                            prevPos = newPos;

                            if (params.phylo_raw)
                                my_fprintf(file_score_raw[file_index], "%.3f", lpr_per_codon[thread_id][xx]);

                            if (params.phylo_smooth)
                                scores.push_back(lpr_per_codon[thread_id][xx]);
                        }
                    }

                    if (strand == '-') {
                        for (size_t xx = 0, bls_pos = aln.length() % 3; xx < lpr_per_codon[thread_id].size(); ++xx, bls_pos += 3)
                        {
                            const float bls_codon_sum = bls_per_bp[thread_id][bls_pos]
                                                        + bls_per_bp[thread_id][bls_pos + 1]
                                                        + bls_per_bp[thread_id][bls_pos + 2];
                            if (bls_codon_sum < params.phylo_threshold * 3)
                            {
                                if (params.phylo_smooth && scores.empty())
                                    startBlockPos = aln.start_pos + ((xx + 1) * 3);
                                continue;
                            }

                            int64_t newPos = aln.start_pos + (xx * 3);

                            if (prevPos + 3 != newPos)
                            {
                                if (params.phylo_raw)
                                    fprintf(file_score_raw[file_index], "fixedStep chrom=%s start=%" PRId64 " step=3 span=3\n", aln.chrom.c_str(), newPos);

                                if (params.phylo_smooth && !scores.empty()) {
                                    if (!params.smooth_wo_wig) {
                                        fprintf(file_score[file_index],
                                                "fixedStep chrom=%s start=%" PRId64 " step=3 span=3\n", aln.chrom.c_str(),
                                                startBlockPos);
                                    }

                                    process_scores(model._hmm, scores, startBlockPos, region, bedregions);

                                    if (!params.smooth_wo_wig) {
                                        for (size_t i = 0; i < region.size(); i++) {
                                            my_fprintf(file_score[file_index], "%.3f", region[i].log_odds_prob);
                                        }
                                    }
                                    if (params.phylo_bed_file) {
                                        for (size_t i = 0; i < bedregions.size(); i++) {
                                            fprintf(file_score_bed[file_index],
                                                    "%s\t%" PRIu32 "\t%" PRIu32 "\t%s:%" PRIu32 "-%" PRIu32 "\t0\t-\t%" PRIu32 "\t%" PRIu32 "\t%" PRIu32 ",%" PRIu32 ",%" PRIu32 "\n",
                                                    aln.chrom.c_str(),
                                                    bedregions[i].region_start, bedregions[i].region_end, aln.chrom.c_str(),
                                                    bedregions[i].region_start + 1,
                                                    bedregions[i].region_end, bedregions[i].region_start,
                                                    bedregions[i].region_end,
                                                    bedregions[i].color, bedregions[i].color, bedregions[i].color);
                                        }
                                    }

                                    scores.clear();
                                    region.clear();
                                    bedregions.clear();
                                    startBlockPos = aln.start_pos + (xx * 3);
                                }
                            }

                            prevPos = newPos;

                            if (params.phylo_raw)
                                my_fprintf(file_score_raw[file_index], "%.3f", lpr_per_codon[thread_id][xx]);

                            if (params.phylo_smooth)
                                scores.push_back(lpr_per_codon[thread_id][xx]);
                        }
                    }

                    if (params.phylo_smooth && !scores.empty())
                    {
                        if (!params.smooth_wo_wig) {
                            fprintf(file_score[file_index],
                                    "fixedStep chrom=%s start=%" PRId64 " step=3 span=3\n", aln.chrom.c_str(),
                                    startBlockPos);
                        }

                        process_scores(model._hmm, scores, startBlockPos, region, bedregions);

                        if (!params.smooth_wo_wig) {
                            for (size_t i = 0; i < region.size(); i++) {
                                my_fprintf(file_score[file_index], "%.3f", region[i].log_odds_prob);
                            }
                        }
                        if (params.phylo_bed_file) {
                            if (strand == '+') {
                                for(size_t i = 0; i < bedregions.size(); i++) {
                                    fprintf(file_score_bed[file_index],
                                            "%s\t%" PRIu32 "\t%" PRIu32 "\t%s:%" PRIu32 "-%" PRIu32 "\t0\t+\t%" PRIu32 "\t%" PRIu32 "\t%" PRIu32 ",%" PRIu32 ",%" PRIu32 "\n",
                                            aln.chrom.c_str(),
                                            bedregions[i].region_start, bedregions[i].region_end, aln.chrom.c_str(),
                                            bedregions[i].region_start + 1,
                                            bedregions[i].region_end, bedregions[i].region_start,
                                            bedregions[i].region_end,
                                            bedregions[i].color, bedregions[i].color, bedregions[i].color);
                                }
                            } else {
                                for (size_t i = 0; i < bedregions.size(); i++) {
                                    fprintf(file_score_bed[file_index],
                                            "%s\t%" PRIu32 "\t%" PRIu32 "\t%s:%" PRIu32 "-%" PRIu32 "\t0\t-\t%" PRIu32 "\t%" PRIu32 "\t%" PRIu32 ",%" PRIu32 ",%" PRIu32 "\n",
                                            aln.chrom.c_str(),
                                            bedregions[i].region_start, bedregions[i].region_end, aln.chrom.c_str(),
                                            bedregions[i].region_start + 1,
                                            bedregions[i].region_end, bedregions[i].region_start,
                                            bedregions[i].region_end,
                                            bedregions[i].color, bedregions[i].color, bedregions[i].color);
                                }
                            }
                        }
                        scores.clear();
                        region.clear();
                        bedregions.clear();

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
        maf_rd.print_progress();

        if (params.phylo_power)
            fclose(file_power);

        for (uint8_t i = 0; i < 6; ++i)
        {
            if (params.phylo_raw)
                fclose(file_score_raw[i]);

            if (params.phylo_smooth){
                fclose(file_score[i]);
                fclose(file_score_bed[i]);
            }
        }
    }

//    printf("\x1b[KMerging temporary output files ...\r");

    // merge power file
    if (params.phylo_power)
        merge_job_output_files(output_folder + "/PhyloCSFpower.wig", jobs, file_id > 1);

    // merge 6 frame files for scores and raw scores
    for (char strand = '+'; strand <= '-'; strand += 2)
    {
        for (unsigned frame = 1; frame <= 3; ++frame)
        {
            if (params.phylo_raw)
                merge_job_output_files(output_folder + "/PhyloCSFRaw" + std::string(1, strand) + std::to_string(frame) + ".wig", jobs, file_id > 1);

            if (params.phylo_smooth) {
                merge_job_output_files(
                        output_folder + "/PhyloCSF" + std::string(1, strand) + std::to_string(frame) + ".wig", jobs,
                        file_id > 1);
                merge_job_output_files(
                        output_folder + "/PhyloCSF" + std::string(1, strand) + std::to_string(frame) + "Regions.bed", jobs,
                        file_id > 1);
                if (!params.phylo_bed_file) {
                    std::string bed_file_name = output_folder + "/PhyloCSF" + std::string(1, strand) + std::to_string(frame) + "Regions.bed";
                    const char *address = bed_file_name.c_str();
                    remove(address);
                }
                if (params.smooth_wo_wig) {
                    std::string wig_file_name = output_folder + "/PhyloCSF" + std::string(1, strand) + std::to_string(frame) + ".wig";
                    const char *address_wig = wig_file_name.c_str();
                    remove(address_wig);
                }
            }
        }
    }
}

int main_build_tracks(int argc, char **argv)
{
    ArgParse args("phylocsf++ build-tracks",
                  "Computes PhyloCSF tracks for each codon and all 6 frames from alignment files in MAF \n"
                  "format as well as the power track containing the branch length scores (confidence of \n"
                  "the PhyloCSF scores). Optionally the PhyloCSF scores can be smoothened and posterior \n"
                  "probabilities can be computed. Outputs wig files.");

    TrackCLIParams params;

    const std::string model_list = get_list_of_models();

    char threshold_default_str[10];
    sprintf(threshold_default_str, "%.1f", params.phylo_threshold);

    args.add_option("output-raw-phylo", ArgParse::Type::BOOL, "Compute all 6 unsmoothened PhyloCSF tracks. Default: " + std::to_string(params.phylo_raw), ArgParse::Level::GENERAL, false);
    args.add_option("output-phylo", ArgParse::Type::BOOL, "Compute all 6 smoothened PhyloCSF tracks. Requires coding exon coordinates and genome length. Default: " + std::to_string(params.phylo_smooth), ArgParse::Level::GENERAL, false);
    args.add_option("output-regions", ArgParse::Type::BOOL, "Generate bed files with coordinates of potential protein coding regions. Requires coding exon coordinates and genome length. Default: " + std::to_string(params.phylo_bed_file), ArgParse::Level::GENERAL,false);
    args.add_option("power-threshold", ArgParse::Type::FLOAT, "Minimum confidence score to output PhyloCSF score. Default: " + std::string(threshold_default_str), ArgParse::Level::GENERAL, false);
    args.add_option("genome-length", ArgParse::Type::INT, "Total genome length (needed for --output-phylo).", ArgParse::Level::GENERAL, false);
    args.add_option("coding-exons", ArgParse::Type::STRING, "BED-like file (chrom name, strand, phase, start, end) with coordinates of coding exons (needed for --output-phylo).", ArgParse::Level::GENERAL, false);
    args.add_option("threads", ArgParse::Type::INT, "Parallelize scoring of multiple MSAs in a file using multiple threads. Default: " + std::to_string(params.threads), ArgParse::Level::GENERAL, false);
    args.add_option("output", ArgParse::Type::STRING, "Directory where tracks in wig format will be written to. If it does not exist, it will be created. Default: output files are stored next to the input files.", ArgParse::Level::GENERAL, false);
    args.add_option("mapping", ArgParse::Type::STRING, "If the MSAs don't use common species names (like Human, Chimp, etc.) you can pass a two-column tsv file with a name mapping.", ArgParse::Level::GENERAL, false);
    args.add_option("species", ArgParse::Type::STRING, "Comma separated list of species to reduce <model> to a subset of species to improve running time, e.g., \"Human,Chimp,Seaturtle\"", ArgParse::Level::GENERAL, false);
    args.add_option("model-info", ArgParse::Type::STRING, "Output the organisms included in a specific model. Included models are: " + model_list + ".", ArgParse::Level::HELP, false);
    // TODO: add dry-run option (to check whether all species names can be mapped and it doesn't fail after 2 days of runtime)
    // or even better: continue where the tool left off!

    args.add_positional_argument("model", ArgParse::Type::STRING, "Path to parameter files, or one of the following predefined models: " + model_list + ".", true);
    args.add_positional_argument("alignments", ArgParse::Type::STRING, "One or more files containing multiple sequence alignments. Formats: MAF and multi FASTA. Multiple MSAs can be stored in a single file separated by empty lines.", true, true);
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
        params.phylo_threshold = args.get_bool("power-threshold");
    if (args.is_set("output-regions"))
        params.phylo_bed_file = args.get_bool("output-regions");

    // when --output-regions is true and --output-phylo is false, we have to give only bed files, not wig files
    if (params.phylo_bed_file && !params.phylo_smooth) {
        params.phylo_smooth = true;
        params.smooth_wo_wig = true;
    }

    // this has to hold true: phylo_smooth => (args.is_set("genome-length") && args.is_set("coding-exons"))
    if (params.phylo_smooth && (!args.is_set("genome-length") || !args.is_set("coding-exons")))
    {
        printf(OUT_ERROR "For smoothened tracks (--output-phylo) you need to provide --genome-length and --coding-exons.\n" OUT_RESET);
        return -1;
    }

    // this has to hold true: output-regions => (args.is_set("genome-length") && args.is_set("coding-exons"))
    if (params.phylo_bed_file && (!args.is_set("genome-length") || !args.is_set("coding-exons")))
    {
        printf(OUT_ERROR "To generate bed file of potential protein coding regions, you need to provide --genome-length and --coding-exons.\n" OUT_RESET);
        return -1;
    }

    if (args.is_set("output"))
    {
        params.output_path = args.get_string("output");
        // create a directory if it doesn't exist yet
        if (create_directory(params.output_path))
            printf(OUT_INFO "Created the output directory.\n" OUT_RESET);
    }

    if (args.is_set("mapping"))
    {
        std::string mapping_file = args.get_string("mapping");
        update_sequence_name_mapping(mapping_file);
    }

    if (!args.is_set_positional("model") || !args.is_set_positional("alignments"))
    {
        printf(OUT_ERROR "No model or alignments provided.\n" OUT_RESET);
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

    std::string species = "";
    if (args.is_set("species"))
        species = args.get_string("species");

    Model model;
    load_model(model, args.get_positional_argument("model"), species,
               params.phylo_smooth, genome_length, coding_exons_path);

    // see whether there are species in the model that did not occur in any
    std::vector<std::vector<bool> > species_seen_in_alignment;
    species_seen_in_alignment.resize(params.threads);
    const uint16_t leaves = (model.phylo_array.size() + 1)/2;
    for (uint32_t i = 0; i < params.threads; ++i)
        species_seen_in_alignment[i].resize(leaves, false);

    // run for every alignment file
    for (uint16_t i = 1; i < args.positional_argument_size(); ++i)
    {
        run_tracks(args.get_positional_argument(i).c_str(), model, params, i, args.positional_argument_size() - 1, species_seen_in_alignment);
    }

    printf("Done!\n");

    // OR all bit vectors and store them in the first (0th)
    for (uint32_t t = 1; t < params.threads; ++t)
    {
        for (uint32_t pos = 0; pos < species_seen_in_alignment[0].size(); ++pos)
        {
            species_seen_in_alignment[0][pos] = species_seen_in_alignment[0][pos] | species_seen_in_alignment[t][pos];
        }
    }

    for (uint32_t pos = 0; pos < species_seen_in_alignment[0].size(); ++pos)
    {
        if (!species_seen_in_alignment[0][pos])
        {
            printf(OUT_INFO "WARNING: %s in the model does not occur in alignment file(s). Check --species to select a subset (this affects the power/confidence track).\n" OUT_RESET, model.phylo_array[pos].label.c_str());
        }
    }

    return 0;
}
