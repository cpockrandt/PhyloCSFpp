#include "common.hpp"

#include "arg_parse.hpp"
#include "models.hpp"

#include "run.hpp"

#include <omp.h>

struct ScoreMSACLIParams
{
    unsigned threads = omp_get_max_threads();
    bool comp_phylo = true;
    bool comp_anc = false;
    bool comp_bls = true;
    std::string output_path = "";
    algorithm_t strategy = MLE;

    // internal flag for annotate-with-msa
    // score every codon separately, flatten and compute the average
    bool avg_codon_score = false;
    bool avg_codon_score_smoothing = false;
    bool avg_codon_score_smoothing_threshold = false;
};

struct scoring_result
{
    std::string seq;
    uint64_t start;
    uint64_t end;
    char strand;
    float phylo;
    float anc;
    float bls;

    scoring_result(const std::string & seq, const uint64_t start, const uint64_t end, const char strand,
                   const float phylo, const float anc, const float bls)
        : seq(seq), start(start), end(end), strand(strand), phylo(phylo), anc(anc), bls(bls) { }
};

void run_scoring_msa(const std::string & alignment_path, const Model & model, const ScoreMSACLIParams & params,
                     const uint32_t file_id, const uint32_t files)
{
    unsigned jobs = (params.threads > 1) ? params.threads * 10 : 1;

    std::vector<Data> data(params.threads);

    // store output next to alignment file if no output_path was specified
    std::string output_file_path = params.output_path;
    if (output_file_path == "")
    {
        output_file_path = alignment_path + ".scores";
    }
    else
    {
        size_t pos = alignment_path.find_last_of('/');
        if (pos != std::string::npos) // alignment_path == ".../aln.maf"
            output_file_path += "/" + alignment_path.substr(pos + 1) + ".scores";
        else                          // alignment_path == "aln.maf"
            output_file_path += "/" + alignment_path + ".scores";
    }

    // open file for writing so the user gets an error message as soon as possible
    FILE *output_file = fopen(output_file_path.c_str(), "w");
    if (output_file == NULL)
    {
        printf("Error creating file!");
        exit(1);
    }

    // write header
    fprintf(output_file, "# PhyloCSF scores computed with PhyloCSF++ %s (%s, %s)\n", ArgParse::version.c_str(), ArgParse::git_hash.c_str(), ArgParse::git_date.c_str());

    fprintf(output_file, "seq\tstart\tend\tstrand");
    if (params.comp_phylo)
        fprintf(output_file, "\tphylocsf-score");
    if (params.comp_anc)
        fprintf(output_file, "\tanc-score");
    if (params.comp_bls)
        fprintf(output_file, "\tblc-score");
    fprintf(output_file, "\n");

    // prepare alignment
    std::vector<alignment_t> alignments;
    const uint16_t leaves = (model.phylo_array.size() + 1)/2;
    for (unsigned i = 0; i < params.threads; ++i)
    {
        alignments.emplace_back(leaves);
    }

    parallel_maf_reader maf_rd(alignment_path.c_str(), jobs, &model.seqid_to_phyloid, false);
    jobs = maf_rd.get_jobs(); // maybe file is too small and a smaller number of jobs is used
    maf_rd.setup_progressbar(file_id, files);

    std::vector<std::vector<scoring_result> > all_results(jobs);

    std::vector<std::vector<double> > lpr_per_codon(params.threads), bls_per_bp(params.threads);

    #pragma omp parallel for num_threads(params.threads) schedule(dynamic, 1)
    for (unsigned job_id = 0; job_id < jobs; ++job_id)
    {
        std::mt19937 gen;
        unsigned thread_id = omp_get_thread_num();
        alignment_t & aln = alignments[thread_id];
        std::vector<scoring_result> & results = all_results[job_id];

        while (maf_rd.get_next_alignment(aln, job_id))
        {
            maf_rd.print_progress();

            aln.translate();

            results.emplace_back(aln.chrom, aln.start_pos, aln.start_pos + aln.length() - 1, aln.strand, NAN, NAN, NAN);

            // TODO: only compute phylocsf score and/or anc score if necessary
            if (!params.avg_codon_score)
            {
                if (params.comp_phylo || params.comp_anc)
                {
                    try {
                        gen.seed(42); // NOTE: make sure that parallelization does not affect resultss
                        std::tuple<float, float> result = run(data[thread_id], model, aln, params.strategy, gen);

                        if (params.comp_phylo)
                            results.back().phylo = std::get<0>(result);

                        if (params.comp_anc)
                            results.back().anc = std::get<1>(result);
                    }
                    catch (const std::runtime_error & e)
                    {
                        // printf("%s\n", e.what());
                    }
                }

                if (params.comp_bls)
                {
                    const float bls_score = compute_bls_score<false>(model.phylo_tree, aln, model, bls_per_bp[thread_id]);
                    results.back().bls = bls_score;
                }
            }
            else
            {
                assert(params.strategy == FIXED);

                bls_per_bp[thread_id].clear();
                const float bls_score = compute_bls_score<true>(model.phylo_tree, aln, model, bls_per_bp[thread_id]);
                results.back().bls = bls_score;

                std::vector<double> scores;
                std::vector<scored_region> region;

                lpr_per_codon[thread_id].clear();
                run_tracks(data[thread_id], model, aln, lpr_per_codon[thread_id]);
                data[thread_id].clear();

                float sum_smoothened_scores = 0.0;
                uint64_t cnt_smoothened_scores = 0;

                float sum_all_lpr = 0.0;

                int64_t prevPos = -4;
                int64_t startBlockPos = aln.start_pos;
                for (size_t xx = 0, bls_pos = 0; xx < lpr_per_codon[thread_id].size(); ++xx, bls_pos += 3)
                {
                    sum_all_lpr += lpr_per_codon[thread_id][xx];
                    const float bls_codon_sum = bls_per_bp[thread_id][bls_pos]
                                              + bls_per_bp[thread_id][bls_pos + 1]
                                              + bls_per_bp[thread_id][bls_pos + 2];
                    if (params.avg_codon_score_smoothing_threshold && bls_codon_sum < 0.1 * 3) // TODO: document this fixed threshold
                    {
                        if (scores.empty())
                            startBlockPos = aln.start_pos + ((xx + 1) * 3);
                        continue;
                    }

                    int64_t newPos = aln.start_pos + (xx * 3);

                    if (prevPos + 3 != newPos)
                    {
                        if (!scores.empty())
                        {
                            // fprintf(file_score[file_index], "fixedStep chrom=%s start=%ld step=3 span=3\n", aln.chrom.c_str(), startBlockPos);
                            process_scores(model._hmm, scores, startBlockPos, region, SCORE_CODON);

                            for(size_t i = 0; i < region.size(); i++)
                                sum_smoothened_scores += region[i].log_odds_prob;
                                //my_fprintf(file_score[file_index], "%.3f", region[i].log_odds_prob);
                            cnt_smoothened_scores += region.size();

                            scores.clear();
                            region.clear();
                            startBlockPos = aln.start_pos + (xx * 3);
                        }
                    }

                    prevPos = newPos;

                    scores.push_back(lpr_per_codon[thread_id][xx]);
                }

                if (!scores.empty())
                {
                    //fprintf(file_score[file_index], "fixedStep chrom=%s start=%ld step=3 span=3\n", aln.chrom.c_str(), startBlockPos);
                    process_scores(model._hmm, scores, startBlockPos, region, SCORE_CODON);

                    for(size_t i = 0; i < region.size(); i++)
                        sum_smoothened_scores += region[i].log_odds_prob;
                        //my_fprintf(file_score[file_index], "%.3f", region[i].log_odds_prob);
                    cnt_smoothened_scores += region.size();

                    scores.clear();
                    region.clear();
                }

                results.back().phylo = sum_smoothened_scores / cnt_smoothened_scores;

                // if we don't want smoothing, just overwrite the value with the unsmoothened avg value
                if (!params.avg_codon_score_smoothing)
                    results.back().phylo = sum_all_lpr / lpr_per_codon[thread_id].size();
            }

            data[thread_id].clear();

            for (auto & seq : aln.seqs)
                seq = "";
        }
    }

    if (files == 1)
        printf("Writing scores to file ...\r");
    else
        printf("File %d of %d: Writing scores to file ...\r", file_id, files);

    // merge and write output
    for (const auto & job_results : all_results)
    {
        for (const auto & result : job_results)
        {
            fprintf(output_file, "%s\t%ld\t%ld\t%c", result.seq.c_str(), result.start, result.end, result.strand);
            if (params.comp_phylo)
                fprintf(output_file, "\t%.6f", result.phylo);
            if (params.comp_anc)
                fprintf(output_file, "\t%.6f", result.anc);
            if (params.comp_bls)
                fprintf(output_file, "\t%.6f", result.bls);
            fprintf(output_file, "\n");
        }
    }
    fclose(output_file);
}

int main_score_msa(int argc, char **argv)
{
    ArgParse args("phylocsf++ score-msa",
                  "Computes PhyloCSF scores for whole alignments in FASTA or MAF files. Output is \n"
                  "written to bed file(s). Other scores such as the ancestral sequence composition \n"
                  "sores and branch length scores can be computed as well. Only one forward frame \n"
                  "is computed, i.e., for other frames reverse the alignments and/or remove the \n"
                  "first one or two bases.");

    ScoreMSACLIParams params;

    const std::string model_list = get_list_of_models();

    std::string default_strategy;
    switch (params.strategy)
    {
        case MLE:
            default_strategy = "MLE";
            break;
        case FIXED:
            default_strategy = "FIXED";
            break;
        case OMEGA:
            default_strategy = "OMEGA";
            break;
    }

    args.add_option("strategy", ArgParse::Type::STRING, "PhyloCSF scoring algorithm: MLE, FIXED or OMEGA. Default: " + default_strategy, ArgParse::Level::GENERAL, false);
    args.add_option("comp-phylo", ArgParse::Type::BOOL, "Compute the PhyloCSF score for each alignment. Default: " + std::to_string(params.comp_phylo), ArgParse::Level::GENERAL, false);
    args.add_option("comp-anc", ArgParse::Type::BOOL, "Compute the ancestral sequence composition score (only in MLE and FIXED mode). Default: " + std::to_string(params.comp_anc), ArgParse::Level::GENERAL, false);
    args.add_option("comp-bls", ArgParse::Type::BOOL, "Compute the branch length score (confidence). Default: " + std::to_string(params.comp_bls), ArgParse::Level::GENERAL, false);

    args.add_option("score-codons", ArgParse::Type::BOOL, "Output confidence score (branch length score). Default: " + std::to_string(params.avg_codon_score), ArgParse::Level::GENERAL, false);
    args.add_option("score-codons-smooth", ArgParse::Type::BOOL, "Output confidence score (branch length score). Default: " + std::to_string(params.avg_codon_score_smoothing), ArgParse::Level::GENERAL, false);
    args.add_option("score-codons-smooth-with-threshold", ArgParse::Type::BOOL, "Output confidence score (branch length score). Default: " + std::to_string(params.avg_codon_score_smoothing_threshold), ArgParse::Level::GENERAL, false);

    args.add_option("threads", ArgParse::Type::INT, "Parallelize scoring of multiple MSAs in a file using multiple threads. Default: " + std::to_string(params.threads), ArgParse::Level::GENERAL, false);
    args.add_option("output", ArgParse::Type::STRING, "Path where tracks in wig format will be written. If it does not exist, it will be created. Default: output files are stored next to the input files.", ArgParse::Level::GENERAL, false);

    args.add_option("mapping", ArgParse::Type::STRING, "If the MSAs don't use common species names (like Human, Chimp, etc.) you can pass a two-column tsv file with a name mapping.", ArgParse::Level::GENERAL, false);
    args.add_option("species", ArgParse::Type::STRING, "Comma separated list of species to reduce <model> to a subset of species to improve running time, e.g., \"Human,Chimp,Seaturtle\"", ArgParse::Level::GENERAL, false);
    args.add_option("model-info", ArgParse::Type::STRING, "Output the organisms included in a specific model. Included models are: " + model_list + ".", ArgParse::Level::GENERAL, false);
    // TODO: ignore sequences not occuring in model (instead of failing)

    args.add_option("genome-length", ArgParse::Type::INT, "Total genome length (needed for --output-phylo).", ArgParse::Level::GENERAL, false);
    args.add_option("coding-exons", ArgParse::Type::STRING, "BED-like file (chrom name, strand, phase, start, end) with coordinates of coding exons (needed for --output-phylo).", ArgParse::Level::GENERAL, false);

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

    if (args.is_set("score-codons") && args.get_bool("score-codons"))
    {
        params.avg_codon_score = true;
        params.avg_codon_score_smoothing = false;
        params.avg_codon_score_smoothing_threshold = false;
    }
    if (args.is_set("score-codons-smooth") && args.get_bool("score-codons-smooth"))
    {
        params.avg_codon_score = true;
        params.avg_codon_score_smoothing = true;
        params.avg_codon_score_smoothing_threshold = false;
    }
    if (args.is_set("score-codons-smooth-with-threshold") && args.get_bool("score-codons-smooth-with-threshold"))
    {
        params.avg_codon_score = true;
        params.avg_codon_score_smoothing = true;
        params.avg_codon_score_smoothing_threshold = true;
    }

    if (args.is_set("strategy"))
    {
        std::string strategy = args.get_string("strategy");
        std::transform(strategy.begin(), strategy.end(), strategy.begin(), toupper);

        if (strategy == "MLE")
            params.strategy = MLE;
        else if (strategy == "FIXED")
            params.strategy = FIXED;
        else if (strategy == "OMEGA")
            params.strategy = OMEGA;
        else
        {
            printf(OUT_ERROR "Please choose a valid strategy (MLE, FIXED or OMEGA)!\n" OUT_RESET);
            return -1;
        }
    }

    if (args.is_set("comp-phylo"))
        params.comp_phylo = args.get_bool("comp-phylo");
    if (args.is_set("comp-anc"))
        params.comp_anc = args.get_bool("comp-anc");
    if (args.is_set("comp-bls"))
        params.comp_bls = args.get_bool("comp-bls");

    if (params.strategy == OMEGA && params.comp_anc)
    {
        printf(OUT_ERROR "The ancestral sequence composition cannot be computed in the Omega mode!\n" OUT_RESET);
        return -1;
    }

    if (!params.comp_phylo && !params.comp_anc && !params.comp_bls)
    {
        printf(OUT_ERROR "At least one score needs to be computed (phylo, anc or bls)!\n" OUT_RESET);
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

    // load and prepare model
    uint64_t genome_length = 0;
    if (args.is_set("genome-length"))
        genome_length = args.get_int("genome-length");

    std::string coding_exons_path = "";
    if (args.is_set("coding-exons"))
        coding_exons_path = args.get_string("coding-exons");

    // TODO: only pass "true" if coding-exons and genome-length are set

    Model model;
    load_model(model, args.get_positional_argument("model"), args.get_string("species"), true/*false*/, genome_length, coding_exons_path);

    // run for every alignment file
    for (uint16_t i = 1; i < args.positional_argument_size(); ++i)
    {
         run_scoring_msa(args.get_positional_argument(i), model, params, i, args.positional_argument_size() - 1);
    }

    printf(OUT_DEL "Done!\n");

    return 0;
}
