#include "arg_parse.hpp"
#include "models.hpp"

#include "run.hpp"

#include <omp.h>

struct RegionCLIParams
{
    unsigned threads = omp_get_max_threads();
    bool comp_phylo = true;
    bool comp_anc = false;
    bool comp_bls = false;
    std::string output_path = "";
    algorithm_t strategy = MLE;
};

struct region_result
{
    std::string seq;
    uint64_t start;
    uint64_t end;
    char strand;
    float phylo;
    float anc;
    float bls;

    region_result(const std::string & seq, const uint64_t start, const uint64_t end, const char strand,
                  const float phylo, const float anc, const float bls)
        : seq(seq), start(start), end(end), strand(strand), phylo(phylo), anc(anc), bls(bls) { }
};

void run_regions(const std::string & alignment_path, const Model & model, const RegionCLIParams & params)
{
    unsigned jobs = (params.threads > 1) ? params.threads * 10 : 1;

    std::vector<Data> data(params.threads);

    // store output next to alignment file if no output_path was specified
    std::string output_file_path = params.output_path;
    if (output_file_path == "")
    {
        // TODO: could overwrite a file created prior by the user
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
    fprintf(output_file, "seq\tstart\tend");
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

    std::vector<std::vector<region_result> > all_results(jobs);

    #pragma omp parallel for num_threads(params.threads) schedule(dynamic, 1)
    for (unsigned job_id = 0; job_id < jobs; ++job_id)
    {
        unsigned thread_id = omp_get_thread_num();
        alignment_t & aln = alignments[thread_id];
        std::vector<region_result> & results = all_results[job_id];

        while (maf_rd.get_next_alignment(aln, job_id))
        {
            aln.translate();

            printf("%s\t%ld\t%ld\t%c", aln.chrom.c_str(), aln.start_pos, aln.start_pos + aln.length() - 1, aln.strand);
            results.emplace_back(aln.chrom, aln.start_pos, aln.start_pos + aln.length() - 1, aln.strand, NAN, NAN, NAN);

            // TODO: only compute phylocsf score and/or anc score if necessary
            if (params.comp_phylo || params.comp_anc)
            {
                try {
                    std::tuple<float, float> result = run(data[thread_id], model, aln, params.strategy);

                    if (params.comp_phylo)
                    {
                        results.back().phylo = std::get<0>(result);
                        printf("\t%.6f", std::get<0>(result));
                    }

                    if (params.comp_anc)
                    {
                        results.back().anc = std::get<1>(result);
                        printf("\t%.6f", std::get<1>(result));
                    }
                }
                catch (const std::runtime_error &)
                {
                    if (params.comp_phylo)
                        printf("\t%.6f", NAN);
                }
            }

            if (params.comp_bls)
            {
                std::vector<double> bls_scores_per_base;
                const float bls_score = compute_bls_score<false>(model.phylo_tree, aln, model, bls_scores_per_base);
                results.back().bls = bls_score;
                printf("\t%.6f", bls_score);
            }

            printf("\n");

            data[thread_id].clear();

            for (auto & seq : aln.seqs)
                seq = "";
        }
    }

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

int main_region(int argc, char **argv)
{
    ArgParse args("phylocsf++ region",
                  "Computes PhyloCSF scores for whole alignments in FASTA or MAF files. Output is \n"
                  "written to bed file(s). Other scores such as the ancestral sequence composition \n"
                  "sores and branch length scores can be computed as well. Only one forward frame \n"
                  "is computed, i.e., for other frames reverse the alignments and/or remove the \n"
                  "first one or two bases."); // For ORF finding and scoring of transcripts in GTF/GFF "
                  // "files, check out PhyloCSF++Anno: https://github.com/cpockrandt/phylocsfppanno");

    RegionCLIParams params;

    std::string model_list = get_list_of_models();

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

    args.add_option("threads", ArgParse::Type::INT, "Parallelize scoring of multiple MSAs in a file using multiple threads. Default: " + std::to_string(params.threads), ArgParse::Level::GENERAL, false);
    args.add_option("output", ArgParse::Type::STRING, "Path where tracks in wig format will be written. If it does not exist, it will be created. Default: output files are stored next to the input files.", ArgParse::Level::GENERAL, false);

    args.add_option("mapping", ArgParse::Type::STRING, "If the MSAs don't use common species names (like Human, Chimp, etc.) you can pass a two-column tsv file with a name mapping.", ArgParse::Level::GENERAL, false);
    args.add_option("species", ArgParse::Type::STRING, "Comma separated list of species to reduce <model> to a subset of species to improve running time, e.g., \"Human,Chimp,Seaturtle\"", ArgParse::Level::GENERAL, false);
    args.add_option("model-info", ArgParse::Type::STRING, "Output the organisms included in a specific model. Included models are: " + model_list + ".", ArgParse::Level::GENERAL, false);
    // TODO: ignore sequences not occuring in model (instead of failing)

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
            print_error_msg("Please choose a valid strategy (MLE, FIXED or OMEGA)");
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
        print_error_msg("The ancestral sequence composition cannot be computed in the Omega mode!");
        return -1;
    }

    if (!params.comp_phylo && !params.comp_anc && !params.comp_bls)
    {
        print_error_msg("At least one score needs to be computed (phylo, anc or bls)!");
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
    Model model;
    load_model(model, args.get_positional_argument("model"), args.get_string("species"), false, 0, "");

    // run for every alignment file
    for (uint16_t i = 1; i < args.positional_argument_size(); ++i)
    {
        // printf("%d, %s\n", i, args.get_positional_argument(i).c_str());
         run_regions(args.get_positional_argument(i), model, params);
    }

    return 0;
}
