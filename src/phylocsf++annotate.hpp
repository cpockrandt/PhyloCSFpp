#include "common.hpp"

#include "arg_parse.hpp"
#include "models.hpp"

#include "run.hpp"

#include <omp.h>


#include <stdio.h>
#include <sys/wait.h>

//void run_annotate(const std::string & alignment_path, const Model & model, const AnnotateCLIParams & params)
//{
//
//}

int main_annotate(int /*argc*/, char** /*argv*/)
{
    ArgParse args("phylocsf++ annotate",
                  "Computes PhyloCSF scores for CDS features in GFF/GTF files and outputs them in \n"
                  "the same format. You can either pass track files (in wig format) or have the \n"
                  "multiple sequence alignment computed for the CDS (this requires mmseqs2 to be \n"
                  "installed).");

    // check whether mmseqs is installed
    int ret = system("mmseqs --help > /dev/null 2>&1");
    if (WEXITSTATUS(ret))
    {
        print_error_msg("ERROR: MMseqs2 seems not to be installed.\n");
        exit(-1);
    }

    ret = system("mmseqs result2dnamsa --help > /dev/null 2>&1");
    if (WEXITSTATUS(ret))
    {
        print_error_msg("ERROR: A version of MMseqs2 seems to be installed that is too old. `mmseqs result2dnamsa --help` does not seem to work.\n");
        exit(-1);
    }

//    AnnotateCLIParams params;
//
//    std::string model_list = get_list_of_models();
//
//    std::string default_strategy;
//    switch (params.strategy)
//    {
//        case MLE:
//            default_strategy = "MLE";
//            break;
//        case FIXED:
//            default_strategy = "FIXED";
//            break;
//        case OMEGA:
//            default_strategy = "OMEGA";
//            break;
//    }
//
//    args.add_option("strategy", ArgParse::Type::STRING, "PhyloCSF scoring algorithm: MLE, FIXED or OMEGA. Default: " + default_strategy, ArgParse::Level::GENERAL, false);
//    args.add_option("comp-phylo", ArgParse::Type::BOOL, "Compute the PhyloCSF score for each alignment. Default: " + std::to_string(params.comp_phylo), ArgParse::Level::GENERAL, false);
//    args.add_option("comp-anc", ArgParse::Type::BOOL, "Compute the ancestral sequence composition score (only in MLE and FIXED mode). Default: " + std::to_string(params.comp_anc), ArgParse::Level::GENERAL, false);
//    args.add_option("comp-bls", ArgParse::Type::BOOL, "Compute the branch length score (confidence). Default: " + std::to_string(params.comp_bls), ArgParse::Level::GENERAL, false);
//
//    args.add_option("threads", ArgParse::Type::INT, "Parallelize scoring of multiple MSAs in a file using multiple threads. Default: " + std::to_string(params.threads), ArgParse::Level::GENERAL, false);
//    args.add_option("output", ArgParse::Type::STRING, "Path where tracks in wig format will be written. If it does not exist, it will be created. Default: output files are stored next to the input files.", ArgParse::Level::GENERAL, false);
//
//    args.add_option("mapping", ArgParse::Type::STRING, "If the MSAs don't use common species names (like Human, Chimp, etc.) you can pass a two-column tsv file with a name mapping.", ArgParse::Level::GENERAL, false);
//    args.add_option("species", ArgParse::Type::STRING, "Comma separated list of species to reduce <model> to a subset of species to improve running time, e.g., \"Human,Chimp,Seaturtle\"", ArgParse::Level::GENERAL, false);
//    args.add_option("model-info", ArgParse::Type::STRING, "Output the organisms included in a specific model. Included models are: " + model_list + ".", ArgParse::Level::GENERAL, false);
//    // TODO: ignore sequences not occuring in model (instead of failing)
//
//    args.add_positional_argument("model", ArgParse::Type::STRING, "Path to parameter files, or one of the following predefined models: " + model_list + ".", true);
//    args.add_positional_argument("alignments", ArgParse::Type::STRING, "One or more files containing multiple sequence alignments. Formats: MAF and multi FASTA. Multiple MSAs can be stored in a single file separated by empty lines.", true, true);
//    args.parse_args(argc, argv);
//
//    // additional help strings
//    if (args.is_set("model-info"))
//    {
//        const std::string & model_name = args.get_string("model-info");
//        return print_model_info(model_name);
//    }
//
//    args.check_args(); // check here because if "--model-info" is set, we don't want to require mandatory arguments
//
//    // retrieve flags/options that were set
//    if (args.is_set("threads"))
//        params.threads = args.get_int("threads");
//
//    if (args.is_set("strategy"))
//    {
//        std::string strategy = args.get_string("strategy");
//        std::transform(strategy.begin(), strategy.end(), strategy.begin(), toupper);
//
//        if (strategy == "MLE")
//            params.strategy = MLE;
//        else if (strategy == "FIXED")
//            params.strategy = FIXED;
//        else if (strategy == "OMEGA")
//            params.strategy = OMEGA;
//        else
//        {
//            print_error_msg("Please choose a valid strategy (MLE, FIXED or OMEGA)");
//            return -1;
//        }
//    }
//
//    if (args.is_set("comp-phylo"))
//        params.comp_phylo = args.get_bool("comp-phylo");
//    if (args.is_set("comp-anc"))
//        params.comp_anc = args.get_bool("comp-anc");
//    if (args.is_set("comp-bls"))
//        params.comp_bls = args.get_bool("comp-bls");
//
//    if (params.strategy == OMEGA && params.comp_anc)
//    {
//        print_error_msg("The ancestral sequence composition cannot be computed in the Omega mode!");
//        return -1;
//    }
//
//    if (!params.comp_phylo && !params.comp_anc && !params.comp_bls)
//    {
//        print_error_msg("At least one score needs to be computed (phylo, anc or bls)!");
//        return -1;
//    }
//
//    if (args.is_set("output"))
//    {
//        params.output_path = args.get_string("output");
//        // create a directory if it doesn't exist yet
//        if (create_directory(params.output_path))
//            print_info_msg("Created the output directory.");
//    }
//
//    if (args.is_set("mapping"))
//    {
//        std::string mapping_file = args.get_string("mapping");
//        update_sequence_name_mapping(mapping_file);
//    }
//
//    if (!args.is_set_positional("model") || !args.is_set_positional("alignments"))
//    {
//        print_error_msg("No model or alignments provided.");
//        return -1;
//    }
//
//    const std::string alignment_path = args.get_positional_argument("alignments");
//
//    // load and prepare model
//    Model model;
//    load_model(model, args.get_positional_argument("model"), args.get_string("species"), false, 0, "");
//
//    // run for every alignment file
//    for (uint16_t i = 1; i < args.positional_argument_size(); ++i)
//    {
//        // printf("%d, %s\n", i, args.get_positional_argument(i).c_str());
//         run_regions(args.get_positional_argument(i), model, params);
//    }

    return 0;
}
