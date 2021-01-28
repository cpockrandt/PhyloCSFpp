#include "arg_parse.hpp"
#include "models.hpp"

#include "run.hpp"

#include <omp.h>

int main_region(int argc, char **argv)
{
    ArgParse args("phylocsf++ region",
                  "Computes PhyloCSF scores for whole alignments in FASTA or MAF files. Output is \n"
                  "written to bed file(s). Other scores such as the ancestral sequence composition \n"
                  "sores and branch length scores can be computed as well. Only one forward frame \n"
                  "is computed, i.e., for other frames reverse the alignments and/or remove the \n"
                  "first one or two bases.");//For ORF finding and scoring of transcripts in GTF/GFF "
                  // "files, check out PhyloCSF++Anno: https://github.com/cpockrandt/phylocsfppanno");

    // default values
    unsigned threads = omp_get_max_threads();
    bool comp_phylo = true;
    bool comp_anc = false;
    bool comp_bls = false;
    std::string strategy = "MLE";

    std::string model_list = "";
    for (const auto & m : models)
        model_list += m.first + ", ";
    assert(model_list.size() >= 2);
    model_list.erase(model_list.size() - 2);

    args.add_option("strategy", ArgParse::Type::STRING, "PhyloCSF scoring algorithm: MLE, FIXED or OMEGA. Default: " + strategy, ArgParse::Level::GENERAL, false);
    args.add_option("comp-phylo", ArgParse::Type::BOOL, "Compute the PhyloCSF score for each alignment. Default: " + std::to_string(comp_phylo), ArgParse::Level::GENERAL, false);
    args.add_option("comp-anc", ArgParse::Type::BOOL, "Compute the ancestral sequence composition score (only in MLE and FIXED mode). Default: " + std::to_string(comp_anc), ArgParse::Level::GENERAL, false);
    args.add_option("comp-bls", ArgParse::Type::BOOL, "Compute the branch length score (confidence). Default: " + std::to_string(comp_bls), ArgParse::Level::GENERAL, false);

    args.add_option("threads", ArgParse::Type::INT, "Parallelize scoring of multiple MSAs in a file using multiple threads. Defaut: " + std::to_string(threads), ArgParse::Level::GENERAL, false);
    args.add_option("output", ArgParse::Type::STRING, "Path where tracks in wig format will be written. If it does not exist, it will be created. Default: output files are stored next to the input files.", ArgParse::Level::GENERAL, false);

    args.add_option("mapping", ArgParse::Type::STRING, "If the MSAs don't use common species names (like Human, Chimp, etc.) you can pass a two-column tsv file with a name mapping.", ArgParse::Level::GENERAL, false);
    args.add_option("species", ArgParse::Type::STRING, "Comma separated list of species to reduce <model> to a subset of species to improve running time, e.g., \"Human,Chimp,Seaturtle\"", ArgParse::Level::GENERAL, false);
    args.add_option("model-info", ArgParse::Type::STRING, "Output the organisms included in a specific model. Included models are: " + model_list + ".", ArgParse::Level::GENERAL, false);
    // TODO: dry-run option (to check whether all species names can be mapped)

    args.add_positional_argument("model", ArgParse::Type::STRING, "Path to parameter files, or one of the following predefined models: " + model_list + ".", false);
    args.add_positional_argument("alignments", ArgParse::Type::STRING, "One or more files containing multiple sequence alignments. Formats: MAF and multi FASTA. Multiple MSAs can be stored in a single file separated by empty lines.", false);
    args.parse_args(argc, argv);

    // retrieve flags/options that were set
    if (args.is_set("threads"))
        threads = args.get_int("threads");

    if (args.is_set("strategy"))
    {
        strategy = args.get_bool("strategy");
        std::transform(strategy.begin(), strategy.end(), strategy.begin(), toupper);
        if (strategy != "MLE" && strategy != "FIXED" && strategy != "OMEGA")
        {
            print_error_msg("Please choose a valid strategy (MLE, FIXED or OMEGA)");
            return -1;
        }
    }

    if (args.is_set("comp-phylo"))
        comp_phylo = args.get_bool("comp-phylo");
    if (args.is_set("comp-anc"))
        comp_anc = args.get_bool("comp-anc");
    if (args.is_set("comp-bls"))
        comp_bls = args.get_bool("comp-bls");

    if (strategy == "OMEGA" && comp_anc)
    {
        print_error_msg("The ancestral sequence composition cannot be computed in the Omega mode!");
        return -1;
    }

    if (!comp_phylo && !comp_anc && !comp_bls)
    {
        print_error_msg("At least one score needs to be computed (phylo, anc or bls)!");
        return -1;
    }

    return 0;
}
