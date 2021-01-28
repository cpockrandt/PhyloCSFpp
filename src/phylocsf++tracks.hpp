#include "arg_parse.hpp"
#include "models.hpp"

#include "run.hpp"

#include <omp.h>

int main_tracks(int argc, char **argv)
{
    ArgParse args("phylocsf++ tracks",
                  "Computes PhyloCSF score tracks for each codon and all 6 frames from alignments \n"
                  "in FASTA or MAF files. Output is written to wig files. Optionally a PhyloCSF \n"
                  "power track containing the branch length scores can written to a wig file \n"
                  "(confidence of the PhyloCSF scores).");

    // default values
    unsigned threads = omp_get_max_threads();
    bool phylo_smooth = false;
    bool phylo_raw = true;
    bool phylo_power = true;
    float phylo_threshold = 0.1;

    std::string model_list = get_list_of_models();

    args.add_option("output-phylo", ArgParse::Type::BOOL, "Compute all 6 smoothened PhyloCSF tracks. Requires coding exon coordinates and genome length. Default: " + std::to_string(phylo_smooth), ArgParse::Level::GENERAL, false);
    args.add_option("output-raw-phylo", ArgParse::Type::BOOL, "Compute all 6 unsmoothened PhyloCSF tracks. Default: " + std::to_string(phylo_raw), ArgParse::Level::GENERAL, false);
    args.add_option("output-power", ArgParse::Type::BOOL, "Output confidence track (branch length score). Default: " + std::to_string(phylo_power), ArgParse::Level::GENERAL, false);
    args.add_option("power-threshold", ArgParse::Type::FLOAT, "Minimum confidence score to output PhyloCSF score. Default: " + std::to_string(phylo_threshold), ArgParse::Level::GENERAL, false);
    args.add_option("genome-length", ArgParse::Type::INT, "Total genome length (needed for --output-phylo).", ArgParse::Level::GENERAL, false);
    args.add_option("coding-exons", ArgParse::Type::STRING, "TODO file with coordinates of coding exons (needed for --output-phylo).", ArgParse::Level::GENERAL, false);
    args.add_option("threads", ArgParse::Type::INT, "Parallelize scoring of multiple MSAs in a file using multiple threads. Defaut: " + std::to_string(threads), ArgParse::Level::GENERAL, false);
    args.add_option("output", ArgParse::Type::STRING, "Path where tracks in wig format will be written. If it does not exist, it will be created. Default: output files are stored next to the input files.", ArgParse::Level::GENERAL, false);
    args.add_option("mapping", ArgParse::Type::STRING, "If the MSAs don't use common species names (like Human, Chimp, etc.) you can pass a two-column tsv file with a name mapping.", ArgParse::Level::GENERAL, false);
    args.add_option("species", ArgParse::Type::STRING, "Comma separated list of species to reduce <model> to a subset of species to improve running time, e.g., \"Human,Chimp,Seaturtle\"", ArgParse::Level::GENERAL, false);
    args.add_option("model-info", ArgParse::Type::STRING, "Output the organisms included in a specific model. Included models are: " + model_list + ".", ArgParse::Level::GENERAL, false);
    // TODO: dry-run option (to check whether all species names can be mapped)

    args.add_positional_argument("model", ArgParse::Type::STRING, "Path to parameter files, or one of the following predefined models: " + model_list + ".", false);
    args.add_positional_argument("alignments", ArgParse::Type::STRING, "One or more files containing multiple sequence alignments. Formats: MAF and multi FASTA. Multiple MSAs can be stored in a single file separated by empty lines.", false);
    args.parse_args(argc, argv);

    // additional help strings
    if (args.is_set("model-info"))
    {
        const std::string & model_name = args.get_string("model-info");
        return print_model_info(model_name);
    }

    // retrieve flags/options that were set
    if (args.is_set("threads"))
        threads = args.get_int("threads");
    if (args.is_set("output-phylo"))
        phylo_smooth = args.get_bool("output-phylo");
    if (args.is_set("output-raw-phylo"))
        phylo_raw = args.get_bool("output-raw-phylo");
    if (args.is_set("power-threshold"))
        phylo_power = args.get_bool("power-threshold");

    // this has to hold true: phylo_smooth => (args.is_set("genome-length") && args.is_set("coding-exons"))
    if (phylo_smooth && (!args.is_set("genome-length") || !args.is_set("coding-exons")))
    {
        print_error_msg("For smoothened tracks (--output-phylo) you need to provide --genome-length and --coding-exons.");
        return -1;
    }

    std::string output_path = "";
    if (args.is_set("output"))
    {
        output_path = args.get_string("output");
        // create a directory if it doesn't exist yet
        if (create_directory(output_path))
            print_info_msg("Created the output directory.");
    }

    uint64_t genome_length;
    if (args.is_set("genome-length"))
        genome_length = args.get_int("genome-length");

    // TODO: open coding-exons

    // TODO: open mapping file

    // read species string (allow common and scientific names)


    // model
    // alignments

    if (!args.is_set_positional("model") || !args.is_set_positional("alignments"))
    {
        print_error_msg("No model or alignments provided.");
        return -1;
    }

    return 0;
}
