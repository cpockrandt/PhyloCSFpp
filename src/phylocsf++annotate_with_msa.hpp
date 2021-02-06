#include "common.hpp"
#include "arg_parse.hpp"
#include "gff_reader.hpp"

#include <omp.h>

#include <stdio.h>
#include <sys/wait.h>

struct AnnotateWithMSACLIParams
{
    std::string output_path = "";
    std::string mmseqs2_bin = "mmseqs2";
    unsigned threads = omp_get_max_threads();
    algorithm_t strategy = MLE;
    bool comp_bls = true;
};

void run_annotate_with_msa(const std::string & gff_path, const AnnotateWithMSACLIParams & params,
                           const uint32_t file_id, const uint32_t files)
{
    // read gff: extract CDS coordinates
    // read reference genome: store CDS sequences in fasta
    // do all the stuff wth mmseqs
    // score alignments
    // read, copy and annotate gff

    (void)gff_path;
    (void)params;
    (void)file_id;
    (void)files;
}

int main_annotate_with_msa(int argc, char** argv)
{
    ArgParse args("phylocsf++ annotate-with-msa",
                  "Computes PhyloCSF scores for CDS features in GFF/GTF files and outputs them in \n"
                  "the same format. Requires MMseqs2 to be installed and reference genomes to compute multiple sequence"
                  "alignments from scratch.");

    AnnotateWithMSACLIParams params;

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

    args.add_option("threads", ArgParse::Type::INT, "Parallelize scoring of multiple MSAs in a file using multiple threads. Default: " + std::to_string(params.threads), ArgParse::Level::GENERAL, false);
    args.add_option("output", ArgParse::Type::STRING, "Path where tracks in wig format will be written. If it does not exist, it will be created. Default: output files are stored next to the input files.", ArgParse::Level::GENERAL, false);

    args.add_option("strategy", ArgParse::Type::STRING, "PhyloCSF scoring algorithm: MLE, FIXED or OMEGA. Default: " + default_strategy, ArgParse::Level::GENERAL, false);
    args.add_option("comp-power", ArgParse::Type::BOOL, "Output confidence score (branch length score). Default: " + std::to_string(params.comp_bls), ArgParse::Level::GENERAL, false);

    args.add_positional_argument("genome-file", ArgParse::Type::STRING, "Two-column text file with species name and path to reference (fasta file).", true);
    args.add_positional_argument("model", ArgParse::Type::STRING, "Path to parameter files, or one of the following predefined models: " + model_list + ".", true);
    args.add_positional_argument("gff-files", ArgParse::Type::STRING, "One or more GFF/GTF files with coding exons to be scored.", true, true);
    args.parse_args(argc, argv);

    args.check_args(); // check here because if "--model-info" is set, we don't want to require mandatory arguments

    // phylocsf++ annotate --genomes-file PATH --strategy STRING --mmseqs-bin PATH <gtf file(s)>
    // both: --output, --confidence

    // check whether MMseqs2 is installed
    if (args.is_set("mmseqs-bin"))
        params.mmseqs2_bin = args.get_int("mmseqs-bin");

    if (system_with_return(params.mmseqs2_bin + " --help > /dev/null 2>&1"))
    {
        printf("ERROR: MMseqs2 seems not to be installed. `%s --help` does not seem to work.\n", params.mmseqs2_bin.c_str());
        return -1;
    }

    if (system_with_return(params.mmseqs2_bin + " mmseqs result2dnamsa --help > /dev/null 2>&1"))
    {
        printf("ERROR: A version of MMseqs2 seems to be installed that is too old. `%s result2dnamsa --help` does not seem to work.\n", params.mmseqs2_bin.c_str());
        return -1;
    }

    // retrieve flags/options that were set
    if (args.is_set("threads"))
        params.threads = args.get_int("threads");

    if (args.is_set("output"))
    {
        params.output_path = args.get_string("output");
        // create a directory if it doesn't exist yet
        if (create_directory(params.output_path))
            printf(OUT_INFO "Created the output directory.\n" OUT_RESET);
    }

    // run for every gff file
    for (uint16_t i = 1; i < args.positional_argument_size(); ++i)
    {
         run_annotate_with_msa(args.get_positional_argument(i), params, i, args.positional_argument_size() - 1);
    }

    printf(OUT_DEL "Done!\n");

    return 0;
}
