#include "common.hpp"

#include "arg_parse.hpp"
#include "models.hpp"

#include "run.hpp"

#include <omp.h>

#include <stdio.h>
#include <sys/wait.h>

#include <bigWig.h>

struct AnnotateCLIParams
{
    unsigned threads = omp_get_max_threads();
    std::string output_path = "";

    bigWigFile_t *bw_files[7] = { NULL }; // in this order: +1, +2, +3, -1, -2, -3, power (confidence)
};

void run_annotate(const std::string & gff_path, const AnnotateCLIParams & params)
{
    (void) gff_path;
    (void) params;
}

int main_annotate(int argc, char** argv)
{
    ArgParse args("phylocsf++ annotate",
                  "Computes PhyloCSF scores for CDS features in GFF/GTF files and outputs them in \n"
                  "the same format. You can either pass track files (in wig format) or have the \n"
                  "multiple sequence alignment computed for the CDS (this requires mmseqs2 to be \n"
                  "installed).");

    // phylocsf++ annotate --tracks PATH <gtf file(s)>
    // phylocsf++ annotate --genomes-file PATH --model 100vertebrates --strategy STRING --mmseqs-bin PATH <gtf file(s)>
    // both: --output, --confidence

    // check whether mmseqs is installed
//    if (system_with_return("mmseqs --help > /dev/null 2>&1"))
//    {
//        print_error_msg("ERROR: MMseqs2 seems not to be installed.\n");
//        return -1;
//    }
//
//    if (system_with_return("mmseqs result2dnamsa --help > /dev/null 2>&1"))
//    {
//        print_error_msg("ERROR: A version of MMseqs2 seems to be installed that is too old. `mmseqs result2dnamsa --help` does not seem to work.\n");
//        return -1;
//    }

    AnnotateCLIParams params;

    args.add_option("tracks", ArgParse::Type::STRING, "Path to PhyloCSF bigWig track for frame +1 (expects the other 5 frames to be in the same directory, optionally the power track)", ArgParse::Level::GENERAL, false);

    args.add_option("threads", ArgParse::Type::INT, "Parallelize scoring of multiple MSAs in a file using multiple threads. Default: " + std::to_string(params.threads), ArgParse::Level::GENERAL, false);
    args.add_option("output", ArgParse::Type::STRING, "Path where tracks in wig format will be written. If it does not exist, it will be created. Default: output files are stored next to the input files.", ArgParse::Level::GENERAL, false);

    args.add_positional_argument("gff-files", ArgParse::Type::STRING, "One or more GFF/GTF files with coding exons to be scored.", true, true);
    args.parse_args(argc, argv);

    args.check_args(); // check here because if "--model-info" is set, we don't want to require mandatory arguments

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

    // load tracks
    if (bwInit(1 << 17) != 0)
    {
        printf(OUT_ERROR "Something failed in bwInit().\n" OUT_RESET);
        return -1;
    }

    std::string bw_path = args.get_string("tracks");
    const size_t bw_path_suffix_pos = bw_path.find("+1");
    if (bw_path_suffix_pos == std::string::npos)
    {
        printf(OUT_ERROR "Could not find '+1' in tracks file name. Expecting a name like 'PhyloCSF+1.bw'.\n" OUT_RESET);
        return -1;
    }
    for (uint16_t i = 0; i < 7; ++i)
    {
        std::string suffix;
        if (i == 6)
            suffix = "power";
        else
            suffix = ((i < 3) ? "+" : "-") + std::to_string((i % 3) + 1);
        bw_path.replace(bw_path_suffix_pos, 2 /* length of "+1" */, suffix);
//        printf("%s\n", bw_path.c_str());
        // TODO: check whether file exists (and print error unless it's the power track)
        params.bw_files[i] = bwOpen(const_cast<char * >(bw_path.c_str()), NULL, "r");
        if (!params.bw_files[i])
        {
            printf(OUT_ERROR "An error occurred while opening the PhyloCSF file '%s'.\n" OUT_RESET, bw_path.c_str());
            return -1;
        }
    }

    // run for every gff file
    const std::string alignment_path = args.get_positional_argument("gff-files");
    for (uint16_t i = 1; i < args.positional_argument_size(); ++i)
    {
         run_annotate(args.get_positional_argument(i), params);
    }

    printf(OUT_DEL "Done!\n");

    return 0;
}
