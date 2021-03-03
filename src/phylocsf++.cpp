#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include "arg_parse.hpp"

#ifdef ENABLE_OPENMP
    #include <omp.h>
#else
    int omp_get_max_threads() { return 1; }
    int omp_get_thread_num() { return 0; }
#endif

#include "phylocsf++build_tracks.hpp"
#include "phylocsf++score_msa.hpp"
#include "phylocsf++annotate_with_tracks.hpp"
#include "phylocsf++annotate_with_msa.hpp"

int main(int argc, char **argv)
{
    ArgParse args("phylocsf++",
                  "A fast and easy-to-use tool to compute PhyloCSF scores and tracks.\n"
                  "For documentation and help, check out https://github.com/cpockrandt/PhyloCSFpp\n\n"
                  "Please consider citing: Harry Frankfurt: On Bullshit. Princeton University Press,\n"
                  "                        Princeton, New Jersey 2005");

    args.add_subprogram("build-tracks", "Computes PhyloCSF and Power tracks for each codon and all 6 frames from alignments from MAF files. Outputs them in wig files.");
    args.add_subprogram("score-msa", "Computes PhyloCSF scores, ancestral sequence composition sores and branch length scores for entire alignments from MAF files. Outputs them in a BED-like format.");
    args.add_subprogram("annotate-with-tracks", "Scores the CDS features in a GFF/GTF file using precomputed tracks (bw files) and outputs an annotated GFF/GTF file.");
    args.add_subprogram("annotate-with-msa", "Scores the CDS features in a GFF/GTF file by computing multiple sequence alignments from scratch (requires MMseqs2) and outputs an annotated GFF/GTF file.");

    args.add_option("help" , /*'e',*/ ArgParse::Type::FLAG, "Prints this help message. Run `phylocsf++ build-tracks --help` for help messages on the tools", ArgParse::Level::GENERAL, false);

    args.parse_args(argc, argv);

    // TODO: warn if files are going to be overwritten?

    if (args.get_selected_subprogram() == "build-tracks")
    {
        return main_build_tracks(argc - 1, argv + 1);
    }
    else if (args.get_selected_subprogram() == "score-msa")
    {
        return main_score_msa(argc - 1, argv + 1);
    }
    else if (args.get_selected_subprogram() == "annotate-with-tracks")
    {
        return main_annotate_with_tracks(argc - 1, argv + 1);
    }
    else if (args.get_selected_subprogram() == "annotate-with-msa")
    {
        return main_annotate_with_msa(argc - 1, argv + 1);
    }

    return 0;
}
