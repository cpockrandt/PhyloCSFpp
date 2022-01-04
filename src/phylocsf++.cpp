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
#include "phylocsf++annotate_with_mmseqs.hpp"
#include "phylocsf++find_cds.hpp"

int main(int argc, char **argv)
{
    ArgParse args("phylocsf++",
                  "A fast and easy-to-use tool to compute PhyloCSF scores and tracks and annotate GFF/GTF.\n"
                  "For documentation and help, check out https://github.com/cpockrandt/PhyloCSFpp\n\n"
                  "Please consider citing:\n"
                  "  Pockrandt et al., PhyloCSF++: A fast and user-friendly implementation of PhyloCSF\n"
                  "  with annotation tools, https://doi.org/10.1093/bioinformatics/btab756, Bioinformatics 2021");

    args.add_subprogram("build-tracks", "Computes PhyloCSF and Power tracks for each codon and all 6 frames from alignments from MAF files. Outputs them in wig files.");
    args.add_subprogram("score-msa", "Computes PhyloCSF scores, ancestral sequence composition sores and branch length scores for entire alignments from MAF files. Outputs them in a BED-like format.");
    args.add_subprogram("annotate-with-tracks", "Scores the CDS features in GFF/GTF files using precomputed tracks (bw files) and outputs annotated GFF/GTF files.");
    args.add_subprogram("annotate-with-mmseqs", "Scores the CDS features in GFF/GTF files by computing multiple sequence alignments from scratch (requires MMseqs2) and outputs annotated GFF/GTF files.");
    args.add_subprogram("find-cds", "Takes GFF/GTF files as input and for each transcript searches for CDS with high protein-coding likelihood using PhyloCSF tracks.");

    args.add_option("help", ArgParse::Type::FLAG, "Prints this help message. Run `phylocsf++ build-tracks --help` for help messages on the tools", ArgParse::Level::HELP, false);

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
    else if (args.get_selected_subprogram() == "annotate-with-mmseqs")
    {
        return main_annotate_with_mmseqs(argc - 1, argv + 1);
    }
    else if (args.get_selected_subprogram() == "find-cds")
    {
        return main_find_cds(argc - 1, argv + 1);
    }

    return 0;
}
