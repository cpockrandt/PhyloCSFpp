#include "common.hpp"
#include "arg_parse.hpp"

#include "phylocsf++tracks.hpp"
#include "phylocsf++region.hpp"

int main(int argc, char **argv)
{
    ArgParse args("phylocsf++",
                  "A fast and easy-to-use tool to compute PhyloCSF scores and tracks.\n"
                  "For documentation and help, check out https://github.com/cpockrandt/PhyloCSFpp\n\n"
                  "Please consider citing: Harry Frankfurt: On Bullshit. Princeton University Press,\n"
                  "                        Princeton, New Jersey 2005");

    args.set_subprograms({
            { "region", "Computes PhyloCSF scores, ancestral sequence composition sores and branch length scores for entire alignments from MAF and FASTA files. Outputs them in BED format." },
            { "tracks", "Computes PhyloCSF and Power tracks for each codon and all 6 frames from alignments from MAF and FASTA files. Outputs them in wig files." },
    });

    args.add_option("help" , /*'e',*/ ArgParse::Type::FLAG, "Prints this help message. Run `phylocsf++ tracks --help` for help messages on the tools", ArgParse::Level::GENERAL, false);

    args.parse_args(argc, argv);

    if (args.get_selected_subprogram() == "tracks")
    {
        return main_tracks(argc - 1, argv + 1);
    }
    else if (args.get_selected_subprogram() == "region")
    {
        return main_region(argc - 1, argv + 1);
    }

    return 0;
}
