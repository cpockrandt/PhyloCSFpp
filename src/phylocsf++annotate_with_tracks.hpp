#include "common.hpp"
#include "arg_parse.hpp"
#include "gff_reader.hpp"

#include <unordered_map>

#include <stdio.h>
#include <sys/wait.h>

#include <bigWig.h>
#include <bwValues.c>
#include <io.c>
#include <bwWrite.c>
#include <bwRead.c>

struct AnnotateWithTracksCLIParams
{
    std::string output_path = "";
    std::string bw_path = ""; // needed for header output

    bigWigFile_t *bw_files[7] = { NULL }; // in this order: +1, +2, +3, -1, -2, -3, power/confidence
    std::unordered_map<std::string, uint64_t> chrom_sizes;

    bool comp_bls = true;
};

void count_scores(float & sum, uint64_t & count, bigWigFile_t *file, const char *chrom,
                  const uint64_t begin, const uint64_t end)
{
    bwOverlappingIntervals_t *intervals = bwGetValues(file, (char*) chrom, begin, end, 0);
    if (intervals != NULL)
    {
        for (uint32_t i = 0; i < (*intervals).l; ++i)
            sum += (*intervals).value[i];
        count = (*intervals).l;
        bwDestroyOverlappingIntervals(intervals);
    }
}

void run_annotate_with_tracks(const std::string & gff_path, const AnnotateWithTracksCLIParams & params,
                              std::unordered_set<std::string> & missing_sequences,
                              const uint32_t file_id, const uint32_t files)
{
    gff_reader reader(gff_path.c_str());
    reader.setup_progressbar(file_id, files);

    // store output next to GFF file if no output_path was specified
    std::string output_file_path = params.output_path;
    if (output_file_path == "")
    {
        output_file_path = gff_path;
    }
    else
    {
        const size_t pos = gff_path.find_last_of('/');
        if (pos != std::string::npos) // gff_path == ".../file.gff"
            output_file_path += "/" + gff_path.substr(pos + 1);
        else                          // gff_path == "file.gff"
            output_file_path += "/" + gff_path;
    }

    // insert/append ".PhyloCSF++" (in)to filename
    const size_t output_file_path_last_dot = output_file_path.find_last_of('.');
    if (output_file_path_last_dot == std::string::npos)
        output_file_path += ".PhyloCSF++"; // NoFileEnding -> NoFileEnding.PhyloCSF++
    else
        output_file_path.insert(output_file_path_last_dot, ".PhyloCSF++"); // Human.Genes.gtf -> Human.Genes.PhyloCSF++.gtf

    FILE *gff_out = fopen(output_file_path.c_str(), "w");
    if (gff_out == NULL)
    {
        printf(OUT_ERROR "Error creating file %s!\n" OUT_RESET, output_file_path.c_str());
        exit(1);
    }

    // TODO: move header line to the end of the already existing header block
    fprintf(gff_out, "# PhyloCSF scores computed with PhyloCSF++ %s (%s, %s) and precomputed tracks %s\n", ArgParse::version.c_str(), ArgParse::git_hash.c_str(), ArgParse::git_date.c_str(), params.bw_path.c_str());

    gff_transcript t;
    while (reader.get_next_transcript<true>(t))
    {
        float transcript_sum = 0.0;
        float transcript_power_sum = 0.0;
        uint64_t transcript_count = 0;
        uint64_t transcript_power_count = 0;

        if (t.CDS.size() > 0)
        {
            // check whether chromosome exists in tracks
            const auto chrom_sizes_it = params.chrom_sizes.find(t.chr);
            if (chrom_sizes_it == params.chrom_sizes.end())
            {
                t.phylo_score = NAN;
                t.phylo_power = NAN;
                // only report a missing sequence once
                if (missing_sequences.find(t.chr) == missing_sequences.end())
                {
                    missing_sequences.insert(t.chr);
                    printf(OUT_DEL OUT_INFO "Sequence %s from the GFF file does not occur in the tracks. Skipping ...\n" OUT_RESET, t.chr.c_str());
                }
            }
            else
            {
                const uint64_t chr_len = chrom_sizes_it->second;

                for (auto & c : t.CDS)
                {
                    uint8_t wig_phase;
                    if (t.strand == '+')
                        wig_phase = 0 + (c.phase + c.begin - 1) % 3;
                    else
                        wig_phase = 3 + (chr_len - c.end - 1 + c.phase + 1) % 3;

                    float cds_sum = 0.0;
                    uint64_t cds_count = 0;

                    float cds_power_sum = 0.0;
                    uint64_t cds_power_count = 0;

                    count_scores(cds_sum, cds_count, params.bw_files[wig_phase], t.chr.c_str(), c.begin - 1, c.end);
                    c.phylo_score = cds_sum / cds_count;

                    if (params.comp_bls)
                    {
                        count_scores(cds_power_sum, cds_power_count, params.bw_files[6], t.chr.c_str(), c.begin - 1, c.end);
                        c.phylo_power = cds_power_sum / cds_power_count;
                    }

                    transcript_sum += cds_sum;
                    transcript_count += cds_count;

                    transcript_power_sum += cds_power_sum;
                    transcript_power_count += cds_power_count;
                }

                t.phylo_score = transcript_sum / transcript_count;
                t.phylo_power = transcript_power_sum / transcript_power_count;
            }
        }

        // output gtf file (with annotation)
        uint16_t CDS_id = 0;
        for (auto line_tuple : t.lines)
        {
            const feature_t f = std::get<0>(line_tuple);
            const std::string & line = std::get<1>(line_tuple);
            if (f == OTHER || t.CDS.size() == 0) // no annotation necessary
            {
                fprintf(gff_out, "%s\n", line.c_str());
            }
            else
            {
                float score;
                float power;
                if (f == TRANSCRIPT)
                {
                    score = t.phylo_score;
                    power = t.phylo_power;
                }
                else
                {
                    score = t.CDS[CDS_id].phylo_score;
                    power = t.CDS[CDS_id].phylo_power;
                    ++CDS_id;
                }

                // TODO: use correct format (gff/gtf)
                if (params.comp_bls)
                    fprintf(gff_out, "%s phylocsf_mean \"%.3f\"; phylocsf_power_mean \"%.3f\";\n", line.c_str(), score, power);
                else
                    fprintf(gff_out, "%s phylocsf_mean \"%.3f\";\n", line.c_str(), score);
            }
        }

        reader.print_progress();
    }

    fclose(gff_out);
}

int main_annotate_with_tracks(int argc, char** argv)
{
    ArgParse args("phylocsf++ annotate-with-tracks",
                  "Computes PhyloCSF scores for CDS features in GFF/GTF files and outputs them in \n"
                  "the same format. Requires PhyloCSF tracks in BigWig (*.bw) format.");

    AnnotateWithTracksCLIParams params;

    args.add_option("output", ArgParse::Type::STRING, "Path where tracks in wig format will be written. If it does not exist, it will be created. Default: output files are stored next to the input files.", ArgParse::Level::GENERAL, false);
    args.add_option("comp-power", ArgParse::Type::BOOL, "Output confidence score (branch length score). Requires the file PhyloCSFpower.bw in the same directory as the tracks. Default: " + std::to_string(params.comp_bls), ArgParse::Level::GENERAL, false);

    args.add_positional_argument("tracks", ArgParse::Type::STRING, "Path to the bigWig file PhyloCSF+1.bw (expects the other 5 frames to be in the same directory, optionally the power track).", true);
    args.add_positional_argument("gff-files", ArgParse::Type::STRING, "One or more GFF/GTF files with coding exons to be scored.", true, true);
    args.parse_args(argc, argv);

    args.check_args();

    // retrieve flags/options that were set
    if (args.is_set("comp-power"))
        params.comp_bls = args.get_bool("comp-power");

    if (args.is_set("output"))
    {
        params.output_path = args.get_string("output");
        // create a directory if it doesn't exist yet
        if (create_directory(params.output_path))
            printf(OUT_INFO "Created the output directory.\n" OUT_RESET);
    }

    // load tracks
    // if (bwInit(1 << 17) != 0) // only needed for remote files (CURL)
    // {
    //     printf(OUT_ERROR "Something failed in bwInit().\n" OUT_RESET);
    //     return -1;
    //}

    params.bw_path = args.get_positional_argument("tracks");
    const size_t bw_path_suffix_pos = params.bw_path.find("+1");
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
        params.bw_path.replace(bw_path_suffix_pos, 2, suffix); // NOTE: length of "+1" is 2

        if (i < 6 || params.comp_bls)
        {
            params.bw_files[i] = bwOpen(const_cast<char * >(params.bw_path.c_str()), NULL, "r");
            if (!params.bw_files[i])
            {
                printf(OUT_ERROR "An error occurred while opening the PhyloCSF file '%s'.\n" OUT_RESET, params.bw_path.c_str());

                // check whether the user has used an unindexed wig file, then print a useful hint
                if (access(params.bw_path.c_str(), F_OK) == 0 &&
                    params.bw_path.size() >= 4 && params.bw_path.compare(params.bw_path.size() - 4, 4, ".wig") == 0
                )
                    printf(OUT_INFO "It seems you provided a *.wig file. You need to simply index them first with wigToBigWig and then use the *.bw files.\n" OUT_RESET);
                else if (i == 6)
                    printf(OUT_INFO "If you set `--comp-power 0` you can compute it without PhyloCSFpower.bw.\n" OUT_RESET);

                return -1;
            }
        }
    }

    const chromList_t *chr_list = bwReadChromList(params.bw_files[0]);
    for (int64_t i = 0; i < chr_list->nKeys; ++i)
    {
        params.chrom_sizes.emplace(chr_list->chrom[i], chr_list->len[i]);
    }

    // run for every gff file
    std::unordered_set<std::string> missing_sequences;
    for (uint16_t i = 1; i < args.positional_argument_size(); ++i)
    {
         run_annotate_with_tracks(args.get_positional_argument(i), params, missing_sequences,
                                  i, args.positional_argument_size() - 1);
    }

    printf(OUT_DEL "Done!\n");

    return 0;
}
