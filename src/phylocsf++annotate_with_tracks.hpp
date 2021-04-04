#include "common.hpp"
#include "arg_parse.hpp"
#include "gff_reader.hpp"

#include <unordered_map>

#include <cmath>
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
};

void count_weighted_scores(float & weighted_score_sum, float & weighted_power_sum, float & power_sum, uint64_t & power_count,
                           bigWigFile_t *file, bigWigFile_t *power_file, const char *chrom, const uint64_t begin, const uint64_t end)
{
    // last parameter has to be 1 - because we need to retrieve missing values (otherwise, we cannot match them from both files)
    // the power file might have values, where the score files doesn't (because the confidence was too low)
    bwOverlappingIntervals_t *intervals = bwGetValues(file, (char*) chrom, begin, end, 1);
    bwOverlappingIntervals_t *intervals_power = bwGetValues(power_file, (char*) chrom, begin, end, 1);
    if (intervals != NULL && intervals_power != NULL)
    {
        for (uint32_t i = 0; i < (*intervals).l; ++i)
        {
            if (!std::isnan((*intervals).value[i]) && !std::isnan((*intervals_power).value[i]))
            {
                weighted_score_sum += (*intervals).value[i] * (*intervals_power).value[i];
                weighted_power_sum += (*intervals_power).value[i];
            }

            // weighted_power_sum only considers confidence values that are above the threshold
            // (if it's below the threshold, the phylocsf score is missing/nan).
            if (!std::isnan((*intervals_power).value[i]))
            {
                power_sum += (*intervals_power).value[i];
            }

            // if there is no confidence available (i.e., no alignment), we still count the codon
            // to reduce the mean confidence
            ++power_count;
        }
    }

    if (intervals)
        bwDestroyOverlappingIntervals(intervals);
    if (intervals_power)
        bwDestroyOverlappingIntervals(intervals_power);
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
        float transcript_all_power_sum = 0.0;
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

                    float cds_power_sum = 0.0;
                    float cds_all_power_sum = 0.0;
                    uint64_t cds_power_count = 0;

                    // cds_sum and cds_power_sum are for the weighted computation (and only consider positions where a score is available)
                    count_weighted_scores(cds_sum, cds_power_sum, cds_all_power_sum, cds_power_count, params.bw_files[wig_phase], params.bw_files[6], t.chr.c_str(), c.begin - 1, c.end);
                    c.phylo_score = cds_sum / cds_power_sum;
                    // here we consider all positions, even if the confidence is below the threshold (when the tracks were computed)
                    c.phylo_power = (cds_power_count == 0) ? 0 : cds_all_power_sum / cds_power_count;

                    transcript_sum += cds_sum;
                    transcript_power_sum += cds_power_sum;
                    transcript_all_power_sum += cds_all_power_sum;
                    transcript_power_count += cds_power_count;
                }

                t.phylo_score = transcript_sum / transcript_power_sum;
                t.phylo_power = (transcript_power_count == 0) ? 0 : transcript_all_power_sum / transcript_power_count;
            }
        }

        // output gtf file (with annotation)
        bool first_processed_line = true;
        bool is_gff = true;
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
                // detect format
                if (first_processed_line)
                {
                    first_processed_line = false;
                    is_gff = is_gff_format(line);
                }

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

                // workaround for "-nan" output (which makes our tests fail)
                if (std::isnan(score))
                    score = NAN;

                if (is_gff)
                {
                    fprintf(gff_out, "%s;phylocsf_score_weighted_mean=%.3f;phylocsf_power_mean=%.3f\n", line.c_str(), score, power);
                }
                else
                {
                    fprintf(gff_out, "%s phylocsf_score_weighted_mean \"%.3f\"; phylocsf_power_mean \"%.3f\";\n", line.c_str(), score, power);
                }
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

    args.add_option("output", ArgParse::Type::STRING, "Path where output GFF/GTF will be written to. If it does not exist, it will be created. Default: output files are stored next to the input files.", ArgParse::Level::GENERAL, false);

    args.add_positional_argument("tracks", ArgParse::Type::STRING, "Path to the bigWig file PhyloCSF+1.bw (expects the other 5 frames to be in the same directory, optionally the power track).", true);
    args.add_positional_argument("gff-files", ArgParse::Type::STRING, "One or more GFF/GTF files with coding exons to be scored.", true, true);
    args.parse_args(argc, argv);

    args.check_args();

    // retrieve flags/options that were set
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

    std::string bw_path = args.get_positional_argument("tracks");
    params.bw_path = bw_path; // wo don't want to change params.bw_path below
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
        bw_path.replace(bw_path_suffix_pos, 2, suffix); // NOTE: length of "+1" is 2

        params.bw_files[i] = bwOpen(const_cast<char * >(bw_path.c_str()), NULL, "r");
        if (!params.bw_files[i])
        {
            // check whether the user has used an unindexed wig file, then print a useful hint
            if (access(bw_path.c_str(), F_OK) == 0 &&
                bw_path.size() >= 4 && bw_path.compare(bw_path.size() - 4, 4, ".wig") == 0
                    )
            {
                printf(OUT_ERROR "An error occurred while opening the PhyloCSF file '%s'.\n" OUT_RESET, bw_path.c_str());
                printf(OUT_INFO "It seems you provided a *.wig file. You need to simply index them first with wigToBigWig and then use the *.bw files.\n" OUT_RESET);
                return -1;
            }
            else
            {
                printf(OUT_ERROR "Could not find PhyloCSF track file '%s'.\n" OUT_RESET, bw_path.c_str());
                return -1;
            }
        }
    }

    chromList_t *chr_list = bwReadChromList(params.bw_files[0]);
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

    destroyChromList(chr_list);
    for (uint16_t i = 0; i < 7; ++i)
        bwClose(params.bw_files[i]);

    printf(OUT_DEL "Done!\n");

    return 0;
}