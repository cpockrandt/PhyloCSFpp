#include "common.hpp"

#include "arg_parse.hpp"

#include <sstream>

#define BEGIN(exon) std::get<0>(exon)
#define END(exon) std::get<1>(exon)

struct FindCDSCLIParams
{
    std::string output_path = "";
    std::string bw_path = ""; // needed for header output

    bigWigFile_t *bw_files[7] = { NULL }; // in this order: +1, +2, +3, -1, -2, -3, power/confidence
    std::unordered_map<std::string, uint64_t> chrom_sizes;

    bool comp_bls = true;
};

void find_all_codons(const std::string & dna, const std::string & codon, std::array<std::vector<uint32_t>, 3> & hits)
{
    size_t pos = dna.find(codon, 0);
    while (pos != std::string::npos)
    {
        hits[pos % 3].push_back(pos);
        pos = dna.find(codon, pos + 1);
    }
}

std::vector<std::tuple<uint32_t, uint32_t> > get_all_ORFs(const std::string & spliced_transcript_seq, const char strand)
{
    std::vector<std::tuple<uint32_t, uint32_t> > closed_reading_frames;

    std::array<std::vector<uint32_t>, 3> startpos;
    std::array<std::vector<uint32_t>, 3> stoppos;

    const uint32_t min_length = 25 * 3;
    const uint32_t max_length = 9999999;

    find_all_codons(spliced_transcript_seq, "ATG", startpos);
    find_all_codons(spliced_transcript_seq, "TAA", stoppos);
    find_all_codons(spliced_transcript_seq, "TAG", stoppos);
    find_all_codons(spliced_transcript_seq, "TGA", stoppos);

    for (uint8_t i = 0; i < 3; ++i)
    {
        sort(stoppos[i].begin(), stoppos[i].end());

        for (auto start : startpos[i])
        {
            // translation cannot continue after the first stop-codon
            auto stop_it = std::find_if(stoppos[i].begin(), stoppos[i].end(), [start] (uint32_t element) { return (start < element); });

            if (stop_it != stoppos[i].end())
            {
                uint32_t stop = *stop_it;

                if (strand == '+')
                {
                    stop += 2; // include stop codon in sequence
                    if (stop - start + 1 >= min_length && stop - start + 1 <= max_length)
                    {
                        closed_reading_frames.emplace_back(start, stop);
                        auto crf = spliced_transcript_seq.substr(start, stop - start + 1);
                        // printf("%d\t%d\t%s\n", start, stop, crf.c_str());
                    }
                }
                else if (strand == '-')
                {
                    uint32_t stop_new = spliced_transcript_seq.size() - start - 3;
                    stop_new += 2; // include stop codon in sequence

                    start = spliced_transcript_seq.size() - stop - 3;

                    if (stop_new - start + 1 >= min_length && stop_new - start + 1 <= max_length)
                    {
                        closed_reading_frames.emplace_back(start, stop_new);
                        auto crf = spliced_transcript_seq.substr(start, stop_new - start + 1);
                        // printf("%d\t%d\t%s\n", start, stop_new, crf.c_str());
                    }
                }
            }
        }
    }

    return closed_reading_frames;
}

template <typename iter_t>
void annotateCDSPhases(iter_t it, iter_t end)
{
    // phase = last CDS has `phase` bases extra
    // c.phase = new CDS has `phase2` bases that we skip
    uint8_t phase = 0;
    for (; it != end; ++it)
    {
        cds_entry & c = *it;
        c.phase = (3 - phase) % 3;
        phase = (phase + c.end - c.begin) % 3;
    }
}

// phylo_scores[exon_id][coord_in_exon]
void compute_PhyloCSF_for_transcript(const gff_transcript & t, bigWigFile_t * const * wigFiles,
                                     std::array<std::vector<std::vector<float> >, 4> & extracted_scores)
{
    for (uint8_t i = 0; i < 4; ++i)
        extracted_scores[i].resize(t.exons.size());

    for (uint16_t exon_id = 0; exon_id < t.exons.size(); ++exon_id)
    {
        const auto & exon = t.exons[exon_id];

        for (uint8_t i = 0; i < 4; ++i)
            extracted_scores[i][exon_id].resize(END(exon) - BEGIN(exon), -999.0f);

        bwOverlappingIntervals_t *intervals = NULL;

        // power
        intervals = bwGetValues(wigFiles[6], (char*)t.chr.c_str(), BEGIN(exon), END(exon), 0);
        if (intervals != NULL)
        {
            for (uint32_t i = 0; i < (*intervals).l; ++i)
                extracted_scores[3][exon_id][(*intervals).start[i] - BEGIN(exon)] = (*intervals).value[i];
            bwDestroyOverlappingIntervals(intervals);
            intervals = NULL;
        }

        if (t.strand == '+')
        {
            for (uint8_t phase = 0; phase < 3; ++phase)
            {
                // (cds_phase + exon.beginPos) % 3 == phase
                intervals = bwGetValues(wigFiles[0 + phase], (char*)t.chr.c_str(), BEGIN(exon), END(exon), 0);
                if (intervals != NULL)
                {
                    for (uint32_t i = 0; i < (*intervals).l; ++i)
                        extracted_scores[phase][exon_id][(*intervals).start[i] - BEGIN(exon)] = (*intervals).value[i];
                    bwDestroyOverlappingIntervals(intervals);
                    intervals = NULL;
                }
            }
        }
        else // negative strand
        {
            for (uint8_t phase = 0; phase < 3; ++phase)
            {
                // (chrLen - c.endPos - 1 + cds_phase + 1) % 3 == phase
                intervals = bwGetValues(wigFiles[3 + phase], (char*)t.chr.c_str(), BEGIN(exon), END(exon), 0);
                if (intervals != NULL)
                {
                    for (int32_t i = (*intervals).l - 1; i >= 0; --i)
                        extracted_scores[phase][exon_id][(END(exon) - 1) - (*intervals).start[i]] = (*intervals).value[i];
                    bwDestroyOverlappingIntervals(intervals);
                    intervals = NULL;
                }
            }
        }
    }
}

template <typename TExons, typename TIter>
std::tuple<float, float>
compute_PhyloCFS(TExons & exons, TIter first, TIter last, char const strand,
                 const std::array<std::vector<std::vector<float> >, 4> & extracted_scores,
                 const uint32_t first_exon_id_in_CDS, const uint32_t last_exon_id_in_CDS, const uint32_t chrLen = 0)
{
    float total_phylo_sum2 = 0;
    uint32_t total_phylo_count2 = 0;
    uint32_t cds_id = 0;

    for (auto it = first; it != last; ++it, ++cds_id)
    {
        cds_entry & c = *it;

        // new_phase = last CDS has `new_phase` bases extra
        // new_phase2 = new CDS has `new_phase2` bases that we skip (because we exclude them since we already retrieved score from CDS entry before)
        // uint32_t new_phase2 = (3 - new_phase) % 3;

        uint32_t cds_phase;
        switch (c.phase)
        {
            case '0':
                cds_phase = 0;
                break;
            case '1':
                cds_phase = 1;
                break;
            case '2':
                cds_phase = 2;
                break;
            default:
                cds_phase = 99;
        }

        uint32_t exon_id;
        uint32_t phylo_start, phylo_end;
        const std::vector<float> * phased_score = NULL;

        // compute phylo score
        if (strand == '+')
        {
            exon_id = first_exon_id_in_CDS + cds_id;
            phased_score = &extracted_scores[(cds_phase + c.begin) % 3][exon_id];
            phylo_start = c.begin - BEGIN(exons[exon_id]);
            phylo_end = END(exons[exon_id]) - c.end;
        }
        else if (strand == '-')
        {
            // formula derived experimentally from wig files. `start` and `end` are 0-based
            // (chrLen - start) % 3
            // use end-coordinates instead because we start with the last CDS (and we use new_phase2 to deal with codons spanning two CDS entries)
            // (chrLen - (end - 3 + 1)) % 3
            // end-coordinates are stored as half-open interval, hence substract one more
            // (chrLen - (end - 3 + 2)) % 3
            // skip `new_phase2` bases since we considered new_phase = 3 - new_phase2 bases of the same codon in the previous CDS entry already
            // (chrLen - (end - 3 + 2 - new_phase2)) % 3
            exon_id = last_exon_id_in_CDS - cds_id;
            phased_score = &extracted_scores[(chrLen - c.end - 1 + cds_phase + 1) % 3][exon_id];
            phylo_start = END(exons[exon_id]) - c.end;
            phylo_end = c.begin - BEGIN(exons[exon_id]);
        }

        float phylo_sum = 0;
        uint32_t phylo_count = (*phased_score).size() - phylo_end - phylo_start;
        total_phylo_count2 += phylo_count;

        uint32_t i = 0;
        for (auto phylo_it = (*phased_score).begin() + phylo_start; phylo_it != (*phased_score).end() - phylo_end; ++phylo_it)
        {
            if (*phylo_it != -999.0f)
            {
                phylo_sum += *phylo_it;
            }
            ++i;
        }

        c.phylo_score = (phylo_count > 0) ? phylo_sum / phylo_count : NAN;
        c.phylo_power = NAN; // TODO

        total_phylo_sum2 += phylo_sum;
    }

    return std::make_tuple(
        (total_phylo_count2 > 0) ? total_phylo_sum2 / total_phylo_count2 : NAN, // total phylo_score
        NAN // total phylo_power
    );
}

void run_find_cds(const std::string & gff_path, const FindCDSCLIParams & params,
                  const std::unordered_map<std::string, std::string> & genome,
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
    unsigned total = 0;
//    unsigned correct = 0;
    while (reader.get_next_transcript<true>(t, true))
    {
//        if (t.CDS.size() == 0)
//            continue;

        const auto chrom_sizes_it = params.chrom_sizes.find(t.chr);
        if (chrom_sizes_it == params.chrom_sizes.end())
        {
            // only report a missing sequence once
            if (missing_sequences.find(t.chr) == missing_sequences.end())
            {
                missing_sequences.insert(t.chr);
                printf(OUT_DEL OUT_INFO "Sequence %s from the GFF file does not occur in the tracks. Skipping ...\n" OUT_RESET, t.chr.c_str());
            }
        }

        const auto genome_it = genome.find(t.chr);
        if (genome_it == genome.end())
        {
            // only report a missing sequence once
            if (missing_sequences.find(t.chr) == missing_sequences.end())
            {
                missing_sequences.insert(t.chr);
                printf(OUT_DEL OUT_INFO "Sequence %s from the GFF file does not occur in the genome. Skipping ...\n" OUT_RESET, t.chr.c_str());
            }
        }

        const std::string & chr = genome_it->second;

//        std::string seq = "";
//        for (const auto & c : t.CDS)
//            seq += chr.substr(c.begin - 1, c.end - c.begin + 1);
//
//        if (t.strand == '+' && t.CDS.back().end + 3 - 1 < chr.size())
//            seq += chr.substr(t.CDS.back().end, 3);
//        else if (t.strand == '-' && t.CDS.front().begin >= 4)
//            seq = chr.substr(t.CDS.front().begin - 4, 3) + seq;
//
//        str_to_upper(seq);
//        if (t.strand == '-')
//        {
//            std::reverse(seq.begin(), seq.end());
//            for (uint64_t j = 0; j < seq.size(); ++j)
//                seq[j] = complement(seq[j]);
//        }

        ++total;
//        if (seq.size() % 3 == 0 &&
//            seq.substr(0, 3) == "ATG" && (
//                seq.substr(seq.size() - 3) == "TAA" ||
//                seq.substr(seq.size() - 3) == "TGA" ||
//                seq.substr(seq.size() - 3) == "TAG"))
//        {
//            ++correct;
//        }
//        else
//        {
//            continue;
//        }

        std::string concatenated_exon_seq = "";
        for (const auto & e : t.exons)
            concatenated_exon_seq += chr.substr(std::get<0>(e) - 1, std::get<1>(e) - std::get<0>(e) + 1);
        // printf("%s\n", concatenated_exon_seq.c_str());

        str_to_upper(concatenated_exon_seq);
        if (t.strand == '-')
        {
            std::reverse(concatenated_exon_seq.begin(), concatenated_exon_seq.end());
            for (uint64_t j = 0; j < concatenated_exon_seq.size(); ++j)
                concatenated_exon_seq[j] = complement(concatenated_exon_seq[j]);
        }
        auto orfs = get_all_ORFs(concatenated_exon_seq, t.strand);

        std::array<std::vector<std::vector<float> >, 4> extracted_scores; // phase 0, phase 1, phase 2, power
        compute_PhyloCSF_for_transcript(t, params.bw_files, extracted_scores);

        for (auto const & orf : orfs)
        {
            // compute CDS
            uint32_t len = 0;

            std::vector<cds_entry> CDS;

            uint32_t first_exon_id_in_CDS = 0;
            uint32_t last_exon_id_in_CDS = 0;
            for (auto const e : t.exons) // don't use t.exons because they already have donor and acceptor annotation
            {
                uint32_t const len_new = len + END(e) - BEGIN(e);

                cds_entry c(BEGIN(e), END(e), 3 /*invalid phase, overwrite later*/);

                if (len < BEGIN(orf))
                    c.begin += BEGIN(orf) - len;

                if (len_new > END(orf)) // TODO: > or >=
                    c.end -= len_new - END(orf) - 1;

                if (len_new >= BEGIN(orf) && len <= END(orf) && !(c.begin == c.end)) // TODO: of-by-one check!
                {
                    // c.type = "CDS";
                    CDS.push_back(std::move(c));
                    ++last_exon_id_in_CDS;
                }
                else if (CDS.size() == 0) // otherwise this will be called after the last CDS, if exons are remaining in the suffix
                {
                    ++first_exon_id_in_CDS;
                    ++last_exon_id_in_CDS;
                }

                len += END(e) - BEGIN(e);
            }
            --last_exon_id_in_CDS;

            // compute phases
            if (t.strand == '+')
                annotateCDSPhases(CDS.begin(), CDS.end());
            else
                annotateCDSPhases(CDS.rbegin(), CDS.rend());

            std::tuple<float, float> phylo_stats;

            if (t.strand == '+')
            {
                phylo_stats = compute_PhyloCFS(t.exons, CDS.begin(), CDS.end(), t.strand, extracted_scores, first_exon_id_in_CDS, last_exon_id_in_CDS, 0);
            }
            else if (t.strand == '-') // leave this if statement in here to skip transcripts without strand information
            {
                phylo_stats = compute_PhyloCFS(t.exons, CDS.rbegin(), CDS.rend(), t.strand, extracted_scores, first_exon_id_in_CDS, last_exon_id_in_CDS, concatenated_exon_seq.size());
            }

            t.phylo_score = std::get<0>(phylo_stats);
            t.phylo_power = std::get<1>(phylo_stats);

            if (t.phylo_score < 0.0)
                continue;

            // output gtf file (with annotation)
            bool first_processed_line = true;
            bool is_gff = true;
            for (auto line_tuple : t.lines)
            {
                const feature_t f = std::get<0>(line_tuple);
                const std::string & line = std::get<1>(line_tuple);

                // detect format
                if (first_processed_line && f == EXON)
                {
                    first_processed_line = false;
                    is_gff = is_gff_format(line);
                }

                if (f == TRANSCRIPT)
                {
                    if (is_gff)
                    {
                        if (params.comp_bls)
                            fprintf(gff_out, "%s;phylocsf_mean=%.3f;phylocsf_power_mean=%.3f\n", line.c_str(), t.phylo_score, t.phylo_power);
                        else
                            fprintf(gff_out, "%s;phylocsf_mean=%.3f\n", line.c_str(), t.phylo_score);
                    }
                    else
                    {
                        if (params.comp_bls)
                            fprintf(gff_out, "%s phylocsf_mean \"%.3f\"; phylocsf_power_mean \"%.3f\";\n", line.c_str(), t.phylo_score, t.phylo_power);
                        else
                            fprintf(gff_out, "%s phylocsf_mean \"%.3f\";\n", line.c_str(), t.phylo_score);
                    }
                }
                else
                {
                    fprintf(gff_out, "%s\n", line.c_str());
                }
            }
            // print computed CDS
            for (const auto c : CDS)
            {
                printf("CDS: %lld, %lld, %.3f\n", c.begin, c.end, c.phylo_score);

                fprintf(gff_out, "%s\tPhyloCSF++\tCDS\t%lld\t%lld\t.\t%c\t%d\t", t.chr.c_str(), c.begin, c.end, t.strand, c.phase);

                if (is_gff)
                {
                    if (params.comp_bls)
                        fprintf(gff_out, "phylocsf_mean=%.3f;phylocsf_power_mean=%.3f", t.phylo_score, t.phylo_power);
                    else
                        fprintf(gff_out, "phylocsf_mean=%.3f", t.phylo_score);
                }
                else
                {
                    if (params.comp_bls)
                        fprintf(gff_out, "phylocsf_mean \"%.3f\"; phylocsf_power_mean \"%.3f\";", t.phylo_score, t.phylo_power);
                    else
                        fprintf(gff_out, "phylocsf_mean \"%.3f\";", t.phylo_score);
                }

                fprintf(gff_out, "\n");
            }
            printf("------------------------\n");
//            const float new_score = (phylo_sum > 0) ? (phylo_sum / phylo_count) : NAN;
//
//            if (CDS.size() > 0 && new_score >= 0.0f)
//            {
//                std::tuple<float,float> phylo_power_stats = get_stats(phylo_scores_power);
//                _scores.emplace_back(new_score, std::get<2>(phylo_stats), CDS, std::get<3>(phylo_stats), std::get<0>(phylo_power_stats), std::get<1>(phylo_power_stats));
//
//                if (new_score > best_score)
//                {
//                    best_score = new_score;
//                }
//            }
        }
        printf("=======================-\n");

        reader.print_progress();
    }

    fclose(gff_out);

    // printf("%d / %d\n", correct, total);
}

int main_find_cds(int argc, char **argv)
{
    ArgParse args("phylocsf++ find-cds",
                  "TODO");

    FindCDSCLIParams params;

    args.add_option("output", ArgParse::Type::STRING, "Path where output GFF/GTF will be written to. If it does not exist, it will be created. Default: output files are stored next to the input files.", ArgParse::Level::GENERAL, false);

    args.add_positional_argument("fasta", ArgParse::Type::STRING, "Path to the genome fasta file.", true);
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

        params.bw_files[i] = bwOpen(const_cast<char * >(params.bw_path.c_str()), NULL, "r");
        if (!params.bw_files[i])
        {
            // check whether the user has used an unindexed wig file, then print a useful hint
            if (access(params.bw_path.c_str(), F_OK) == 0 &&
                params.bw_path.size() >= 4 && params.bw_path.compare(params.bw_path.size() - 4, 4, ".wig") == 0
            )
            {
                printf(OUT_ERROR "An error occurred while opening the PhyloCSF file '%s'.\n" OUT_RESET, params.bw_path.c_str());
                printf(OUT_INFO "It seems you provided a *.wig file. You need to simply index them first with wigToBigWig and then use the *.bw files.\n" OUT_RESET);
                return -1;
            }
            else if (i == 6)
            {
                params.comp_bls = false;
                printf(OUT_INFO "PhyloCSFpower.bw not found. Annotation will not have confidence scores.\n" OUT_RESET);
            }
        }
    }

    std::unordered_map<std::string, std::string> genome;
    load_fasta_file(args.get_positional_argument("fasta"), genome);
//    for (auto & chr : genome)
//    {
//        str_to_upper(chr.second);
//    }

    const chromList_t *chr_list = bwReadChromList(params.bw_files[0]);
    for (int64_t i = 0; i < chr_list->nKeys; ++i)
    {
        params.chrom_sizes.emplace(chr_list->chrom[i], chr_list->len[i]);
    }

    // run for every gff file
    std::unordered_set<std::string> missing_sequences;
    for (uint16_t i = 2; i < args.positional_argument_size(); ++i)
    {
        run_find_cds(args.get_positional_argument(i), params, genome, missing_sequences, i, args.positional_argument_size() - 2);
    }

    printf(OUT_DEL "Done!\n");

    return 0;
}
