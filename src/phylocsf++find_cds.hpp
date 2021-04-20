#include "common.hpp"
#include "arg_parse.hpp"

#define BEGIN(exon) std::get<0>(exon)
#define END(exon) std::get<1>(exon)

enum find_mode_t
{
    ALL,
    LONGEST,
    BEST_SCORE
};

struct FindCDSCLIParams
{
    std::string output_path = "";
    std::string bw_path = ""; // needed for header output

    bigWigFile_t *bw_files[7] = { NULL }; // in this order: +1, +2, +3, -1, -2, -3, power/confidence
    std::unordered_map<std::string, uint64_t> chrom_sizes;

    find_mode_t mode = LONGEST;
    float min_score = 0.0;
    uint32_t min_codons = 25;

    bool evaluate = false;
};

std::string mode_to_str(const find_mode_t mode) noexcept
{
    std::string mode_str;
    switch (mode)
    {
        case ALL:
            mode_str = "ALL";
            break;
        case LONGEST:
            mode_str = "LONGEST";
            break;
        case BEST_SCORE:
            mode_str = "BEST_SCORE";
            break;
    }
    return mode_str;
}

void find_all_codons(const std::string & dna, const std::string & codon, std::array<std::vector<uint32_t>, 3> & hits) noexcept
{
    size_t pos = dna.find(codon);
    while (pos != std::string::npos)
    {
        hits[pos % 3].push_back(pos);
        pos = dna.find(codon, pos + 1);
    }
}

std::vector<std::tuple<uint32_t, uint32_t> > get_all_ORFs(const std::string & spliced_transcript_seq_orig,
                                                          const char strand, const FindCDSCLIParams & params)
{
    std::vector<std::tuple<uint32_t, uint32_t> > orfs;
    std::array<std::vector<uint32_t>, 3> startpos;
    std::array<std::vector<uint32_t>, 3> stoppos;

    std::string spliced_transcript_seq = spliced_transcript_seq_orig;

    if (strand == '-')
    {
        std::reverse(spliced_transcript_seq.begin(), spliced_transcript_seq.end());
        for (uint64_t j = 0; j < spliced_transcript_seq.size(); ++j)
            spliced_transcript_seq[j] = complement(spliced_transcript_seq[j]);

        // printf("on pos strand: %s\n", spliced_transcript_seq_orig.c_str());
        // printf("on neg strand: %s\n", spliced_transcript_seq.c_str());
    }

    // seq is already in reverse complement form when we are on the neg. strand
    find_all_codons(spliced_transcript_seq, "ATG", startpos);
    find_all_codons(spliced_transcript_seq, "TAA", stoppos);
    find_all_codons(spliced_transcript_seq, "TAG", stoppos);
    find_all_codons(spliced_transcript_seq, "TGA", stoppos);

//    for (uint8_t i = 0; i < 3; ++i)
//    {
//        printf("Start[%d]: ", i);
//        for (auto x : startpos[i])
//            printf("%d, ", x);
//        printf("\n");
//
//        printf("Stop[%d]: ", i);
//        for (auto x : stoppos[i])
//            printf("%d, ", x);
//        printf("\n");
//    }

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

                    uint32_t orf_length = stop - start + 1;
                    assert(orf_length % 3 == 0);
                    if (3 * params.min_codons <= orf_length)
                    {
                        orfs.emplace_back(start, stop);
                        // auto orf_seq = spliced_transcript_seq.substr(start, orf_length);
                        // printf("%d\t%d\t%s\n", start, stop, orf_seq.c_str());
                    }
                }
                else if (strand == '-')
                {
                    uint32_t stop_rev = spliced_transcript_seq.size() - start - 3;
                    stop_rev += 2; // include stop codon in sequence
                    uint32_t start_rev = spliced_transcript_seq.size() - stop - 3;

                    uint32_t orf_length = stop_rev - start_rev + 1;
                    assert(orf_length % 3 == 0);

                    if (3 * params.min_codons <= orf_length)
                    {
                        orfs.emplace_back(start_rev, stop_rev);
                        // auto orf_seq_pos = spliced_transcript_seq.substr(start, orf_length);
                        // auto orf_seq_neg = spliced_transcript_seq_orig.substr(start_rev, orf_length);
                        // printf("on pos: %s\n", orf_seq_pos.c_str());
                        // printf("on neg: %s\n", orf_seq_neg.c_str());
                    }
                }
            }
        }
    }

    return orfs;
}

template <typename iter_t>
void annotateCDSPhases(iter_t it, iter_t end) noexcept
{
    // phase = last CDS has `phase` bases extra
    // c.phase = new CDS has `phase2` bases that we skip
    uint8_t phase = 0;
    for (; it != end; ++it)
    {
        cds_entry & c = *it;
        c.phase = ((3 - phase) % 3);
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
        if (t.strand == '-')
            std::reverse(extracted_scores[3][exon_id].begin(), extracted_scores[3][exon_id].end());

        if (t.strand == '+')
        {
            for (uint8_t phase = 0; phase < 3; ++phase)
            {
                // (cds_phase + exon.beginPos) % 3 == phase
                intervals = bwGetValues(wigFiles[0 + phase], (char*)t.chr.c_str(), BEGIN(exon), END(exon), 0);
                if (intervals != NULL)
                {
                    // weighten the scores with the power/confidence
                    for (uint32_t i = 0; i < (*intervals).l; ++i)
                    {
                        if (extracted_scores[3][exon_id][(*intervals).start[i] - BEGIN(exon)] == -999.0f)
                        {
                            extracted_scores[phase][exon_id][(*intervals).start[i] - BEGIN(exon)] = -999.0f;
                        }
                        else
                        {
                            extracted_scores[phase][exon_id][(*intervals).start[i] - BEGIN(exon)] = (*intervals).value[i] * extracted_scores[3][exon_id][(*intervals).start[i] - BEGIN(exon)];
                        }
                    }
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
                    // weighten the scores with the power/confidence
                    for (int32_t i = (*intervals).l - 1; i >= 0; --i)
                    {
                        if (extracted_scores[3][exon_id][(END(exon) - 1) - (*intervals).start[i]] == -999.0f)
                        {
                            extracted_scores[phase][exon_id][(END(exon) - 1) - (*intervals).start[i]] = -999.0f;
                        }
                        else
                        {
                            extracted_scores[phase][exon_id][(END(exon) - 1) - (*intervals).start[i]] = (*intervals).value[i] * extracted_scores[3][exon_id][(END(exon) - 1) - (*intervals).start[i]];
                        }
                    }
                    bwDestroyOverlappingIntervals(intervals);
                    intervals = NULL;
                }
            }
        }
    }
}

template <typename TExons, typename TIter>
std::tuple<float, float>
compute_PhyloCSF(TExons & exons, TIter first, TIter last, char const strand,
                 const std::array<std::vector<std::vector<float> >, 4> & extracted_scores,
                 const uint32_t first_exon_id_in_CDS, const uint32_t last_exon_id_in_CDS, const uint32_t chrLen)
{
    float total_phylo_sum = 0.0;
    float total_power_sum = 0.0;
    uint32_t total_phylo_count = 0;
    uint32_t total_power_count = 0;

    uint32_t cds_id = 0;

    for (auto it = first; it != last; ++it, ++cds_id)
    {
        cds_entry & c = *it;

        // new_phase = last CDS has `new_phase` bases extra
        // new_phase2 = new CDS has `new_phase2` bases that we skip (because we exclude them since we already retrieved score from CDS entry before)
        // uint32_t new_phase2 = (3 - new_phase) % 3;
        assert(c.phase < 3);

        uint32_t exon_id;
        uint32_t phylo_start, phylo_end;
        const std::vector<float> * phased_score = NULL;

        // compute phylo score
        if (strand == '+')
        {
            exon_id = first_exon_id_in_CDS + cds_id;
            phased_score = &extracted_scores[(c.phase + c.begin) % 3][exon_id];
            phylo_start = c.begin - BEGIN(exons[exon_id]);
            phylo_end = END(exons[exon_id]) - c.end;
        }
        else // if (strand == '-')
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
            phased_score = &extracted_scores[(chrLen - c.end - 1 + c.phase + 1) % 3][exon_id];
            phylo_start = END(exons[exon_id]) - c.end;
            phylo_end = c.begin - BEGIN(exons[exon_id]);
        }

        float phylo_sum = 0.0;
        uint32_t phylo_count = 0;
        for (auto phylo_it = (*phased_score).begin() + phylo_start; phylo_it != (*phased_score).end() - phylo_end; ++phylo_it)
        {
            if (*phylo_it != -999.0f)
            {
                phylo_sum += *phylo_it;
                ++phylo_count;
            }
        }
        total_phylo_count += phylo_count;

        auto & power_scores = extracted_scores[3][exon_id];
        float power_sum = 0.0;
        uint32_t power_count = power_scores.size() - phylo_end - phylo_start;
        total_power_count += power_count;
        for (auto phylo_it = power_scores.begin() + phylo_start; phylo_it != power_scores.end() - phylo_end; ++phylo_it)
        {
            if (*phylo_it != -999.0f)
                power_sum += *phylo_it;
        }

        c.phylo_score = (phylo_count > 0) ? phylo_sum / power_sum : NAN; // weighten the scores with the power/confidence (or not)
        c.phylo_power = (power_count > 0) ? power_sum / power_count : NAN;

        total_phylo_sum += phylo_sum;
        total_power_sum += power_sum;
    }

    return std::make_tuple(
        (total_phylo_count > 0) ? total_phylo_sum / total_power_sum : NAN, // weighten the scores with the power/confidence (or not)
        (total_power_count > 0) ? total_power_sum / total_power_count : NAN // total phylo_power
    );
}

void output_transcript(const gff_transcript & t, const std::vector<cds_entry> & CDS, FILE * gff_out)
{
    bool first_processed_line = true;
    bool is_gff = true;
    for (auto line_tuple : t.lines)
    {
        const feature_t f = std::get<0>(line_tuple);
        const std::string & line = std::get<1>(line_tuple);

        // detect format
        if (first_processed_line && f == TRANSCRIPT)
        {
            first_processed_line = false;
            is_gff = is_gff_format(line);
        }

        if (f == TRANSCRIPT)
        {
            if (is_gff)
            {
                fprintf(gff_out, "%s;phylocsf_weighted_mean=%.3f", line.c_str(), t.phylo_score);
                fprintf(gff_out, ";phylocsf_power_mean=%.3f", t.phylo_power);
            }
            else
            {
                fprintf(gff_out, "%s phylocsf_weighted_mean \"%.3f\";", line.c_str(), t.phylo_score);
                fprintf(gff_out, " phylocsf_power_mean \"%.3f\";", t.phylo_power);
            }
            fprintf(gff_out, "\n");
        }
        else if (f != feature_t::CDS) // do not output any original CDS entries from the input
        {
            fprintf(gff_out, "%s\n", line.c_str());
        }
    }
    // print computed CDS
    for (const auto c : CDS)
    {
        // transform 0-based indices from seqan back to 1-based indices of gff
        fprintf(gff_out, "%s\tPhyloCSF++\tCDS\t%" PRIu64 "\t%" PRIu64 "\t.\t%c\t%d\t", t.chr.c_str(), c.begin + 1, c.end, t.strand, c.phase);

        if (is_gff)
        {
            fprintf(gff_out, "phylocsf_weighted_mean=%.3f", c.phylo_score);
            fprintf(gff_out, ";phylocsf_power_mean=%.3f", c.phylo_power);
        }
        else
        {
            fprintf(gff_out, "phylocsf_weighted_mean \"%.3f\";", c.phylo_score);
            fprintf(gff_out, " phylocsf_power_mean \"%.3f\";", c.phylo_power);
        }

        fprintf(gff_out, "\n");
    }
}

bool is_stop_coding(const std::string & codon) noexcept
{
    assert(codon.size() == 3);
    return codon == "TAA" || codon == "TGA" || codon == "TAG";
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

    fprintf(gff_out, "# CDS predicted with PhyloCSF++ %sles (%s, %s) and precomputed tracks %s "
                     "(options: --mode %s --min-score %f --min-codons %" PRIu32 ")\n",
                     ArgParse::version.c_str(), ArgParse::git_hash.c_str(), ArgParse::git_date.c_str(),
                     params.bw_path.c_str(),
                     mode_to_str(params.mode).c_str(),
                     params.min_score, params.min_codons);

    std::string concatenated_exon_seq;

    uint64_t nbr_transcripts = 0;
    uint64_t nbr_transcripts_with_annotated_orf = 0;
    uint64_t nbr_transcripts_with_annotated_orf_satisfies_criteria = 0;
    uint64_t nbr_transcripts_with_annotated_orf_matches_predicted_stop = 0;
    uint64_t nbr_transcripts_with_annotated_orf_matches_predicted_start_and_stop = 0;
    uint64_t nbr_transcripts_with_annotated_invalid_orf = 0;
    uint64_t nbr_transcripts_without_annotated_orf_with_orf_predicted = 0;

    std::string annotated_cds_seq, computed_cds_seq;
    std::vector<std::tuple<float, std::string, std::vector<cds_entry> > > hits;

    gff_transcript t;
    while (reader.get_next_transcript<true>(t, true, true))
    {
        ++nbr_transcripts;

        // change 1-based gff-coordinates to 0-based gff-coordinates (like seqan)
        for (auto & e : t.exons)
        {
            std::get<0>(e)--;
        }

        // check whether chromosome exists in tracks
        const auto chrom_sizes_it = params.chrom_sizes.find(t.chr);
        if (chrom_sizes_it == params.chrom_sizes.end())
        {
            // only report a missing sequence once
            if (missing_sequences.find(t.chr) == missing_sequences.end())
            {
                missing_sequences.insert(t.chr);
                printf(OUT_DEL OUT_INFO "Sequence '%s' from the GFF file does not occur in the tracks. Skipping ...\n" OUT_RESET, t.chr.c_str());
            }
            continue;
        }
        const uint64_t chr_len = chrom_sizes_it->second;

        const auto genome_it = genome.find(t.chr);
        if (genome_it == genome.end())
        {
            // only report a missing sequence once
            if (missing_sequences.find(t.chr) == missing_sequences.end())
            {
                missing_sequences.insert(t.chr);
                printf(OUT_DEL OUT_INFO "Sequence '%s' from the GFF file does not occur in the genome. Skipping ...\n" OUT_RESET, t.chr.c_str());
            }
            continue;
        }
        const std::string & chr = genome_it->second;

        // TODO: run it with both strands
        if (t.strand != '+' && t.strand != '-')
        {
            printf(OUT_DEL OUT_INFO "Skipped a transcript with no strand information.\n" OUT_RESET);
            continue;
        }

        annotated_cds_seq = "";
        if (params.evaluate && t.CDS.size() > 0)
        {
            // NOTE: CDS coordinates from input file are 1-based, we have to decrement the begin position
            for (uint16_t i = 0; i < t.CDS.size(); ++i)
            {
                const auto & c = t.CDS[i];
                // the CDS can have 1-2 bases before that START codon (i.e., phase != 0): we need to remove those bases
                if (i == 0 && t.strand == '+' && c.phase != 0)
                    annotated_cds_seq += chr.substr(c.begin + c.phase - 1, c.end - (c.begin + c.phase) + 1);
                else if (i == t.CDS.size() - 1 && t.strand == '-' && c.phase != 0)
                    annotated_cds_seq += chr.substr(c.begin - 1, c.end - c.phase - c.begin + 1);
                else
                    annotated_cds_seq += chr.substr(c.begin - 1, c.end - c.begin + 1);
            }

            // compute reverse complement for transcripts on the negative strand
            str_to_upper(annotated_cds_seq);
            if (t.strand == '-')
            {
                std::reverse(annotated_cds_seq.begin(), annotated_cds_seq.end());
                for (uint64_t j = 0; j < annotated_cds_seq.size(); ++j)
                    annotated_cds_seq[j] = complement(annotated_cds_seq[j]);
            }

            // if STOP codon is not part of the last resp. first CDS, append it
            const std::string last_codon = annotated_cds_seq.substr(annotated_cds_seq.size() - 3);
            if (!is_stop_coding(last_codon))
            {
                std::string extended_codon = "";
                if (t.strand == '+')
                {
                    if (t.CDS.back().end + 3 - 1 < chr.size())
                    {
                        extended_codon = chr.substr(t.CDS.back().end, 3);
                        str_to_upper(extended_codon);
                    }
                }
                else // t.strand == '-'
                {
                    if (t.CDS.front().begin >= 4)
                    {
                        extended_codon = chr.substr(t.CDS.front().begin - 4, 3);
                        str_to_upper(extended_codon);
                        std::reverse(extended_codon.begin(), extended_codon.end());
                        for (uint64_t j = 0; j < extended_codon.size(); ++j)
                            extended_codon[j] = complement(extended_codon[j]); // TODO: test this!
                    }
                }

                if (is_stop_coding(extended_codon))
                    annotated_cds_seq += extended_codon;
            }

            ++nbr_transcripts_with_annotated_orf;
            if (!(annotated_cds_seq.size() % 3 == 0 && annotated_cds_seq.size() >= 3
                && annotated_cds_seq.substr(0, 3) == "ATG"
                && is_stop_coding(annotated_cds_seq.substr(annotated_cds_seq.size() - 3))))
            {
                ++nbr_transcripts_with_annotated_invalid_orf;
                // printf("\nWARNING: %s, %s, %ld, %c, %s, %lld, %lld, %ld\n", annotated_cds_seq.substr(0, 3).c_str(), annotated_cds_seq.substr(annotated_cds_seq.size() - 3).c_str(), annotated_cds_seq.size() % 3, t.strand,
                //        t.chr.c_str(), t.begin, t.end, t.CDS.size());
                // continue;
            }
        }

        concatenated_exon_seq = "";
        for (const auto & e : t.exons)
        {
            // printf("Exon: %s\n", chr.substr(std::get<0>(e) - 1, std::get<1>(e) - std::get<0>(e) + 1).c_str());
            concatenated_exon_seq += chr.substr(std::get<0>(e), std::get<1>(e) - std::get<0>(e));
        }
        str_to_upper(concatenated_exon_seq);

        auto orfs = get_all_ORFs(concatenated_exon_seq, t.strand, params);

        // sort by length/score descending. hence, we can abort as first as we see the first matching ORF
        if (params.mode == LONGEST)
        {
            std::sort(orfs.begin(), orfs.end(), [](const std::tuple<uint32_t, uint32_t> & x, const std::tuple<uint32_t, uint32_t> & y) {
                return std::get<1>(x) - std::get<0>(x) > std::get<1>(y) - std::get<0>(y);
            });
        }

        std::array<std::vector<std::vector<float> >, 4> extracted_scores; // phase 0, phase 1, phase 2, power
        compute_PhyloCSF_for_transcript(t, params.bw_files, extracted_scores);

        hits.clear();

        std::string longest_or_best_seq = "";
        bool CDS_outputted_or_found = false;
        bool found_an_orf_satisfying_criteria = false;
        float best_CDS_score = -999.0f;
        std::vector<cds_entry> longest_or_best_CDS;
        std::tuple<float, float> longest_or_best_phylo_stats;
        for (auto const & orf : orfs)
        {
            // orf: 0-based closed intervals. E.g., orf = [0..10] is of length 11
            // gtf: 1-based closed intervals. E.g. first triplet is [1,3]

            // compute CDS
            uint32_t len = 0;
            uint32_t first_exon_id_in_CDS = 0;
            uint32_t last_exon_id_in_CDS = 0;
            std::vector<cds_entry> CDS;
            for (auto const e : t.exons)
            {
                uint32_t const len_new = len + END(e) - BEGIN(e);

                cds_entry c(BEGIN(e), END(e), 3 /*invalid phase, overwrite later*/);

                // beginning of exon is not part of the ORF
                if (len < BEGIN(orf))
                    c.begin += BEGIN(orf) - len;

                // end of exon is not part of the ORF
                if (len_new > END(orf))
                    c.end -= len_new - END(orf) - 1;

                // does orf overlap with exon? check for empty exons (I have seen empty exons before somewhere ...)
                if (BEGIN(orf) <= len_new && len <= END(orf) && c.begin < c.end)
                {
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

            std::tuple<float, float> phylo_stats;
            if (t.strand == '+')
            {
                annotateCDSPhases(CDS.begin(), CDS.end());
                phylo_stats = compute_PhyloCSF(t.exons, CDS.begin(), CDS.end(), t.strand, extracted_scores, first_exon_id_in_CDS, last_exon_id_in_CDS, 0);
            }
            else
            {
                annotateCDSPhases(CDS.rbegin(), CDS.rend());
                phylo_stats = compute_PhyloCSF(t.exons, CDS.rbegin(), CDS.rend(), t.strand, extracted_scores, first_exon_id_in_CDS, last_exon_id_in_CDS, chr_len);
            }

            // ignore transcripts with neg. score or that contain a negative CDS
            // bool consider_transcript_not_positive = false;
            // if (t.phylo_score <= 0.0 || isnan(t.phylo_score))
            //     consider_transcript_not_positive = true;
            // for (const auto & c : CDS)
            // {
            //     if (c.phylo_score <= 0.0 || isnan(c.phylo_score))
            //         consider_transcript_not_positive = true;
            // }

            // output gtf file (with annotation)
            const float phylo_score = std::get<0>(phylo_stats);
            if (phylo_score >= params.min_score) // ORF meets criteria
            {
                found_an_orf_satisfying_criteria = true;

                if (params.evaluate)
                {
                    computed_cds_seq = "";
                    // NOTE: CDS stem from exons, and thus have been decremented to 0-based coordinates
                    for (const auto & c : CDS)
                        computed_cds_seq += chr.substr(c.begin, c.end - c.begin);
                    str_to_upper(computed_cds_seq);
                    if (t.strand == '-')
                    {
                        std::reverse(computed_cds_seq.begin(), computed_cds_seq.end());
                        for (uint64_t j = 0; j < computed_cds_seq.size(); ++j)
                            computed_cds_seq[j] = complement(computed_cds_seq[j]);
                    }
                    hits.emplace_back(phylo_score, computed_cds_seq, CDS);
                }

                if (params.mode == BEST_SCORE) // remember what ORF has to best score
                {
                    CDS_outputted_or_found = true;
                    if (phylo_score > best_CDS_score)
                    {
                        best_CDS_score = phylo_score;
                        longest_or_best_CDS = CDS;
                        longest_or_best_phylo_stats = phylo_stats;
                        if (params.evaluate)
                            longest_or_best_seq = computed_cds_seq;
                    }
                }
                else if (params.mode == ALL)
                {
                    CDS_outputted_or_found = true;
                    t.phylo_score = std::get<0>(phylo_stats);
                    t.phylo_power = std::get<1>(phylo_stats);
                    output_transcript(t, CDS, gff_out);
                }
                else if (params.mode == LONGEST && !CDS_outputted_or_found)
                {
                    CDS_outputted_or_found = true;
                    longest_or_best_CDS = CDS;
                    longest_or_best_phylo_stats = phylo_stats;
                    if (params.evaluate)
                        longest_or_best_seq = computed_cds_seq;

                    if (params.mode == LONGEST && !params.evaluate) // if we evaluate we need to see all ORFs to check whether the annotated one is found (can be made faster!)
                        break; // we are done (ORFs are sorted by length in descending order)
                }
            }
        }

        if (CDS_outputted_or_found && (params.mode == LONGEST || params.mode == BEST_SCORE))
        {
            t.phylo_score = std::get<0>(longest_or_best_phylo_stats);
            t.phylo_power = std::get<1>(longest_or_best_phylo_stats);
            output_transcript(t, longest_or_best_CDS, gff_out);

            if (longest_or_best_seq == annotated_cds_seq)
                ++nbr_transcripts_with_annotated_orf_matches_predicted_start_and_stop;

            size_t suffix_length = std::min(annotated_cds_seq.size(), longest_or_best_seq.size());
            if (annotated_cds_seq != "" && annotated_cds_seq.substr(std::max<size_t>(0, annotated_cds_seq.size() - suffix_length)) == longest_or_best_seq.substr(std::max<size_t>(0, longest_or_best_seq.size() - suffix_length)))
                ++nbr_transcripts_with_annotated_orf_matches_predicted_stop;
        }
        else if (!CDS_outputted_or_found)
        {
            output_transcript(t, std::vector<cds_entry>{}, gff_out);
        }

        if (found_an_orf_satisfying_criteria && annotated_cds_seq == "")
        {
            ++nbr_transcripts_without_annotated_orf_with_orf_predicted;
        }

        if (std::find_if(hits.begin(), hits.end(), [&annotated_cds_seq](const std::tuple<float, std::string, std::vector<cds_entry> > & x) { return std::get<1>(x) == annotated_cds_seq; }) != hits.end())
        {
            ++nbr_transcripts_with_annotated_orf_satisfies_criteria;
        }

        reader.print_progress();
    }

    fclose(gff_out);

    if (params.evaluate)
    {
        printf("%-73s %6" PRIu64 "\n", "Transcripts in total:", nbr_transcripts);
        printf("--------------------------------------------------------\n");
        printf("%-73s %6" PRIu64 "\n", "Transcripts with annotated CDS:", nbr_transcripts_with_annotated_orf);
        printf("%-73s %6" PRIu64 "\n", "- annotated CDS is a proper CDS (start/stop codon, no frame shifts,...):", nbr_transcripts_with_annotated_orf - nbr_transcripts_with_annotated_invalid_orf);
        printf("%-73s %6" PRIu64 "\n", "- annotated CDS satisfies criteria (min-codons and min-score):", nbr_transcripts_with_annotated_orf_satisfies_criteria);
        if (params.mode != ALL)
        {
            printf("%-73s %6" PRIu64 "\n", "- predicted ORF matches annotated ORF (same frame and stop position):", nbr_transcripts_with_annotated_orf_matches_predicted_stop);
            printf("%-73s %6" PRIu64 "\n", "- predicted ORF matches annotated ORF (same start and stop position):", nbr_transcripts_with_annotated_orf_matches_predicted_start_and_stop);
        }
        printf("--------------------------------------------------------------------------------\n");
        const uint64_t nbr_transcripts_without_annotated_cds = nbr_transcripts - nbr_transcripts_with_annotated_orf;
        printf("%-73s %6" PRIu64 "\n", "Transcripts without annotated CDS:", nbr_transcripts_without_annotated_cds);
        printf("%-73s %6" PRIu64 "\n", "- Transcripts with an ORF predicted:", nbr_transcripts_without_annotated_orf_with_orf_predicted);
        //printf("%-73s %6" PRIu64 "\n", "- Transcripts with no ORF predicted:", nbr_transcripts_without_annotated_orf - nbr_transcripts_without_annotated_orf_with_orf_predicted);
    }
}

int main_find_cds(int argc, char **argv)
{
    ArgParse args("phylocsf++ find-cds",
                  "Takes a GFF file with transcripts and scores every possible ORF to find novel \n"
                  "protein-coding genes and outputs a GFF file with annotated CDS features.\n"
                  "Either all CDS meeting the criteria can be outputted, or the longest or \n"
                  "best-scoring one.\n"
                  "'find-cds' and especially '--evaluate' is experimental. It is recommended to \n"
                  "only have 'transcript', 'exon', and optionally 'CDS' features (see --evaluate) \n"
                  "in the GFF file. Other features can mess up the output.");

    FindCDSCLIParams params;

    std::string default_mode = mode_to_str(params.mode);

    args.add_option("output", ArgParse::Type::STRING, "Path where output GFF/GTF will be written to. If it does not exist, it will be created. Default: output files are stored next to the input files.", ArgParse::Level::GENERAL, false);
    args.add_option("mode", ArgParse::Type::STRING, "Which CDS to report: ALL, LONGEST, BEST_SCORE. Default: " + default_mode, ArgParse::Level::GENERAL, false);
    args.add_option("min-score", ArgParse::Type::FLOAT, "Only consider ORFs with a minimum weighted PhyloCSF mean score (range from -15 to +15, >0 more likely to be protein-coding). Default: " + std::to_string(params.min_score), ArgParse::Level::GENERAL, false);
    args.add_option("min-codons", ArgParse::Type::INT, "Only consider ORFs with a minimum codon length. Default: " + std::to_string(params.min_codons), ArgParse::Level::GENERAL, false);
    args.add_option("evaluate", ArgParse::Type::BOOL, "Additionally reads CDS features in input GFF/GTF and outputs statistics how well PhyloCSF++ has predicted CDS. This does not affect the output GFF file, i.e., CDS from the input file are not included, only predicted CDS are outputted (regardless whether they matched the CDS in the input or not). Default: " + std::to_string(params.evaluate), ArgParse::Level::GENERAL, false);

    args.add_positional_argument("fasta", ArgParse::Type::STRING, "Path to the genome fasta file.", true);
    args.add_positional_argument("tracks", ArgParse::Type::STRING, "Path to the bigWig file PhyloCSF+1.bw (expects the other 5 frames to be in the same directory, optionally the power track).", true);
    args.add_positional_argument("gff-files", ArgParse::Type::STRING, "One or more GFF/GTF files with coding exons to be scored.", true, true);
    args.parse_args(argc, argv);

    args.check_args();

    if (args.is_set("mode"))
    {
        std::string mode = args.get_string("mode");
        std::transform(mode.begin(), mode.end(), mode.begin(), toupper);

        if (mode == "ALL")
            params.mode = ALL;
        else if (mode == "LONGEST")
            params.mode = LONGEST;
        else if (mode == "BEST_SCORE")
            params.mode = BEST_SCORE;
        else
        {
            printf(OUT_ERROR "Please choose a valid mode (ALL, LONGEST, BEST_SCORE)!\n" OUT_RESET);
            return -1;
        }
    }

    // retrieve flags/options that were set
    if (args.is_set("output"))
    {
        params.output_path = args.get_string("output");
        // create a directory if it doesn't exist yet
        if (create_directory(params.output_path))
            printf(OUT_INFO "Created the output directory.\n" OUT_RESET);
    }

    if (args.is_set("min-codons"))
        params.min_codons = args.get_int("min-codons");

    if (args.is_set("min-score"))
        params.min_score = args.get_float("min-score");

    if (args.is_set("evaluate"))
        params.evaluate = args.get_bool("evaluate");

    if (params.evaluate)
        printf(OUT_INFO "Warning: the --evaluate mode is experimental.\n" OUT_RESET);

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

    std::unordered_map<std::string, std::string> genome;
    load_fasta_file(args.get_positional_argument("fasta"), genome);

    const chromList_t *chr_list = bwReadChromList(params.bw_files[0]);
    for (int64_t i = 0; i < chr_list->nKeys; ++i)
    {
        params.chrom_sizes.emplace(chr_list->chrom[i], chr_list->len[i]);
    }

    // run for every gff file
    std::unordered_set<std::string> missing_sequences;
    for (uint16_t i = 2; i < args.positional_argument_size(); ++i)
    {
        run_find_cds(args.get_positional_argument(i), params, genome, missing_sequences, i - 1, args.positional_argument_size() - 2);
    }

    printf(OUT_DEL "Done!\n");

    return 0;
}
