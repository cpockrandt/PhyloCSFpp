#include "run.hpp"

#include "src/estimate_hmm_parameter.hpp"
#include "src/create_tracks.hpp"

//#ifdef OPENMP
#include <omp.h>
//#endif

void my_fprintf(FILE* f, char *format_str, const float d)
{
    char buf[12];
    sprintf(buf, format_str, d);
    for (int8_t i = strlen(buf); i >= 0; --i)
    {
        if (isdigit(buf[i]))
        {
            if (buf[i] != '0')
                break;
            else
                buf[i] = 0; // cut off the 0
        }
    }
    fprintf(f, "%s\n", buf);
}

int main(int /*argc*/, char ** /*argv*/)
{
    unsigned threads = 1;
    unsigned jobs = 1;

//    uint32_t genome_length = 3252208893;
//    char model_str[] = "58mammals";
//    char selected_species_str[] = ""; // "Human,Chimp,Gorilla,Orangutan,Gibbon,Rhesus,Crab_eating_macaque,Baboon,Green_monkey,Marmoset,Squirrel_monkey,Bushbaby,Chinese_tree_shrew,Squirrel,Lesser_Egyptian_jerboa,Prairie_vole,Chinese_hamster,Golden_hamster,Mouse,Rat,Naked_mole_rat,Guinea_pig,Chinchilla,Brush_tailed_rat,Rabbit,Pika,Pig,Alpaca,Wild_bactrian_camel,Dolphin,Killer_whale,Tibetan_antelope,Cow,Sheep,Domestic_goat,Horse,White_rhinoceros,Cat,Dog,Ferret,Panda,Pacific_walrus,Weddell_seal,Black_flying_fox,Megabat,Big_brown_bat,Davids_myotis,Microbat,Hedgehog,Shrew,Star_nosed_mole,Elephant,Cape_elephant_shrew,Manatee,Cape_golden_mole,Tenrec,Aardvark,Armadillo"; // "Dog,Cow,Horse,Human,Mouse,Rat";
//    char aln_path[] = "/home/chris/dev-uni/PhyloCSF++/chr4_GL000008v2_random.maf"; // chr4_GL000008v2_random.maf"; chr4_small // "/home/chris/Downloads/chr22.516alignments.maf"; // chr22.head.maf";
////    char aln_path[] = "/home/chris/dev-uni/PhyloCSF++/ALDH2.exon5.maf";
////    char scores_path[] = "/home/chris/dev-uni/PhyloCSF++/chr22.head.orig.results.formatted";
//    char hmm_data_path[] = "/home/chris/dev-uni/PhyloCSF++/HumanCodingExonsV29.txt";

    uint32_t genome_length = 1065365434;
    char model_str[] = "53birds";
    char selected_species_str[] = "";
    char aln_path[] = "/home/chris/dev-uni/PhyloCSF++/phylo_galgal/galGal6_chrM.maf";
    char hmm_data_path[] = "/home/chris/dev-uni/PhyloCSF++/phylo_galgal/ChickenCodingExonsV2.txt";
    const std::string output_folder = "/home/chris/dev-uni/PhyloCSF++/galgal_out";

// NOTE: FINAL YEAST
//    uint32_t genome_length = 12157105;
//    char model_str[] = "7yeast";
//    char selected_species_str[] = "";
//    char aln_path[] = "/home/chris/dev-uni/PhyloCSF++/phylo_yeast/chrIV.maf";
//    char hmm_data_path[] = "/home/chris/dev-uni/PhyloCSF++/phylo_yeast/YeastCodingExons.txt";
//    const std::string output_folder = "/home/chris/dev-uni/PhyloCSF++/yeast_out";

    const hmm_parameter hmm_param = estimate_hmm_params_for_genome(hmm_data_path, genome_length);
    const hmm hmm = get_coding_hmm(hmm_param);

    std::vector<Data> data_fixed_mle(threads);
    char delim[] = ",";

    char *ptr = strtok(selected_species_str, delim);
    while (ptr != NULL)
    {
        for (size_t thread_id = 0; thread_id < threads; ++thread_id)
        {
            data_fixed_mle[thread_id].selected_species.emplace(ptr);
        }
        ptr = strtok(NULL, delim);
    }

    for (size_t thread_id = 0; thread_id < threads; ++thread_id)
    {
        // TODO: check whether selected species exist!
        data_fixed_mle[thread_id].load_model(model_str, false);
    }

    std::unordered_map<std::string, uint16_t> fastaid_to_alnid;
    // we use the order of data.phylo_array
    {
        FILE *file = fopen("/home/chris/dev-uni/PhyloCSF++/commonNames_assemblies.txt", "r");
        if (file != NULL)
        {
            char line[BUFSIZ];
            char c1[BUFSIZ];
            char c2[BUFSIZ];
            while (fgets(line, sizeof line, file) != NULL)
            {
                sscanf(line, "%s\t%s", c1, c2); // PhyloCSF identifier => maf identifier

//                bool found = false;
                for (uint16_t i = 0; i < data_fixed_mle[0].phylo_array.size(); ++i)
                {
                    if (strcmp(data_fixed_mle[0].phylo_array[i].label.c_str(), c1) == 0)
                    {
                        // remove leading digits (e.g., galGal6 -> galGal)
                        for (size_t c2_i = 0; c2_i < strlen(c2); ++c2_i)
                        {
                            if (isdigit(c2[c2_i]))
                            {
                                c2[c2_i] = 0;
                                break;
                            }
                        }

                        fastaid_to_alnid.emplace(c2, i); // maf identifier => id in vector for ids and seqs
                        fastaid_to_alnid.emplace(c1, i); // maf identifier => id in vector for ids and seqs
//                        found = true;
                        break;
                    }
                }
//                if (!found)
//                    printf("ERROR: %20s\t%10s\tmapping missing!\n", c1, c2);
            }
        }
        fclose(file);

        // TODO: check whether there are mappings missing (not sure) or superfluous mappings (also not sure)
    }

    // prepare alignment
    std::vector<alignment_t> alignments;
    for (unsigned i = 0; i < threads; ++i)
    {
        alignments.emplace_back(data_fixed_mle[0].phylo_array); // TODO: remove id's outside of aln struct because they are read-only
    }

    std::vector<std::vector<double> > lpr_per_codon(threads), bls_per_bp(threads);

    parallel_maf_reader maf_rd(aln_path, jobs, &fastaid_to_alnid);
    jobs = maf_rd.get_jobs(); // maybe file is too small and a smaller number of jobs is used

    FILE ** file_power = (FILE**) malloc(jobs * sizeof(FILE*));
    for (size_t job_id = 0; job_id < jobs; ++job_id)
    {
        const std::string filename_power = output_folder + "/PhyloCSFpower.wig." + std::to_string(job_id);
        file_power[job_id] = fopen(filename_power.c_str(), "w");
        if (file_power[job_id] == NULL)
        {
            printf("Error!");
            exit(1);
        }
    }

    FILE *** file_score_raw = (FILE***) malloc(jobs * sizeof(FILE**));
    FILE *** file_score = (FILE***) malloc(jobs * sizeof(FILE**));
    for (size_t job_id = 0; job_id < jobs; ++job_id)
    {
        file_score_raw[job_id] = (FILE**) malloc(6 * sizeof(FILE*));
        file_score[job_id] = (FILE**) malloc(6 * sizeof(FILE*));

        for (uint8_t i = 0; i < 6; ++i)
        {
            char strand = (i < 3) ? '+' : '-';
            unsigned frame = (i % 3) + 1;
            const std::string filename_score_raw = output_folder + "/PhyloCSFRaw" + std::string(1, strand) + std::to_string(frame) + ".wig." + std::to_string(job_id);
            const std::string filename_score     = output_folder + "/PhyloCSF" + std::string(1, strand) + std::to_string(frame) + ".wig." + std::to_string(job_id);
            file_score_raw[job_id][i] = fopen(filename_score_raw.c_str(), "w");
            file_score[job_id][i]= fopen(filename_score.c_str(), "w");

            if (file_score_raw[job_id][i] == NULL || file_score[job_id][i] == NULL)
            {
                printf("Error creating file!");
                exit(1);
            }
        }
    }

//    #pragma omp parallel for num_threads(threads) default(none) shared(jobs, alignments, maf_rd, data_fixed_mle, data_omega, frame, stdout, computed_results, lpr_per_codon, bls_per_codon)
    for (unsigned job_id = 0; job_id < jobs; ++job_id) // TODO: split it in more parts than there are threads
    {
        unsigned thread_id = 0;//omp_get_thread_num();
        auto & aln = alignments[thread_id];
        // TODO: merge arrays file_range_pos and _end are thread-safe for cache locality. should be thread-safe

        std::vector<double> scores;
        std::vector<scored_region> region;

        while (maf_rd.get_next_alignment(aln, job_id))
        {
//            if (aln.chrom == "chrM")
//                aln.chrom = "chrMT";

            // on first iteration, compute bls scores (used by all 6 frames then!)
            bls_per_bp[thread_id].clear(); // TODO: should all be thread_id not job_id
            compute_bls_score(data_fixed_mle[thread_id].phylo_tree, aln, bls_per_bp[thread_id]);

            const uint64_t orig_start_pos = aln.start_pos;

            for (char strand = '+'; strand <= '-'; strand += 2)
            {
                for (unsigned frame = 1; frame <= 3; ++frame)
                {
                    const uint8_t file_index = frame - 1 + (strand == '+' ? 0 : 3);
                    aln.update_seqs(orig_start_pos, strand, frame);

                    std::tuple<double, double, double> results_fixed;

                    results_fixed = run(data_fixed_mle[thread_id], aln, algorithm_t::FIXED, lpr_per_codon[thread_id], bls_per_bp[thread_id]);
                    data_fixed_mle[thread_id].clear();

                    if (strand == '-')
                    {
                        std::reverse(lpr_per_codon[thread_id].begin(), lpr_per_codon[thread_id].end());

                        const size_t excess_basepairs = aln.length() % 3;
                        aln.start_pos += excess_basepairs;
                    }

                    if (frame == 3 && strand == '+')
                    {
                        // since we iterate over a codon array, there must be 3 bp for each codon
                        // the last 0-2 remaining basepairs in the bls array do not have a codon entry
                        assert(lpr_per_codon[thread_id].size() * 3 <= bls_per_bp[thread_id].size());

                        fprintf(file_power[job_id], "fixedStep chrom=%s start=%ld step=3 span=3\n", aln.chrom.c_str(), aln.start_pos);
                        for (uint32_t pos = frame - 1; pos < lpr_per_codon[thread_id].size() * 3; pos += 3)
                        {
                            const float bls_codon_avg = (bls_per_bp[thread_id][pos]
                                                      +  bls_per_bp[thread_id][pos + 1]
                                                      +  bls_per_bp[thread_id][pos + 2]) / 3.0;

                            my_fprintf(file_power[job_id], "%.4f", bls_codon_avg);
                        }
                    }

                    size_t bls_pos = 0;
                    // compute offset
                    if (strand == '+')
                        bls_pos += (frame - 1);
                    else
                        bls_pos += aln.length() % 3;

                    int64_t prevPos = -4;
                    int64_t startBlockPos = aln.start_pos;
                    for (uint32_t xx = 0; xx < lpr_per_codon[thread_id].size(); ++xx, bls_pos += 3)
                    {
                        const float bls_codon_sum = bls_per_bp[thread_id][bls_pos]
                                                  + bls_per_bp[thread_id][bls_pos + 1]
                                                  + bls_per_bp[thread_id][bls_pos + 2];
                        if (bls_codon_sum < 0.1f * 3)
                        {
                            if (scores.empty())
                                startBlockPos = aln.start_pos + ((xx + 1) * 3);
                            continue;
                        }

                        int64_t newPos = aln.start_pos + (xx * 3);

                        if (prevPos + 3 != newPos)
                        {
                            fprintf(file_score_raw[job_id][file_index], "fixedStep chrom=%s start=%ld step=3 span=3\n", aln.chrom.c_str(), newPos);
                            if (!scores.empty())
                            {
                                fprintf(file_score[job_id][file_index], "fixedStep chrom=%s start=%ld step=3 span=3\n", aln.chrom.c_str(), startBlockPos);

                                process_scores(hmm, scores, startBlockPos, region, SCORE_CODON);

                                for(size_t i = 0; i < region.size(); i++)
                                {
                                    my_fprintf(file_score[job_id][file_index], "%.3f", region[i].log_odds_prob);
                                }

                                scores.clear();
                                region.clear();
                                startBlockPos = aln.start_pos + (xx * 3);
                            }
                        }

                        prevPos = newPos;
                        my_fprintf(file_score_raw[job_id][file_index], "%.3f", lpr_per_codon[thread_id][xx]);
                        scores.push_back(lpr_per_codon[thread_id][xx]);
                    }

                    if (!scores.empty())
                    {
                        fprintf(file_score[job_id][file_index], "fixedStep chrom=%s start=%ld step=3 span=3\n", aln.chrom.c_str(), startBlockPos);
                        process_scores(hmm, scores, startBlockPos, region, SCORE_CODON);

                        for(size_t i = 0; i < region.size(); i++)
                        {
                            my_fprintf(file_score[job_id][file_index], "%.3f", region[i].log_odds_prob);
                        }

                        scores.clear();
                        region.clear();
                    }
                }

                // compute reverse complement for neg. strand
                for (auto & seq : aln.seqs)
                {
                    std::reverse(seq.begin(), seq.end());
                    for (uint64_t j = 0; j < seq.size(); ++j)
                    {
                        seq[j] = complement(seq[j]);
                    }
                }
            }

            for (auto & seq : aln.seqs)
                seq = "";
        }
    }

    for (size_t job_id = 0; job_id < jobs; ++job_id)
    {
        fclose(file_power[job_id]);

        for (uint8_t i = 0; i < 6; ++i)
        {
            fclose(file_score_raw[job_id][i]);
            fclose(file_score[job_id][i]);
        }
        free(file_score_raw[job_id]);
        file_score_raw[job_id] = NULL;
        free(file_score[job_id]);
        file_score[job_id] = NULL;
    }

    free(file_power);
    file_power = NULL;

    free(file_score_raw);
    file_score_raw = NULL;

    free(file_score);
    file_score = NULL;

    // merge files
//    std::ofstream of_a("a.txt", std::ios_base::binary | std::ios_base::app);
//    std::ifstream if_b("b.txt", std::ios_base::binary);
//
//    of_a.seekp(0, std::ios_base::end);
//    of_a << if_b.rdbuf();

    return 0;
}
