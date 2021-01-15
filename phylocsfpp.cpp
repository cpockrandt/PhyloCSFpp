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

    const std::string filename_power = output_folder + "/PhyloCSFpower-chrM.wig";

    const hmm_parameter hmm_param = estimate_hmm_params_for_genome(hmm_data_path, genome_length);
    const hmm hmm = get_coding_hmm(hmm_param);

    std::vector<Data> data_omega(threads), data_fixed_mle(threads);
    char delim[] = ",";

    char *ptr = strtok(selected_species_str, delim);
    while (ptr != NULL)
    {
        for (size_t thread_id = 0; thread_id < threads; ++thread_id)
        {
            data_omega[thread_id].selected_species.emplace(ptr);
            data_fixed_mle[thread_id].selected_species.emplace(ptr);
        }
        ptr = strtok(NULL, delim);
    }

    for (size_t thread_id = 0; thread_id < threads; ++thread_id)
    {
        // TODO: check whether selected species exist!
        data_omega[thread_id].load_model(model_str, true);
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
                for (uint16_t i = 0; i < data_omega[0].phylo_array.size(); ++i)
                {
                    if (strcmp(data_omega[0].phylo_array[i].label.c_str(), c1) == 0)
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
        alignments.emplace_back(data_omega[0].phylo_array); // TODO: remove id's outside of aln struct because they are read-only
    }

    // open correct values from orig PhyloCSF++
//    // awk 'BEGIN{ OFS="\t"; print "ALN-ID", "FIXED", "FIXED-ANC", "MLE", "MLE-ANC", "OMEGA", "BLS" } ($0 != "" && !($0 ~ /^\/tmp/)) { if ($3 ~ /^\/tmp/) { $3 = "nan" } if ($2 == "Fixed:") { bls[$1] = $4 } anc[$1][$2] = $5; f[$1][$2] = $3 } END { for (id in f) print id, f[id]["Fixed:"], anc[id]["Fixed:"], f[id]["MLE:"], anc[id]["MLE:"], f[id]["Omega:"], bls[id] }' chr22.head.orig.results > chr22.head.orig.results.formatted
//    std::vector<comp_result> correct_results;
//    {
//        FILE *file = fopen(scores_path, "r");
//        if (file != NULL)
//        {
//            char line[BUFSIZ];
//            size_t id;
//            double fixed_score, fixed_anc, mle_score, mle_anc, omega_score, bls;
//            fgets(line, sizeof line, file); // skip header line
//            while (fgets(line, sizeof line, file) != NULL)
//            {
//                sscanf(line, "%lu\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", &id, &fixed_score, &fixed_anc, &mle_score, &mle_anc, &omega_score, &bls);
//                assert(id == correct_results.size());
//                correct_results.emplace_back(fixed_score, fixed_anc, mle_score, mle_anc, omega_score, bls);
//            }
//        }
//        fclose(file);
//    }
////    printf("%ld\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n\n", 0,  results[0].fixed_score, results[0].fixed_anc, results[0].mle_score, results[0].mle_anc, results[0].omega_score, results[0].bls);

    parallel_maf_reader tmp_maf_rd(aln_path, jobs, &fastaid_to_alnid);
    jobs = tmp_maf_rd.get_jobs();

//    printf("ALN-ID\tFIXED\tFIXED-ANC\tMLE\tMLE-ANC\tOMEGA\tBLS\n");

    std::vector<std::vector<comp_result> > computed_results(jobs);

    std::vector<std::vector<double> > lpr_per_codon(jobs), bls_per_codon(jobs);

    FILE *file_power = fopen(filename_power.c_str(), "w");
    if (file_power == NULL)
    {
        printf("Error!");
        exit(1);
    }

    for (char strand = '+'; strand <= '-'; strand += 2)
    {
        for (unsigned frame = 1; frame <= 3; ++frame)
        {
            parallel_maf_reader maf_rd(aln_path, jobs, &fastaid_to_alnid);

            const std::string filename_score_raw = output_folder + "/PhyloCSFRaw" + std::string(1, strand) + std::to_string(frame) + "-chrM.wig";
            const std::string filename_score     = output_folder + "/PhyloCSF" + std::string(1, strand) + std::to_string(frame) + "-chrM.wig";
            FILE *file_score_raw = fopen(filename_score_raw.c_str(), "w");
            FILE *file_score = fopen(filename_score.c_str(), "w");
            if (file_score_raw == NULL || file_score == NULL)
            {
                printf("Error!");
                exit(1);
            }

        //    #pragma omp parallel for num_threads(threads) default(none) shared(jobs, alignments, maf_rd, data_fixed_mle, data_omega, frame, stdout, computed_results, lpr_per_codon, bls_per_codon)
            for (unsigned job_id = 0; job_id < jobs; ++job_id) // TODO: split it in more parts than there are threads
            {
                unsigned thread_id = 0;//omp_get_thread_num();
                auto & aln = alignments[thread_id];
                size_t a_id = 0;
                // TODO: merge arrays file_range_pos and _end are thread-safe for cache locality. should be thread-safe

                std::vector<double> scores;
                std::vector<scored_region> region;

                while (/*a_id < correct_results.size() && */maf_rd.get_next_alignment(aln, job_id, frame, strand))
                {
                    if (aln.chrom == "chrM")
                        aln.chrom = "chrMT";
//                    printf("%20ld, %20ld\n", aln.start_pos, aln.seqs[0].size());
        //            auto & r = correct_results[0];
        //            printf("%d\n", aln.start_pos);

                    std::tuple<double, double, double> results_fixed, results_mle, results_omega;

                    {
                        results_fixed = run(data_fixed_mle[thread_id], aln, algorithm_t::FIXED, lpr_per_codon[job_id], bls_per_codon[job_id]);
                        data_fixed_mle[thread_id].clear();

        //                const double score = fabs(std::get<0>(results_fixed) - r.fixed_score);
        //                const double anc = fabs(std::get<2>(results_fixed) - r.fixed_anc);
        //                const double bls = fabs(std::get<1>(results_fixed) - r.bls);
        //                if (score >= 0.000001 || anc >= 0.000001 || bls >= 0.000001)
        //                    printf("%ld\t%s:\t%f\t%f\t%f\t*****\n", a_id, "Fixed", score, anc, bls);
                    }

        //            {
        //                results_mle = run(data_fixed_mle[thread_id], aln, algorithm_t::MLE);
        //                data_fixed_mle[thread_id].clear();
        ////                const double score = fabs(std::get<0>(results_mle) - r.mle_score);
        ////                const double anc = fabs(std::get<2>(results_mle) - r.mle_anc);
        ////                if (score >= 0.000001 || anc >= 0.000001)
        ////                    printf("%ld\t%s:\t%f\t%f\t--------\t*****\n", a_id, "MLE", score, anc);
        //            }

                    if (strand == '-')
                    {
                        std::reverse(lpr_per_codon[job_id].begin(), lpr_per_codon[job_id].end());
                        std::reverse(bls_per_codon[job_id].begin(), bls_per_codon[job_id].end());

                        aln.start_pos += (aln.seqs[0].size() % 3);
                        for (uint8_t i = 0; i < aln.seqs[0].size() % 3; ++i)
                            bls_per_codon[job_id].erase(bls_per_codon[job_id].begin());
                    }

                    if (frame == 3 && strand == '+')
                    {
                        fprintf(file_power, "fixedStep chrom=%s start=%ld step=3 span=3\n", aln.chrom.c_str(), aln.start_pos);
                        for (uint32_t xx = 0; xx < lpr_per_codon[job_id].size(); ++xx)
                        {
                            size_t pos = xx * 3;

                            float bls_codon_avg = bls_per_codon[job_id][pos];
                            if (pos + 2 < bls_per_codon[job_id].size())
                            {
                                bls_codon_avg += bls_per_codon[job_id][pos + 1];
                                bls_codon_avg += bls_per_codon[job_id][pos + 2];
                                bls_codon_avg /= 3;
                            } else if (pos + 1 < bls_per_codon[job_id].size()) {
                                assert(false);
                                bls_codon_avg += bls_per_codon[job_id][pos + 1];
                                bls_codon_avg /= 2;
                            }
                            else
                            {
                                assert(false);
                            }
\
                            my_fprintf(file_power, "%.4f", bls_codon_avg);
                        }
                    }

                    int64_t prevPos = -4;
                    int64_t startBlockPos = aln.start_pos;
                    for (uint32_t xx = 0; xx < lpr_per_codon[job_id].size(); ++xx)
                    {
                        const float bls_codon_sum = bls_per_codon[job_id][xx * 3] + bls_per_codon[job_id][xx * 3 + 1] + bls_per_codon[job_id][xx * 3 + 2];
                        if (bls_codon_sum < 0.1f * 3)
                        {
                            if (scores.empty())
                                startBlockPos = aln.start_pos + ((xx + 1) * 3);
                            continue;
                        }

                        int64_t newPos = aln.start_pos + (xx * 3);

                        if (prevPos + 3 != newPos)
                        {
                            fprintf(file_score_raw, "fixedStep chrom=%s start=%ld step=3 span=3\n", aln.chrom.c_str(), newPos);
                            if (!scores.empty())
                            {
                                fprintf(file_score, "fixedStep chrom=%s start=%ld step=3 span=3\n", aln.chrom.c_str(), startBlockPos);

                                process_scores(hmm, scores, startBlockPos, region, SCORE_CODON);

                                for(size_t i = 0; i < region.size(); i++)
                                {
                                    my_fprintf(file_score, "%.3f", region[i].log_odds_prob);
                                }

                                scores.clear();
                                region.clear();
                                startBlockPos = aln.start_pos + (xx * 3);
                            }
                        }

                        prevPos = newPos;
                        my_fprintf(file_score_raw, "%.3f", lpr_per_codon[job_id][xx]);
                        scores.push_back(lpr_per_codon[job_id][xx]);
        //                printf("\t\t\t\tCodon %ld\t%f\t%f\n", aln.start_pos + (xx * 3), lpr_per_codon[job_id][xx], bls_per_codon[job_id][xx]);
                    }

                    if (!scores.empty())
                    {
                        fprintf(file_score, "fixedStep chrom=%s start=%ld step=3 span=3\n", aln.chrom.c_str(), startBlockPos);
                        process_scores(hmm, scores, startBlockPos, region, SCORE_CODON);

                        for(size_t i = 0; i < region.size(); i++)
                        {
                            my_fprintf(file_score, "%.3f", region[i].log_odds_prob);
                        }

                        scores.clear();
                        region.clear();
                    }

                    computed_results[job_id].emplace_back(
                            std::get<0>(results_fixed), 0/*std::get<2>(results_fixed)*/,
                            0/*std::get<0>(results_mle)*/, 0/*std::get<2>(results_mle)*/,
                            0,
                            std::get<1>(results_fixed));

        //            {
        //                try
        //                {
        //                    results_omega = run(data_fixed_mle[thread_id], aln, algorithm_t::OMEGA);
        ////                    const double score = fabs(std::get<0>(results_omega) - r.omega_score);
        ////                    if (score >= 0.000001)
        ////                        printf("%ld\t%s:\t%f\t--------\t--------\t*****\n", a_id, "Omega", score);
        //                }
        //                catch (const std::runtime_error& error)
        //                {
        //                    std::get<0>(results_omega) = NAN;
        //                    printf("%ld\t%s:\t%s\t--------\t--------\t*****\n", a_id, "Omega", error.what());
        //                }
        //
        //                data_omega[thread_id].clear();
        //            }

        //            printf("%ld\t%f\t%f\t%f\t%f\t%f\t%f\n",
        //                       a_id,
        //                       std::get<0>(results_fixed), std::get<2>(results_fixed),
        //                       std::get<0>(results_mle), std::get<2>(results_mle),
        //                       std::get<0>(results_omega), std::get<1>(results_fixed)
        //            );
        //            fflush(stdout);

                    for (auto & seq : aln.seqs)
                        seq = "";
                    ++a_id;
                }
            }
            fclose(file_score_raw);
            fclose(file_score);
        }
    }

    fclose(file_power);

      // compute differences
//    size_t j = 0;
//    for (unsigned job_id = 0; job_id < jobs; ++job_id) // TODO: split it in more parts than there are threads
//    {
//        for (size_t i = 0; i < computed_results[job_id].size(); ++i)
//        {
//            const auto & comp = computed_results[job_id][i];
//            printf("%f\t%f\n", comp.fixed_score, comp.bls);
////            const auto & cor = correct_results[j];
////            if (fabs(comp.fixed_score - cor.fixed_score) > 1e-6 ||
////                fabs(comp.fixed_anc   - cor.fixed_anc  ) > 1e-6 ||
////                fabs(comp.bls         - cor.bls        ) > 1e-6
////            )
////                printf("%f\t%f\t\t%f\t%f\t\t%f\t%f\n", comp.fixed_score, cor.fixed_score, comp.fixed_anc, cor.fixed_anc, comp.bls, cor.bls);
//            ++j;
//        }
//    }

    // create wig file for track
//    std::vector<double> scores = {-13.6462, -8.2680, -7.7692, -7.1405, -3.1386, -3.6159, -9.1372, -14.4608, 2.5891};
//    std::vector<scored_region> region;
//    process_scores(hmm, scores, /*block_start_pos*/204316, region, SCORE_CODON);
//    scores.clear();
//
////    std::string phyloCSFregionDir = "/home/chris/dev-uni/PhyloCSF++/chr4_GL000008v2_random.Strand-.Frame1.fixed.out";
////    std::string chrom = "chr4_GL000008v2_random";
////    std::string strand = "-";
////    std::string frame = "1";
////
////    std::vector<scored_region> region = create_track(hmm_param, phyloCSFregionDir);
//    for(size_t i = 0; i < region.size(); i++)
//    {
//        std::cout << "chr4_GL000008v2_random" << "\t" << "-" << "\t" << "1" << "\t"
//                  << region[i].region_start <<"\t" << region[i].region_end << "\t"
//                  << region[i].log_odds_prob << std::endl;
//    }

    return 0;
}

//    char selected_species_str[] = "Dog,Cow,Horse,Human,Mouse,Rat";
//    auto new_results = run("/home/chris/dev-uni/chr4_100vert_alignment.fa", "100vertebrates", "", algorithm_t::FIXED);
//    auto new_results = run("/tmp/test", "100vertebrates", "Dog,Cow,Horse,Human,Mouse,Rat", algorithm_t::MLE);
//    char aln_path[] = "/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/tal-AA-tiny3.fa";
//    char model_str[] = "23flies";
//    "Dog,Cow,Horse,Human,Mouse,Rat";
//    run("/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/tal-AA-tiny3.fa", "23flies", "", algorithm_t::OMEGA);
//    run("/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/tal-AA.fa", "12flies", "", algorithm_t::MLE);
//    run("/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/ALDH2.exon5.fa", "100vertebrates", "", algorithm_t::MLE);
//    printf("\nSOLL:\n"
//           "361.687566\t0.731391\t48.024101\n"
//           "297.623468\t0.731391\t48.257442\n"
//           "-176.807976\t0.058448\t-33.030802\n"   // FIXED
//           "-179.110842\t0.058448\t-32.878297\n"); // MLE
