#pragma once

#include <cfloat>
#include <stdio.h>
#include <random>       // std::default_random_engine
#include <cstdint>
#include <string.h>

#include <string>
#include <map>
#include <list>
#include <algorithm>

struct hmm_parameter{
    double coding_prior;
    double coding_length_in_codons;
    double non_coding_weights[3];
    double non_coding_lengths_in_codons[3];
    hmm_parameter(double coding_prior, double coding_length_in_codons,
                  double non_coding_weights[3], double non_coding_lengths_in_codons[3]) :
            coding_prior(coding_prior), coding_length_in_codons(coding_length_in_codons)
    {
        this->non_coding_weights[0] = non_coding_weights[0];
        this->non_coding_weights[1] = non_coding_weights[1];
        this->non_coding_weights[2] = non_coding_weights[2];
        this->non_coding_lengths_in_codons[0] = non_coding_lengths_in_codons[0];
        this->non_coding_lengths_in_codons[1] = non_coding_lengths_in_codons[1];
        this->non_coding_lengths_in_codons[2] = non_coding_lengths_in_codons[2];
    }
};

enum simplex_type{
    Reflection,
    Expansion,
    Contraction,
    Reduction
};

double fjj(const std::vector<uint32_t> &points, const double * class_probs, const double simplex_i){
    double f_of_params = 0.0;
    for(size_t p = 0; p < points.size(); p++){
        double tau = pow(10.0, simplex_i);
        double log_likelihood = (-static_cast<double>(points[p]) / tau - log(tau));
        f_of_params -= class_probs[p] * log_likelihood;
    }
    return f_of_params;
}

simplex_type simplex_one_step(const std::vector<uint32_t> &points, const double * class_probs, std::vector<std::pair<double, double> > & simplex){
    size_t n = simplex.size() - 1;
    double centroid = 0.0;
    for(size_t i = 0; i < n; i ++){
        centroid += simplex[i].first;
    }
    centroid /= n;
    double reflection = centroid + (centroid - simplex[n].first);
    simplex_type stepType;
    double freflection = fjj(points, class_probs, reflection);
    if(simplex[0].second <= freflection && freflection < simplex[n - 1].second){
        stepType = Reflection;
        simplex[n] = std::make_pair(reflection, freflection);
    }else if (freflection < simplex[0].second) {
        double expansion = centroid + 2 * (centroid - simplex[n].first);
        double fexpansion = fjj(points, class_probs,expansion);
        if(fexpansion < freflection){
            stepType = Expansion;
            simplex[n] = std::make_pair(expansion, fexpansion);
        }else {
            stepType = Reflection;
            simplex[n] = std::make_pair(reflection, freflection);
        }
    } else {
        double contraction = centroid - .5 * (centroid - simplex[n].first);
        double fcontraction = fjj(points, class_probs, contraction);
        if(fcontraction < simplex[n].second){
            stepType = Contraction;
            simplex[n] = std::make_pair(contraction, fcontraction);
        }else{
            stepType = Reduction;
            for(size_t i = 1; i < n + 1; i++){
                double newx = simplex[0].first + 0.5 * (simplex[i].first- simplex[0].first);
                simplex[i].first  = newx;
                simplex[i].second = fjj(points, class_probs, newx);
            }
        }
    }
    std::sort(simplex.begin(), simplex.end(), [](const std::pair<double,double> &left, const std::pair<double,double> &right) {
        return left.second < right.second;
    });
    return stepType;
}

std::pair<double, double> nelder_mead(const std::vector<uint32_t> &points, const double * class_probs,
                                      std::vector<double> & initial_simplex, const double xtol_vec, const uint32_t max_steps){
    std::vector<std::pair<double, double> > simplex;
    simplex.resize(initial_simplex.size());
    for(size_t i = 0; i < initial_simplex.size(); i++){
        double initial_param = fjj(points, class_probs, initial_simplex[i]);
        simplex[i].first = initial_simplex[i];
        simplex[i].second = initial_param;
    }
    std::sort(simplex.begin(), simplex.end(),  [](const std::pair<double,double> &left, const std::pair<double,double> &right) {
        return left.second < right.second;
    });
    bool prev_step_expansion_or_reduction = true;
    for(size_t i = 0; i <= max_steps; i++){
        if (prev_step_expansion_or_reduction == false){
            double min_x = DBL_MAX;
            double max_x = DBL_MIN;
            for(size_t s = 0; s < simplex.size(); s++){
                min_x = std::min(min_x,simplex[s].first);
                max_x = std::max(max_x,simplex[s].first);
            }
            if(max_x - min_x < xtol_vec){
                return simplex[0];
            }
        }
        simplex_type step_type = simplex_one_step(points, class_probs, simplex);
        prev_step_expansion_or_reduction = (step_type == Expansion || step_type == Reduction);
    }

    printf("ERROR: nelder_mead did not converged in %d steps.\n", max_steps);
    exit(-1);
}


double minimize(const std::vector<uint32_t> &points, double * class_probs, double guess, double xscales, double relxtol){
    const size_t max_num_steps = 30;
    std::vector<double> initial_simplex;
    initial_simplex.push_back(guess);
    double vertex = guess + xscales;
    initial_simplex.push_back(vertex);
    return nelder_mead(points, class_probs, initial_simplex, relxtol * xscales, max_num_steps).first;
}

struct mixture{
    double param_seqs[3];
    double class_prior[3];
    mixture(double param_seqs[3], double class_prior[3])
    {
        this->param_seqs[0] = param_seqs[0];
        this->param_seqs[1] = param_seqs[1];
        this->param_seqs[2] = param_seqs[2];
        this->class_prior[0] = class_prior[0];
        this->class_prior[1] = class_prior[1];
        this->class_prior[2] = class_prior[2];
    }
};


mixture infer_mixture(const std::vector<uint32_t> & points, const double param_guess[3],
                      const double guess_prior[3], const uint32_t num_steps,
                      const double relxtol){
    size_t num_points = points.size();
    size_t num_dist = 3;
    double param_seqs[3];
    memcpy(param_seqs, param_guess, sizeof(double) * 3);
    double class_priors[3];
    memcpy(class_priors, guess_prior,sizeof(double) * 3);
    double ** class_probes = new double*[3];
    for(size_t i = 0; i < 3; i++){
        class_probes[i]=new double[num_points];
        memset(class_probes[i], 0, sizeof(double ) * num_points);
    }
    double likelihoods[3];
    for(size_t it = 0; it < num_steps; it++){
        for(size_t i = 0; i < num_points; i++) {
            for(size_t j = 0; j < num_dist; j++){
                double tau = pow(10.0, param_seqs[j]);
                double log_density_and_grad = (-static_cast<double>(points[i]) / tau - log(tau));
                likelihoods[j] = class_priors[j] * exp(log_density_and_grad);
            }
            double total = likelihoods[0] + likelihoods[1] + likelihoods[2];
            for(size_t j = 0; j < num_dist; j++){
                class_probes[j][i]  = (total != 0.0) ? likelihoods[j] / total : 1.0 / num_dist;
            }
        }
        for (size_t j = 0; j < num_dist; j++) {
            double sum = 0.0;
            for(size_t i = 0; i < num_points; i++) {
                sum += class_probes[j][i];
            }
            class_priors[j] = sum / static_cast<double>(num_points);
        }
        for(size_t j = 0; j < num_dist; j++) {
            if(param_seqs[j] == 0){
                continue;
            }
            double xcales = 0.1;
            param_seqs[j] = minimize(points, class_probes[j], param_seqs[j], xcales, relxtol);
        }
    }
    delete [] class_probes[0];
    delete [] class_probes[1];
    delete [] class_probes[2];
    delete [] class_probes;
    return mixture(param_seqs, class_priors);
}

struct gap_mixture{
    double gap_mean_lengths_codon[3];
    double class_prior[3];
    gap_mixture(double gap_mean_lengths_codon[3], double class_prior[3]) {
        this->gap_mean_lengths_codon[0] = gap_mean_lengths_codon[0];
        this->gap_mean_lengths_codon[1] = gap_mean_lengths_codon[1];
        this->gap_mean_lengths_codon[2] = gap_mean_lengths_codon[2];
        this->class_prior[0] = class_prior[0];
        this->class_prior[1] = class_prior[1];
        this->class_prior[2] = class_prior[2];
    }
};

gap_mixture estimate_gap_mixture_model(std::vector<uint32_t> &gaps_nt,
                                       const uint32_t num_steps,
                                       const double relxtol){

    const uint32_t MAX_NUM_GAPS = 20000;
    if(gaps_nt.size() > MAX_NUM_GAPS){
        std::default_random_engine rng(0);
        std::shuffle(gaps_nt.begin(), gaps_nt.end(), rng);
        gaps_nt.resize(MAX_NUM_GAPS);
    }
//    double tmp[1000] ={765,22056,105433,2503,600,66755,864,23876,285045,12318,264,2651,5418,427547,9582,2500,78743,9844,3649,1305,5481,19008,10368,1563,3159,82050,4740,188670,75,5829,15441,6057,3,25092,20814,2193,167139,8526,1698,6231,330206,453742,19097,7191,602,666172,22080,9348,10826,1590197,1290,1335,3045,8041,32556,3198,2812,3189,819,2025,76514,8184,56172,1921,6964,12792,2535,569388,6594,1578,774,2913,4818,4717,2022,7844,1444785,111550,1548,3651,4184,9288,21504,13705,141186,3848,1878,498,252,26872,17615,471,31927,1833,15613,17956,17196,10747,442,196566,3,49401,519,1413,930,21051,219685,2544,6691,36189,1443,1230,471,143849,13176,231730,2436,27396,42237,18370,5350,249,22498,2564,144866,9746,2543,130776,6627,73152,8427,8145,10092,72,3855,87,15543,7509,126,531,768854,3348,5622,348346,323232,10258,866,794359,1297671,2332,7857,15063,342,84,8382,18104,24047,4506,153,696,5972,13450,20834,2850,1044,1205,373722,1155,2310,479368,3948,99403,1503,4245,14526,270,57777,122154,11697,6925,702,32968,30729,869589,208203,21679,62614,4984,11297,189378,3954,1731,7118,1860,5008,129,159,4632,600,645099,108,24841,615,13914,93231,55100,7302,78,1019,4287,17055,7480,95769,401112,4035,2301,560,7196,5615,816,387,831,147,8711,3970,5022,1668,3810,6360,11706,669,14645,234,29541,2358,7001,2350,369,22100,81560,16370,3683,3897,5706,89202,7377,5578,319236,363,71600,942,9939,15397,5064,109491,5936,1266,1116,13316,385,22279,238907,3813,168,3599,15297,4008,14193,172737,3375,1138,37167,16355,4158,659706,4470,522,1194,306,285,6733,51460,8636,34497,8394,909933,432,9888,4460,10966,1740,337158,723,1364210,10206,891,81,945144,12486,10984,99,687,1787,7557,429035,25041,192,1395,360,15935,3402,846,1536,10709,176159,221809,4065,129030,39901,10894,285,507,6668,92410,267,4908,2607,48097,21293,12585,101591,386480,1404,20606,5202,2307229,450620,27382,7375,324803,1362,90753,134170,1507,759,7935,433648,1884,3153,11613,822843,410127,834,5717,5508,39954,42993,3211,33740,25129,68757,9511,396504,233447,1905,48953,31125,210,4544,4401,77043,7411,2673,1311,6905,67517,1721,2556,11368,4316,1458,5105,6286,4828,80718,35071,286,1665,6131,38748,72519,618,74678,9350,33098,5438,65229,21942,4216,126,108437,705,153,282858,132410,5892,24718,67978,26794,8397,7644,4458,1077,19185,2268,5664,1797,5108,174823,7320,674142,34564,367,1542,5539,2505,7518,25995,6813,13620,42480,3471,29140,2892,9288,130614,4626,3553,19298,7386,1824,6703,2193,255,28934,10253,2028,729,26090,10328,2475,219,99,6133,2571,5956,185020,708,12983,23994,1020,2992,17871,4449,3501,165378,9785,5719,832,3353,12911,648,14617,1848,507,92798,453184,5631,387493,11771,4797,126226,84,2407,7805,80634,1352,20854,102,168,429,20176,1689,20898,1434,10578,126,3080,105,226,13702,8238,88851,7464,146040,5278,321,906,23590,638,13313,52458,5397,1137874,4014,169029,7888,1776,3271,609,2061,22380,21865,448,2289,4780,6727,405916,15501,14526,225297,273,2433,4353,3303,300,11136,34748,1476,7419,10317,810,3552,55442,753,55274,7496,53928,13747,3510,27204,7676,4875,164510,63657,3965,30705,2766,6130,3996,8175,3246,36483,229255,42480,514,5575,18073,3243,1923,46255,3777,7375,52473,6180,16453,358389,26536,3097,10034,10945,2802,5391,4743,52032,9786,1135,602259,19459,22128,7323,3664,12807,3120,5993,2722,7244,777,1908,3909,2203,2001,995557,2121,1185,1563,864,420,17850,3396,20610,7478,4518,80059,180,3184,14283,4332,5820,48710,252,16064,2579,1317,1032,46055,17022,7176,1372,25025,67513,6938,8469,44123,4101,15435,38438,2538,8664,10878,258,71309,7405,14856,7434,3527,5695,3989,1545,2094,14633,5685,2052,2673,8094,16478,174,8182,13962,67398,4106,126883,149702,996,17991,5319,6009,93,6651,3402,339,84,22420,11226,2022,14219,9659,1569,7305,1351,5917,1176,4527,1738,14571,3564,18679,2217,2689,933,5811,16754,198,777,8308,12331,65858,231,348,4887,450,206,2331,582,257013,2271,10483,122018,540,3828,1548,195430,174,1514,363,2055,807,132921,1085,11097,39071,9364,447,35232,7589,2020,10952,5760,1405,540,768,17162,4848,12476,165,27163,247,2961,2045,537157,17743,11946,168124,5889,359,2022,1594,23379,20462,4735,2988,1813537,752,2560,1830,18606,9827,41676,3009,74166,97035,74943,1173,1981,11861,4786,4161,3072,6276,566658,41958,40030,3354,2007,38339,13564,900,114,377,681,261100,2724,8286,474,1176913,217175,2691,3069,20414,2238,6222,6130,12114,14800,2838,24654,2285,1683,10221,301456,7859,229486,36432,1555,323,16342,75,125119,3043,13914,12969,2919,9207,50766,112876,1833,2020,704,20232,49143,6691,276,6852,771,1569,330,11751,2889,903,31919,70369,7026,222,10344,11649,291,2139,1478,22702,5502,1500,615,14130,4437,777,20521,4257,4425,7260,14880,212289,5709,90,144,8468,9327,438,63120,6249,393,7887,5107,3174,2217,24061,37320,3039,76230,101790,1446,7947,891,48003,28,172693,17911,1116,1050,1056506,252,2037,1287,3025,1653,990,1830,1042,3699,14775,399,6279,6567,2406,210,153,3153,825,26603,379620,4742,43866,6481,2757,15351,155130,8744,471,10335,1512,1887,26293,10503,14108,6639,6424,246,3636,49692,984,76523,3051,4116,249402,1151649,6124,9822,312,6052,17328,20088,1200,1029,37273,42570,4581,11185,245273,1327,10064,2070,2857,121549,210,33704,4899,3831,4190,780,87,15920,8299,175335,79623,534,4680,3204,4587,762,12606,427,3118,831,279285,779,15606,3673,24993,632,1923,5387,1897,141,7337,2822,296828,117477,5874,255,3212,738,171,15947,1769,882,298974,1323,3722,81000,2235,1074,867,423,1638,27784,12593,280202,739,9455,288,26601};
//    gaps_nt.clear();
//    for(size_t i =0; i < MAX_NUM_GAPS; i++){
//        gaps_nt.push_back(tmp[i]);
//    }

    uint32_t guess_lengths[3] = {3000, 80000, 100};
    double guess_priors[3] = {30, 10, 1};
    double guess_priors_sum = guess_priors[0] + guess_priors[1] + guess_priors[2];
    double param_guesses[3];
    for(size_t i = 0; i < 3; i++){
        guess_priors[i] = guess_priors[i] / guess_priors_sum;
        param_guesses[i] = log10(guess_lengths[i]);
    }

    mixture mix = infer_mixture(gaps_nt, param_guesses,
                                guess_priors,
                                num_steps, relxtol);
    double gap_mean_lengths_codon[3];
    for(size_t i = 0; i < 3; i++){
        gap_mean_lengths_codon[i] = pow(10, mix.param_seqs[i]) / 3;
    }
    return gap_mixture(gap_mean_lengths_codon, mix.class_prior);
}

hmm_parameter estimate_hmm_params_for_genome(const char * path_exon_list, const uint32_t genome_length){
    struct range{
        uint32_t start;
        uint32_t end;
        range( uint32_t start, uint32_t end)
                : start(start), end(end)  {}

        static bool compare(const range &first, const range &second){
            if (first.start < second.start)
                return true;
            if (second.start < first.start)
                return false;
            if (first.end < second.end)
                return true;
            if (second.end < first.end)
                return false;
            return false;
        }
    };

    FILE* exon_file = fopen(path_exon_list, "r");
    if(exon_file == NULL){
        printf("ERROR: could not open %s.\n", path_exon_list);
        exit(-1);
    }

    char line[1024];
    std::map<std::string, std::list<range> > exons_by_chr_map;
    while (fgets(line, sizeof(line), exon_file)) {
        const char * delim = " \t";
        char *token = strtok(line, delim);
        std::string chrom(token);
        token = strtok(NULL, delim);
        std::string strand(token);
        token = strtok(NULL, delim);
        std::string frame(token);
        token = strtok(NULL, delim);
        uint32_t start=atoi(token);
        token = strtok(NULL, delim);
        uint32_t end=atoi(token);
        std::string key = chrom + ":" + strand + ":" + frame;
        exons_by_chr_map[key].emplace_back(start, end);
    }
    fclose(exon_file);
    for(auto it = exons_by_chr_map.begin(); it != exons_by_chr_map.end(); it++){
        it->second.sort(range::compare);
    }
    uint64_t num_exons = 0;
    size_t total_coding_length_nt = 0;
    std::vector<uint32_t> gaps_nt;
    for(auto map_it = exons_by_chr_map.begin(); map_it != exons_by_chr_map.end(); map_it++){
        std::list<range> & pair_list = map_it->second;
        size_t i = 0;
        auto tail = pair_list.begin();
        auto head = pair_list.begin();
        head++;
        while( i < pair_list.size() - 1){
            uint32_t start1 = tail->start;
            uint32_t end1 = tail->end;
            uint32_t start2 = head->start;
            uint32_t end2 = head->end;
            if(start2 <= end1){
                if(end1 - start1 >= end2 - start2){
                    head=pair_list.erase(head);
                } else{
                    tail=pair_list.erase(tail);
                    head=tail;
                    head++;
                }
            }else{
                // non overlapping gap
                head++;
                tail++;
                i++;
            }
        }
        for(auto it = pair_list.begin(); it != pair_list.end(); it++){
            uint32_t end1 = it->end;
            it++;
            if(it == pair_list.end())
                break;
            uint32_t start2 = it->start;

            if(start2 > end1 + 1){
                gaps_nt.push_back(start2 - end1 - 1);
            }
        }
        num_exons += pair_list.size();
        for(auto it = pair_list.begin(); it != pair_list.end(); it++){
            total_coding_length_nt += it->end - it->start + 1;
        }
    }
    gap_mixture gap_mix = estimate_gap_mixture_model(gaps_nt, 20, 0.001);

    double codingPrior = static_cast<double>(total_coding_length_nt) / static_cast<double>(genome_length) / 6.0; // Prior for being coding in a particular frame
    double codingLengthInCodons = static_cast<double>(total_coding_length_nt) / static_cast<double>(num_exons) / 3.0; // Mean length in codons
    return hmm_parameter(codingPrior, codingLengthInCodons, gap_mix.class_prior, gap_mix.gap_mean_lengths_codon);
}
