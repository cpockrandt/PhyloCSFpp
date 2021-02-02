#pragma once

#include <cfloat>
#include <stdio.h>
#include <random>
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
    std::map<std::string, std::list<range> > exons_by_chr_map;  // TODO:  can't we make this an unordered_map?
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
