#pragma once

#include <string>
#include <map>
#include <list>
#include <stdio.h>
#include <algorithm>
#include <cstdint>
#include <string.h>

struct hmm{
    uint32_t num_states;
    double init_probs[4];
    double trans_probs[4][4];

    hmm() {}

    hmm(uint32_t num_states, double init_probs[4],
        double trans_probs[4][4]){
        this->num_states = num_states;
        memcpy(this->trans_probs, trans_probs, 4 * 4 * sizeof(double));
        memcpy(this->init_probs, init_probs, 4 *  sizeof(double));
    }

    static double emit_prob(uint32_t state, double score){
        return (state == 0 ) ? pow(10,(score / 10)) : 1;
    }

    static std::vector<uint32_t> get_best_path_by_viterbi(const hmm & hmm, const std::vector<double> &observation_seq) {
        double init_obs = observation_seq[0];
        std::vector<std::vector<uint32_t>> best_states;

        double prev_best_probs[4];
        for(size_t state = 0; state < hmm.num_states; state++ ){
            prev_best_probs[state] = hmm.init_probs[state] * emit_prob(state, init_obs);
        }

        for(size_t pos = 1; pos < observation_seq.size(); pos++) {
            std::vector<uint32_t> cur_best_states;
            double cur_obs = observation_seq[pos];
            double max_prob = 0.0;
            double cur_best_probs[4];
            for (uint32_t cur_state = 0; cur_state < hmm.num_states; cur_state++) {
                double cur_max_Prob = 0.0;
                uint32_t max_prev_state = 0;
                for (uint32_t prev_state = 0; prev_state < hmm.num_states; prev_state++) {
                    double currProb = prev_best_probs[prev_state] * hmm.trans_probs[prev_state][cur_state];
                    if(currProb > cur_max_Prob){
                        max_prev_state = prev_state;
                        cur_max_Prob = currProb;
                    }
                }
                cur_max_Prob *= emit_prob(cur_state, cur_obs);
                cur_best_states.push_back(max_prev_state);
                cur_best_probs[cur_state]= cur_max_Prob;
                max_prob = std::max(max_prob, cur_max_Prob);
            }
            best_states.push_back(cur_best_states);

            // Scale to avoid floating underflow
            for (size_t state = 0; state < hmm.num_states; state++) {
                cur_best_probs[state] = cur_best_probs[state] / max_prob;
                prev_best_probs[state] = cur_best_probs[state];
            }
        }
        std::vector<uint32_t> result;
        result.push_back(indexmax(prev_best_probs));
        std::reverse(best_states.begin(), best_states.end());
        for(size_t i = 0; i < best_states.size(); i++){
            uint32_t prevState = result[result.size()-1];
            result.push_back(best_states[i][prevState]);
        }
        std::reverse(result.begin(), result.end());
        return result;
    }

    static uint32_t indexmax(double probs[4]) {
        // Return the index of the maximum element in seq
        uint32_t best_ind = UINT32_MAX;
        double max_elm = 0.0;
        for (size_t i = 0; i < 4; i++) {
            if (best_ind == UINT32_MAX || probs[i] > max_elm){
                max_elm = probs[i];
                best_ind = i;
            }
        }
        return best_ind;
    }

    static double ** state_posterior_probabilities(const hmm &hmm, const std::vector<double> &observation_seq) {
        size_t num_obs = observation_seq.size();

        double ** forward = new double*[num_obs];
        forward[0] = new double[hmm.num_states * num_obs];
        for(size_t i = 1; i < num_obs; i++){
            forward[i] = forward[0] + i * hmm.num_states;
        }
        for(size_t state = 0; state < hmm.num_states; state++){
            forward[0][state] = hmm.init_probs[state] * emit_prob(state, observation_seq[0]);
        }
        for(size_t pos = 1; pos < num_obs; pos++){
            double obs = observation_seq[pos];
            double maxf = 0.0;
            for(size_t state = 0; state < hmm.num_states; state++){
                double sum=0.0;
                for(size_t prevState = 0; prevState < hmm.num_states; prevState++){
                    sum += forward[pos - 1][prevState] * hmm.trans_probs[prevState][state];
                }
                forward[pos][state] = sum * emit_prob(state, obs);
                maxf = std::max(maxf, forward[pos][state]);
            }
            for(size_t state = 0; state < hmm.num_states; state++) {
                forward[pos][state] /= maxf;
            }
        }

        double ** backward = new double*[num_obs];
        backward[0] = new double[hmm.num_states * num_obs];
        for(size_t i = 1; i < num_obs; i++){
            backward[i] = backward[0] + i * hmm.num_states;
        }
        for(size_t state = 0; state < hmm.num_states; state++) {
            backward[num_obs-1][state] = 1.0;
        }

        for(int32_t pos = static_cast<int32_t>(num_obs)-2; pos > -1; pos--){
            double obs = observation_seq[pos + 1];
            double maxb = 0.0;
            for(size_t state = 0; state < hmm.num_states; state++) {
                double sum = 0.0;
                for(size_t nextState = 0; nextState < hmm.num_states; nextState++) {
                    sum += (hmm.trans_probs[state][nextState] * emit_prob(nextState, obs) *
                            backward[pos + 1][nextState]);
                }
                backward[pos][state] = sum;
                maxb = std::max(maxb, backward[pos][state]);
            }
            for(size_t state = 0; state < hmm.num_states; state++) {
                backward[pos][state] /= maxb;
            }
        }

        double ** result = forward;
        for(size_t i = 1; i < num_obs; i++){
            result[i] = result[0] + i * hmm.num_states;
        }
        for(size_t pos = 0; pos < num_obs; pos++){
            double total = 0.0;
            for(size_t state = 0; state < hmm.num_states; state++) {
                total += (forward[pos][state] * backward[pos][state]);
            }
            for(size_t state = 0; state < hmm.num_states; state++) {
                result[pos][state] = (forward[pos][state] * backward[pos][state]) / total;
            }
        }
        delete [] backward[0];
        delete [] backward;
        return result;
    }
};


hmm get_coding_hmm(const hmm_parameter & param){
    size_t numNonCoding = 3;
//    assert abs(1 - sum(nonCodingWeights)) < 1e-9;
    double unnormalized_noncoding_priors[3];
    double c_to_nc_probs[3];
    double nc_to_c_probs[3];

    for(size_t i = 0; i < 3; i++){
        unnormalized_noncoding_priors[i] = param.non_coding_weights[i] * param.non_coding_lengths_in_codons[i];
        c_to_nc_probs[i]= param.non_coding_weights[i] / param.coding_length_in_codons;
        nc_to_c_probs[i] = 1.0 / param.non_coding_lengths_in_codons[i];

    }
    double init_probs[4];

    init_probs[0] = param.coding_prior;
    double unnormalized_noncoding_priors_sum = unnormalized_noncoding_priors[0] + unnormalized_noncoding_priors[1] + unnormalized_noncoding_priors[2];
    for(size_t i = 0; i < 3; i++){
        init_probs[i+1] = (1 - param.coding_prior) * unnormalized_noncoding_priors[i] / unnormalized_noncoding_priors_sum;
    }
    double trans_probs[4][4];
    trans_probs[0][0] = (1.0 - (c_to_nc_probs[0] + c_to_nc_probs[1] + c_to_nc_probs[2]));
    for (size_t j = 0; j < 3; j++) {
        trans_probs[0][j+1] = c_to_nc_probs[j];
    }
    for(size_t i = 1; i < 4; i++) {
        trans_probs[i][0] = nc_to_c_probs[i - 1];
        for (size_t j = 1; j < 4; j++) {
            if(i == j){
                trans_probs[i][j] = 1.0 - nc_to_c_probs[i - 1];
            }else{
                trans_probs[i][j] = 0.0;
            }
        }
    }

    return hmm(1 + numNonCoding, init_probs, trans_probs);
}

struct scored_region{
    const uint32_t region_start;
    const uint32_t region_end;
    const double log_odds_prob;
    scored_region(uint32_t region_start, uint32_t region_end,
                  double log_odds_prob) :
                      region_start(region_start),
                      region_end(region_end),
                      log_odds_prob(log_odds_prob){}
};

struct scored_bed_region{
    const uint32_t region_start;
    const uint32_t region_end;
    const double log_odds_prob;
    const uint32_t color;
    scored_bed_region(uint32_t region_start, uint32_t region_end,
                  double log_odds_prob, uint32_t color) :
            region_start(region_start),
            region_end(region_end),
            log_odds_prob(log_odds_prob),
            color(color){}
};

double compute_log_odds(double prob) {
    const double MAX_LOG_ODDS = 15.0;
    if (prob < pow(10, -MAX_LOG_ODDS)) {
        return -MAX_LOG_ODDS;
    } else if (prob > 1 - pow(10, -MAX_LOG_ODDS)) {
        return MAX_LOG_ODDS;
    } else {
        return log10(prob / (1 - prob));
    }
}

uint32_t computing_color_code(double log_odd) {
    if (log_odd < 1) {
        return 210;
    } else if (log_odd > 7) {
        return 0;
    } else {
        return 210 - 30 * int(floor(log_odd));
    }
}
/*
 * Processes a set of contiguous scores and returns coding regions with its respective maximal log-odds score
 */
void process_scores(hmm const & hmm, std::vector<double> &scores,
                    uint32_t blockStartPos, std::vector<scored_region> & result, std::vector<scored_bed_region> & bedresult){
    double ** state_probabilities = hmm::state_posterior_probabilities(hmm, scores);
    std::vector<uint32_t> path = hmm::get_best_path_by_viterbi(hmm, scores);
    for (size_t cur_codon_count = 0; cur_codon_count < scores.size(); cur_codon_count++) {
        uint32_t chunk_start_pos = blockStartPos + 3 * cur_codon_count;
        uint32_t chunk_end_pos = chunk_start_pos + 2;
        double log_odds = compute_log_odds(state_probabilities[cur_codon_count][0]);
        result.emplace_back(chunk_start_pos, chunk_end_pos, log_odds);
    }
    uint32_t starting_position=0;
    uint32_t end_position;
    uint32_t starting_count=0;
    uint32_t end_count;
    for (size_t i=0; i<path.size()-1; i++) {
        double prob=0;
        uint32_t color=0;
        if (i==0 && path[i]==0) {
            starting_position=blockStartPos-1;
            starting_count=0;
            if (path[i+1]!=0) {
                prob = state_probabilities[starting_count][0];
                double log_odd = compute_log_odds(prob);
                computing_color_code(log_odd);
                bedresult.emplace_back(starting_position, starting_position+3, prob, color);
            }
        } else if (path[i+1]==0 && path[i]!=0) {
            if (i!=path.size()-2) {
                starting_position = blockStartPos + 3 * i + 2;
                starting_count = i + 1;
            }
            else {
                end_position=blockStartPos+3*i+5;
                end_count=i+1;
                prob = state_probabilities[end_count][0];
                double log_odd = compute_log_odds(prob);
                computing_color_code(log_odd);
                bedresult.emplace_back(end_position-3, end_position, prob, color);
            }
        } else if (path[i+1]!=0 && path[i]==0) {
            end_position=blockStartPos+3*i+2;
            end_count=i;
            for (uint32_t codon=starting_count; codon<=end_count; codon++) {
                if (prob < state_probabilities[codon][0]) {
                    prob = state_probabilities[codon][0];
                }
            }
            double log_odd = compute_log_odds(prob);
            computing_color_code(log_odd);
            bedresult.emplace_back(starting_position, end_position, prob, color);
        } else if (i==path.size()-2 && path[i+1]==0 && path[i]==0) {
            end_position=blockStartPos+3*i+5;
            end_count=i+1;
            for (uint32_t codon=starting_count; codon<=end_count; codon++) {
                if (prob < state_probabilities[codon][0]) {
                    prob = state_probabilities[codon][0];
                }
            }
            double log_odd = compute_log_odds(prob);
            computing_color_code(log_odd);
            bedresult.emplace_back(starting_position, end_position, prob, color);
        }
    }
    delete [] state_probabilities[0];
    delete [] state_probabilities;
}
