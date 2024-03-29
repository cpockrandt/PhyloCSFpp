#pragma once

// TODO: represent subset as bitvector, store bit vector for each node and then do bitwise and check whether it is non-zero
// also: can we store these nodes in an array instead of a tree structure? still need to consider "arrived_at_lca"
double newick_sum_branch_lengths(const newick_node* node, const std::unordered_set<std::string> & subset, bool arrived_at_lca = false, int16_t overlap_parent = -1)
{
    if (node->left == NULL) // leaf
    {
        assert(subset.find(node->label) != subset.end());
        return node->branch_length;
    }
    else
    {
        if (overlap_parent == -1)
        {
            overlap_parent = newick_overlap_size(node, subset);
        }

        const uint16_t overlap_left_child = newick_overlap_size(node->left, subset);
        const uint16_t overlap_right_child = overlap_parent - overlap_left_child;

        double bl = 0.0;

        // start counting branch lengths only when we arrived at a node where some selected species are in the left and some in the right child.
        // this is called the least common ancestor (lca). For that node and all parent nodes we don't want to add up the branch lengths.
        if (arrived_at_lca)
        {
            bl = node->branch_length;
        }

        if (overlap_left_child > 0 && overlap_right_child > 0)
            arrived_at_lca = true;

        if (overlap_left_child > 0)
            bl += newick_sum_branch_lengths(node->left, subset, arrived_at_lca, overlap_left_child);
        if (overlap_right_child > 0)
            bl += newick_sum_branch_lengths(node->right, subset, arrived_at_lca, overlap_right_child);

        return bl;
    }
}

template <bool score_per_base>
double compute_bls_score(const newick_node* node, const alignment_t & alignment, const Model & model, std::vector<double> & scores_per_base)
{
    const uint64_t lo = 0;
    const uint64_t hi = alignment.seqs[0].size();
    double bl_total = 0.0;

    std::unordered_set<std::string> all_species;

    for (uint16_t i = 0; i < (model.phylo_array.size() + 1)/2; ++i)
        all_species.insert(model.phylo_array[i].label);

    const double all_species_branch_length = newick_sum_branch_lengths(node, all_species); // TODO: not necessary to compare species names when counting all branch lengths

    for (uint64_t i = lo; i < hi; ++i)
    {
        // determine subspecies set that has A,C,G or T at position i
        std::unordered_set<std::string> subset;
        for (uint16_t species_id = 0; species_id < alignment.seqs.size(); ++species_id)
        {
            if (alignment.seqs[species_id].size() > 0 && get_dna_id(alignment.seqs[species_id][i]) <= 3)
                subset.insert(model.phylo_array[species_id].label/*alignment.ids[species_id]*/);
        }

        if (subset.size() >= 2) // NOTE: if only one sequence has a DNA4 base, Ocaml produces an (empty?) subtree, and we seem to produce a tree with some branch length! that's why we have this if statement here!
        {
            const double bl = newick_sum_branch_lengths(node, subset);
            bl_total += bl;

            if (score_per_base)
                scores_per_base.push_back(bl / all_species_branch_length);
        }
        else if (score_per_base)
        {
            scores_per_base.push_back(0.0);
        }
    }

    const double divisor = all_species_branch_length * (hi - lo);

    return bl_total / divisor;
}