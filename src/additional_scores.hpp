#pragma once

double newick_sum_branch_lengths(newick_node* node, const std::unordered_set<std::string> & subset, bool arrived_at_lca = false) {
    if (node->left == NULL) // leaf
    {
        assert(subset.find(node->label) != subset.end());
//        std::cout << "leaf: " << node->branch_length << '\n';
        return node->branch_length;
    }
    else
    {
        const uint16_t overlap_left_child = newick_overlap_size(node->left, subset);
        const uint16_t overlap_right_child = newick_overlap_size(node->right, subset); // store nbr of hits in parent node to avoid two passes!
//        std::cout << overlap_left_child << " ... " << overlap_right_child << '\n';

        double bl = 0.0;

        // start counting branch lengths only when we arrived at a node where some selected species are in the left and some in the right child.
        // this is called the least common ancestor (lca). For that node and all parent nodes we don't want to add up the branch lengths.
        if (arrived_at_lca)
        {
            bl = node->branch_length;
//            std::cout << "rek: " << node->branch_length << '\n';
        }

        if (overlap_left_child > 0 && overlap_right_child > 0)
            arrived_at_lca = true;

        if (overlap_left_child > 0)
            bl += newick_sum_branch_lengths(node->left, subset, arrived_at_lca);
        if (overlap_right_child > 0)
            bl += newick_sum_branch_lengths(node->right, subset, arrived_at_lca);

        return bl;
    }
}

//newick_node * newick_go_to_lca(newick_node * n, const std::unordered_set<std::string> & subset)
//{
//    uint16_t overlap_left_child = newick_overlap_size(n->left, subset);
//    uint16_t overlap_right_child = subset.size() - overlap_left_child;
//
//    while (n->left != NULL && (overlap_left_child == 0 || overlap_right_child == 0))
//    {
//        if (overlap_left_child == 0)
//        {
//            n = n->right;
//            if (n->left != NULL)
//            {
//                overlap_left_child = newick_overlap_size(n->left, subset);
//                overlap_right_child -= overlap_left_child;
//            }
//        }
//        else
//        {
//            n = n->left;
//            if (n->left != NULL)
//            {
//                overlap_right_child = newick_overlap_size(n->right, subset);
//                overlap_left_child -= overlap_right_child;
//            }
//        }
//    }
//
//    return n;
//}

double compute_bls_score(newick_node* node, const alignment_t & alignment)
{
    const uint64_t lo = 0;
    const uint64_t hi = alignment.seqs[0].size();
    double bl_total = 0.0;
    for (uint64_t i = lo; i < hi; ++i)
    {
        // determine subspecies set that has A,C,G or T at position i
        std::unordered_set<std::string> subset;
        for (uint16_t species_id = 0; species_id < alignment.seqs.size(); ++species_id)
        {
            if (get_dna_id(alignment.seqs[species_id][i]) <= 3)
                subset.insert(alignment.ids[species_id]);
        }
//        std::cout << subset.size() << ' ' << newick_sum_branch_lengths(node, subset) << '\n';
        bl_total += newick_sum_branch_lengths(node, subset);
    }

    std::unordered_set<std::string> all_species(alignment.ids.begin(), alignment.ids.end());
    const double divisor = newick_sum_branch_lengths(node, all_species) * (hi-lo);
//    std::cout << "sum: " << bl_total << '\n';
//    std::cout << "divisor: " << newick_sum_branch_lengths(node, all_species) << '\n';

    return bl_total/divisor;
}