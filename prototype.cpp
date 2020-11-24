#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>

#include "src/newick.hpp"
#include "src/models.hpp"

std::ostream& operator<<(std::ostream & os, const newick_elem & e)
{
    os  << e.id << '\t'
        << e.label << '\t'
        << e.parent_id << '\t'
        << '(' << e.child1_id << ", " << e.child2_id << ")\t"
        << e.sibling_id << '\t'
        << e.branch_length;
    return os;
}

std::ostream& operator<<(std::ostream & os, const std::vector<newick_elem> & v)
{
    for (uint32_t i = 0; i < v.size(); ++i)
    {
        os << v[i] << '\n';
    }
    return os;
}

std::ostream& operator<<(std::ostream & os, const std::vector<double> & v)
{
    for (uint32_t i = 0; i < v.size(); ++i)
    {
        char buf[10];
        sprintf (buf, "%f\t", v[i]);
        os << buf;
    }
    std::cout << '\n';
    return os;
}

std::ostream& operator<<(std::ostream & os, const std::vector<std::vector<double> > & v)
{
    for (uint32_t i = 0; i < v.size(); ++i)
    {
        for (uint32_t j = 0; j < v[i].size(); ++j)
        {
            char buf[10];
            sprintf (buf, "%.3f ", v[i][j]);
            os << buf << ' ';
        }
        std::cout << '\n';
    }
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream & os, const std::vector<T> & v)
{
    std::cout << '[';
    for (uint32_t i = 0; i < v.size(); ++i)
    {
        os << v[i] << ' ';
    }
    std::cout << ']' << '\n';
    return os;
}

#include "src/alignment_reader.hpp"
#include "src/ecm.hpp"
#include "src/instance.hpp"
#include "src/fixed_lik.hpp"

#include <gsl/gsl_blas.h>

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

int main(int argc, char ** argv)
{
    std::unordered_set<std::string> selected_species;
    std::string aln_path, model_path, model_str; //, selected_species;

    model_str = "23flies";
    aln_path = "/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/tal-AA-tiny3.fa";
    model_path = "/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Parameters/23flies";
    selected_species = {"dmel", "dsim", "dsec"};

    model_str = "100vertebrates";
    aln_path = "/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/ALDH2.exon5.fa";
    model_path = "/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Parameters/100vertebrates";
    selected_species = {/*"Human", "Mouse", "Rat", "Cow", "Horse", "Dog"*/};

    empirical_codon_model ecm_coding, ecm_noncoding;
    if (false)
    {
        ecm_coding.open(std::string(model_path + "_coding.ECM").c_str());
        ecm_noncoding.open(std::string(model_path + "_noncoding.ECM").c_str());
    }
    else
    {
        if (!ecm_coding.load(model_str.c_str(), models, 1))
        {
            printf("ERROR: The model \"23flies\" does not exist. Please choose one from the following set or give a path to your model:\n");
            for (const auto & m : models)
                printf("\t- %s\n", m.first.c_str());
            return -1;
        }
        ecm_noncoding.load(model_str.c_str(), models, 0);
    }

    // open newick (and ignore spaces)
    // TODO: make sure not to delete species, that are NOT in the alignment, but that are passed as "selected_species". relevant for branch length score
    newick_node* root = newick_open(std::string(model_path + ".nh").c_str());
    assert(root->branch_length == 0.0);

//    std::cout << newick_print(root) << '\n';
    if (selected_species.size() > 0)
    {
        std::unordered_set<std::string> missing_species = selected_species; // merge into one set with bool
        newick_check_missing_species(root, missing_species);
        if (missing_species.size() > 0)
        {
            printf("ERROR: The following selected species are missing in the phylogenetic tree (TODO: output filename):\n");
            for (const std::string & sp : selected_species)
                printf("\t- %s\n", sp.c_str());
            newick_free(root);
            return -1;
        }

        newick_reduce(root, selected_species, true);
        std::cout << newick_print(root) << '\n';
        assert(root->branch_length == 0.0);
    }

    // annotate tree (ids, parents, siblings, etc)
    int16_t leaf_id = 0;
    int16_t inner_node_id = newick_count_leaves(root); // (nbr_nodes / 2) + 1;
    // merge both traversals into one: during parsing remember how many nodes the tree has. then we now how many leaves. use two counters, one for leaves and for internal nodes
    newick_annotate_nodes(root, leaf_id, inner_node_id);

    // transform tree into array notation
    std::vector<newick_elem> newick_flattened;
    newick_flatten(root, newick_flattened);
//    for (const auto & elem : newick_flattened)
//        std::cout << elem << '\n';

    // read alignment
    alignment_t alignment;
    const uint16_t nbr_leaves = (newick_flattened.size() + 1)/2;
    alignment.ids.resize(nbr_leaves);
    alignment.seqs.resize(nbr_leaves, "");
    alignment.peptides.resize(nbr_leaves, {});
    for (uint16_t i = 0; i < nbr_leaves; ++i)
    {
        alignment.ids[i] = newick_flattened[i].label;
    }

    read_alignment(aln_path.c_str(), alignment);
//    for (uint16_t i = 0; i < alignment.seqs.size(); ++i)
//        std::cout << alignment.ids[i] << '\t' << alignment.seqs[i] << '\t' << print_peptide(alignment.peptides[i]) << '\n';

    // initialize model
    const uint16_t nbr_leaves_in_tree = newick_count_leaves(root);
    instance_t instance_coding, instance_noncoding;
    PhyloCSFModel_make(instance_coding, ecm_coding, newick_flattened);
    PhyloCSFModel_make(instance_noncoding, ecm_noncoding, newick_flattened);

//    std::cout << "Coding instance:\n" << instance_coding;
//    std::cout << "NonCoding instance:\n" << instance_noncoding;

    double lpr_c, lpr_nc;
    double elpr_anc_c, elpr_anc_nc;
    lpr_leaves(instance_coding, alignment, 1.0, nbr_leaves_in_tree, lpr_c, elpr_anc_c);
    lpr_leaves(instance_noncoding, alignment, 1.0, nbr_leaves_in_tree, lpr_nc, elpr_anc_nc);

    const double decibans_score = 10.0 * (lpr_c - lpr_nc) / log(10.0);
    const double anchestral_score = 10.0 * (elpr_anc_c - elpr_anc_nc) / log(10.0); // TODO
    const double bls_score = compute_bls_score(root, alignment);
    printf("%f\t%f\t%f", decibans_score, bls_score, anchestral_score);

    newick_free(root);

    return 0;
}
