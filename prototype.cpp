#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>

#include "src/newick.hpp"

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

double compute_bls_score(const uint64_t lo, const uint64_t hi)
{
    assert(hi >= lo);
//    let bls nt aln which_row lo hi =
//        let pres = function '-' | '.' | 'N' -> false | _ -> true
    double bl_total = 0.0;
//        for i = lo to hi do
//            let subtree =
//                Newick.subtree
//                    fun sp -> try pres aln.(SMap.find sp which_row).[i] with Not_found -> false
//                    nt
//            bl_total := !bl_total +. Option.map_default Newick.total_length 0. subtree
//        !bl_total /. (Newick.total_length nt *. float (hi-lo+1))
    return 0.0;
}

int main(int argc, char ** argv)
{
    std::string aln_path, model_path; //, selected_species;

    aln_path = "/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/tal-AA-tiny.fa";
    model_path = "/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Parameters/12flies";

    aln_path = "/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/ALDH2.exon5.fa";
    model_path = "/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Parameters/100vertebrates";

    std::unordered_set<std::string> selected_species = {"Human", "Mouse", "Rat", "Cow", "Horse", "Dog"};

    empirical_codon_model ecm_coding, ecm_noncoding;
    ecm_coding.open(std::string(model_path + "_coding.ECM").c_str());
    ecm_noncoding.open(std::string(model_path + "_noncoding.ECM").c_str());

    // open newick (and ignore spaces)
    newick_node* root = newick_open(std::string(model_path + ".nh").c_str());
    if (selected_species.size() > 0)
    {
        root = newick_reduce(root, selected_species);
//        std::cout << newick_print(root) << '\n';
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

//    std::cout << "Coding instance:\n";
//    std::cout << instance_coding;
//    std::cout << "NonCoding instance:\n";
//    std::cout << instance_noncoding;

    const double lpr_c = lpr_leaves(instance_coding, alignment, 1.0, nbr_leaves_in_tree);
    const double lpr_nc = lpr_leaves(instance_noncoding, alignment, 1.0, nbr_leaves_in_tree);

    const double decibans_score = (10.0 * (lpr_c - lpr_nc) / log(10.0));
    const double bls_score = 0.0; // TODO: compute_bls_score();
    printf("%f\t%f\t%f", decibans_score, bls_score);

    newick_free(root);

    return 0;
}
