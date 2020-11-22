#include <iostream>
#include <vector>
#include <string>

#include "src/alignment_reader.hpp"
#include "src/newick.hpp"
#include "src/ecm.hpp"
#include "src/instance.hpp"
#include "src/fixed_lik.hpp"

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

int main(int argc, char ** argv)
{
    std::string aln_path = "/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/tal-AA-tiny.fa";
    std::string model_path = "/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Parameters/23flies";

//    std::string aln_path = "/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/ALDH2.exon5.fa";
//    std::string model_path = "/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Parameters/100vertebrates";

    empirical_codon_model ecm_coding, ecm_noncoding;
    ecm_coding.open(std::string(model_path + "_coding.ECM").c_str());
    ecm_noncoding.open(std::string(model_path + "_noncoding.ECM").c_str());
//    ecm_coding.print_model();

    // read alignment
    alignment_t alignment;
    read_alignment(aln_path.c_str(), alignment);
    for (uint16_t i = 0; i < alignment.seqs.size(); ++i)
        std::cout << alignment.ids[i] << '\t' << alignment.seqs[i] << '\t' << print_peptide(alignment.peptides[i]) << '\n';

    // open newick (and ignore spaces)
    newick_node* root = newick_open(std::string(model_path + ".nh").c_str());
//    std::cout << newick_print(root) << '\n';
    // transform tree into array notation
    std::vector<newick_elem> newick_flattened;
    newick_flatten(root, newick_flattened);
    for (const auto & elem : newick_flattened)
        std::cout << elem << '\n';

    // initialize model
    const uint16_t nbr_leaves_in_tree = newick_count_leaves(root);
    instance_t instance_coding, instance_noncoding;
    PhyloCSFModel_make(instance_coding, ecm_coding, newick_flattened);
    PhyloCSFModel_make(instance_noncoding, ecm_noncoding, newick_flattened);

    const double lpr_c = lpr_leaves(instance_coding, alignment, 1.0, nbr_leaves_in_tree);
    const double lpr_nc = lpr_leaves(instance_noncoding, alignment, 1.0, nbr_leaves_in_tree);

    const double decibans_score = (10.0 * (lpr_c - lpr_nc) / log(10.0));
    printf("%.5f\t%.5f\t%.5f\n", lpr_c, lpr_nc, decibans_score);

    newick_free(root);

    return 0;
}
