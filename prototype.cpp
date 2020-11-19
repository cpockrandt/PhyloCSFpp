#include <iostream>
#include <vector>
#include <string>

#include "src/alignment_reader.hpp"
#include "src/newick.hpp"
#include "src/ecm.hpp"
#include "src/instance.hpp"

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

    empirical_codon_model ecm_coding, ecm_noncoding;
    ecm_coding.open(std::string(model_path + "_coding.ECM").c_str());
    ecm_noncoding.open(std::string(model_path + "_noncoding.ECM").c_str());
//    ecm_coding.print_model();

    // read alignment
    std::vector<std::string> ids;
    std::vector<std::string> seqs;
    std::vector<std::string> peptides;
    read_alignment(aln_path.c_str(), ids, seqs);

    for (uint16_t i = 0; i < seqs.size(); ++i)
    {
        peptides.push_back(translate(seqs[i]));
        std::cout << ids[i] << '\t' << seqs[i] << '\t' << peptides[i] << '\n';
    }

    // open newick (and ignore spaces)
    newick_node* root = newick_open(std::string(model_path + ".nh").c_str());
//    std::cout << newick_print(root) << '\n';
    // transform tree into array notation
    std::vector<newick_elem> newick_flattened;
    newick_flatten(root, newick_flattened);
    for (const auto & elem : newick_flattened)
        std::cout << elem << '\n';

    // TODO: initialize model
    instance_t instance_coding, instance_noncoding;
    make(instance_coding, ecm_coding, newick_flattened);
//    make(instance_noncoding, ecm_noncoding, newick_flattened);

    return 0;
}
