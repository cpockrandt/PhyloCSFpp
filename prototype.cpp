#include <iostream>
#include <vector>
#include <string>
#include <regex>

#include "src/newick.hpp"

// #include <seqan/sequence.h>
//
// using namespace seqan;

std::regex reg1("0\\.");
std::regex reg2("0+([\\),])");

bool test_newick(std::string str)
{
    newick_node* root;
    std::string str_input = str;

    str = regex_replace(str, reg1, ".");
    str = regex_replace(str, reg2, "$1");

    std::cout << "Before: " << str << "\n";
    root = newick_parse(str_input);
    std::string after = newick_print(root);

    after = regex_replace(after, reg1, ".");
    after = regex_replace(after, reg2, "$1");

    std::cout << "After : " << after << "\n" << std::endl;
    newick_free(root);

    if (str != after)
    {
        std::cout << "ERROR!\n";
        exit(-1);
    }

    return str == after;
}

struct newick_elem
{
//    uint16_t id; // same as index in vector1
    uint16_t child1_id;
    uint16_t child2_id;
    uint16_t sibling_id;
    uint16_t parent_id;
    float branch_length;
    std::string label;
};

int main(int argc, char ** argv)
{
    // std::vector<std::string> ids {
    //     "dmel",
    //     "dana",
    //     "dpse",
    //     "dwil",
    //     "dvir"
    // };
    //
    // std::vector<std::string> seq {
    //     "ATGCTGGATCCCACTGGAACATACCGGCGACCACGCGACACGCAGGACTCCCGCCAAAAGAGGCGACAGGACTGCCTGGATCCAACCGGGCAGTAC",
    //     "ATGCTGGATCCCACAGGAACGTACAGGCGACCCCGCGACTCGCAGGACACACGCCAGAAGCGGCGCCAGGATTGCCTGGATCCAACCGGGCAGTAC",
    //     "ATGCTGGATCCCACAGGAACCTACCGTCGCCCACGCGACACAAAAGACACACGCCAGAAACGGCGTCAGGACTGCCTAGATCCCACCGGGCAGTAC",
    //     "ATGCTGGATCCAACTGGTACATATCGCCGACCAAGGGATACACAGGACACACGCCATAAACGGC---------GCCTGGATCCAACTGGACAATAT",
    //     "ATGCTGGATCCAACGGGCACCTATCGGCGGCCGCGTGAGTCGCGTGACACGCGCCACAAGCAGCGGCAGCTGGCGCTGGATCCTACCGGGCAGTAC"
    // };

    std::string newick_str = "((((((Scer:.0836,Spar:.0617):.0432,Smik:.1154):.0165,Skud:.1731):.0275,Sbay:.1594):.0770,Scas:.1099):.0550,Sklu:.2034)";

    test_newick("(A:.420000,BB:0.111)");
    test_newick("((A:.123,BB:.456):0.111,CCC:0.789)");
    test_newick("(DDDD:0.1,((A:0.2,BB:0.3):0.4,CCC:0.5):0.6)");

    test_newick(newick_str);
    test_newick("((((((((((((((((((Human:0.00655,Chimp:0.00684):0.00422,Gorilla:0.008964):0.009693,Orangutan:0.01894):0.003471,Gibbon:0.02227):0.01204,(((Rhesus:0.004991,Crab_eating_macaque:0.004991):0.003,Baboon:0.008042):0.01061,Green_monkey:0.027):0.025):0.02183,(Marmoset:0.03,Squirrel_monkey:0.01035):0.01965):0.07261,Bushbaby:0.13992):0.013494,Chinese_tree_shrew:0.174937):0.002,(((Squirrel:0.125468,(Lesser_Egyptian_jerboa:0.1,((Prairie_vole:0.08,(Chinese_hamster:0.04,Golden_hamster:0.04):0.04):0.06,(Mouse:0.084509,Rat:0.091589):0.047773):0.06015):0.1):0.022992,(Naked_mole_rat:0.1,(Guinea_pig:0.065629,(Chinchilla:0.06,Brush_tailed_rat:0.1):0.06):0.05):0.06015):0.025746,(Rabbit:0.114227,Pika:0.201069):0.101463):0.015313):0.020593,(((Pig:0.12,((Alpaca:0.047275,Wild_bactrian_camel:0.04):0.04,((Dolphin:0.034688,Killer_whale:0.039688):0.03,(Tibetan_antelope:0.1,(Cow:0.1,(Sheep:0.05,Domestic_goat:0.05):0.05):0.01):0.013592):0.025153):0.020335):0.02,(((Horse:0.059397,White_rhinoceros:0.025):0.05,(Cat:0.098612,(Dog:0.052458,(Ferret:0.05,(Panda:0.02,(Pacific_walrus:0.02,Weddell_seal:0.02):0.02):0.03):0.03):0.02):0.049845):0.006219,((Black_flying_fox:0.05,Megabat:0.063399):0.05,(Big_brown_bat:0.02,(Davids_myotis:0.04,Microbat:0.04254):0.05):0.06):0.033706):0.004508):0.011671,(Hedgehog:0.221785,(Shrew:0.169562,Star_nosed_mole:0.1):0.1):0.056393):0.021227):0.023664,(((((Elephant:0.002242,Cape_elephant_shrew:0.05):0.04699,Manatee:0.1):0.049697,(Cape_golden_mole:0.03,Tenrec:0.235936):0.01):0.03,Aardvark:0.03):0.02,Armadillo:0.169809):0.006717):0.234728,(Opossum:0.125686,(Tasmanian_devil:0.1,Wallaby:0.072008):0.05):0.2151):0.071664,Platypus:0.456592):0.109504,(((((Rock_pigeon:0.1,((Saker_falcon:0.1,Peregrine_falcon:0.1):0.03,(((Collared_flycatcher:0.04,((White_throated_sparrow:0.034457,Medium_ground_finch:0.041261):0.015,Zebra_finch:0.06):0.01):0.052066,Tibetan_ground_jay:0.06):0.025,(Budgerigar:0.046985,(Puerto_Rican_parrot:0.026,Scarlet_macaw:0.026):0.01):0.04):0.064703):0.06):0.05,(Mallard_duck:0.1,(Chicken:0.041254,Turkey:0.085718):0.031045):0.09):0.22,American_alligator:0.25):0.045143,((Green_seaturtle:0.1,Painted_turtle:0.1):0.05,(Chinese_softshell_turtle:0.1,Spiny_softshell_turtle:0.1):0.05):0.04):0.01,Lizard:0.447):0.122):0.05,Frog_X._tropicalis:0.977944):0.1,Coelacanth:0.977944):0.111354,(((((((Tetraodon:0.124159,(Fugu:0.103847,Yellowbelly_pufferfish:0.1):0.1):0.09759,(Nile_tilapia:0.1,(Princess_of_Burundi:0.05,(Burtons_mouthbreeder:0.05,(Zebra_mbuna:0.05,Pundamilia_nyererei:0.05):0.05):0.05):0.1):0.1):0.09759,(Medaka:0.38197,Southern_platyfish:0.4):0.1):0.015,Stickleback:0.246413):0.045,Atlantic_cod:0.25):0.22564,(Zebrafish:0.430752,Mexican_tetra:0.4):0.3):0.143632,Spotted_gar:0.4):0.326688):0.2,Lamprey:0.975747)");

    // TODO: transform this into array notation
    std::vector<newick_elem>

    return 0;
}
