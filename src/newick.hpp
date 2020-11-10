#include <iostream>
#include <string>

#include <cassert>

struct newick_node
{
    std::string label = "";
    float branch_length = .0f;
    newick_node* left_child = NULL;
    newick_node* right_child = NULL;
    newick_node* parent = NULL;
};

void newick_parse(std::string& str, newick_node* node, int i, const bool debug)
{
    if (debug)
        std::cout << std::string(2*i, '-') << "Current string: " << str << '\n';

    if (str[0] != '(' && str[0] != ',' && str[0] != ')' && str[0] != ':')
    {
        while (str[0] != '(' && str[0] != ',' && str[0] != ')' && str[0] != ':')
        {
            node->label += str[0];
            str = str.substr(1);
        }
        if (debug)
            std::cout << std::string(2*i, '-') << "Set label " << node->label << '\n';

        str = str.substr(1); // remove ':'

        std::string len = "";
        while (std::isdigit(str[0]) || str[0] == '.')
        {
            len += str[0];
            str = str.substr(1);
        }
        node->branch_length = std::stof(len);
        if (debug)
            std::cout << std::string(2*i, '-') << "Set branch length " << std::to_string(node->branch_length) << '\n';
    }

    if (str[0] == '(')
    {
        newick_node* left = new newick_node;
        left->parent = node;
        node->left_child = left;
        //
        str = str.substr(1); // remove '('
        if (debug)
            std::cout << std::string(2*i, '-') << "Created left child\n";
        newick_parse(str, left, i + 1, debug);

        newick_node* right = new newick_node;
        right->parent = node;
        node->right_child = right;
        if (debug)
            std::cout << std::string(2*i, '-') << "Created right child\n";
        newick_parse(str, right, i + 1, debug);
    }

    if (str[0] == ',')
    {
        str = str.substr(1); // remove ','
        if (debug)
            std::cout << std::string(2*i, '-') << "Go up one parent to go to the right next!\n";
        return; // go up one
    }

    if (str[0] == ')')
    {
        str = str.substr(1); // remove ')'

        if (str.size() == 0 || str[0] == ';')
            return;

        assert(str[0] == ':');
        str = str.substr(1); // remove ':'

        std::string len = "";
        while (std::isdigit(str[0]) || str[0] == '.')
        {
            len += str[0];
            str = str.substr(1);
        }
        node->parent->branch_length = std::stof(len);
        if (debug)
            std::cout << std::string(2*i, '-') << "Set branch length " << std::to_string(node->parent->branch_length) << '\n';
    }
}

inline newick_node* newick_parse(std::string& str, const bool debug)
{
    // TODO: remove trailing semicolon
    newick_node* root = new newick_node;
    root->parent = NULL;

    newick_parse(str, root, 0, debug);

    return root;
}

inline std::string newick_print(newick_node* n)
{
    assert(n != NULL);

    if (n->left_child == NULL && n->right_child == NULL)
    {
        return n->label + ":" + std::to_string(n->branch_length);
    }
    else
    {
        std::string bl = "";
        if (n->parent != NULL)
            bl = ":" + std::to_string(n->branch_length);
        return "(" + newick_print(n->left_child) + "," + newick_print(n->right_child) + ")" + bl;
    }
}

inline void newick_free(newick_node* n)
{
    if (n->left_child != NULL)
        newick_free(n->left_child);
    if (n->right_child != NULL)
        newick_free(n->right_child);
    free(n);
}