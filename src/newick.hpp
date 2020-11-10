#include <iostream>
#include <string>

#include <cassert>

struct newick_node
{
    std::string label = "";
    float branch_length = .0f;
    newick_node* left = NULL;
    newick_node* right = NULL;
    newick_node* parent = NULL; // TODO: I think we can easily eliminate this field
};

void newick_parse(std::string& str, newick_node* node)
{
    if (str[0] != '(' && str[0] != ',' && str[0] != ')' && str[0] != ':')
    {
        while (str[0] != '(' && str[0] != ',' && str[0] != ')' && str[0] != ':')
        {
            node->label += str[0];
            str = str.substr(1);
        }

        str = str.substr(1); // remove ':'

        std::string len = "";
        while (std::isdigit(str[0]) || str[0] == '.')
        {
            len += str[0];
            str = str.substr(1);
        }
        node->branch_length = std::stof(len);
    }

    if (str[0] == '(')
    {
        newick_node* left = new newick_node;
        left->parent = node;
        node->left = left;
        //
        str = str.substr(1); // remove '('
        newick_parse(str, left);

        newick_node* right = new newick_node;
        right->parent = node;
        node->right = right;
        newick_parse(str, right);
    }

    if (str[0] == ',')
    {
        str = str.substr(1); // remove ','
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
    }
}

inline newick_node* newick_parse(std::string& str)
{
    // TODO: remove trailing semicolon
    newick_node* root = new newick_node;
    root->parent = NULL;

    newick_parse(str, root);

    return root;
}

inline std::string newick_print(newick_node* n)
{
    assert(n != NULL);

    if (n->left == NULL && n->right == NULL)
    {
        return n->label + ":" + std::to_string(n->branch_length);
    }
    else
    {
        std::string bl = "";
        if (n->parent != NULL)
            bl = ":" + std::to_string(n->branch_length);
        return "(" + newick_print(n->left) + "," + newick_print(n->right) + ")" + bl;
    }
}

inline void newick_free(newick_node* n)
{
    if (n->left != NULL)
        newick_free(n->left);
    if (n->right != NULL)
        newick_free(n->right);
    free(n);
}