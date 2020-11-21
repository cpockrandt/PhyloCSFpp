#include <iostream>
#include <string>
#include <algorithm>

#include <cassert>

struct newick_node
{
    int16_t id = -999; // not necessary here, but helps with flattening
    std::string label = "";
    double branch_length = .0f;
    newick_node* left = NULL;
    newick_node* right = NULL;
    newick_node* parent = NULL; // TODO: I think we can easily eliminate this field
    newick_node* sibling = NULL; // not necessary here, but helps with flattening
};

struct newick_elem
{
    int16_t id; // same as index in vector1
    int16_t child1_id;
    int16_t child2_id;
    int16_t sibling_id;
    int16_t parent_id;
    float branch_length;
    std::string label;
};

void newick_parse(std::string& str, newick_node* node, int16_t & nbr_nodes)
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
        node->branch_length = std::stod(len);
    }

    if (str[0] == '(')
    {
        nbr_nodes += 2;

        newick_node* left = new newick_node;
        left->parent = node;
        node->left = left;
        //
        str = str.substr(1); // remove '('
        newick_parse(str, left, nbr_nodes);

        newick_node* right = new newick_node;
        right->parent = node;
        node->right = right;
        newick_parse(str, right, nbr_nodes);
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
        node->parent->branch_length = std::stod(len);
    }
}

inline void newick_annotate_nodes(newick_node * node, int16_t & leaf_id, int16_t & inner_node_id)
{
    assert((node->left == NULL) == (node->right == NULL));
    if (node->left) // inner node
    {
        node->left->sibling = node->right;
        node->right->sibling = node->left;
        newick_annotate_nodes(node->left, leaf_id, inner_node_id);
        newick_annotate_nodes(node->right, leaf_id, inner_node_id);

        node->id = inner_node_id;
        ++inner_node_id;
    }
    else // leaf
    {
        node->id = leaf_id;
        ++leaf_id;
    }
}

uint16_t newick_count_leaves(const newick_node * const node)
{
    assert((node->left == NULL) == (node->right == NULL));
    if (node->left) // inner node
    {
        return newick_count_leaves(node->left) + newick_count_leaves(node->right);
    }
    else // leaf
    {
        return 1;
    }
}

inline newick_node* newick_parse(std::string& str)
{
    // TODO: remove trailing semicolon
    newick_node* root = new newick_node;
    root->parent = NULL;

    int16_t nbr_nodes = 1;
    newick_parse(str, root, nbr_nodes);

    int16_t leaf_id = 0;
    int16_t inner_node_id = (nbr_nodes / 2) + 1;
    // merge both traversals into one: during parsing remember how many nodes the tree has. then we now how many leaves. use two counters, one for leaves and for internal nodes
    newick_annotate_nodes(root, leaf_id, inner_node_id);

    return root;
}

inline newick_node* newick_open(const char * const file_path)
{
    char * linebuf;
    size_t linesiz = 0;

    FILE *fptr;
    if ((fptr = fopen(file_path, "r")) == NULL)
    {
        printf("Error! opening file");
        return NULL;
    }

    std::string str = "";
    while (getline(&linebuf, &linesiz, fptr) > 0)
    {
        str += std::string(linebuf);
    }

    free(linebuf);
    linebuf = NULL;

    fclose(fptr);

    // remove any newlines, spaces, etc.
    str.erase(std::remove_if(str.begin(), str.end(), isspace), str.end());
//    std::cout << str << '\n';
    return newick_parse(str);
}

void newick_flatten_add_leaves(newick_node * node, std::vector<newick_elem> & newick_flattened)
{
    // store leaves first
    assert((node->left == NULL) == (node->right == NULL));
    if (node->left)
    {
        newick_flatten_add_leaves(node->left, newick_flattened);
        newick_flatten_add_leaves(node->right, newick_flattened);
    }
    else
    {
        // encountered leaf
        newick_elem elem;
        elem.id = node->id; // id
        elem.child1_id = -1; // child2_id
        elem.child2_id = -1; // child1_id
        elem.sibling_id = node->sibling->id; // sibling_id
        elem.parent_id = node->parent->id; // parent_id
        elem.branch_length = node->branch_length; // branch_lengt
        elem.label = node->label; // label

        newick_flattened.push_back(std::move(elem));
    }
}

void newick_flatten_add_inner_nodes(newick_node * node, std::vector<newick_elem> & newick_flattened)
{
    // store leaves first
    assert((node->left == NULL) == (node->right == NULL));
    if (node->left)
    {
        newick_flatten_add_inner_nodes(node->left, newick_flattened);
        newick_flatten_add_inner_nodes(node->right, newick_flattened);

        // encountered inner node
        newick_elem elem;
        elem.id = node->id; // id
        elem.child1_id = (node->left != NULL) ? node->left->id : -1; // child1_id
        elem.child2_id = (node->right != NULL) ? node->right->id : -1; // child2_id
        elem.sibling_id = (node->parent != NULL) ? node->sibling->id : -1; // sibling_id
        elem.parent_id = (node->parent != NULL) ? node->parent->id : -1; // parent_id
        elem.branch_length = node->branch_length; // branch_lengt
        elem.label = node->label; // label

        newick_flattened.push_back(std::move(elem));
    }
}

inline void newick_flatten(newick_node * root, std::vector<newick_elem> & newick_flattened)
{
    // TODO: merge this into only one traversal
    newick_flatten_add_leaves(root, newick_flattened);
    newick_flatten_add_inner_nodes(root, newick_flattened);
}

inline std::string newick_print(newick_node * n)
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