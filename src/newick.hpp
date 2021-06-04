#pragma once

#include <string>
#include <algorithm>
#include <unordered_set>

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

void _newick_parse(std::string& str, newick_node* node)
{
    if (str[0] != '(' && str[0] != ',' && str[0] != ')' && str[0] != ':')
    {
        while (str[0] != '(' && str[0] != ',' && str[0] != ')' && str[0] != ':')
        {
            node->label += tolower(str[0]);
            str = str.substr(1); // TODO: inefficient!
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
        newick_node* left = new newick_node;
        left->parent = node;
        node->left = left;
        //
        str = str.substr(1); // remove '(' // TODO: inefficient
        _newick_parse(str, left);

        newick_node* right = new newick_node;
        right->parent = node;
        node->right = right;
        _newick_parse(str, right);
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

// temporary wrapper to make a copy of the string
void newick_parse(std::string str, newick_node* node)
{
    _newick_parse(str, node);
    assert(node->branch_length == 0.0);
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

inline newick_node* newick_open(const char * const file_path)
{
    char * linebuf = NULL;
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

    // TODO: remove trailing semicolon
    newick_node* root = new newick_node;
    root->parent = NULL;

    newick_parse(str, root);

    assert(root->branch_length == 0.0);

    return root;
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
        elem.branch_length = node->branch_length; // branch_length
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
    // annotate tree (ids, parents, siblings, etc)
    int16_t leaf_id = 0;
    int16_t inner_node_id = newick_count_leaves(root); // (nbr_nodes / 2) + 1;
    // merge both traversals into one: during parsing remember how many nodes the tree has. then we now how many leaves. use two counters, one for leaves and for internal nodes
    newick_annotate_nodes(root, leaf_id, inner_node_id);

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
    delete n;
}

bool newick_is_leaf(const newick_node* node)
{
    assert((node->left == NULL) == (node->right == NULL));
    return node->left == NULL;
}

// cannot be called on NULL!
uint16_t newick_overlap_size(const newick_node* node, const std::unordered_set<std::string> & subset)
{
    if (newick_is_leaf(node))
        return subset.find(node->label) != subset.end();
    else
        return newick_overlap_size(node->left, subset) + newick_overlap_size(node->right, subset);
}

void newick_check_missing_species(newick_node* node, std::unordered_set<std::string> & selected_species)
{
    if (node->left == NULL)
    {
        // leaf
        selected_species.erase(node->label);
    }
    else
    {
        newick_check_missing_species(node->left, selected_species);
        newick_check_missing_species(node->right, selected_species);
    }
}

void newick_reduce(newick_node* node, const std::unordered_set<std::string> & subset)
{
    if (node->left == NULL)
        return;

    const uint16_t overlap_left_child = newick_overlap_size(node->left, subset);
    const uint16_t overlap_right_child = newick_overlap_size(node->right, subset); // store nbr of hits in parent node to avoid two passes!

    if (overlap_left_child == 0)
    {
        newick_node* left_node_old = node->left;
        newick_node* right_node_old = node->right;

        // merge current node and right node, because left node is about to get deleted
        node->left = right_node_old->left;
        node->right = right_node_old->right;

        // update parents and add branch lengths when we merge two nodes
        if (node->left != NULL)
        {
            assert(node->right != NULL);
            node->left->parent = node;
            node->right->parent = node;
        }
        else
        {
            assert(node->right == NULL);
            node->label = right_node_old->label;
        }

        if (node->parent != NULL) // parent node does not have branch length
            node->branch_length += right_node_old->branch_length;

        // delete left subtree
        delete right_node_old;
        newick_free(left_node_old);

        // continue on right subtree (which is now "node")
        newick_reduce(node, subset);
    }
    else if (overlap_right_child == 0)
    {
        newick_node* left_node_old = node->left;
        newick_node* right_node_old = node->right;

        // merge current node and left node, because right node is about to get deleted
        node->left = left_node_old->left;
        node->right = left_node_old->right;

        // update parents and add branch lengths when we merge two nodes
        if (node->left != NULL)
        {
            assert(node->right != NULL);
            node->left->parent = node;
            node->right->parent = node;
        }
        else
        {
            assert(node->right == NULL);
            node->label = left_node_old->label;
        }

        if (node->parent != NULL) // parent node does not have branch length
            node->branch_length += left_node_old->branch_length;

        // delete right subtree
        delete left_node_old;
        newick_free(right_node_old);

        // continue on left subtree (which is now "node")
        newick_reduce(node, subset);
    }
    else
    {
        newick_reduce(node->left, subset);
        newick_reduce(node->right, subset);
    }
}