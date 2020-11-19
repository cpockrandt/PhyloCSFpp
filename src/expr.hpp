#include <iostream>
#include <string>
#include <algorithm>

#include <cassert>

enum expr_op : uint8_t {
    ADD,
    SUB,
    MUL,
    DIV,
    VAL,
    VAR
};

struct expr_node // TODO: really bad implementation, try to get rid of expressions after reimplementation works!
{
    expr_op op;
    expr_node * left = NULL;
    expr_node * right = NULL;

    uint64_t variable_int; // for op == VAR
    double value_double; // for op == VAL

    expr_node(expr_op op, expr_node * left, expr_node * right, uint64_t variable_int, double value_double)
    {
        this->op = op;
        this->left = left;
        this->right = right;
        this->variable_int = variable_int;
        this->value_double = value_double;
    }
};

std::string expr_to_string(const expr_node * const n)
{
    std::string res = "";
    switch (n->op)
    {
        case ADD:
            res = "(" + expr_to_string(n->left) + " + " + expr_to_string(n->right) + ")";
            break;
        case SUB:
            res = "(" + expr_to_string(n->left) + " - " + expr_to_string(n->right) + ")";
            break;
        case MUL:
            res = "(" + expr_to_string(n->left) + " * " + expr_to_string(n->right) + ")";
            break;
        case DIV:
            res = "(" + expr_to_string(n->left) + " / " + expr_to_string(n->right) + ")";
            break;
        case VAL:
            res = std::to_string(n->value_double);
            break;
        case VAR:
            res = "Var " + std::to_string(n->variable_int);
            break;
        default:
            std::cerr << "Invalid expression type " << (unsigned)n->op << "!\n";
    }
    return res;
}

double expr_eval(expr_node * const n, const std::vector<double> & variables)
{
    double result;
    switch (n->op)
    {
        case ADD:
            result = expr_eval(n->left, variables) + expr_eval(n->right, variables);
            break;
        case SUB:
            result = expr_eval(n->left, variables) - expr_eval(n->right, variables);
            break;
        case MUL:
            result = expr_eval(n->left, variables) * expr_eval(n->right, variables);
            break;
        case DIV:
            result = expr_eval(n->left, variables) / expr_eval(n->right, variables);
            break;
        case VAL:
            result = n->value_double;
            break;
        case VAR:
            assert(n->variable_int < variables.size());
            result = variables[n->variable_int];
            break;
        default:
            std::cerr << "Invalid expression type " << (unsigned)n->op << "!\n";
    }
    return result;
}

void expr_free(expr_node * n)
{
    if (n->left != NULL)
        expr_free(n->left);
    if (n->right != NULL)
        expr_free(n->right);
    free(n);
}