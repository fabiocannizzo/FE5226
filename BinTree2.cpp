/*
 Author: Fabio Cannizzo
 NUS FE5226 Example Program

 Construct a binary tree and demonstrate a method
 to visit it in pre-order using callbacks
*/

#include <iostream>

using namespace std;

struct Node
{
    int value;   // the value of the node
    Node *left;  // the node to the left
    Node *right; // the node to the right
};

// create a node and assign to it the given value and pointers
Node *createNode(int value, Node *left, Node *right)
{
    Node *n = new Node;
    n->value = value;
    n->left = left;
    n->right = right;
    return n;
}

Node *buildTree()
{
    Node *n1, *n2, *n3;

    n1 = createNode(4, NULL, NULL);
    n1 = createNode(2, n1, NULL);

    n2 = createNode(7, NULL, NULL);
    n2 = createNode(6, NULL, n2);

    n3 = createNode(5, NULL, NULL);
    n3 = createNode(3, n3, n2);

    n1 = createNode(1, n1, n3);

    return n1;
}

void visitTreeInPreOrder(Node *node, void(*action)(Node*, void*), void *args )
{
    if (node != NULL) {
        action(node, args);                             // performs action
        visitTreeInPreOrder(node->left, action, args);  // visit left
        visitTreeInPreOrder(node->right, action, args); // visit right 
    }
}

void printNode(Node *n, void *)
{
    cout << n->value << " ";  // no end of line
}

void incrNode(Node *n, void *a)
{
    int *p = reinterpret_cast<int*>(a);
    n->value += *p;
}

int main()
{
    Node *treeTop = buildTree();
    
    visitTreeInPreOrder( treeTop, &printNode, NULL );
    cout << endl;

    int x = 3;
    visitTreeInPreOrder( treeTop, &incrNode, &x  );
    
    visitTreeInPreOrder( treeTop, &printNode, NULL );
    cout << endl;
    
    return 0;
}
