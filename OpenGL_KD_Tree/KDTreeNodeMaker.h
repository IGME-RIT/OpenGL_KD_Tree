#ifndef _KDTREENODEMAKER_H
#define _KDTREENODEMAKER_H

#include "all_Fixed_Data.h"

#pragma region KD_Tree and Node maker

enum cutType {
	VERTICAL = 100,
	HORIZONTAL,
	LEAF,
	INITIAL_CUT = VERTICAL
};

// Node in the k-d tree
class KD_Node
{
public:
	//constructor and destructor
	KD_Node();
	KD_Node(Tree_Point t_point);
	~KD_Node();

	bool isleaf() { return (type == LEAF); }
	bool isRoot() { return root; }

	//set and get functions
	Tree_Point Get_point() { return this->p_tree; }
	KD_Node* Get_leftTree() { return this->left_Tree; };
	KD_Node* Get_rightTree() { return this->right_Tree; };
	int Get_type() { return this->type; }
	int Get_depth() { return this->depth; }
	void set_Type(int T) { this->type = T; }
	void set_Root(bool bool_r) { this->root = bool_r; }
	void set_Depth(int D) { this->depth = D; }
	void set_Point(Tree_Point t_point) { this->p_tree = t_point; }
	void set_Left(KD_Node* L) { this->left_Tree = L; }
	void set_Right(KD_Node* R) { this->right_Tree = R; }

	// print functions
	void print_Info();
	void print_Type();
	void print_Point(Tree_Point p) { this->p1 = p; };

private:
	int type;
	bool root;
	int depth;
	Tree_Point p_tree; // if node is a leaf, p_tree is contained in the region; else point is the median.
	KD_Node* left_Tree;
	KD_Node* right_Tree;
	Tree_Point p1;


};

class KD_tree
{
public:
	//constructor and destructor
	KD_tree();
	KD_tree(std::vector<Tree_Point> all_points);
	~KD_tree();

	//helper functions
	void initial_tree(std::vector<Tree_Point> all_points);
	KD_Node* KD_Maker(std::vector<Tree_Point> X_point, std::vector<Tree_Point> Y_point, int depth);
	void remaker_KD(std::vector<Tree_Point> all_points);
	int compute_KD_Height();
	void levelOrderPts();
	void add_Level(KD_Node* node, int depth);
	void deallocate_tree(KD_Node* node);


	// check for even points
	bool Even_pts() { return (num_kd_nodes % 2 == 0); }

	//get functions
	int get_treeHeight() { return this->height_Tree; }
	int get_num_KDNodes() { return this->num_kd_nodes; }
	KD_Node * get_Root() { return this->root; }
	std::vector<KD_Node*> get_Points() { return this->ordered_pts_set; }

	//set functions
	void set_TreeRoot(KD_Node * r) { this->root = r; }
	void set_TreeHeight(int h) { this->height_Tree = h; }
	void init_Height() { this->height_Tree = 1; }
	void set_numKDNodes(int num) { this->num_kd_nodes = num; }
	void set_Pts(std::vector<Tree_Point> all_points) { this->pts_set = all_points; }

	//print functions
	void printTree();
	void printInfo();
	void printNumNodes();


private:
	KD_Node * root;
	int height_Tree;
	int num_kd_nodes;
	std::vector<Tree_Point> pts_set;
	std::vector<KD_Node*> ordered_pts_set;


	//these structs sort the point based on coordinate.
	struct sort_X_coords
	{
		bool operator() (const Tree_Point &P, const Tree_Point &Q)
		{
			if (fabs(P.xpos - Q.xpos) < EPSILON)
			{
				return (P.ypos < Q.ypos);
			}
			else
			{
				return (P.xpos < Q.xpos);
			}

		}

	}sort_X_coords;

	struct _sort_Y_coords
	{
		bool operator() (const Tree_Point &P, const Tree_Point &Q)
		{
			if (fabs(P.ypos - Q.ypos) < EPSILON)
			{
				return (P.xpos < Q.xpos);
			}
			else
			{
				return (P.ypos < Q.ypos);
			}

		}

	}sort_Y_coords;
};

#pragma endregion

#pragma region Kd tree

//Constructor
KD_Node::KD_Node()
{
	Tree_Point empty_point;
	empty_point.xpos = 0.0f;
	empty_point.ypos = 0.0f;
	set_Point(empty_point);
	set_Type(LEAF);
	set_Root(false);
	set_Depth(0);
	set_Left(NULL);
	set_Right(NULL);
}


KD_Node::KD_Node(Tree_Point t_point)
{
	set_Point(t_point);

	//sets each point as a leaf when the kd-tree is built
	set_Type(LEAF);
	set_Root(false);
	set_Depth(1);
	set_Left(NULL);
	set_Right(NULL);
}


//print functions
void KD_Node::print_Info()
{
	print_Point(Get_point());
	print_Type();
	if (!isRoot()) {
		printf("Depth: %d\n", Get_depth());
	}

}


//Print functions
void KD_Node::print_Type()
{
	string type;
	if (isRoot()) return; // tree will state if it's printing the root
	if (Get_type() == VERTICAL) {
		type = "VERTICAL";
	}
	else if (Get_type() == HORIZONTAL) {
		type = "HORIZONTAL";
	}
	else if (Get_type() == LEAF) {
		type = "LEAF";
	}
	printf("Type: %s\n", type.c_str());
}

//Destructor
KD_Node::~KD_Node() {

}

//Default Constructor
KD_tree::KD_tree()
{
	set_TreeRoot(NULL);
	set_TreeHeight(0);
	set_numKDNodes(0);

}

//Constructor
KD_tree::KD_tree(std::vector<Tree_Point> all_points)
{
	sort(all_points.begin(), all_points.end(), sort_X_coords); //sort input vector by x-coordinates
	set_Pts(all_points); // saves the points as initialized to class
	set_numKDNodes(all_points.size());
	init_Height(); //sets the height of k-d tree

	initial_tree(all_points);
}

void KD_tree::remaker_KD(std::vector<Tree_Point> all_points)
{
	//clears vectors and deallocate tree
	pts_set.clear();
	ordered_pts_set.clear();
	deallocate_tree(get_Root());
	set_TreeRoot(NULL);


	//build tree again with same constructor
	sort(all_points.begin(), all_points.end(), sort_X_coords); //sort input vector by x-coordinates
	set_Pts(all_points); // saves the points as initialized to class
	set_numKDNodes(all_points.size());
	init_Height();

	initial_tree(all_points);
}

void KD_tree::initial_tree(std::vector<Tree_Point> all_points)
{

	if (num_kd_nodes == 1) //if there is only one coord in all_points
	{
		KD_Node * root = new KD_Node(pts_set[0]);
		root->set_Depth(get_treeHeight());
		root->set_Root(true);
		set_TreeRoot(root);


	}
	else if (num_kd_nodes == 2) //if there is two coords in all_points
	{
		KD_Node * root = new KD_Node(pts_set[1]);
		root->set_Depth(get_treeHeight() + 1);
		root->set_Type(INITIAL_CUT);
		root->set_Root(true);
		set_TreeRoot(root);

		KD_Node * left_Tree = new KD_Node(pts_set[0]);
		root->set_Depth(get_treeHeight() + 1);
		root->set_Root(true);
		set_TreeRoot(root);

		set_TreeHeight(compute_KD_Height());

	}
	else
	{
		// Initialize root node
		int median_idx = num_kd_nodes / 2;

		KD_Node * root = new KD_Node(all_points[median_idx]);
		root->set_Type(INITIAL_CUT);
		root->set_Depth(get_treeHeight());
		root->set_Root(true);
		set_TreeRoot(root);


		// making vectors for sorted list of positions
		std::vector<Tree_Point> x_left_arr, x_right_arr, y_left_arr, y_right_arr;

		// initialize by copying data from pts_set
		// for loops done separately because left/right arrays could
		// be of different sizes 
		// if the number of points is even, the right array is one 
		// index smaller than the left
		// Example:
		// 0 1 2 3 4, median index 2 --> left [0,1] right [3,4]
		// 0 1 2 3, median 2 --> left [0,1] right [3]

		for (int i = 0; i < median_idx; i++)
		{
			x_left_arr.push_back(pts_set[i]);
			y_left_arr.push_back(pts_set[i]);
		}
		for (int i = median_idx + 1; i < num_kd_nodes; i++)
		{
			x_right_arr.push_back(pts_set[i]);
			y_right_arr.push_back(pts_set[i]);
		}

		sort(x_left_arr.begin(), x_left_arr.end(), sort_X_coords);
		sort(x_right_arr.begin(), x_right_arr.end(), sort_X_coords);
		sort(y_left_arr.begin(), y_left_arr.end(), sort_Y_coords);
		sort(y_right_arr.begin(), y_right_arr.end(), sort_Y_coords);

		//builds tree with recursion
		//use the depth of the root (distance of the root node and its leaf) for
		//recursive call because the height is updated during the call
		root->set_Left(KD_Maker(x_left_arr, y_left_arr, root->Get_depth() + 1));
		root->set_Right(KD_Maker(x_right_arr, y_right_arr, root->Get_depth() + 1));

		//a sanity check
		set_TreeHeight(compute_KD_Height());
	}
	//clear points for a clean start
	//traverse tree and save the points in level order
	pts_set.clear();
	levelOrderPts();
}

/*
Using recursion to build K-D Tree; takes vectors of points sorted by X and Y also the current depth(height)
of the tree and returns a pointer to a node, which is stored as the left or right node of the parent.

Base:
If vector of points only contains one point then this point is a leaf, and the children of this node are NULL.
If vector of points is empty then all conditions will be checked and failed.

If and Else block:
If the vector only has one point (base case)
If the vector has more than one point and the depth is odd, then the cut type is VERTICAL.
If the vector has more than one point and the depth is even, then the cut type is HORIZONTAL.

Parameters:
Vectors of the same n points sorted by their x-coordinates, and sorted by their y-coordinates. The size of these is approximately half (-1 or +1) of the size of the vectors at the previous level.
The value of depth is the height of the AFTER the new nodes have been created, so that the height of the tree does not change if the vector is empty.
*/
KD_Node* KD_tree::KD_Maker(std::vector<Tree_Point> X_point, std::vector<Tree_Point> Y_point, int depth)
{
	int num = X_point.size();
	int middle_position = num / 2;

	//Base Case: if only one point, return leaf containing point If the height of the tree if ODD at the new node, then the cut type is VERTICAL,
	//The root the height is 1 and the initial cut is VERTICAL,
	//Thus if the height is even then the cut type is HORIZONTAL 
	
	if (num == 1)
	{
		KD_Node* node = new KD_Node(X_point[0]);
		node->set_Depth(depth);//default type is LEAF
		node->set_Left(NULL);//a sanity check
		node->set_Right(NULL);

		set_TreeHeight(depth);
		return node;
	}

	else if (depth % 2 != 0 && num > 1) //Height is odd at new node, with cut type VERTICAL
	{
		//Create new node and add median point to the tree
		KD_Node* odd_node = new KD_Node(X_point[middle_position]);
		odd_node->set_Type(VERTICAL);
		odd_node->set_Depth(depth);

		//Split x-coordinates into two parts and sort the parts by y-coordinates then next cut is HORIZONTAL
		std::vector<Tree_Point> x_left_odd, x_right_odd, y_left_odd, y_right_odd;

		//Copy data
		for (int i = 0; i < middle_position; i++)
		{
			x_left_odd.push_back(X_point[i]);
			y_left_odd.push_back(X_point[i]);
		}
		for (int i = middle_position + 1; i < num; i++)
		{
			x_right_odd.push_back(X_point[i]);
			y_right_odd.push_back(X_point[i]);
		}

		//Sort y-coordinates
		sort(y_left_odd.begin(), y_left_odd.end(), sort_Y_coords);
		sort(y_right_odd.begin(), y_right_odd.end(), sort_Y_coords);

		//Recursive calls
		odd_node->set_Left(KD_Maker(x_left_odd, y_left_odd, depth + 1));
		odd_node->set_Right(KD_Maker(x_right_odd, y_right_odd, depth + 1));

		set_TreeHeight(depth);
		return odd_node;
	}
	else if (depth % 2 == 0 && num > 1) //Depth is even and num > 1, cut is HORIZONTAL
	{
		//Add node to the tree
		KD_Node * even_node = new KD_Node(Y_point[middle_position]);
		even_node->set_Type(HORIZONTAL);
		even_node->set_Depth(depth);

		//Split x-coordinates into two parts and sort the parts by y-coordinates
		//LEFT refers to the part above the horizontal line, RIGHT refers to the points below it.
		std::vector<Tree_Point> x_left_even, x_right_even, y_left_even, y_right_even;
		x_right_even;
		for (int i = 0; i < middle_position; i++)
		{
			x_left_even.push_back(Y_point[i]);
			y_left_even.push_back(Y_point[i]);
		}
		for (int i = middle_position + 1; i < num; i++)
		{
			x_right_even.push_back(Y_point[i]);
			y_right_even.push_back(Y_point[i]);
		}

		//Sort x-coordinates
		sort(x_left_even.begin(), x_left_even.end(), sort_X_coords);
		sort(x_right_even.begin(), x_right_even.end(), sort_X_coords);

		//Recursive calls
		even_node->set_Left(KD_Maker(x_left_even, y_left_even, depth + 1));
		even_node->set_Right(KD_Maker(x_right_even, y_right_even, depth + 1));

		set_TreeHeight(depth);
		return even_node;
	}

	//ELSE
	// if the arrays are empty (ex: there is no right leaf for Node n), 
	// return NULL to set n->right to NULL
	return NULL;
}

/* Returns an int whose value is the height of the tree.
For checking correctness: HEIGHT = CEIL ( lg(NUM_NODES) )
traverses tree from root to leftmost leaf, increments a counter, and returns that value
If n is even, then the "left side" array is larger than the right – thus,
the longest path from root to leaf is the leftmost one.
Height is computed on the fly; this method was originally written to make sure that those
values were correct. It is redundant to leave it in, but it serves as a sanity check.
*/
int KD_tree::compute_KD_Height()
{
	KD_Node* temp_node = root;
	int height = 1;

	while (!temp_node->isleaf())
	{
		height++;
		temp_node = temp_node->Get_leftTree();
	}
	return height;
}

/* LEVEL ORDER POINTS
adds all nodes in KD-Tree to the pts vector using a simple
breadth first traversal
Iterates through the "levels" of the tree – calls addLevel to traverse
tree to find appropriate nodes.
*/

void KD_tree::levelOrderPts()
{
	for (int depth = 1; depth < height_Tree + 1; depth++)
	{
		add_Level(get_Root(), depth);
	}
}

/* ADD LEVEL adds all nodes at depth x to the vector, from left->right
called by levelOrderPts()
BASE CASES: node is NULL (parent only has one child node, or invalid input);
node is at desired depth, so add to the vector and exit
RECURSIVE CALLS: traverse subtrees and print all nodes of desired depth
*/

void KD_tree::add_Level(KD_Node* tree_node, int depth)
{
	if (!tree_node)
	{
		return;
	}
	if (tree_node->Get_depth() == depth)
	{
		ordered_pts_set.push_back(tree_node);
		return;
	}
	add_Level(tree_node->Get_leftTree(), depth);
	add_Level(tree_node->Get_rightTree(), depth);

}

void KD_tree::printInfo()
{

	printf("TREE INFO\n");
	printNumNodes();
	printf("Height: %d\n", get_treeHeight());
	//get_Points().at(0);
	printf("Root: ");

}

/* PRINTS NODES IN TREE -- LEVEL ORDER TREE TRAVERSAL?
*/
void KD_tree::printTree()
{
	if (!get_Root())
	{
		printf("Root is NULL; tree is empty.\n");
	}
	printf("\n\nPRINTING TREE:\n");
	for (int i = 0; i < ordered_pts_set.size(); i++)
	{
		ordered_pts_set[i]->print_Info();
		printf("\n");
	}

}

void KD_tree::printNumNodes() {
	printf("Number of nodes: %d\n", get_num_KDNodes());
}

/* Helper function for the destructor
Recurses through entire tree and call Node destructor to deallocate
(delete) every node in the tree.
*/
void KD_tree::deallocate_tree(KD_Node* tree_node)
{
	if (tree_node->isleaf() || !tree_node)
	{
		return;
	}
	deallocate_tree(tree_node->Get_leftTree());
	deallocate_tree(tree_node->Get_rightTree());
	delete(tree_node);
}


/* DESTRUCTOR
delete nodes from leaves to the root's left/right children,
then delete root
*/
KD_tree::~KD_tree()
{
	deallocate_tree(get_Root());
}

#pragma endregion

#endif _KDTREENODEMAKER_H
