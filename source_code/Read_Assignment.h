
#ifndef READ_ALIGNMENT_H
#define READ_ALIGNMENT_H

#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <list>
#include <stack>
#include <string>
#include <iterator>
#include <algorithm>
#include <functional>
#include <ctime>
#include <time.h>
#include <random>
#include <assert.h>
//#include "Entropy.h"

using namespace std;

#define _max(x, y) x>y?x:y 
#define _min(x, y) x<y?x:y  
#define Lim1 75
#define Lim2 300
#define r_Type double
#define m_floor 0.000000001
#define MAX_LINE 500
#define DEFAULT_MEAN 120
#define DEFAULT_VAR 30
typedef unsigned char uchar;
typedef pair<int, float> PAIR;
typedef pair<int, int> PAIR_INT;
typedef vector<int> VEC_INT;

// Read
struct Rd{
	string name;					// read name
	int Left1, Right2;				// left end of left read, right end of right read
	int Left2, Right1;				// right end of left read, left end of right read
	int insert_Len;					// fragment length
	short l1,l2;					// segment the left read falls in
	short r1,r2;					// segment the left read falls in
	float match_Score;				// match score of the read
	vector<int> left_readIdx;		// exons in which the boundaries of left read fall
	vector<int> right_readIdx;		// exons in which the boundaries of left read fall
	vector<vector<int>> bi_pos;		// splitting points of the exons which may fall in the introns
	vector<int> left_pos;			// possible splitting points of left read
	vector<int> right_pos;			// possible splitting points of right read
	
};

// Interval for identifying exons
struct Interval{
	vector<vector<int>> site_New1, site_New2;		// starting site and stopping site of the exon
	vector<vector<int>> cnt_New1, cnt_New2;			// number of reads in the segment
	vector<vector<int>> gap1, gap2;					// distance of splice sites
	vector<int> serial;								// segment type
	vector<int> ch_type;							// exon type
};

// Path(possible transcript)
struct Path{
	int length;						// length of path
	int exon_Number;				// number of exons
	int* exon_Node;					// exons contained in the path
	int* bound;						// boundaries of the exons
	double entro_proAccu;			// accumulated probabilities for calculating H(X|Y)
	double entro_Abun1;				// entropy, for calculating H(X|Y)
	double entro_Abun2;				// entropy, for calculating H(Y|X)
	double assemble_Score;			// generation probability of the reads by the given set of isoforms
	double link_Cost1;				// cost of the first and last nodes
	double link_Cost2;				// cost of path length
    double entro_Den;				// sum of the generation probabilities
	int max_Id1, max_Id2;			// auxiliary indices
	int min_Id1, min_Id2;			// auxiliary indices
	double connect_FirstWeight;		// connectivity cost to constrain the first exon
	double connect_LastWeight;		// connectivity cost to constrain the last exon
};

// Exon
struct Exon_Node{
	int bound1, bound2;				// left boundary and right boundary of exon
	int len;						// length of exon
};

// Transcript
struct Trans{
	string trsName;					// transcript name
	double abundance;				// abundance of transcript
	double FPKM;					// abundance measured in FPKM
	double FPKM1;					// abundance measured in FPKM
	int length;						// length of transcript
	int start, stop;				// starting site and stopping site of transcript
	char orientation;				// orientation of splice site
	vector<int> subexon_Serial;		// serial of subexon contained in the transcript
	vector<bool> cds;				// flag indicating the coding region
	vector<int> exon_len;			// exon length
	vector<int> exon_start;			// starting site of exon
	vector<int> exon_stop;			// stopping site of exon
};


// Gene Structure
struct Gene{
	string gene_Name;				// gene name
	int start, stop;				// starting site and stopping site of gene
	int subexon_Num;				// number of exons
	char orientation;				// orientation of splice site
	vector<Exon_Node> subexon;		// subexons of gene
	vector<int> exon_Type;			// exon type
	vector<Trans> trans;			// transcripts of the gene
	vector<PAIR_INT> non_Link;      // exons without linkage
};


////////////////////////////////////////////////////////
//// Structure of Graph
class Graph_Trans{
public:
	int vtx_Number;							// number of vertex(node)
	int cvtx_Number;						// number of vertex including the introns
	float** link_Mtx;						// linkage matrix of the directed graph
	float** clink_Mtx;						// linkage matrix of the directed graph including the introns as vertex
	vector<int> junc_Link1;					// linkage from the junctions
	vector<int> junc_Link2;					// linkage from the junctions
	vector<Exon_Node> exon_Vec;				// vector of exons
	vector<int> exon_Type;					// exon type
	vector<vector<int>> paths;				// possible paths
	vector<int> in_Degree;					// in-degree of node
	vector<int> out_Degree;					// out-deree of node
	vector<VEC_INT> pre_Nodes;				// previous node
	vector<VEC_INT> pos_Nodes;				// next node
	vector<vector<bool>> pre_Linked;		// flag for linkage with the previous node
	vector<vector<bool>> pos_Linked;		// flag for linkage with the next node
	vector<float> in_Flow;					// in-flows of node
	vector<float> out_Flow;					// out-flow of node
	double aver_in;							// average in-flow
	double aver_out;						// average out-flow
	Gene cur_Gene;							// current considered gene
	vector<int> annotated_Idx;				// index of annotated transcipt

	int junc_Num;							// number of junctions
	int** junc_Rank;						// matrix of ranks of junctions
	double* first_NodeCost;					// cost for being the first node in the path
	double* last_NodeCost;					// cost for being the first node in the path
	vector<vector<Exon_Node>> intron_Add;	// vector of introns
	vector<vector<int>> serial_Add;			// serials of introns
	int add_Num;							// number of added exons
	vector<PAIR_INT> non_Link;				// exons without linkage
		
public:
	~Graph_Trans();			// Destructor function
	void clear();			// Reset function
    void clear_JuncRank();	// Clear the linkage
	void load_Graph(char* filename, int vex_Num);		// Load the graph structure	
	void load_Trans(char* filename);					// Load transcript information
	void graph_Initial(int exon_Num);					// Build directed graph containing both exons and introns based on the read data
	void graph_Initial_1(int exon_Num);					// Build directed graph based on the read data
	void graph_Initial(vector<int>& add, int add_Num);	// Build directed graph containing both exons and introns based on the read dat
	void graph_Build(int exon_Num);						// Build directed graph
	void junc_Index();									// Build indices of junctions
	void cost_FirstLastNode();							// Calculate the cost of the exon for being the first node or the last node of the path

	// Modify the linkage in the matrix and clear low-confidence linkage
	void graph_Modify();

	// Constraints on the paths with the costs of the first node and the last node
	bool graph_Constraint_Initial();

	// Constraints on the paths when annotations are available
	bool graph_Constraint_Annotated();

	// Enhance the constraints for being the first node or the last node 
	bool graph_Constraint();
	
	// Enhance the constraints for being the first node or the last node 
	bool graph_Constraint_1(vector<float>& in_Flow, vector<float>& out_Flow);

    // Search path by DFS(Depth First Search)
	bool path_DFS(bool* mark);

    // Search path by DFS(Depth First Search): Two-Stack method
	bool path_DFS1();

	//  Path search with constraints when there are too many possible paths
	bool path_DFS1_constraint();

	// Path search with adaptive constraints
	bool path_DFS1_Adaptive();

    // Obtain information of the paths including path structures and path entropies
	void path_Speci(vector<Path>& path);
	
	// Output the graph structure
	void print_Graph();
	
	// Output the paths
	void print_Path();
};

//////////////////////////////////////////////////////////
//// Public function
// Compare the values of elements with the type PAIR
int cmp(PAIR& x, PAIR& y);
//  Compare the values of elements with the type PAIR_INT
int cmp_1(PAIR_INT& x, PAIR_INT& y);

//  Load the gene annotations
void load_GTF(char* path, char* filename, char* ext);

// Convert to valid filename
string valid_filename(string filename);

// Check whether the given gene exists
int isExisting(vector<Gene>& gene, string gene_Name);

// Check whether the given transcript exists
int isExisting_Trans(vector<Trans>& trs, string trs_Name);

#endif