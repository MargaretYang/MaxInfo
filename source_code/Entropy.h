
#ifndef ENTROPY_H
#define ENTROPY_H

#include "Read_Assignment.h"
typedef vector<Exon_Node> VEC_NODE;

/////////////////////////////////////////////////////////
//// Read Set
class Rd_Set{
public:
	string chromosome_Name;  // chromosom name
	string gene_Name;        // gene name
	int locate_RdNum;        // number of aligned reads
	int rd_Num;              // number of reads
	int rd_Len;              // read length
	int exon_Num;            // number of exons
	float p1, p2;            // auxiliary variable
	r_Type** rd_Pro;         // generation probability
	r_Type** rd_Pro1;        // posterior probability
	vector<int> idx;         // index of valid transcripts
	vector<int> m_idx;       // index of valid transcripts after re-selection
	vector<int> sel_Idx;     // selected transcripts
	vector<vector<int>> rd_Left;  // left reads in each segment
	vector<vector<int>> rd_Right; // right reads in each segment
	float* fragLen_Pro;           // probility distribution of fragment length
	float* fragLen_Interval;      // probility distribution of fragment length by intervals
	int lim1;                     // lower limit of fragment length
	int lim2;                     // upper limit of fragment length
	int locate_start;             // starting site of a read sampling region
	int locate_stop;              // stopping site of a read sampling region
	vector<Rd> rd_Vec;            // vector of reads
	vector<Rd> singlerd_Vec;      // vector of single-end reads
	vector<int> annotated_Idx;    // index of annotated transcript

	double mean;                  // mean of fragment length
	double sigma;                 // variance of fragment length
	vector<Path> path;            // set of feasible paths
	double px_Likelihood;         // likelihood
	double cxy;                   // relative entropy
	float** junc_Mtx;             // cost of junction
	int junc_Num;                 // number of junctions
	char orientation;			  // orientation of splice site
	
public:

	// Destructor function
	~Rd_Set(){ 
		clear_RdVec();
	}

#pragma region File processing function
///////////////////////////////////////////////////////////////////////////////
////// File processing function
	// Load gene annotations
	void load_Annotation(char* file1, char* file2, char* ext);

	// Generate file of GTF format according to the annotations
	void gen_transBed(char* file1, char* file2);

	// Generate the abundances of transcripts
	void gen_Exply(char* file1, char* file2);

	// Load annotations of transcripts
	void load_GTF(char* file1, char* file2);

	// Determine the exons contained by each transcript
	void locate_transSub(vector<Gene>& gene, char* path, char* filename, char* ext);

	// Output the gene information to files
	void output_transBed(vector<Gene>& gene, char* filename);

	// Identify the index of the transcript
	int indexTrans(Gene& single_Gene, string trs_Name);

	// Load gene annotations
	void load_exonAnnotation(char* filename, char* path);
	
	// Output the gene annotations to files to provide information of exons
	void output_GeneAnnotation(char* filename, vector<Gene>& Gene_Set);

	// Select the genes according to their transcripts and subexons
	void transcript_Select(char* root_path, char* filename);

	// Downsampling of sequencing data
	void subSampling(char* filename1, char* filename2, int subFactor);
	
	// Merge the sequencing data
	void merge(char* filename1, char* filename2, vector<string>& gene_Name);

	// Evaluate the accuracy of exon assembly
	float exon_assemblyEvaluation(char* filename1, char* filename2, float& precision, float& recall);

	// Separted original sequencing data into the form of paired-end reads
	void separate_pairedData(char* filename1, char* filename2, char* filename3);

	// Merge the data of paired-end reads into one file
	void Merge_pairedData(char* filename1, char* filename2, char* filenam3);

	// Extract the read data in the specific gene
	void extract_pairedData(char* filename1, char* filename2, int site1, int site2, string chrom_Name);

///////////////////////////////////////////////////////////////////////////////
#pragma endregion

#pragma region Result Analysis
///////////////////////////////////////////////////////////////////////////////
	// Analysis of the isoform identification results
	// Input£ºfilename1 annotated transcripts; filename2 transcritp abundances; filename3 reconstrcuted transcripts; filename4 comparison results; gene_Name
	void quantification_Compare(char* filename1, char* filename2, char* filename3, char* filename4, string gene_Name);
#pragma endregion

	// Load the results of read alignment
	void load_Assignment(char* filename, Graph_Trans& gph, int exon_Num);
	// Load read information
	void load_Read(char* filename);
	// Load transcript information
	void load_Trans(char* filename);

	// Estimate the probability distribution of fragment length based on the sequencing data using discrete intervals
	bool fragLen_Estimation(int& read_Number, float len_Mean, float len_Var=DEFAULT_VAR);
	int  fragLen_Estimation_Initial(vector<int>& read_Length);

	// Estimate the probability distribution of fragment length
	void fragLen_Dis(double mean, double sigma, int flim1, int flim2);

	// Estimate the probability distribution of fragment length based on the sequencing data using discrete intervals
	void fragLen_DisInterval(double mean, double sigma, int flim1, int film2);

	// Calculate the cost of the first node and the last node
	void cost_Linkage(Graph_Trans& gph);

	// Build cost matrix for the junctions
	void junc_MtxBuild_1(Graph_Trans& gph);

	// Output the read information to the file
	void output_Read(char* filename);

	// Select possible transcripts after calculating P(x_i|y_k)
	void select_ValidTrans();

	// Select possible transcripts
	void select_ValidTransScore();

	// Select transcripts when annotations are available
	void select_ValidTransAnnoted(Graph_Trans& gph);

	// Output the transcript information to the file
	void output_Trans(char* filename);

	// Print transcript information
	void print_Trans();

	// Initialize the generation probability of read (p(x_i|t_k)
	void initial_RdPro(int trans_Num);

	// Connectivity cost of the path, used to constrain the first exon and the last exon
	void connect_Weight(Graph_Trans& gph);

	// Clear the vectors
	void clear();

	// Clear the generation probability array of the reads
	void clear_RdPro(int rd_num);

	// Clear the vector
	void clear_RdVec();

	// Calcuate the probablity that the read is generated from a given transcript
	void path_Entro_Locate(Graph_Trans& gph, Rd& cur_Rd, r_Type* rd_pro);

	// Calcuate the entropy of a single read
	void path_Entro1(Graph_Trans& gph, string filename);

	// Calculate the generation probability of the read
	double rd_Entro(Rd &rd, Path &path, vector<Exon_Node> &exon_Vec, int id1, int id2, int id3, int id4);

	// Normalize p(x|y)
	void rd_Normalize();

	// EM algorithm: calculate the expectation of transcript abundance
	float abun_Estimation(float* abun_Est, float* abun, int interval);

	// EM algorithm: calculate the expectation of transcript abundance
	void abun_Estimation_Re(float* abun_Est, float* abun_Est_Re, float* abun, int interval);

	// Calculate H(X)
	double Hx1(float* abundance);
	double Hx1(float* abundance, vector<int>& sel);
	
	// Calculate P(X)
	double Px(float* abundance, int interval = 1);

	// Calculate H(Y|X)
	double Hy(float* abundance);

	// Calculate H(Y|X)
	double Hyx(float* abundance, int interval = 1);

	// Calculate I(X,Y) = H(X) - H(X|Y)
	double Ixy(float* abundance);
	double Ixy1(float* abun, vector<int>& sel, vector<int>& nsel);
	
	// Calculate I(Y,X) = H(Y) - H(Y|X)
	double Iyx1(float* abundance, int interval = 1);

	// Calculate I(Y,X)=H(Y)-H(Y|X)
	double Iyx2(float* abundance);

private:
	float UnifoRand(){ return float(rand())/float(RAND_MAX);};

};


#endif