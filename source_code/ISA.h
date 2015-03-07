#ifndef _ISA_H_
#define _ISA_H_

//#include <Eigen/Eigen>
//using namespace Eigen;

#include "Read_Assignment.h"
#include "Entropy.h"

class ISA
{
public:
	ISA() : m_samples(nullptr), m_energy(nullptr) {;};
	ISA(int dataDim, float sigma = 1, int numSample = 200, int iteration = 20);
	~ISA()
	{
		if(m_samples) delete m_samples;
		if(m_energy) delete m_energy;
	}

public:
	void    initialization(float* x, int initSampleNum /* = 1 */, int iters /* = 20 */);

	// initial solution od the distribution of abundance
	void    initial_Sampling(Rd_Set& rd_Set, float* ini_Sample, int initSampleNum);
	
	// parameter estimation by EM algorithm
	float   Optimize_EM(float* abun, Rd_Set& rd_Set, Graph_Trans& gph, vector<int>& sel_Idx, vector<int>& nsel_Idx, int iter1=20, int iter_Ctr=7);

	// parameter estimation by Variational Bayesian
	void    Optimize_VB(Rd_Set& rd_Set, Graph_Trans& gph, vector<int>& sel_Idx, vector<int>& nsel_Idx, int iter1=20, int iter_Ctr=7);

	// parameter estimation by MDL
	void    Optimize_MDL(Rd_Set& rd_Set, Graph_Trans& gph, vector<int>& sel_Idx, vector<int>& nsel_Idx, int iter1=20, int iter_Ctr=7);

	// EM algorithm with mutual information constraints
	float   Optimize_EM_Re(float* abun, Rd_Set& rd_Set, Graph_Trans& gph, vector<int>& sel_Idx, vector<int>& nsel_Idx, int iter_Ctr=7, float gamma=0.01, int interval=1);

	// interator control
	bool	Optimize_IteCtr(float* x, Rd_Set& rd_Set, Graph_Trans& gph, char* root_path, int initSampleNum=1, int iters=20);
    
	// output optimum
	void    output_Solution(float* x);

	// output sample
	void	output_Sample(char* filename, float* sample);

	// calculate objective function
	float   obj_Func(Rd_Set& rd_Set, Graph_Trans& gph, float* sol);

	// calculate objective function
	float   obj_FuncMI(Rd_Set& rd_Set, Graph_Trans& gph, float* sol);

	// calculate objective function with sparse constraints
	float   obj_FuncRe(Rd_Set& rd_Set, Graph_Trans& gph, float* sol);

	// calculate objective function with mutual information and graph cost constraints
	float   cost_Graph(Rd_Set& rd_Set, Graph_Trans& gph, float* abun);

	// graph cost 1: junction cost
	float   cost_Junction(Rd_Set& rd_Set, Graph_Trans& gph, float* abun);

	// graph cost 2: connection cost
	float   cost_Linkage(Rd_Set& rd_Set, Graph_Trans& gph, float* abun);

	// swap some parts of two solution
	// input£ºindex£¬start point of swap£¬length of swap
	void    exchange_Intra(int k1, int k2, int i1, int i2, int len);
	
	// swap two dimensions in one solution
	// input£ºstart point of swap
	void    exchange_Inter(int k, int i1, int i2);
	// exchange solution
	void    exchange();

	// select path
	void    path_SFMS(float* abun, Rd_Set& rd_Set, Graph_Trans& gph, vector<int>& sel, vector<int>& nsel, vector<PAIR>& pair_Vec, int select_pathNum);  // forward selection of path
	void    path_SFMS_Auto(float* abun, Rd_Set& rd_Set, Graph_Trans& gph, vector<int>& sel, vector<int>& nsel, vector<PAIR>& pair_Vec, int select_pathNum);  // forward selection of path
	void    path_SFMS_Auto1(float* abun, Rd_Set& rd_Set, Graph_Trans& gph, vector<int>& sel, vector<int>& nsel, vector<PAIR>& pair_Vec, int interval);  // forward selection of path
	bool	path_Compare_Annotated(float* abun, Rd_Set& rd_Set, Graph_Trans& gph, vector<int>& sel, vector<int>& nsel); // add annotated isoforms
	bool	path_Limit_Annotated(float* abun, Rd_Set& rd_Set, Graph_Trans& gph, vector<int>& sel, vector<int>& nsel);	// limited to annotated isoforms
	
	// check whether the paths overlap
	float	isSubset(Path& path1, Path& path2);
	float	isoform_Similarity(Rd_Set& rd_Set, vector<int>& sel, int sel_idx);

	void	setMode(int mode){method_Mode = mode;};
	
private:
	// sparse sample
	void    rank_Dim(vector<PAIR>& pair_Vec, float* sol);
	float	UnifoRand() { return float(rand())/float(RAND_MAX);};
	void    sparse_Normalize(float* abun, vector<int>& sel_Idx, vector<int>& nsel_Idx);  // normalize on the sparse array

	// output abundance
	void    report(Rd_Set& rd_Set, Graph_Trans& gph, char* filename, vector<int>& sel_Idx);	
  
private:
	int         trans_Num;   // number of isoform
	int			m_dataDim;   // dimension of valid isoform
	int			m_numSample;
	int			m_iteration;
	float		m_dynamicFactor;
	float		m_annealFactor;
	float       m_curValue;  // optimum
	int			method_Mode; // mode

	//related to the parameters
	float       floor_Const; // offset

	float*		m_samples;
	float*		m_energy;
	float*      m_weights;   // weight of solution
	float*      m_CurOpt;    // current optimum
	float*      m_CurAver;   // average of optimum

	vector<int>  v_Idx;      // selected dimension
	vector<int>  n_Idx;      // unselected dimension
}; 

// check whether the array is normalized
float is_Normalized(float* abun, int data_Dim);

#endif
