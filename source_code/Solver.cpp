#include <vector>
#include <algorithm>
#include <numeric>
#include "time.h"

#include "ISA.h"

using namespace std;

// Initialize the solver
ISA::ISA(int dataDim, float sigma /* = 1 */, int numSample /* = 200 */, int iteration /* = 20 */)
{
	srand(time(NULL));

	m_dataDim = dataDim;
	m_numSample = numSample;
	m_iteration = iteration;

	m_dynamicFactor = 1.0f;
	m_annealFactor = .7f;

	m_samples = new float[m_dataDim*m_numSample];
	m_energy = new float[m_numSample];

	m_curValue = 100000;
	floor_Const = 0;
}

// Initialization of the solver
void ISA::initialization(float* x,int initSampleNum /* = 1 */, int iters /* = 20 */)
{
	if(iters>m_iteration||iters<=0) 
		iters = m_iteration;

	// 1) initialization
	if(initSampleNum<1) {
		printf("you should input one samples at least, i.e., initSampleNum should be large than 1!\n");
		return;
	}
	else if(initSampleNum==m_numSample)
		memcpy(m_samples,x,sizeof(float)*initSampleNum*m_dataDim);
	else {
		// randomly select the initialization
		for(int i=0;i<m_numSample;i++) {
			int idx = int(UnifoRand()*(initSampleNum-1));
			memcpy(m_samples+i*m_dataDim,x+idx*m_dataDim,sizeof(float)*m_dataDim);
		}
	}
}
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
//////////////////////////////////////////////////////////////////////
// Calculate the value of the objective function
float ISA::obj_Func(Rd_Set& rd_Set, Graph_Trans& gph, float* cur_Abun)
{
	float obj_Graph = cost_Graph(rd_Set,gph,cur_Abun);

	float* aver_Abun = new float[m_dataDim];
	float s = 0;
	for(int i = 0; i<m_dataDim; i++)
	{
		s += cur_Abun[i];
		aver_Abun[i] = cur_Abun[i];
	}
	for(int i = 0; i<m_dataDim; i++)
	{
		aver_Abun[i] /= (s+1e-12);
	}

	float rate = 0;
	float obj1 = - rd_Set.Px(aver_Abun)/(rd_Set.rd_Num);  // likelihood

	float obj = obj1;

	delete[] aver_Abun;

	float lambda1 = 1, lambda2 = 0;
	double obj_Joint = lambda1*obj + lambda2*obj_Graph;
	obj_Joint = obj_Joint - floor_Const;

	return obj_Joint;
}

//////////////////////////////////////////////////////////////////////
// Calculate the value of the objective function
float ISA::obj_FuncMI(Rd_Set& rd_Set, Graph_Trans& gph, float* cur_Abun)
{
	// abun_Infer(cur_Sol,temp_Abun);
	float obj_Graph = cost_Graph(rd_Set,gph,cur_Abun);
	// float obj1 = 0;

	float* aver_Abun = new float[m_dataDim];
	float s = 0;
	for(int i = 0; i<m_dataDim; i++)
	{
		s += cur_Abun[i];
		aver_Abun[i] = cur_Abun[i];
	}
	for(int i = 0; i<m_dataDim; i++)
		aver_Abun[i] /= (s+1e-12);

	float rate = 0;
	float obj1 = - rd_Set.Px(aver_Abun)/(rd_Set.rd_Num);
	float obj2 = - rd_Set.Iyx1(aver_Abun); 

	float obj = 1*obj2 + 0*obj1;

	delete[] aver_Abun;

	float lambda1 = 1, lambda2 = 0;
	double obj_Joint = lambda1*obj + lambda2*obj_Graph;

	obj_Joint = obj_Joint - floor_Const;

	return obj_Joint;
}

// Calculate objective function
float ISA::obj_FuncRe(Rd_Set& rd_Set, Graph_Trans& gph, float* cur_Abun)
{
	float obj_Graph = cost_Graph(rd_Set,gph,cur_Abun);
	// float obj1 = 0;
	
	float* aver_Abun = new float[m_dataDim];
	float s = 0;
	for(int i = 0; i<m_dataDim; i++)
	{
		s += cur_Abun[i];
		aver_Abun[i] = cur_Abun[i];
	}
	for(int i = 0; i<m_dataDim; i++)
	{
		aver_Abun[i] /= (s+1e-12);
		cur_Abun[i] /= (s+1e-12);
	}

	float rate = 0;
	float obj1 = - rd_Set.Px(aver_Abun)/(rd_Set.rd_Num);  // likelihood 
	float obj2 = 0;
	float obj3 = - rd_Set.Iyx1(aver_Abun);

	float gamma1 = 1, gamma2 = 0, gamma3 = 0, gamma_S = 0;
	float obj = gamma1*obj1 + gamma2*obj2 - gamma_S*s;   // joint objective of likelihood and mutual information

	delete[] aver_Abun;

	float lambda1 = 1, lambda2 = 0;
	double obj_Joint = lambda1*obj + lambda2*obj_Graph;
	obj_Joint = obj_Joint - floor_Const;

	return obj_Joint;
}

// Initial sampling of the abundances
void ISA::initial_Sampling(Rd_Set& rd_Set, float* ini_Sample, int initSampleNum)
{	
	std::tr1::mt19937 eng;       // a core engine class£ºMersenne Twister generator
	std::tr1::uniform_real<float>  unif1(0,1);    // respoinding to low abundance

	double b1 = 10, b2 = 20;
	std::tr1::uniform_real<float>  unif2(b1,b2);  // respoinding to high abundance
	eng.seed((unsigned int)time(NULL));

	int speci = 0, sum = 0;
	int N = m_numSample;

#pragma region generate initial abundances
	float* samples = new float[m_dataDim];
	float* ini_Abun = new float[m_dataDim];
	for(int k = 0; k<m_dataDim; k++)
	{
		int id = rd_Set.idx[k];
		float temp = rd_Set.path[id].assemble_Score;
		if(temp<0.5)
			ini_Abun[k] = 2*temp;
		else
			ini_Abun[k] = (temp-0.5)*2*(b2-b1) + b1;
		sum += ini_Abun[k];
	}
	for(int i = speci; i<N; i++)
	{	
		float sum = 0;
		for(int k = 0; k<m_dataDim; k++)
		{
			int id = rd_Set.idx[k];
			float temp = rd_Set.path[id].assemble_Score;
			if(temp<0.5)
			    samples[k] = (0.5*unif1(eng)+ini_Abun[k])/2;
			else
			    samples[k] = (0.5*unif2(eng)+ini_Abun[k])/2;
			sum += samples[k];
		}
		for(int k = 0; k<m_dataDim; k++)
			samples[k] /= sum;

		memcpy(ini_Sample+i*m_dataDim,samples,sizeof(float)*m_dataDim);
	}
#pragma endregion

	delete[] samples;
	delete[] ini_Abun;

}

// Output the samples
void ISA::output_Sample(char* filename, float* sample)
{
	ofstream outfile;
	outfile.open(filename,ios::out);

	for(int i = 0; i<m_numSample; i++)
	{
		for(int k = 0; k<m_dataDim; k++)
		{
			outfile<<sample[i*m_dataDim+k]<<" ";
		}
		outfile<<endl;
	}

	outfile.close();
}

// Choose sparse dimensions for optimization
void ISA::sparse_Normalize(float* abun, vector<int>& sel_Idx, vector<int>& nsel_Idx)
{
	int id_Num = sel_Idx.size();  // valid dimensions
	int n_Num = nsel_Idx.size();  // invalid dimensions

	double sum = 0;
	for(int j = 0; j<id_Num; j++)
	{
		int k = sel_Idx[j];  // number of selected dimensions
		sum += abun[k];
	}
	for(int j = 0; j<id_Num; j++)
	{
		int k = sel_Idx[j];  // number of selected dimensions
		abun[k] /= (sum+1e-09);  // normalize on the selected dimensions
	}
	for(int j = 0; j<n_Num; j++)
	{
		int k = nsel_Idx[j];
		abun[k] = 0;  // abundances on the other dimensions are set to zero
	}

	float sum1 = is_Normalized(abun,m_dataDim);
	if(abs(sum1-1)>1e-04)
		cout<<"abundance is not normalized! "<<sum1<<endl;
}

// Examine whether the array is normalized
float is_Normalized(float* abun, int data_Dim)
{
	float sum = 0;
	for(int i = 0; i<data_Dim; i++)
		sum += abun[i];

	return sum;
}