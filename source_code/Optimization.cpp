#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include "time.h"

#include "ISA.h"


// EM-algorithm with mutual information constraints
float ISA::Optimize_EM_Re(float* abun, Rd_Set& rd_Set, Graph_Trans& gph, vector<int>& sel_Idx, vector<int>& nsel_Idx, int iter_Ctr, float gamma, int interval)
{
	// normalize sparse array

	float* abun_Est = new float[m_dataDim];   // expectation of abundance
	float* abun_Est_Reg = new float[m_dataDim];   // entropy of abundance

	float L1 = 0, L2 = 1000;     // function value
	int rd_Num = rd_Set.rd_Num;  // number of samples

	int iter = 0;
	int iter_Limit = 5000;
	float thresh = 1e-03;     // convergence
	float lamda = 1;  // mutual information constraints
	float obj = 0, obj_pre = 0;
	int id_Num = sel_Idx.size();
	rd_Set.sel_Idx = sel_Idx;  // selected path

	//while(abs(1-L1/(L2+1e-09))>thresh&&iter<iter_Limit)
	float delta = 1000;
	float* abun_pre = new float[m_dataDim];

	//int interval = (int)rd_Set.rd_Num*1.0/sub_num;
	//int interval = 1;

	// while(abs(1-obj_pre/(obj+1e-09))>thresh&&iter<iter_Limit)
	while(delta>thresh&&iter<iter_Limit)
	{

		obj_pre = obj;
		memcpy(abun_pre,abun,sizeof(float)*m_dataDim); // last step solution

		for(int j = 0; j<m_dataDim; j++){
			abun_Est[j] = 0;
			abun_Est_Reg[j] = 0;
		}

		// E step
		rd_Set.abun_Estimation_Re(abun_Est,abun_Est_Reg,abun,interval);  // expectation of abundance
		
		// M step
		int num = sel_Idx.size();   // number of valid path
		float den = 0;
		for(int j = 0; j<num; j++){
			int k = sel_Idx[j];     // index of valid path
			abun[k] = lamda*abun_Est[k] + gamma*abun_Est_Reg[k];
			den += abun[k];
		}
		// calculate abundance
		if(abs(den)>0)
		{
			for(int j = 0; j<num; j++)
			{
				int k = sel_Idx[j];    // index of valid path
				abun[k] /= den;
			}
		}
		else
			cout<<"regularized EM algorithm: denominator is zero!"<<endl;

		/*for(int j=0;j<m_dataDim;j++)
		printf("%.4f ", abun[j]);
		printf("\n");*/

		L2 = rd_Set.Px(abun,interval);    // calculate function value
		float obj2 = - rd_Set.Hyx(abun,interval);  // entropy
		obj = lamda*L2 + gamma*obj2;     // objective function

		for(int j = 0; j<m_dataDim; j++)
			abun_Est[j] = 0;   // reset to 0

		delta = 0;
		for(int j = 0; j<m_dataDim; j++)
			delta += (abun_pre[j]-abun[j])*(abun_pre[j]-abun[j]);
		delta = sqrt(delta);
		
		// cout<<iter<<" "<<L2<<" "<<L2/rd_Num<<" "<<obj<<endl;
		iter++;
	} // end while

	ofstream outfile;
	// outfile.open("D:\\Splicing Graphensembl\\data-sdhd\\v4\\em_Sol.txt",ios::out);
	outfile.open("G:\\My documents\\FragmentStreaming\\D Melanogaster\\chr3R\\em_Sol.txt",ios::out);

	memcpy(m_CurOpt,abun,sizeof(float)*m_dataDim);  // copy to optimum

	// write to file
	// cout<<"likelihood: "<<iter<<" "<<L2<<endl;
	for(int k = 0; k<m_dataDim; k++)
	{
		outfile<<m_CurOpt[k]<<" ";    // save solution
		
	}

	//printf("\n");
	outfile<<endl;	
	outfile.close();

	return obj;
}

// parameter estimation by EM algorithm
// sel_Idx£ºselected isoform£¬
float ISA::Optimize_EM(float* abun, Rd_Set& rd_Set, Graph_Trans& gph, vector<int>& sel_Idx, vector<int>& nsel_Idx, int iter1, int iter_Ctr)
{
	float* abun_Est = new float[m_dataDim];   // expectation of abundance

	float L1 = 0, L2 = 1000;  // function value
	float thresh = 1e-05;     // convergence
	int rd_Num = rd_Set.rd_Num;   // number of samples

	int iter = 0;
	int iter_Limit = iter1;
	int id_Num = sel_Idx.size();
	rd_Set.sel_Idx = sel_Idx;  // selected path
	int interval = 1;

	while(abs(1-L1/(L2+1e-09))>thresh&&iter<iter_Limit)
	{
		L1 = L2;  // last function value
		for(int j = 0; j<m_dataDim; j++)
			abun_Est[j] = 0;   // reset to 0

		// E step
		float valid_den = rd_Set.abun_Estimation(abun_Est,abun,interval);  // expectation of abundance

		// M step
		int num = sel_Idx.size();  // number of valid path
		float sum = 0;
		for(int j = 0; j<num; j++)
		{
			int k = sel_Idx[j];    // index of valid path
			abun[k] = abun_Est[k]/valid_den;
			sum += abun[k];
		}

		/*for(int j=0;j<m_dataDim;j++)
		printf("%.4f ", abun[j]);
		printf("\n");*/

		L2 = rd_Set.Px(abun,interval)/rd_Num;    // calculate function value

		/*for(int j=0;j<m_dataDim;j++)
		printf("%.4f ", abun[j]);
		printf("\n");

		cout<<iter<<" "<<L2*rd_Num<<" "<<L2<<endl;*/
		iter++;

	} // end while

	ofstream outfile;
	// outfile.open("D:\\Splicing Graphensembl\\data-sdhd\\v4\\em_Sol.txt",ios::out);
	outfile.open("D:\\Splicing Graph\\Simulator\\em_Sol.txt",ios::out);


	for(int k = 0; k<m_dataDim; k++)
	{
		outfile<<abun[k]<<" ";

	}
	outfile<<endl;
	// printf("\n");
	outfile.close();
	// cout<<"iteration: "<<iter<<endl;

	return L2;  // return function value
}

// optimization and iteration control
bool ISA::Optimize_IteCtr(float* x, Rd_Set& rd_Set, Graph_Trans& gph, char* root_path, int initSampleNum /* = 1 */, int iters /* = 20 */)
{
	m_CurOpt = new float[m_dataDim];   // current optimum
	int sub_num = 1000;   // number of subseet of samples
	int interval = _max(1,(int)rd_Set.rd_Num*1.0/sub_num);  // sample interval
	int thresh = 5000;

	char filename[200];  
	
	vector<int> sel_Idx;
	initialization(x,initSampleNum,iters); // initialization of solution

	bool ctr_Flag = true;

	if(m_dataDim==1)  // 1 valid isoform
	{
		sel_Idx.clear();
		sel_Idx.push_back(0);
		m_CurOpt[0] = 1;   // select the isoform
 	}
	else if(m_dataDim<5000){
	int test_Run = 3;
	if(m_dataDim>500)   // if the number of paths is too large, run once
		test_Run = 1;

	vector<float> mI(test_Run);         // mutual information
	vector<float> con_mI(test_Run);     // entropy of mutual information
	vector<float> mI_Ixy(test_Run);     // calculate mutual information
	vector<float> likelihood(test_Run); // likelihood function
	vector<float> re_ml(test_Run);      // likelihood function with mutual information constraints
	
	rd_Set.sel_Idx = rd_Set.idx;     // initialize sel_Idx
	floor_Const = 0;

	int iter = 0;
	int iter1 = 50, iter2 = 20, iter3 = 500;
	vector<int> sel_tIdx, nsel_Idx, nsel_tIdx;

	int iter_Ctr = 50;
	float function_Value = 0;
	float* abun = new float[m_dataDim];   // estimated optimum abundance
	m_CurAver = new float[m_dataDim];     // average
	float* cnt = new float[m_dataDim];    // count the frequency of each isoform
	for(int i = 0; i<m_dataDim; i++) 
		cnt[i] = 0;

	float fun_Max = -1000000;        // optimum function value
	vector<int> sel_Opt, nsel_Opt;   // selected path of the optimum
	float* sol_Opt = new float[m_dataDim];  // optimum
	float* m_abun = new float[test_Run];  // maximum abundance

	clock_t start_Time, finish_Time, start_Time1, finish_Time1;  // runtime
	double cal_Time = 0;

	start_Time1 = clock();  // timing

	float sum_re = 0;   // normalization constant

	for(int iter = 0; iter<test_Run; iter++)
	{
		memcpy(abun,m_samples+iter*m_dataDim,sizeof(float)*m_dataDim);  // copy the initial solution to current solution

		//cout<<"The first pass:"<<endl;
#pragma region first step iteration
#pragma region select path

		sel_Idx.clear();  
		nsel_Idx.clear();
		rd_Set.sel_Idx.clear();
		vector<int> idx = rd_Set.idx;     // select all possible isoform
		int id_Num = rd_Set.idx.size();   // number of valid path
		int v_Num = idx.size();   // number of valid path

		bool* flag = new bool[m_dataDim];
		for(int i = 0; i<m_dataDim; i++)
			flag[i] = true;

		for(int i = 0; i<m_dataDim; i++){
			if(flag[i]){
				sel_Idx.push_back(i);    // add selected isoform
				sel_tIdx.push_back(idx[i]);
			}
			else
			{
				nsel_tIdx.push_back(idx[i]);
				nsel_Idx.push_back(i);
			}
		}
		rd_Set.sel_Idx = sel_Idx;

#pragma endregion 

		// start_Time = clock();

		sparse_Normalize(abun,sel_Idx, nsel_Idx);  // normalize the sparse array

		float gamma = 0.006;  // adjust parameter
		function_Value = Optimize_EM_Re(abun, rd_Set, gph, sel_Idx, nsel_Idx, iter_Ctr, gamma, interval);   //  parameter estimation by EM algorithm with constraints

		memcpy(m_samples+iter*m_dataDim,abun,sizeof(float)*m_dataDim);  // copy to optimum
	
		vector<PAIR> pair_Vec1;
		rank_Dim(pair_Vec1,abun);   // sort the isoform by abundance
#pragma region write abundance to files
		int num1 = pair_Vec1.size();     // number of isoform
		ofstream pair_file;

        for(int k = 0; k<num1; k++)
		{
			int k1 = pair_Vec1[k].first;
			int index = rd_Set.idx[k1];
			float link_Cost = rd_Set.path[index].link_Cost1;
			// pair_file<<pair_Vec1[k].first<<" "<<index<<" "<<pair_Vec1[k].second<<" "<<link_Cost<<endl;
		}
#pragma endregion write abundance to files

		mI[iter] = rd_Set.Hyx(abun,interval);     // H(Y|X) of optimum
		con_mI[iter] = rd_Set.Iyx1(abun,interval);  // I(X,Y) of optimum
		likelihood[iter] = rd_Set.Px(abun,interval)/rd_Set.rd_Num;  // function value of optimum
		re_ml[iter] = function_Value;      // function value of optimum with constraints
		
		m_abun[iter] = 0;
		for(int i = 0; i<m_dataDim; i++){
			if(abun[i]>m_abun[iter])
				m_abun[iter] = abun[i];
		}
		
	    sum_re += re_ml[iter];
#pragma endregion
	}

#pragma region second step optimization

	float maxabun = 0;
	float* rate = new float[test_Run];
	for(int k = 0; k<test_Run; k++)
		rate[k] = re_ml[k]/sum_re;   // weight of each solution

	for(int k = 0; k<test_Run; k++){
		//m_abun[k] = 0;
		for(int i = 0; i<m_dataDim; i++)
		{
			float temp_abun = m_samples[k*m_dataDim+i];
			/*if(temp_abun>m_abun[k])
				m_abun[k] = temp_abun;*/
			cnt[i] += rate[k]*temp_abun;
		}
		maxabun += rate[k]*m_abun[k];
	}

	delete[] rate;

	memcpy(abun,cnt,sizeof(float)*m_dataDim);
	vector<PAIR> pair_Vec;
	rank_Dim(pair_Vec,abun);   // sort the isoform by abundance

	vector<int> sel, nsel;   // candidate path , unselected path
	start_Time = clock();
	int K = 3;  // number of selected path
	// path_SFMS(abun, rd_Set, gph, sel, nsel, pair_Vec,K);		
	path_SFMS_Auto1(abun, rd_Set, gph, sel, nsel, pair_Vec, interval);

	// compare with annotated isoform and select path
	if(method_Mode==1)
		path_Compare_Annotated(abun, rd_Set, gph, sel, nsel);
	else if(method_Mode==2)
		path_Limit_Annotated(abun, rd_Set, gph, sel, nsel);
	
	finish_Time = clock();
	cal_Time = 1000*double(finish_Time - start_Time)/CLOCKS_PER_SEC;  // runtime
	// printf("Running time: %.4f ms\n",cal_Time);		

	rd_Set.sel_Idx = sel;   // reset rd_Set index
	sel_Idx = sel;   // selected path by forward selection
	nsel_Idx = nsel; // unselected path by forward selection

	int sel_trsNum = sel_Idx.size();


	//start_Time = clock();
	sparse_Normalize(abun, sel_Idx, nsel_Idx);  // normalize the sparse array
	function_Value = Optimize_EM(abun,rd_Set, gph, sel_Idx, nsel_Idx, iter3, iter_Ctr);  // optimization by EM algorithm

	memcpy(sol_Opt,abun,sizeof(float)*m_dataDim);

#pragma endregion second step optimization
	
	finish_Time1 = clock();  
	cal_Time = finish_Time1 - start_Time1;
	/*printf("Running time: %.4f ms\n",cal_Time);*/

	memcpy(m_CurOpt,sol_Opt,sizeof(float)*m_dataDim);  // optimum
	delete[] abun;  // release array
	}
	else
	{
		//cout<<"Too many paths: "<<m_dataDim<<endl;
		ctr_Flag = false;
		return ctr_Flag;
    }

	sprintf(filename,"%s\\exply.txt",root_path);
	report(rd_Set, gph, filename, sel_Idx);  // output abundance
	return ctr_Flag;
}

// parameter estimation by MDL
void ISA::Optimize_MDL(Rd_Set& rd_Set, Graph_Trans& gph, vector<int>& sel_Idx, vector<int>& nsel_Idx, int iter1, int iter_Ctr)
{
	int knz = rd_Set.idx.size();  // number of valid isoform
	int kmin = 4;    // minimum number of valid isoform
	int kmax = knz;	
	float Lmin = 10000;     // optimum function value
	int iter = 0;    // iterater count
	int rd_Num = rd_Set.rd_Num;   // number of isoform

	float* abun = new float[kmax];  // dimension of the solution
	float* abun_Opt = new float[kmax];  // optimum

	float L1 = 0, L2 = 0, L_min = 10000;

	while(knz>kmin)
	{
		do{
			L1 = L2;
			float abun_Sum = 0;
			for(int i = 0; i<kmax; i++)
			{
				float* alpha = new float[kmax];  // estimated probability on each dimention

#pragma region calculate post probability
				for(int j = 0; j<rd_Num; j++)
				{
					float sum = 0;
					for(int l = 0; l<kmax; l++)
					{						
						int m = sel_Idx[l];  // selected isoform
						rd_Set.rd_Pro1[j][m] = rd_Set.rd_Pro[j][m]*abun[m];
						sum += rd_Set.rd_Pro1[j][m];
					}

					for(int l = 0; l<kmax; l++)
					{
						int m = sel_Idx[l];  // selected isoform
						rd_Set.rd_Pro1[j][m] /= sum;       // probability normalization
						alpha[l] += rd_Set.rd_Pro1[j][m];  // increase the probability
					}
				}
#pragma endregion calculate post probability

#pragma region probability estimatation of each model
				float lambda  = 0.5;     // parameter of sparse array
				int N = 1;
				float thresh1 = N*lambda;   // N for number of parameters
				float sum = 0;
				for(int l = 0; l<kmax; i++)
				{
					abun[l] = _max(0,alpha[l]-thresh1);  // set the negtive value to 0
					sum += alpha[l];
				}
				for(int l = 0; l<kmax; l++)
					abun[l] = abun[l]/sum;

#pragma endregion probability estimatation of each model

				if(abun[i]>0)  // update parameter
				{
#pragma region calculate post probability
					for(int j = 0; j<rd_Num; j++)
					{
						for(int l = 0; l<kmax; l++)
						{
							float sum = 0;
							int m = sel_Idx[l];  // selected isoform
							rd_Set.rd_Pro1[j][m] = rd_Set.rd_Pro[j][m]*abun[m];
							sum += rd_Set.rd_Pro1[j][m];
						}

						for(int l = 0; l<kmax; l++)
						{
							int m = sel_Idx[l];  // selected isoform
							rd_Set.rd_Pro1[j][m] /= sum;       // probability normalization
							alpha[l] += rd_Set.rd_Pro1[j][m];  // increase the probability
						}
					}
#pragma endregion calculate post probability
					abun[i] = alpha[i]/rd_Num;  // optimal estimation
					abun_Sum += abun[i];
				}
				else
					knz--;
			} // end for
#pragma region calculate objective function
			float temp_Sum = 0;
			for(int l = 0; l<kmax; l++)
			{
				abun[l] /= abun_Sum;     // normalization
					temp_Sum += log(rd_Num*abun[l]/12);
			}

			float obj1 = 1.0/2*temp_Sum + knz/2*log(double(rd_Num/12.0)) + knz;
			float obj2 = - rd_Set.Px(abun);
			L2 = obj1 + obj2;  
#pragma endregion calculate objective function
			iter++;

		}while(knz>kmin&&abs(1-L1/(L2+1e-09))<1e-04);	

#pragma region set the minimum abundance to 0
		float min1 = 10, sum1 = 0;
		int m_Idx = -1;
		for(int l = 0; l<kmax; l++){
			if(abun[l]>0&&abun[l]<min1){
				min1 = abun[l];
				m_Idx = l;
				sum1 += abun[l];
			}
		}
		abun[m_Idx] = 0;  knz--;
		// normalization
		for(int l = 0; l<kmax; l++)
			abun[l] /= sum1; 

		if(L2<L_min)
		{
			L_min = L2;  // optimal function value
			memcpy(abun_Opt,abun,sizeof(float)*kmax);
		}

#pragma endregion

	} // end while
}

// sort by abundance
void ISA::rank_Dim(vector<PAIR>& pair_Vec, float* sol)
{
	v_Idx.clear(); // clear array
	pair_Vec.clear();

	for(int i = 0; i<m_dataDim; i++)
	{
		pair_Vec.push_back(make_pair(i,sol[i]));
	}

	sort(pair_Vec.begin(),pair_Vec.end(),cmp);

	int i = 0;
	float thresh = 0.5/m_dataDim;  // threshold
	int v_Num = pair_Vec.size();

	vector<PAIR>::iterator iter = pair_Vec.begin();  // pointer of the start point
	for(vector<PAIR>::iterator iter=pair_Vec.begin(); iter!=pair_Vec.end();++iter)
	{
		if(iter->second>thresh)
		{
			v_Idx.push_back((iter->first));   // valid dimension
			i++;
		}
	}
}

// output abundance estiamtino
void ISA::report(Rd_Set& rd_Set, Graph_Trans& gph, char* filename, vector<int>& sel_Idx)
{
	ofstream outfile;
	outfile.open(filename,ios::app);
	if(!outfile)
		cout<<"Open file error!(report)"<<endl;

	int num_path = sel_Idx.size();   // number of paths
	string chr_Name = rd_Set.chromosome_Name;
	string method_Name = "MaxInfo";
	char gene_Id[100],trs_Id[100];
	int gene_id = 1, trs_id = 1;
	
	double mapped_FragNum = 2000000;

	for(int j = 0; j<num_path; j++)
	{
		trs_id = j+1;
		int k1 = sel_Idx[j];
		int i = rd_Set.idx[k1];
		int exon_Num = rd_Set.path[i].exon_Number;  //  number of exons in the paths
#pragma region output isoform
		int idx0 = rd_Set.path[i].exon_Node[0];     // index of first exons in the paths
		int idx1 = rd_Set.path[i].exon_Node[exon_Num-1];  // index of last exons in the paths
		int start = gph.exon_Vec[idx0].bound1;     // start point of isoform
		int stop = gph.exon_Vec[idx1].bound2;      // end point of isoform
		double fpkm = 1e09*1.0*rd_Set.locate_RdNum*m_CurOpt[k1]/(mapped_FragNum*rd_Set.path[i].length);
		sprintf(trs_Id,"%s.%d",rd_Set.gene_Name.c_str(),trs_id);
		outfile<<chr_Name<<"\t"<<method_Name<<"\t"<<"transcript"<<"\t"<<start<<"\t"<<stop<<"\t"
			<<"."<<"\t"<<rd_Set.orientation<<"\t"<<"."<<"\t"<<"gene_id"<<"\t"<<rd_Set.gene_Name<<"\t"<<"transcript_id"<<"\t"<<trs_Id<<"\t"
			<<"exon_number"<<"\t"<<exon_Num<<"\t"<<"FPKM"<<"\t"<<fpkm<<"\t"<<"frac"<<"\t"<<m_CurOpt[k1]<<"\t"
			<<"cov"<<"\t"<<rd_Set.locate_RdNum<<endl;
#pragma endregion output isoform
#pragma region output exon

		for(int k = 0; k<exon_Num; k++)
		{
			int idx = rd_Set.path[i].exon_Node[k];     // index of exons in the paths
			int s1 = gph.exon_Vec[idx].bound1;   // start point of exon
			int s2 = gph.exon_Vec[idx].bound2;   // end point of exon
			outfile<<chr_Name<<"\t"<<method_Name<<"\t"<<"exon"<<"\t"<<s1<<"\t"<<s2<<"\t"
				<<"."<<"\t"<<rd_Set.orientation<<"\t"<<"."<<"\t"<<"gene_id"<<"\t"<<rd_Set.gene_Name<<"\t"
				<<"transcript_id"<<"\t"<<trs_Id<<"\t"
				<<"exon_number"<<"\t"<<k+1<<"\t"<<"FPKM"<<"\t"<<fpkm<<"\t"<<"frac"<<"\t"<<m_CurOpt[k1]<<"\t"
				<<"cov"<<"\t"<<rd_Set.rd_Num<<endl;
		}
#pragma endregion output exon
	}
	outfile.close();
}

// output optimum
void ISA::output_Solution(float* x)
{
	float beta = pow(float(m_iteration), m_annealFactor);
	vector<float> weights(m_numSample);
	float weightSum = 0.f;
	for(int i=0;i<m_numSample;i++) {
		weights[i] = exp(-m_energy[i]*beta);
		// weights[i] = weights[i]>1.f?0:weights[i];
		weightSum += weights[i];
	}
	if(weightSum==0) {
		printf("can not find the solution!\n");
		return;
	}

	memset(x,0,sizeof(float)*m_dataDim);
	for(int i=0;i<m_numSample;i++) {
		weights[i] /= weightSum;
		for(int j=0;j<m_dataDim;j++) {
			x[j] += (m_samples+i*m_dataDim)[j]*weights[i];
		}
	}

	// output optimum
	for(int j=0;j<m_dataDim;j++) 
		printf("%.4f ",m_CurOpt[j]);
	printf("\n");

	// output optimum by weight
	for(int j=0;j<m_dataDim;j++) 
		printf("%.4f ",x[j]);
	printf("\n");

	weights.swap(vector<float>());
}
	