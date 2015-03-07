#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include "time.h"

#include "ISA.h"

// Calculate the joint cost
float ISA::cost_Graph(Rd_Set& rd_Set, Graph_Trans& gph, float* abun)
{
	// junction ´ú¼Û
	double cost1 = cost_Junction(rd_Set, gph, abun);

	int cnt = 0;
	int e_Num = gph.vtx_Number-2;
	double cost2 = 0;	
	double rate = 0;
	for(int i = 0; i<m_dataDim; i++)
	{
		int k = rd_Set.idx[i];
		double temp = rd_Set.path[k].link_Cost1 + rd_Set.path[k].link_Cost2*rate;
		cost2 += abun[i]*temp;    // cost of path connectivity and path length
	}

	float lambda1 = 1, lambda2 = 0.5;

	return lambda1*cost1 + lambda2*cost2;
}

// Graph cost: cost of junctions
float ISA::cost_Junction(Rd_Set& rd_Set, Graph_Trans& gph, float* abun)
{
	int id_Num = rd_Set.idx.size();  // number of valid paths
	int e_Num = gph.vtx_Number - 2;
	double cost = 0, dis = 0, est_Num = 0;
	int rank = 0;

	for(int i = 1; i<=e_Num; i++)
	{
		for(int j = 1; j<=e_Num; j++)
		{
			if(gph.link_Mtx[i][j]>0)
			{
				for(int k = 0; k<id_Num; k++)
				{
					est_Num = rd_Set.junc_Mtx[rank][k]*abun[k];
					dis = gph.link_Mtx[i][j]-est_Num;
					// cost += dis*dis;
					cost += abs(dis);
				}
				rank++;
			}
		}
	}

	double cost1 = sqrt(cost)/(0.5*id_Num*rank);

	return cost1;
}

// Forward selection of the path
void  ISA::path_SFMS(float* abun, Rd_Set& rd_Set, Graph_Trans& gph, vector<int>& sel, vector<int>& nsel, vector<PAIR>& pair_Vec, int select_pathNum)
{
	float* ini_Abun = new float[m_dataDim];
	int id_num = rd_Set.idx.size();     // number of selected paths
	int K = _min(select_pathNum,id_num); 

	float link_Thresh = 1.5;            // threshold for the connectivity cost
	float abun_Thresh = 0.005;          // threshold of the abundance
	bool* flag = new bool[m_dataDim];   // flag indicating whether the path can be selected
	for(int i = 0; i<m_dataDim; i++)
		flag[i] = false;

	int cnt = 0;
	int i = 0;
	for(vector<PAIR>::iterator iter1=pair_Vec.begin(); iter1!=pair_Vec.end();++iter1)
	{
		int idx = iter1->first;
		int idx1 = rd_Set.idx[idx];
		float link_Cost = rd_Set.path[idx1].link_Cost1;  // connectivity cost of the path
		float cur_Abun = iter1->second;

		if(link_Cost<link_Thresh&&cur_Abun>abun_Thresh)
			flag[idx] = true;   // possible to be selected
		else
			flag[idx] = false;  // not to be selected

		if(i==0)  // first node
		{
			if(flag[idx]==true){
				sel.push_back(idx); cnt++;
			}
			else{
				float p = UnifoRand();
				if(p>0.5){
					sel.push_back(idx); cnt++;
				}
				else  nsel.push_back(idx);
			}
		}
		else{
			if(flag[idx]==false)
				nsel.push_back(idx);
			else if(cnt<1){
				sel.push_back(idx);
				cnt++;
			}
			else{
				nsel.push_back(idx);
			}
		}
		i++;
	}

	float func_Max = -1000000; // optimal function value
	vector<int>::iterator iter;
	vector<int> nsel1 = nsel;
	vector<int> sel1 = sel;
	int sel_Idx = -1, idx_Serial = -1;
	int sel_Num = sel.size();  // number of currently selected paths
	int iter3 = 30, iter_Ctr = 10;
	while(sel_Num<K)
	{
		int n1 = nsel.size();
		func_Max = -1000000;
		for(int k = 0; k<n1; k++)
		{
			int l = nsel[k];
			if(flag[l]==true)
			{
				sel1 = sel; 
				nsel1 = nsel;
				sel1.push_back(l);
				iter = nsel1.begin()+k;
				nsel1.erase(iter); 
				memcpy(ini_Abun, abun, sizeof(float)*m_dataDim);
				sparse_Normalize(ini_Abun,sel1, nsel1);  // normalize on the sparse dimensions
				float function_Value =  Optimize_EM(ini_Abun,rd_Set, gph, sel1, nsel1, iter3, iter_Ctr);        // optimization with EM algorithms
				// float function_Value =  Optimize_EM_Re(ini_Abun,rd_Set, gph, sel1, nsel1, iter3, iter_Ctr);  // optimization with regularized EM algorithms
				if(function_Value>func_Max){
					func_Max = function_Value;
					sel_Idx = l;
					idx_Serial = k;
				}
			}
		}
		sel.push_back(sel_Idx);  // add the path and update sel
		iter = nsel.begin()+idx_Serial;
		nsel.erase(iter);        // clear the path and update nsel
		sel_Num = sel.size();    // current selected paths
	}
}

// Forward selection of the path
void  ISA::path_SFMS_Auto(float* abun, Rd_Set& rd_Set, Graph_Trans& gph, vector<int>& sel, vector<int>& nsel, vector<PAIR>& pair_Vec, int select_pathNum)
{
	float* ini_Abun = new float[m_dataDim];
	int id_num = rd_Set.idx.size();    // number of selected paths 
	int K = _min(select_pathNum,id_num);   

	float link_Thresh = 1.5;           // threshold for the connectivity cost
	float abun_Thresh = 0.005;         // threshold of the abundance
	bool* flag = new bool[m_dataDim];  // flag indicating whether the path can be selected
	for(int i = 0; i<m_dataDim; i++)
		flag[i] = false;
	int path_cnt = 0;
	int i = 0;
#pragma region
	int valid_trs = 0;   // number of valid transcripts 
	for(vector<PAIR>::iterator iter1=pair_Vec.begin(); iter1!=pair_Vec.end();++iter1)
	{
		int idx = iter1->first;
		int idx1 = rd_Set.idx[idx];
		float link_Cost = rd_Set.path[idx1].link_Cost1;  // connectivity cost of the path
		float cur_Abun = iter1->second;
		if(link_Cost<link_Thresh&&cur_Abun>abun_Thresh){
			flag[idx] = true;   // possible to be selected
			valid_trs++;
		}
		else
			flag[idx] = false;  // not to be selected

		if(i==0)  // first node
		{
			if(flag[idx]==true){
				sel.push_back(idx); path_cnt++;
			}
			else{    // if the constraint isn't satisfied by the node, the node is selected by a probability.
				float p = UnifoRand();
				if(p>0.5){
					sel.push_back(idx); path_cnt++;
				}
				else  nsel.push_back(idx);
			}
		}
		else{
			if(flag[idx]==false)
				nsel.push_back(idx);
			else if(path_cnt<1){
				sel.push_back(idx);
				path_cnt++;
			}
			else{
				nsel.push_back(idx);
			}
		}
		i++;
	}
#pragma endregion

	float func_Max = -1000000;   // optimal function value
	vector<int>::iterator iter;
	vector<int> nsel1 = nsel;
	vector<int> sel1 = sel;
	int sel_Idx = -1, idx_Serial = -1;
	int sel_Num = sel.size();    // number of currently selected paths
	int iter3 = 30, iter_Ctr = 10;
	float mi_threshold = -1;
	float mi = 0, mi_pre = 0, delta = 10;

	while(path_cnt<valid_trs&&delta>mi_threshold)
	{
		mi_pre = mi;
		int n1 = nsel.size();
		func_Max = -1000000;
		for(int k = 0; k<n1; k++)
		{
			int l = nsel[k];
			if(flag[l]==true)
			{
				sel1 = sel; 
				nsel1 = nsel;
				sel1.push_back(l);
				iter = nsel1.begin()+k;
				nsel1.erase(iter); 
				memcpy(ini_Abun, abun, sizeof(float)*m_dataDim);
				sparse_Normalize(ini_Abun,sel1, nsel1);  // normalize on the sparse dimensions
				float function_Value =  Optimize_EM(ini_Abun,rd_Set, gph, sel1, nsel1, iter3, iter_Ctr);  // optimization with EM algorithm
				if(function_Value>func_Max){
					func_Max = function_Value;
					sel_Idx = l;
					idx_Serial = k;
					memcpy(m_CurOpt,ini_Abun,m_dataDim*sizeof(float));
				}
			}
		}
		sel.push_back(sel_Idx);  
		iter = nsel.begin()+idx_Serial;
		nsel.erase(iter);        
		sel_Num = sel.size();    

		rd_Set.sel_Idx = sel;
		mi = rd_Set.Hyx(m_CurOpt);
		delta = abs(1-mi_pre/(mi+1e-09));   // change of the mutual information
		path_cnt++;
	}
}

// forward selection of path
void  ISA::path_SFMS_Auto1(float* abun, Rd_Set& rd_Set, Graph_Trans& gph, vector<int>& sel, vector<int>& nsel, vector<PAIR>& pair_Vec, int select_pathNum)
{
	float* ini_Abun = new float[m_dataDim];
	int id_num = rd_Set.idx.size();  // candidate path
	int K = _min(select_pathNum,id_num);   // selected path

	vector<PAIR>::iterator iter1=pair_Vec.begin(); 
	float link_Thresh = 1.5;    // threshold of coneection cost
	float abun_Thresh = _max(5*1e-05,iter1->second*0.005);  // threshold of abundance

	bool* flag = new bool[m_dataDim];  // flag of whether the path can be selected
	for(int i = 0; i<m_dataDim; i++)
		flag[i] = false;
	int path_cnt = 0;
	int interval = 5;
#pragma region analysis the first node and the last node
	int valid_trs = 0;   // count valis isoform
	for(iter1=pair_Vec.begin(); iter1!=pair_Vec.end();++iter1)
	{
		int idx = iter1->first;
		int idx1 = rd_Set.idx[idx];   // selected path
		float link_Cost = rd_Set.path[idx1].link_Cost1;  // estimated path cost of the first node and last node
		float cur_Abun = iter1->second;
		if(link_Cost<link_Thresh&&cur_Abun>abun_Thresh){
			flag[idx] = true;   // candidate isoform
			valid_trs++;        // valid isoform
		}
		else
			flag[idx] = false;  // not a isoform

		// cout<<idx<<"\t";

		nsel.push_back(idx);
	}
	cout<<endl;

#pragma endregion selected path

	float func_Max = -1000000;   // optimal value
	vector<int>::iterator iter;
	vector<int> nsel1 = nsel;
	vector<int> sel1 = sel;
	int sel_Idx = -1, idx_Serial = -1;
	int sel_Num = sel.size();  // selected path
	int iter3 = 30, iter_Ctr = 10;
	float mi = 0, mi_pre = 0, delta = 10;  // mutual information

#pragma region write probability to files
	// ofstream outfile;
	// outfile.open("G:\\My documents\\FragmentStreaming\\D Melanogaster\\chr3R\\path_selection1.txt",ios::app);
	//if(!outfile){
	//	cout<<"Open file error!"<<endl;
	//	// exit(-1);
	//}
	//cout<<"record path_SFMS_Auto"<<endl;
	//ofstream outfile1;
	//outfile1.open("G:\\My documents\\FragmentStreaming\\D Melanogaster\\chr3R\\likelihood_selection1.txt",ios::app);
#pragma endregion 

	float temp_abun = 0;
	float temp_thresh = 0.01;
	float re_ml = -1000, re_mlPre = 0;    // regularized likelihood
	float gamma = 0.7;        // regularized parameter
	float mi_threshold = -1;  
	float re_threshold = 0;   
	vector<int> nsel_bak;     // unselected path
	float* t_abun = new float[m_dataDim];
	bool flag_Sel = false;

#pragma region path selection 
	while(delta>re_threshold)  // path available, variation of mutual information exceed the threshold
	{
		re_mlPre = re_ml;
		int n1 = nsel.size();  // number of unselected paths
		func_Max = -1000000;
		flag_Sel = false;      // flag of selected isoform
		sel_Idx = -1;
		for(int k = 0; k<n1; k++)
		{
			int l = nsel[k];   // path index
			if(flag[l]==true)
			{
				sel1 = sel; 
				nsel1 = nsel;
				float simi = isoform_Similarity(rd_Set, sel, rd_Set.idx[l]);  // isoform similarity
				if(simi>0.99){
					flag[l] = false;   // unselected isoform
				}
				else
				{
					sel1.push_back(l);  // add to the set
					iter = nsel1.begin()+k;
					nsel1.erase(iter);  // delete  element

					memcpy(ini_Abun, abun, sizeof(float)*m_dataDim);
					sparse_Normalize(ini_Abun,sel1, nsel1);  // normalize the sparse array
					float function_Value =  Optimize_EM(ini_Abun,rd_Set, gph, sel1, nsel1, iter3, iter_Ctr);  // optimize the function by EM
					if(function_Value>func_Max){	
						flag_Sel = true;
						float add_abun = ini_Abun[l];   // increase the abundance of the isoform
						if(add_abun<0.1)  // if the abundance of the isoform is small, check whether the isoform exists
						{
							if(simi>0.90){
								flag_Sel = false;  // flag of whether the abundance is updated
								flag[l] = false;    // candidate isoform
							}
						}
						if(flag_Sel){
							func_Max = function_Value;
							sel_Idx = l;
							idx_Serial = k;
							memcpy(t_abun,ini_Abun,m_dataDim*sizeof(float));  // copy current isoform to optimum
						}
					}
				}
			}
		}

		if(sel_Idx>-1){
			sel.push_back(sel_Idx);  // add current path index, update sel
			iter = nsel.begin()+idx_Serial;     // delete path

			nsel.erase(iter);        // delete, update nsel
			sel_Num = sel.size();    // selected path

			rd_Set.sel_Idx = sel;
			float ml = rd_Set.Px(t_abun,interval)*interval/rd_Set.rd_Num;

			float con_mi = rd_Set.Hyx(t_abun,interval);

			re_ml = ml - gamma*con_mi;
			delta = re_ml - re_mlPre;

			int num = sel.size();  // selected path
			float add_Abun = t_abun[sel_Idx];

			if(delta>0&&add_Abun>0.005)
			{
				memcpy(m_CurOpt,t_abun,m_dataDim*sizeof(float));
			} 
			else
			{
				nsel.push_back(sel[num-1]);
				sel.pop_back();			
			}

			path_cnt++;
			if(path_cnt>6)
				break;   // number of path
		}
		else
			break;  // no isoform qualified
	}// end while
#pragma endregion path selection

	// outfile.close();
	// outfile1.close();
}

// add annotated isoforms
bool  ISA::path_Compare_Annotated(float* abun, Rd_Set& rd_Set, Graph_Trans& gph, vector<int>& sel, vector<int>& nsel)
{
	int anno_Num = gph.annotated_Idx.size();  // number of indexes
	int id_Num = rd_Set.idx.size();  // candidate path
	int sel_Num = sel.size();        // selected path
	int cnt = 0;
	vector<bool> Flag(anno_Num,false);
	for(int j=0; j<sel_Num; j++)
	{
		int k = sel[j];
		for(int i=0; i<anno_Num; i++)
		{
			if(rd_Set.idx[k]==gph.annotated_Idx[i])
			{
				Flag[i] = true;
				break;
			}
		}
	}
	for(int i=0; i<anno_Num; i++)
	{
		if(Flag[i]==false)   // unselected path
		{
			cnt++;
			int serial = 0;
			for(int j=0; j<id_Num; j++)
			{
				if(rd_Set.idx[j]==gph.annotated_Idx[i])
				{
					serial = j;
					cnt--;
					break;
				}
			}
			sel.push_back(serial);  // add selected path
			int nsel_Num = nsel.size();
			for(int l=0; l<nsel_Num; l++)
			{
				int k = nsel[l];
				if(rd_Set.idx[k]==gph.annotated_Idx[i])
				{
					nsel.erase(nsel.begin()+l);  // delete the index
					break;
				}
			}
		}
	}

	if(cnt==0)
		return true;

	return false;
}

// limited to annotated isoforms
bool  ISA::path_Limit_Annotated(float* abun, Rd_Set& rd_Set, Graph_Trans& gph, vector<int>& sel, vector<int>& nsel)
{
	int anno_Num = gph.annotated_Idx.size();  // number of indexes
	int id_Num = rd_Set.idx.size();  // candidate path
	int sel_Num = sel.size();        // selected path
	int cnt = 0;
	vector<bool> Flag(anno_Num,false), sel_Flag(sel_Num,false);
	for(int j=0; j<sel_Num; j++)
	{
		int k = sel[j];
		for(int i=0; i<anno_Num; i++)
		{
			if(rd_Set.idx[k]==gph.annotated_Idx[i])
			{
				Flag[i] = true;
				sel_Flag[i] = true;
				break;
			}
		}
	}
	for(int j=0; j<sel_Num; j++)
	{
		if(sel_Flag[j]==false)   // unmatched path
		{
			int k = sel[j];
			sel.erase(sel.begin()+j);  // delete the index
			nsel.push_back(k);
			break;
		}
	}
	for(int i=0; i<anno_Num; i++)
	{
		if(Flag[i]==false)   // unselected path
		{
			cnt++;
			int serial = 0;
			for(int j=0; j<id_Num; j++)
			{
				if(rd_Set.idx[j]==gph.annotated_Idx[i])
				{
					serial = j;
					cnt--;
					break;
				}
			}
			sel.push_back(serial);  // add selected path
			int nsel_Num = nsel.size();
			for(int l=0; l<nsel_Num; l++)
			{
				int k = nsel[l];
				if(rd_Set.idx[k]==gph.annotated_Idx[i])
				{
					nsel.erase(nsel.begin()+l);  // delete the index
					break;
				}
			}
		}
	}

	if(cnt==0)
		return true;

	return false;
}

// Examine the overlap of the paths
// Examine whether path 1 is the subset of path 2
float  ISA::isSubset(Path& path1, Path& path2)
{
	int len1 = path1.exon_Number;  // number of exons in path 1
	int len2 = path2.exon_Number;  // number of exons in path 2

	float simi_Cnt = 0;
	int i = 0;
	while(i<len1&&path1.exon_Node[i]!=path2.exon_Node[0])
		i++;
	if(i<len1){
		for(int j = 0; j<len2; j++){
			if(i<len1&&path1.exon_Node[i]==path2.exon_Node[j]){
				simi_Cnt++;
				i++;
			}
			else
				break;
		}
	}

	float simi = simi_Cnt*1.0/len2, similarity = simi;
	if(i==len1)  // isoform 1 is the subset of isoform 2
	{
		if(simi<0.9)
			similarity = 1;
	}
	else if(i>=len2)  // isoform 2 is the subset of isoform 1
	{
		similarity = 0.5;
	}
	else
	{
		similarity = simi_Cnt*1.0/(_max(len1,len2));
	}

	return similarity;
}

// Measure the similarity of one path to the others
float  ISA::isoform_Similarity(Rd_Set& rd_Set, vector<int>& sel, int sel_idx)
{
	int len2 = rd_Set.path[sel_idx].length;
	int num = sel.size();   // number of selected isoforms
	float similarity = 0;
	float simi = 0;
	for(int k = 0; k<num; k++)
	{
		int k1 = sel[k];
		int idx = rd_Set.idx[k1];   // index of the isoform
		simi = isSubset(rd_Set.path[sel_idx],rd_Set.path[idx]);  // similarity of the path to the other paths
		if(simi>0.99)
			return simi;
		if(simi>similarity)
			similarity = simi;
	}

	return similarity;
}