
#include "Read_Assignment.h"
using namespace std;

// Destructor function
Graph_Trans::~Graph_Trans()
{
	// cout<<"destructed"<<endl;
}

// Reset function
void Graph_Trans::clear()
{ 
	vector<Exon_Node>().swap(exon_Vec);
	vector<vector<int>>().swap(paths);
	vector<int>().swap(annotated_Idx);
	exon_Vec.clear();
	paths.clear();
}

// Clear the linkage
void Graph_Trans::clear_JuncRank()
{ 
	int exon_Num = exon_Vec.size();
	for(int i = 0; i<exon_Num; i++)
	{
		if(junc_Rank[i]!=NULL)
			delete[] junc_Rank[i];
	}
	if(junc_Rank!=NULL)
		delete[] junc_Rank;
	junc_Rank = NULL;

	delete[] first_NodeCost;
	delete[] last_NodeCost;
}

// Load the graph structure
void Graph_Trans::load_Graph(char* filename, int vex_Num)
{
	ifstream myfile;
	myfile.open(filename,ios::in);
	if(!myfile)
	{
		std::cout<<"can not open the file."<<endl;
		exit(-1);
	}

	vtx_Number = vex_Num;
	link_Mtx = new float*[vex_Num*vex_Num];

	string word, line;
	int temp = 0, i = 0;
	while(std::getline(myfile,line))
	{
		stringstream stream(line);
		link_Mtx[i] = new float[vex_Num];

		for(int j = 0; j<vex_Num; j++)
		{
			stream>>word;  // matrix element
			std::sscanf(word.c_str(),"%f",&(link_Mtx[i][j]));
		}
		i++;
	}
}

// Output the graph structure
void Graph_Trans::print_Graph()
{
	cout<<"number of vertex: "<<vtx_Number<<endl;
	ofstream outfile;
	char filename[200];
	sprintf(filename,"E:\\Splicing Graph\\data\\SRR585572\\chr2\\graph\\graph_%d.txt",vtx_Number);
	outfile.open(filename,ios::out);
	if(!outfile)
		cout<<"Open file error!(graph)"<<endl;

	for(int i = 0;i<vtx_Number; i++)
	{
		for(int j = 0 ; j<vtx_Number; j++)
		{
			cout<<link_Mtx[i][j]<<" ";
			outfile<<link_Mtx[i][j]<<" ";
		}
		cout<<endl;
		outfile<<endl;
	}
	for(int i = 1;i<vtx_Number-1; i++)
	{
		for(int j = 1 ; j<vtx_Number-1; j++)
		{
			if(link_Mtx[i][j]>0)
				outfile<<i<<" "<<j<<" "<<link_Mtx[i][j]<<endl;
		}
	}
	outfile.close();
}

////////////////////////////////////////////////////////////////////////////
// Modify the linkage in the matrix and clear low-confidence linkage
void Graph_Trans::graph_Modify()
{
	int e_Num = vtx_Number - 2; 
	float cnt = 0;
	int junc_cnt = 0;
	for(int i = 1; i<=e_Num; i++)
	{
		for(int j = 1; j<i; j++)
		{
			if(link_Mtx[j][i]>0)
			{
				cnt += link_Mtx[j][i];
				junc_cnt++;
			}
		}
	}
	
	float mean_Cnt = cnt*1.0/junc_cnt;  
	float thresh = mean_Cnt*0.1;        // threshold for a resonable linkage

#pragma region clear the low-confidence linkage
	//for(int i = 1; i<=e_Num; i++)
	//{
	//	for(int j = 1; j<i; j++)
	//	{
	//		if(link_Mtx[j][i]<thresh)
	//			link_Mtx[j][i] = 0;
	//	}
	//}
#pragma endregion

}

// Build directed graph
void Graph_Trans::graph_Build(int exon_Num)
{
	int vtx_Num = exon_Num + 2;
	vtx_Number = vtx_Num;  
	link_Mtx = new float*[vtx_Num];
	for(int i = 0; i< vtx_Num; i++)
	{
		link_Mtx[i] = new float[vtx_Num];
		for(int j = 0; j<vtx_Num; j++)
			link_Mtx[i][j] = 0;
	}

	// initialize the matrix
	for(int i = 1;i<vtx_Num-1;i++)
		link_Mtx[0][i] = 1;
	for(int i = 1;i<vtx_Num-1;i++)
		link_Mtx[i][vtx_Num-1] = 1;
	for(int i = 1;i<vtx_Num-2;i++)
	{
		link_Mtx[i][i+1] = 1;
		link_Mtx[i][i+2] = 1;
	}
	link_Mtx[2][4] = 0;
}

// Build directed graph based on the read data
void Graph_Trans::graph_Initial(int exon_Num)
{
	if(exon_Num<=0)
		return;

	vtx_Number = exon_Num + 2;
	link_Mtx = new float*[vtx_Number];
	junc_Num = 0;

	// error control
	if(exon_Vec[exon_Num-1].bound2<exon_Vec[exon_Num-1].bound1)
		exon_Vec[exon_Num-1].bound2 = exon_Vec[exon_Num-1].bound1 + 200;	

	for(int i = 0; i<vtx_Number; i++)
	{
		link_Mtx[i] = new float[vtx_Number];
		for(int j = 0; j<vtx_Number; j++)
			link_Mtx[i][j] = 0;         // initialize the linkage
	}

	for(int i = 1; i<vtx_Number-1; i++)
		link_Mtx[0][i] = 1;             // all the exons are initialized as possible starting nodes
	for(int i = 1; i<vtx_Number-1; i++)
		link_Mtx[i][vtx_Number-1] = 1;  // all the exons are initialized as possible stopping nodes

	int mark = 0;
	for(int i=0; i<exon_Num; i++)
	{
		junc_Link1.push_back(0);
		junc_Link2.push_back(0);
		if(exon_Type[i]==4&&mark%2==0&&i<exon_Num-1)
		{
			link_Mtx[i][i+1]++;
			mark++;
		}
	}
}

// Build directed graph based on the read data
void Graph_Trans::graph_Initial_1(int exon_Num)
{
	if(exon_Num<=0)
		return;

	vtx_Number = exon_Num + 2;
	link_Mtx = new float*[vtx_Number];
	junc_Num = 0;

	// error control
	if(exon_Vec[exon_Num-1].bound2<exon_Vec[exon_Num-1].bound1)
		exon_Vec[exon_Num-1].bound2 = exon_Vec[exon_Num-1].bound1 + 200;	

	for(int i = 0; i<vtx_Number; i++)
	{
		link_Mtx[i] = new float[vtx_Number];
		for(int j = 0; j<vtx_Number; j++)
			link_Mtx[i][j] = 0;         // initialize the linkage
	}

	for(int i = 1; i<vtx_Number-1; i++)
		link_Mtx[0][i] = 1;             // all the exons are initialized as possible starting nodes
	for(int i = 1; i<vtx_Number-1; i++)
		link_Mtx[i][vtx_Number-1] = 1;  // all the exons are initialized as possible stopping nodes

	int mark = 0;
	for(int i=0; i<exon_Num; i++)
	{
		junc_Link1.push_back(0);
		junc_Link2.push_back(0);
	}
}

// Build directed graph containing both exons and introns based on the read data
void Graph_Trans::graph_Initial(vector<int>& add, int add_Num)
{
	cvtx_Number = exon_Vec.size() + add_Num + 2;
	clink_Mtx = new float*[cvtx_Number];
	
	int junc_Num1 = 0;
	for(int i = 0; i<cvtx_Number; i++)
	{
		clink_Mtx[i] = new float[cvtx_Number];
		for(int j = 0; j<vtx_Number; j++)
			clink_Mtx[i][j] = 0;         // initialize the linkage
	}
	for(int i = 1; i<cvtx_Number-1; i++)
		clink_Mtx[0][i] = 1;             // all the exons are initialized as possible starting nodes
	for(int i = 1; i<vtx_Number-1; i++)
		clink_Mtx[i][vtx_Number-1] = 1;  // all the exons are initialized as possible stopping nodes

	for(int i = 0;i<add.size();i++)
	{
		int temp = add[i];     // intron added
		clink_Mtx[0][temp+i+2] = 0;      // introns are initialized as without linkage
		clink_Mtx[temp+i+2][0] = 0;
	}
}

// Build indices of junctions
void Graph_Trans::junc_Index()
{
	int e_Num = vtx_Number - 2;
	junc_Rank = new int*[e_Num];
	for(int i = 0; i<e_Num; i++)
		junc_Rank[i] = new int[e_Num];
	int rank = 1;
	for(int i = 1; i<=e_Num; i++)
	{
		for(int j = 1; j<=e_Num; j++)
		{
			if(link_Mtx[i][j]>0)
			{
				junc_Rank[i-1][j-1] = rank;
				rank++;
			}
			else
				junc_Rank[i-1][j-1] = 0;
		}
	}
	junc_Num = rank-1;  // number of junctions
}

// Constraints on the paths with the costs of the first node and the last node
bool Graph_Trans::graph_Constraint_Initial()
{
	int exon_Num = vtx_Number - 2;
	aver_in = 0;    // average in-flow
	aver_out = 0;   // average out-flow
	int in_degree = 0, out_degree = 0;
	float cnt1 = 0, cnt2 = 0, cnt = 0, thresh = 0;

	for(int i = 1; i<=exon_Num; i++)
	{
		for(int j = 1; j<i; j++)
			cnt += link_Mtx[j][i];
		in_Flow.push_back(0);
		out_Flow.push_back(0);
	}
	float mean_Cnt = cnt*1.0/exon_Num;

#pragma region in-degree and out-degree of the nodes
	for(int i=0; i<vtx_Number; i++)
	{
		in_Degree.push_back(0);
		out_Degree.push_back(0);
		VEC_INT temp1, temp2;
		pre_Nodes.push_back(temp1); pos_Nodes.push_back(temp2);
	}

	for(int i=1; i<vtx_Number-1; i++)  
	{
		in_degree = 0; out_degree = 0;
		float in_flow = 0, out_flow = 0;
		for(int j=1; j<i; j++)
		{
			if(link_Mtx[j][i]>0)
			{
				in_flow += link_Mtx[j][i];
				in_degree++;
				pre_Nodes[i].push_back(j); 
			}
		}
		for(int j=i+1; j<vtx_Number-1; j++)
		{
			if(link_Mtx[i][j]>0)
			{
				out_flow += link_Mtx[i][j];
				out_degree++;
				pos_Nodes[i].push_back(j);
			}
		}
		in_Degree[i] = in_degree;
		out_Degree[i] = out_degree;
		aver_in += in_flow;
		aver_out += out_flow;
		in_Flow[i-1] = in_flow; out_Flow[i-1] = out_flow;
	}
#pragma endregion

#pragma region in-flow and out-flow of the nodes
	aver_in /= vtx_Number-2;
	aver_out /= vtx_Number-2;
	float thresh_1 = aver_in*0.2, thresh_2 = aver_out*0.2;

	link_Mtx[0][1] = 1;
	link_Mtx[vtx_Number-2][vtx_Number-1] = 1;

	for(int j = 1; j<exon_Num; j++)
	{
		if(exon_Type[j]==0)  // exon at the boundary of the gene
		{
			if(out_Flow[j]>0)
				link_Mtx[j+1][vtx_Number-1] = 0;  // the first node
			if(in_Flow[j]>0)
				link_Mtx[0][j+1] = 0;  // the last node
		}
		else if(exon_Type[j]==1)  // only the right end of the exon is identified; the left end isn't identified
		{
			if(out_Flow[j]>thresh_1&&out_Flow[j]>0.5*in_Flow[j])   // out-flow is significantly bigger than in-flow
				link_Mtx[j+1][vtx_Number-1] = 0;
		}
		else if(exon_Type[j]==2)  // only the right end of the exon is identified; the right end isn't identified
		{
			if(in_Flow[j]>thresh_2||in_Flow[j]>0.5*out_Flow[j])    // out-flow is significantly bigger than in-flow
				link_Mtx[0][j+1] = 0;
		}
		else if(exon_Type[j]==3||exon_Type[j]==4)
		{
			if(j<0.5*exon_Num)
			{
				if(in_Flow[j]>thresh_1){
					if(out_Flow[j]<0.5*in_Flow[j])
						link_Mtx[0][j+1] = 0;
					if(out_Flow[j]>thresh_2){
						link_Mtx[j+1][vtx_Number-1] = 0;
						if(abs(1-in_Flow[j]/(0.1+out_Flow[j]))<0.2&&exon_Type[j]!=5/*||in_Flow[j]>out_Flow[j]/3.0*/) // 流量较为平衡
						{
							link_Mtx[0][j+1] = 0;
						}
					}
				}
			}
			else 
			{
				if(out_Flow[j]>thresh_2){
					if(in_Flow[j]<0.5*out_Flow[j])
						link_Mtx[j+1][vtx_Number-1] = 0;
					if(in_Flow[j]>thresh_1){
						link_Mtx[0][j+1] = 0;
						if(abs(1-in_Flow[j]/(0.1+out_Flow[j]))<0.2&&exon_Type[j]!=6/*||out_Flow[j]>in_Flow[j]/3.0*/) // 流量较为平衡或前后
						{
							link_Mtx[j+1][vtx_Number-1] = 0;							
						}
					}
				}
			}
		}
	}

#pragma endregion

#pragma region
	int out_start = 0, in_end = 0;
	for(int i=1; i<vtx_Number-1; i++)
	{
		int temp = link_Mtx[0][i];
		out_start += temp;
		in_Degree[i] += temp;
		if(temp>0)
		{
			pos_Nodes[0].push_back(i);
			pre_Nodes[i].insert(pre_Nodes[i].begin(),0);
		}
	}
	out_Degree[0] = out_start;

	for(int i=1; i<vtx_Number-1; i++)
	{
		int temp = link_Mtx[i][vtx_Number-1];
		in_end += temp;
		out_Degree[i] += temp;
		if(temp>0)
		{
			pre_Nodes[vtx_Number-1].push_back(i);
			pos_Nodes[i].push_back(vtx_Number-1);
		}
	}
	in_Degree[vtx_Number-1] = in_end;
#pragma endregion

	return true;
}

// Constraints on the paths when annotations are available
bool Graph_Trans::graph_Constraint_Annotated()
{
	Gene& gene = cur_Gene;
	int exon_Num = vtx_Number - 2;
	aver_in = 0;    // average in-flow
	aver_out = 0;   // average out-flow
	int in_degree = 0, out_degree = 0;
	float cnt1 = 0, cnt2 = 0, cnt = 0, thresh = 0;

	for(int i = 1; i<=exon_Num; i++)
	{
		for(int j = 1; j<i; j++)
			cnt += link_Mtx[j][i];
		in_Flow.push_back(0);
		out_Flow.push_back(0);
	}

	float mean_Cnt = cnt*1.0/exon_Num;

#pragma region in-degree and out-degree of the nodes
	for(int i=0; i<vtx_Number; i++)
	{
		in_Degree.push_back(0);
		out_Degree.push_back(0);
		VEC_INT temp1, temp2;
		pre_Nodes.push_back(temp1); pos_Nodes.push_back(temp2);
	}

	for(int i=1; i<vtx_Number-1; i++)
	{
		in_degree = 0; out_degree = 0;
		float in_flow = 0, out_flow = 0;
		for(int j=1; j<i; j++)
		{
			if(link_Mtx[j][i]>0)
			{
				in_flow += link_Mtx[j][i];
				in_degree++;
				pre_Nodes[i].push_back(j);
			}
		}
		for(int j=i+1; j<vtx_Number-1; j++)
		{
			if(link_Mtx[i][j]>0)
			{
				out_flow += link_Mtx[i][j];
				out_degree++;
				pos_Nodes[i].push_back(j);
			}
		}
		in_Degree[i] = in_degree;
		out_Degree[i] = out_degree;
		aver_in += in_flow;
		aver_out += out_flow;
		in_Flow[i-1] = in_flow; out_Flow[i-1] = out_flow;
	}
#pragma endregion

#pragma region constraints on the nodes
	aver_in /= vtx_Number-2;
	aver_out /= vtx_Number-2;
	float thresh_1 = aver_in*0.2, thresh_2 = aver_out*0.2;

	for(int j = 1; j<exon_Num; j++)
	{
		if(j<0.5*exon_Num)
		{
			if(in_Flow[j]>thresh_1){
				if(out_Flow[j]<0.5*in_Flow[j])
					link_Mtx[0][j+1] = 0;   // low possibility of being the first node
				if(out_Flow[j]>thresh_2){					
					if(abs(1-in_Flow[j]/(0.1+out_Flow[j]))<0.25)   // the flow is balanced
					{
						link_Mtx[j+1][vtx_Number-1] = 0;
						link_Mtx[0][j+1] = 0;
					}
				}
			}
		}
		else
		{
			if(out_Flow[j]>thresh_2){
				if(in_Flow[j]<0.5*out_Flow[j])
					link_Mtx[j+1][vtx_Number-1] = 0;   // low possibility of being the first node
				if(in_Flow[j]>thresh_1){					
					if(abs(1-in_Flow[j]/(0.1+out_Flow[j]))<0.25)   // the flow is balanced
					{
						link_Mtx[0][j+1] = 0;
						link_Mtx[j+1][vtx_Number-1] = 0;
					}
				}
			}
		}
	}

	link_Mtx[0][1] = 1;
	link_Mtx[vtx_Number-2][vtx_Number-1] = 1;

#pragma endregion

#pragma region judgenment of exon type according to the annotations
	vector<short> Flag(exon_Num,0);  // judgenment of exon type according to the annotations
	int num = gene.trans.size();     // number of annotated transcripts
	for(int i=1; i<exon_Num-1; i++)
	{
		int s1 = exon_Vec[i].bound1;   // left end of the exon
		int s2 = exon_Vec[i].bound2;   // right end of the exon
		if(Flag[i]>=0)
		{
			for(int j=0; j<num; j++)
			{
				int exon_num1 = gene.trans[j].exon_start.size();  // number of exons of the annotated transcript
				int t1 = gene.trans[j].exon_start[exon_num1-1];   // left boundary of the last exon
				int t2 = gene.trans[j].exon_stop[0];  // right boundary of the first exon
				if(abs(s2-t2)<5)
					Flag[i] = -1;  //  possible first node
				if(abs(s1-t1)<5)
					Flag[i] = -2;  //  possible last node

				int l=1;
				while(l<exon_num1-1)
				{
					int p1 = gene.trans[j].exon_start[l];
					int p2 = gene.trans[j].exon_stop[l];
					if(abs(s1-p1)<5&&abs(s2-p2)<5)
					{
						Flag[i]++; 
						break;
					}
					l++;
				}
			}
		}
	}
#pragma endregion 

#pragma region constraints on the nodes according to the annotations
	for(int i = 1; i<exon_Num-1; i++)
	{
		if(Flag[i]>0)
		{
			link_Mtx[0][i+1] = 0;
			link_Mtx[i+1][vtx_Number-1] = 0;
		}
		else if(Flag[i]==-1) 
			link_Mtx[0][i+1] = 1;  
		else if(Flag[i]==-2) 
			link_Mtx[i+1][vtx_Number-1] = 1;
	}
#pragma endregion

#pragma region
	int out_start = 0, in_end = 0;
	for(int i=1; i<vtx_Number-1; i++)
	{
		int temp = link_Mtx[0][i];
		out_start += temp;
		in_Degree[i] += temp;
		if(temp>0)
		{
			pos_Nodes[0].push_back(i);
			pre_Nodes[i].insert(pre_Nodes[i].begin(),0);
		}
	}
	out_Degree[0] = out_start;

	for(int i=1; i<vtx_Number-1; i++)
	{
		int temp = link_Mtx[i][vtx_Number-1];
		in_end += temp;
		out_Degree[i] += temp;
		if(temp>0)
		{
			pre_Nodes[vtx_Number-1].push_back(i);
			pos_Nodes[i].push_back(vtx_Number-1);
		}
	}
	in_Degree[vtx_Number-1] = in_end;
#pragma endregion

	return true;
}

// Enhance the constraints for being the first node or the last node 
bool Graph_Trans::graph_Constraint()
{
	int exon_Num = vtx_Number - 2;
	int in_degree = 0, out_degree = 0;

	for(int j = 0; j<exon_Num; j++)
	{
		if(in_Flow[j]>aver_in&&out_Flow[j]<0.5*in_Flow[j])   // in-flow is significantly larger than the out-flow
			link_Mtx[0][j+1] = 0;    // low possiblity of being the first node
		if(out_Flow[j]>aver_out&&in_Flow[j]<0.5*out_Flow[j]) // out-flow is significantly larger than the in-flow
			link_Mtx[j+1][vtx_Number-1] = 0;  // low possiblity of being the last node
		if(out_Flow[j]>0.8*aver_out&&in_Flow[j]>0.8*aver_in&&abs(1-in_Flow[j]/(0.1+out_Flow[j]))<0.25) // the flow is balanced
		{
			link_Mtx[0][j+1] = 0;
			link_Mtx[j+1][vtx_Number-1] = 0;
		}
	}

#pragma region update linkages to the first node and the last node
	for(int j = 1; j<vtx_Number-1; j++)  // no constraints on the first three nodes
	{
		int num = pre_Nodes[j].size();
		if(num>0&&link_Mtx[0][j]==0&&pre_Nodes[j][0]==0)  // delete linkage to the first node 
		{
			pre_Nodes[j].erase(pre_Nodes[j].begin());     // delete the first node
			in_Degree[j]--;
		}
		num = pos_Nodes[j].size();
		if(num>0&&link_Mtx[j][vtx_Number-1]==0&&pos_Nodes[j][num-1]==vtx_Number-1)  // delete linkage to the last node
		{
			pos_Nodes[j].pop_back();  // delete the last node
			out_Degree[j]--;
		}
	}

	int out_start = 0, in_end = 0;	
	pos_Nodes[0].clear();   // reconsider the first node and the last node
	pre_Nodes[vtx_Number-1].clear();
	for(int i=1; i<vtx_Number-1; i++)
	{
		int temp = link_Mtx[0][i];
		out_start += temp;
		if(temp>0)
		{
			pos_Nodes[0].push_back(i);
		}
	}
	out_Degree[0] = out_start;
	for(int i=1; i<vtx_Number-1; i++)
	{
		int temp = link_Mtx[i][vtx_Number-1];
		in_end += temp;
		if(temp>0)
		{
			pre_Nodes[vtx_Number-1].push_back(i);
		}
	}
	in_Degree[vtx_Number-1] = in_end;
#pragma endregion

	return true;
}

// Enhance the constraints for being the first node or the last node 
bool Graph_Trans::graph_Constraint_1(vector<float>& in_Flow, vector<float>& out_Flow)
{
	int exon_Num = vtx_Number - 2;
	float aver_in = 0, aver_out = 0;
	int in_degree = 0, out_degree = 0;
#pragma region in-degree and out-degree of the nodes
	for(int i=1; i<vtx_Number-1; i++)
	{
		in_degree = 0; out_degree = 0;
		float in_flow = 0, out_flow = 0;
		for(int j=1; j<i; j++)
		{
			if(link_Mtx[j][i]>0)
			{
				in_flow += link_Mtx[j][i];
				in_degree++;
				pre_Nodes[i].push_back(j);
			}
		}
		for(int j=i+1; j<vtx_Number-1; j++)
		{
			if(link_Mtx[i][j]>0)
			{
				out_flow += link_Mtx[i][j];
				out_degree++;
				pos_Nodes[i].push_back(j);
			}
		}
		in_Degree[i] = in_degree;
		out_Degree[i] = out_degree;
		aver_in += in_flow;
		aver_out += out_flow;
		in_Flow[i-1] = in_flow; out_Flow[i-1] = out_flow;
	}
#pragma endregion

	aver_in /= vtx_Number-2;
	aver_out /= vtx_Number-2;

	for(int j = 0; j<exon_Num; j++)
	{
		if(in_Flow[j]>aver_in&&out_Flow[j]<0.5*in_Flow[j])   // in-flow is significantly larger than the out-flow
			link_Mtx[0][j+1] = 0;   // low possiblity of being the first node
		if(out_Flow[j]>aver_out&&in_Flow[j]<0.5*out_Flow[j]) // out-flow is significantly larger than the in-flow
			link_Mtx[j+1][vtx_Number-1] = 0;  // low possiblity of being the last node
		if(out_Flow[j]>0.8*aver_out&&in_Flow[j]>0.8*aver_in&&abs(1-in_Flow[j]/(0.1+out_Flow[j]))<0.25) // the flow is balanced
		{
			link_Mtx[0][j+1] = 0;
			link_Mtx[j+1][vtx_Number-1] = 0;
		}
	}

	link_Mtx[0][1] = 1;
	link_Mtx[vtx_Number-2][vtx_Number-1] = 1;

	int out1 = 0, in_end = 0;
	for(int i=1; i<vtx_Number-1; i++)
	{
		int temp = link_Mtx[0][i];
		out1 += temp;
		in_Degree[i] += temp;
		if(temp>0)
		{
			pos_Nodes[0].push_back(i);
			pre_Nodes[i].insert(pre_Nodes[i].begin(),0);
		}
	}
	out_Degree[0] = out1;

	for(int i=1; i<vtx_Number-1; i++)
	{
		int temp = link_Mtx[i][vtx_Number-1];
		in_end += temp;
		out_Degree[i] += temp;
		if(temp>0)
		{
			pre_Nodes[vtx_Number-1].push_back(i);
			pos_Nodes[i].push_back(vtx_Number-1);
		}
	}
	in_Degree[vtx_Number-1] = in_end;

	
	return true;
}





