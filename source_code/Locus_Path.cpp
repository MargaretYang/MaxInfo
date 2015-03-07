
#include "Read_Assignment.h"
using namespace std;

// Calculate the cost of the exon for being the first node or the last node of the path
void Graph_Trans::cost_FirstLastNode()
{
	int e_Num = vtx_Number-2;
	float cnt1 = 0, cnt2 = 0, cnt = 0;
	first_NodeCost = new double[e_Num];
	last_NodeCost = new double[e_Num];

	for(int i = 1; i<=e_Num; i++)
	{
		for(int j = 1; j<i; j++)
			cnt += link_Mtx[j][i];
	}
	double mean_Cnt = cnt*1.0/e_Num;
	// cout<<"Average input flow on each node: "<<mean_Cnt<<endl;
	double eps = 1e-09;

	for(int i = 1; i<=e_Num; i++)
	{
		cnt1 = 0, cnt2 = 0;
		for(int j = 1; j<i; j++)
			cnt1 += link_Mtx[j][i];  // in-flow
		for(int j = i+1; j<=e_Num; j++)  
			cnt2 += link_Mtx[i][j];  // out-flow

		double rate = 1.0*(cnt1-cnt2)/(mean_Cnt+eps);
		// when in-low is smaller than out-flow, the cost for being the first node is reduced
		first_NodeCost[i-1] = exp(rate);
		// when in-low is larger than out-flow, the cost for being the first node is reduced
		last_NodeCost[i-1] = exp(-rate);
	}
}

// Search path by DFS(Depth First Search)
bool Graph_Trans::path_DFS(bool* mark)
{
	int index = 0;
	vector<int> p_Stack;  // stack for saving nodes of the path

	p_Stack.push_back(0); // save the source node
	mark[0] = true;
	bool flag = true;
	int cnt = 0;

	while(!p_Stack.empty())
	{
		int top = p_Stack[p_Stack.size()-1];
		if(top==vtx_Number-1)   // the top element of the stack is the sink node
		{
			int num = p_Stack.size();  // length of the path
			vector<int> cur_Path;      // new path
			for(int i = 1;i<num-1; i++)
			{
				cur_Path.push_back(p_Stack[i]);
			}
			paths.push_back(cur_Path);
			p_Stack.pop_back();    // delete the top element
			// mark[top] = false;  // set the status to being unvisited
		}
		else
		{
			bool neigh = true;    // flag indicating whether there are unvisited nodes
			int i = top+1;			
			while(i<vtx_Number&&(link_Mtx[top][i]==0||mark[i]==true)) // find the next unvisited linked node
				i++;
			if(i<vtx_Number)
			{
				neigh = false;
				p_Stack.push_back(i);   
				mark[i] = true;				
			}

			if(neigh)  // all the neighbor nodes have been visited
			{
				p_Stack.pop_back();  // delete the top element
				mark[top] = true;    // set the status to being visited
				for(int j = top+1; j<vtx_Number; j++) 
				{
					if(link_Mtx[top][j]>0)
						mark[j] = false; // set the status to the neighbor node to being unvisited
				}
			}
		}
	}

	return true;
}

// Search path by DFS(Depth First Search): Two-Stack method
bool Graph_Trans::path_DFS1()
{
	int index = 0;
	paths.clear();
	vector<int> p_Stack;   // stack for saving nodes of the path
	vector<int> idx_Stack; // stack for saving the indices of the nodes

	p_Stack.push_back(0);  // saving the source node
	idx_Stack.push_back(1);
	bool flag = true;
	int cur_Idx = 0;
	int path_Num = 0;
	int thresh = 1500;

	while(!p_Stack.empty())
	{
		int top = p_Stack[p_Stack.size()-1];

		if(top==vtx_Number-1)   // the top element of the stack is the sink node
		{
			int num = p_Stack.size();
			vector<int> cur_Path;
			for(int i = 1;i<num-1; i++)
			{
				cur_Path.push_back(p_Stack[i]);
			}
			paths.push_back(cur_Path);
			path_Num++;
			if(path_Num>thresh)
			{
				flag = false;
				break;
			}
			p_Stack.pop_back();     // delete the top element
			idx_Stack.pop_back();   // delete index of the sink node
			// mark[top] = false;   // set the status to being visited
		}
		else
		{
			bool neigh = true;      // flag indicating whether there are unvisited nodes
			cur_Idx = idx_Stack[idx_Stack.size()-1];

			int i = cur_Idx;			
			while(i<vtx_Number&&(link_Mtx[top][i]==0))  // find the next unvisited linked node
				i++;
			if(i<vtx_Number)  // the neighbor node hasn't visited
			{
				neigh = false;
				cur_Idx = i+1;  // change the index
				idx_Stack[idx_Stack.size()-1] = cur_Idx;
				p_Stack.push_back(i);
				idx_Stack.push_back(i+1);				
			}

			if(neigh)  // all the neighbor nodes have been visited
			{
				p_Stack.pop_back();    // delete the top element of the node stack
				idx_Stack.pop_back();  // delete the top element of the index stack		
			}
		}
	}

	return flag;
}

// Path search with adaptive constraints
bool Graph_Trans::path_DFS1_Adaptive()
{
	graph_Constraint();
	bool search_Flag = path_DFS1_constraint();
	int thresh = 2;

	if(search_Flag==false)   // there are too many paths
	{
		cout<<"Linkage to be modified"<<endl;
		int path_Num = paths.size();
		for(int i=0; i<path_Num; i++)
			vector<int>().swap(paths[i]);
		vector<VEC_INT>().swap(paths);

#pragma region Modify the constraints for being the first and last node 
		for(int j = 4; j<vtx_Number-1; j++)  // no constraints on the first three nodes
		{
			if(link_Mtx[0][j]==1&&(pre_Nodes[j].size()>=thresh))
			{
				link_Mtx[0][j] = 0;
				pre_Nodes[j].erase(pre_Nodes[j].begin());  // delete the first node
				in_Degree[j]--;
			}
		}
		for(int j = 1; j<vtx_Number-4; j++)  // no constraints on the first three nodes
		{
			if(link_Mtx[j][vtx_Number-1]==1&&(pos_Nodes[j].size()>=thresh))
			{
				link_Mtx[j][vtx_Number-1] = 0;
				pos_Nodes[j].pop_back();  // delete the first node
				out_Degree[j]--;
			}
		}

		cout<<"Linkage modified"<<endl;

		search_Flag = path_DFS1_constraint(); // search the path again

		cout<<"Path re-searched"<<endl;
		if(search_Flag==false)  // select 100 paths
		{
			int path_Num = paths.size();
			int interval = int(path_Num/50)-2;
			vector<VEC_INT> temp_paths;
			for(int i=0; i<50; i++)  // select the first 50 paths
			{ 
				temp_paths.push_back(paths[i]);
			}
			for(int i=0; i<50; i++)  // randomly select 50 paths
			{ 
				temp_paths.push_back(paths[i*interval+interval]);
			}
			for(int i=0; i<path_Num; i++)
				vector<int>().swap(paths[i]);
			vector<VEC_INT>().swap(paths);
			paths = temp_paths;
			vector<VEC_INT>().swap(temp_paths);
			search_Flag = true;
		}		
#pragma endregion
	}
	// clear the arrays
	vector<VEC_INT>().swap(pre_Nodes);
	vector<VEC_INT>().swap(pos_Nodes);
	vector<float>().swap(in_Flow);
	vector<float>().swap(out_Flow);

	return search_Flag;
}

// Path search with constraints when there are too many possible paths
bool Graph_Trans::path_DFS1_constraint()
{
	int path_num = 0;      // number of paths
	int index = 0;
	paths.clear();
	vector<int> p_Stack;   // stack for saving the nodes of the path
	vector<int> idx_Stack; // stack for saving the indices of the nodes
	vector<vector<VEC_INT>> node_Paths;  // all the paths from the current node to the source node

	p_Stack.push_back(0);  // save the source node
	idx_Stack.push_back(1);
#pragma region initialize the paths
	for(int i=0; i<vtx_Number; i++)
	{
		vector<bool> mark1(pre_Nodes[i].size(),false);
		vector<bool> mark2(pos_Nodes[i].size(),false);
		pre_Linked.push_back(mark1); pos_Linked.push_back(mark2);
	}

	bool flag = true;
	int cur_Idx = 0;
	int exon_Num = vtx_Number - 2;
	vector<bool> visited_Flag(vtx_Number);  // flag of visit
	vector<int> temp_outDegree = out_Degree;
	vector<int> temp_inDegree = in_Degree;
	for(int i=0; i<vtx_Number; i++)
	{
		visited_Flag[i] = false;
		vector<VEC_INT> temp;
		node_Paths.push_back(temp);
	}

	for(int i=1; i<vtx_Number; i++)
	{
		if(out_Degree[i]==1&&link_Mtx[i][vtx_Number-1]==1)
		{
			vector<int> temp; temp.push_back(i);
			node_Paths[i].push_back(temp);
			visited_Flag[i] = true;
		}
	}

	vector<int> temp;
	node_Paths[vtx_Number-1].push_back(temp);
	visited_Flag[vtx_Number-1] = true;
#pragma endregion

	vector<int> node_Bak;
	int thresh = 1500;
	while(!p_Stack.empty()&&flag==true)
	{
		int top = p_Stack[p_Stack.size()-1];
		if(top==vtx_Number-1)   // top element is the sink node
		{
			int num = p_Stack.size();   // path length
			vector<int> cur_Path; // new path
			for(int i=1; i<num-1; i++)
			{
				cur_Path.push_back(p_Stack[i]);
			}
			paths.push_back(cur_Path);
			p_Stack.pop_back();     // delete the top element of the node stack
			idx_Stack.pop_back();   // delete the top element of the index stack
			path_num++;
			if(path_num>thresh)    
			{
				flag = false;
				break;
			}
		}
		else if(top==-1)  // reach a visited node and stop the search
		{
#pragma region
			int visited_idx = node_Bak[node_Bak.size()-1];
			int num = p_Stack.size();
			vector<int> cur_Path;
			for(int i=1; i<num-1; i++)
			{
				cur_Path.push_back(p_Stack[i]);
			}
			vector<VEC_INT>& path_Set = node_Paths[visited_idx];
			int p_Num = path_Set.size();
			for(int j=0; j<p_Num; j++)
			{
				vector<int> cur_Path1 = cur_Path;
				for(int k=0; k<path_Set[j].size(); k++)
				{
					cur_Path1.push_back(path_Set[j][k]);
				}
				paths.push_back(cur_Path1);
				path_num++;
				if(path_num>thresh)
				{
					flag = false;
					break;
				}
			}

#pragma region
			int idx = p_Stack[num-2];
			int l1=0, l2=0;
			while(pos_Nodes[idx][l2]!=visited_idx)
				l2++;
			while(pre_Nodes[visited_idx][l1]!=idx)
				l1++;
			if(pos_Linked[idx][l2]==false)
			{
				for(int j=0; j<p_Num; j++)
				{
					vector<int> cur_Path2 = path_Set[j];
					cur_Path2.insert(cur_Path2.begin(),idx);
					node_Paths[idx].push_back(cur_Path2);
				}		
				temp_outDegree[idx]--;
				pos_Linked[idx][l2] = true;
			}
			if(pre_Linked[visited_idx][l1]==false)
			{
				temp_inDegree[visited_idx]--;
				pre_Linked[visited_idx][l1] = true;
			}
#pragma endregion				
			if(temp_inDegree[visited_idx]==0)  // all the previous nodes have been visited
				vector<VEC_INT>().swap(node_Paths[visited_idx]);  // delete records of the paths
			if(temp_outDegree[idx]==0)         // all the next nodes have been visited
				visited_Flag[idx] = true;
			node_Bak.pop_back();
			p_Stack.pop_back(); 
			idx_Stack.pop_back();
#pragma endregion
		} 
		else
		{
			bool neigh = true;
			cur_Idx = idx_Stack[idx_Stack.size()-1];

			int i = cur_Idx;			
			while(i<vtx_Number&&(link_Mtx[top][i]==0))
				i++;

			if(i<vtx_Number)
			{
				neigh = false;
				if(visited_Flag[i]==true)  // the node is already visited
				{
					idx_Stack[idx_Stack.size()-1] = i+1;
					p_Stack.push_back(-1);
					node_Bak.push_back(i);
					idx_Stack.push_back(i+1);  // the next node to visit				
				}
				else  // the next node is unvisited
				{
					cur_Idx = i+1;     // change the index
					idx_Stack[idx_Stack.size()-1] = cur_Idx;
					idx_Stack.push_back(i+1);
					p_Stack.push_back(i);
				}				
			}
			if(neigh==true)  // all the neighbor nodes are already visited
			{
				p_Stack.pop_back();    // delete the top element of the node stack
				idx_Stack.pop_back();  // delete the top element of the index stack
			}
		}
	}	

#pragma region clear the arrays
	vector<vector<bool>>().swap(pre_Linked);
	vector<vector<bool>>().swap(pos_Linked);
	vector<bool>().swap(visited_Flag);
	vector<int>().swap(p_Stack);
	vector<int>().swap(idx_Stack);
	vector<int>().swap(node_Bak);
	int p_Size = node_Paths.size();  // number of paths
	for(int i=0; i<p_Size; i++)
	{
		vector<VEC_INT>().swap(node_Paths[i]);
	}
	vector<vector<VEC_INT>>().swap(node_Paths);
	vector<int>().swap(temp_inDegree);   // clear the vector
	vector<int>().swap(temp_outDegree);  // clear the vector
#pragma endregion

	return flag;
}

// Obtain information of the paths including path structures and path entropies
void Graph_Trans::path_Speci(vector<Path> &path)
{
	int p_Num = paths.size();  // number of paths
	int e_Serial = 0;
	float thresh_firstnode = 1, thresh_lastnode = 1;  // costs of the first node and the last node
	float thresh_Score = 0.8;
	int exon_Num = exon_Vec.size();
	if(exon_Num<=3)   // modify the costs of the first node and the last node when there are too few exons
	{
		thresh_firstnode = 1.2;  thresh_lastnode = 1.2;
	}

	for(int i = 0; i<p_Num; i++)
	{
		int subNum = paths[i].size();    // number of exons in a certain path
		int id1 = paths[i][0]-1;         // first node of the path
		int id2 = paths[i][subNum-1]-1;  // last node of the path
		double cost1 = first_NodeCost[id1];
		double cost2 = last_NodeCost[id2];

		// select transcripts according to the costs of the first node and the last node
		// and the generation probability of the reads
		if(cost1<thresh_firstnode&&cost2<thresh_lastnode)
		{
			int length = 0;
			Path cur_Path;
			cur_Path.exon_Number = subNum;
			cur_Path.exon_Node = new int[subNum];
			cur_Path.bound = new int[subNum*2];
			for(int k = 0; k<subNum; k++)
			{
				e_Serial = paths[i][k]-1;          // serial of the exon
				cur_Path.exon_Node[k] = e_Serial;  // exons contained in the path
				cur_Path.bound[2*k] = exon_Vec[e_Serial].bound1;    // left boundary of the exon
				cur_Path.bound[2*k+1] = exon_Vec[e_Serial].bound2;  // right boundary of the exon
				length += exon_Vec[e_Serial].len;
			}
			cur_Path.length = length;
			cur_Path.entro_Abun1 = 0; 
			path.push_back(cur_Path);  // save information of the transcript
		}		
	}

#pragma region add annotated transcripts
	int num = cur_Gene.trans.size();  // number of annotated transcripts
	int num_1 = path.size();          // number of current estimated transcripts
	exon_Num = exon_Vec.size();
	annotated_Idx.clear();
	for(int i=0; i<num; i++)
	{
		Trans& cur_Trans = cur_Gene.trans[i];
		int exon_num = cur_Trans.exon_start.size();
		int k1 = 0;
		while(k1<exon_num-1&&abs(cur_Trans.exon_start[k1+1]-cur_Trans.exon_stop[k1])<2)
			k1++;
		int t1 = 0, t2 = cur_Trans.exon_stop[k1];
		int s1 = 0, s2 = 0;
		bool match_Flag = false;
		for(int j=0; j<num_1; j++)
		{
			int exon_num_1 = path[j].exon_Number;
			int k2 = 0;
			while(k2<exon_num_1-1&&abs(path[j].bound[2*k2+2]-path[j].bound[2*k2+1])<2)
				k2++;
			int p1 = 0, p2 = path[j].bound[2*k2+1]; 
			if(abs(p2-t2)>3)  // right boundaries of the first exon don't match
				continue;
			int k = k1+1, l = k2+1;
			while(k<exon_num-1)
			{
				s1 = cur_Trans.exon_start[k];
				s2 = cur_Trans.exon_stop[k];
				if(l<exon_num_1-1)
				{
					p1 = path[j].bound[2*l];    // left boundary
					p2 = path[j].bound[2*l+1];  // right boundary
					if(abs(p1-s1)>=3||abs(p2-s2)>=3)
						break;
					l++;
				}
				k++;
			}
			if(k<exon_num-1||l<exon_num_1-1)  // the transcript doesn't meet the requirement
				continue;  // consider the next transcript
			else
			{
				p1 = path[j].bound[2*exon_num_1-2];   // left boundary of the last exon
				s1 = cur_Trans.exon_start[exon_num-1];
				if(abs(p1-s1)<=3)
				{
					match_Flag = true;           // the estimated transcript and the annotated transcript are matched
					annotated_Idx.push_back(j);  // serial of the annotated transcript
					break;
				}
			}
		}
		if(match_Flag==false)
		{
			int k=0;
			int s2 = cur_Trans.exon_stop[0];
			int s1 = cur_Trans.exon_start[0];

#pragma region serial of exon in the path
			vector<int> path_Serial;  // serial of exon in the path
			while(k<exon_Num)
			{
				if(abs(s2-exon_Vec[k].bound2)<3)
					break;
				k++;
			}
			int temp = 0;
			if(k<exon_Num)
			{
				for(int l=0; l<k; l++)
				{
					if(abs(s1-exon_Vec[l].bound1)<3)  // the left boundaries coincide
					{
						temp = l;
						break;
					}
				}
				for(int l=temp; l<k; l++)
					path_Serial.push_back(l);
				path_Serial.push_back(k);

				int l1 = 1;
				while(l1<exon_num-1)
				{
					for(int l=k+1; l<exon_Num-1; l++)
					{
						s1 = cur_Trans.exon_start[l1];    
						s2 = cur_Trans.exon_stop[l1];     
						if(abs(s2-exon_Vec[l].bound2)<3)  // the right boundaries coincide
						{
							for(int l2=k+1; l2<=l; l2++)
							{
								if(abs(s1-exon_Vec[l2].bound1)<3)  // the left boundaries coincide
									path_Serial.push_back(l2);
							}
							k = l;
						}
						else if(abs(s1-exon_Vec[l].bound1)<3)      // the left boundaries coincide
						{
							int l3 = l;
							for(l3=l; l3<exon_Num-1; l3++)
							{
								if(abs(s2-exon_Vec[l3].bound2)<3)  // the right boundaries coincide
									path_Serial.push_back(l3);
							}
							k = l3;
						}
					}
					l1++;
				}
				if(l1==exon_num-1)
				{
					s1 = cur_Trans.exon_start[l1];  // left boundary of the first exon
					k++;
					while(k<exon_Num)
					{
						if(abs(s1-exon_Vec[k].bound1)<3)
							break;
						k++;
					}
					if(k<exon_Num)
						path_Serial.push_back(k);
					for(int l=k+1; l<exon_Num; l++)
					{
						if(abs(s2-exon_Vec[l].bound2)<3)
						{
							temp = l;
							break;
						}
					}
					for(int l=k+1; l<=temp; l++)
						path_Serial.push_back(l);
				}
			}
#pragma endregion

			Path cur_Path;
			int length = 0;
			int temp_Num = path_Serial.size();
			cur_Path.exon_Number = temp_Num;
			cur_Path.exon_Node = new int[temp_Num];
			cur_Path.bound = new int[temp_Num*2];
			for(int k = 0; k<temp_Num; k++)
			{
				int e_Serial_1 = path_Serial[k];
				cur_Path.exon_Node[k] = path_Serial[k];   // serial of the exon contained in the path
				cur_Path.bound[2*k] = exon_Vec[e_Serial_1].bound1;   
				cur_Path.bound[2*k+1] = exon_Vec[e_Serial_1].bound2;
				length += exon_Vec[e_Serial_1].len;
			}
			cur_Path.length = length;
			cur_Path.entro_Abun1 = 0;
			path.push_back(cur_Path);  // save information of the transcript
			annotated_Idx.push_back(num_1);
			num_1++;
		}
	}
#pragma endregion

}

// Output the paths
void Graph_Trans::print_Path()
{
	int num = paths.size();  // number of paths
	std::cout<<"number of paths: "<<num<<endl;
	for(int i = 0; i<num; i++)
	{
		int subNum = paths[i].size();  // exon contained by the path
		for(int k = 0; k<subNum; k++)
		{
			std::cout<<paths[i][k]<<" ";
		}
		std::cout<<endl;
	}
}
