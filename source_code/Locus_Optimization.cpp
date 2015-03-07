//////////////////////////////////////////////////////////////////////////
//// TRANSCRIPT ASSEMBLY AND ABUNDANCE ESTIMATION BASED ON RNA-SEQ DATA

#include "Entropy.h"
#include "Batch_Mode.h"
#include "ISA.h"

// Exon Assembly
bool Gene_Batch::exon_AssembleBatch(ifstream& infile, ifstream& infile2, vector<int>& exon_Start,vector<int>& exon_Stop)
{
#pragma region Declaration of variables
	exon_Start.clear(); exon_Stop.clear();
	junction.clear();
	vector<int> exon_Len;   // length of exon

	string line, word, seg[15];
	string::size_type pos, pre;	
	char tPar = ',';	
	int jump_Thresh = 3200;      // threshold of splice site distance
	int jump_Thresh_1 = 2800;    //threshold of splice site distance in the reverse orientation
	int line_Cnt = 0;         
	int exon_start = 0, exon_stop = 1e09, exon_len = 0, temp = 0;
	int junc_start = 0, junc_stop = 1e09, junc_num = 0;
	int p1 = 0, p2 = 0;
	int m_exonStart = 0, m_exonStop = 0;
	string orientation = "";     // orientation of read alignment
	int reverse_cnt = 0;         // counter for the change of orientation 
	bool locate_Flag = true;     // flag for read location
	int reverse_Flag = 2;        // flag for change of orientation
	bool changed_Flag = false;   // flag for change of chromosome
	int intra_Flag = -1;         // flag for intra-overlap of gene loci 
	int overlap_Flag = -1;       // flag for overlop of gene loci  marker for intra-overlap 
	std::streampos pos1 = infile.tellg();   // current file pointer 
	vector<VEC_INT> reverse_Pair;           // junction-end pairs in the reverse orientation
	vector<int> reverse_Left, reverse_Right;   
	reverse_Pair.push_back(reverse_Left);
	reverse_Pair.push_back(reverse_Right);

	intra_activeFlag = false;    // flag for overlap of gene loci
	vector<int>().swap(rev_breakPoints);
	start_sites.clear();   
#pragma endregion
	stringstream stream(pre_line);
	for(int i=0; i<12; i++)
		stream>>seg[i];          // read the fields of each line of the alignment file

	sscanf(seg[1].c_str(),"%d",&junc_start);  // starting point of junction
	sscanf(seg[2].c_str(),"%d",&junc_stop);   // stopping point of junction
	sscanf(seg[4].c_str(),"%d",&junc_num);    // number of junction-spanning reads
	exon_Start.push_back(0);   

#pragma region load in the ends of exons in the stack
	bool stack_Flag0 = false, stack_Flag1 = false;
	if(stack_exonBoundary.size()>0)
	{
		int num1 = stack_exonBoundary[1].size();
		if(num1>0)
		{
			stack_Flag1 = true;
			exon_Start = stack_exonBoundary[1];
			m_exonStart = exon_Start[num1-1];
		}
		int num2 = stack_exonBoundary[0].size();
		if(num2>0)
		{
			stack_Flag0 = true;
			exon_Stop = stack_exonBoundary[0]; 
			m_exonStop = exon_Stop[num2-1];
		}
		vector<int>().swap(stack_exonBoundary[0]);
		vector<int>().swap(stack_exonBoundary[1]);
		stack_exonBoundary.clear();
	}
	stack_exonBoundary.clear();
	vector<int> exon_stop1, exon_start1;
	stack_exonBoundary.push_back(exon_stop1);   // save stopping point of exon
	stack_exonBoundary.push_back(exon_start1);  // save starting point of exon
#pragma endregion

	ending_site1 = 0;
	ending_site2 = 0;
	double spr_Thresh = 200000;
	double spr_Thresh_high = 100000;
	bool spr_Flag = true;

	while(!infile.eof())
	{
		if((line_Cnt==0||intra_Flag==1))
		{
			line = pre_line;
		}

		intra_Flag = -1;
		stringstream stream(line);
		for(int i=0; i<12; i++)
			stream>>seg[i];   // read the 12 fileds of each line of the BED file for junction recognition 

		if(seg[0]!=chromosome_Name)
		{
			next_chromoName = seg[0];  // record change of chromosome
			pre_line = line;
			record_file<<"Chromosome changed: "<<chromosome_Name<<endl;
			changed_Flag = true;
			break;
		}

		sscanf(seg[1].c_str(),"%d",&junc_start);  // starting point of junction 
		sscanf(seg[2].c_str(),"%d",&junc_stop);   // stopping point of junction
		sscanf(seg[4].c_str(),"%d",&junc_num);    // number of junction-spanning reads

#pragma region delimit the boundaries of exons according to the annotations
		pre = 0;
		if((pos=seg[10].find(tPar,pre))!=string::npos)  // search the splitting charater ','
		{
			string str = seg[10].substr(pre,pos-pre);
			sscanf(str.c_str(),"%d",&temp);
			exon_stop = junc_start + temp;        // stopping point of exon
		}

		pre = 0;
		if((pos=seg[11].find(tPar,pre))!=string::npos)  //  search the splitting charater ','
		{
			int len = seg[11].length();
			string str = seg[11].substr(pos+1,len-pos-1);
			sscanf(str.c_str(),"%d",&temp);       // exon length
			exon_start = junc_start + temp + 1;
		}
#pragma endregion

		//////////////////////////////////////////////////////////////////////////
		// error control
		// spr_Thresh = 500000;
		int jump_Distance = abs(exon_start-exon_stop);
		if(jump_Distance>spr_Thresh)  // too large distance between adjacent splice sites
		{
			sus_Junc.push_back(make_pair(exon_stop,exon_start));  // save suspicious error exons 
			getline(infile,line);
			line_Cnt++;
			continue;
		}
		else if(jump_Distance>spr_Thresh_high)  // too large distance between adjacent splice sites
		{
			if(junc_num<3)
			{
				sus_Junc.push_back(make_pair(exon_stop,exon_start));			
				getline(infile,line);
				line_Cnt++;
				continue;
			}			
		}

		if(line_Cnt==0)
		{
			if(pre_orientation!="")
			{
				orientation = pre_orientation;
				pre_orientation = "";
			}
			else
				orientation = seg[5];
			if(stack_Flag1==false)
			{
				m_exonStart = exon_start;
			}
		}
		line_Cnt++;		

		// record the junctions
		junction.push_back(make_pair(exon_stop,exon_start));

		if(seg[5]!=orientation)    // change of orientation
		{
#pragma region cases of change of orientation
			reverse_Flag = true;
			if(exon_stop-m_exonStart>jump_Thresh_1)  // change of gene locus
			{
				pre_line = line;   // record the current line
				break;
			}
			else if(exon_stop-m_exonStart<0)  // overlap of gene loci
			{
				int temp = -1;
				if(intra_activeFlag==false)   // finishing the record of the last splice site
				{
#pragma region record possible starting points of exons
					int n1 = exon_Stop.size(), n2 = exon_Start.size(); // the number of exons
					int temp1 = exon_stop-500;

					if(exon_stop-m_exonStop>50)
						temp = _max(temp1,m_exonStop);
					else if(n2>1)
					{
						int t2  = m_exonStop;
						int l = 0;
						for(l=n2-2; l>=0; l--)  // locate the previous exon
						{
							if(exon_Start[l]<m_exonStop)
								break;
						}
						if(l>=0&&exon_Start[l]>0)
							t2 = _max(temp1,(int)(m_exonStop+exon_Start[l])/2);  // locate the middle point of the previous exon
						temp = _min(t2,exon_stop-100);						
					}
					else
						temp = _max(temp1,0);
					start_sites.push_back(temp);  // record possible starting points
#pragma endregion
				}				
				intra_Flag = -1;
				intra_activeFlag = true;      // marked as active to identify intra-overlap of gene loci
				if(rev_breakPoints.size()>0)  // there are splice sites which were already used for computation
				{
					vector<VEC_INT> reverse_Pair1; 
					vector<int> reverse_Left, reverse_Right;
					reverse_Pair1.push_back(reverse_Left);
					reverse_Pair1.push_back(reverse_Right);
					reverse_Pair1[0].push_back(exon_stop);  // save the exon-end pairs
					if(temp!=-1)
					{
						reverse_Pair1[1].push_back(temp);
					}
					else
						reverse_Pair1[1].push_back(exon_stop-200);
					reverse_Pair1[1].push_back(exon_start);
					intra_Flag = intra_Overlap(infile,orientation,m_exonStart,m_exonStop,exon_stop,exon_start,reverse_Pair1,exon_Start,exon_Stop);
					rev_intraPairs.push_back(reverse_Pair1);
				}
				else
				{
					reverse_Pair[0].push_back(exon_stop);  // save the exon-end pairs
					reverse_Pair[1].push_back(exon_start);
					intra_Flag = intra_Overlap(infile,orientation,m_exonStart,m_exonStop,exon_stop,exon_start,reverse_Pair,exon_Start,exon_Stop);
				}
				if(intra_Flag==2)    // change to the next gene or the next chromosome
				{
					break;
				}
				else if(intra_Flag==1)
				{
					continue;
				}
				else if(intra_Flag==0)
				{
					changed_Flag = true;
					break;
				}
			}
			else  // there is possibly gene overlap
			{
				overlap_Flag = -1;
				pre_orientation = seg[5];  // record the opposite orientation
				vector<VEC_INT>& reverse_Pair1 = reverse_Pair;
				int n1 = rev_intraPairs.size(); 

				if(n1>0)
					reverse_Pair1 = rev_intraPairs[n1-1];
				overlap_Flag = gene_Overlap(infile,orientation,m_exonStart,m_exonStop,exon_stop,exon_start,reverse_Pair1, exon_Start,exon_Stop);

				if(overlap_Flag==2)  // there is no gene overlap, but instead change of gene or chromosome
				{
					break;  // there is no gene overlap, but instead change of gene or chromosome
				}
				else if(overlap_Flag==1) // there is gene overlap
				{
					int s1 = m_exonStart, s2 = m_exonStop;
					int num = reverse_Pair[0].size();  // number of stopping points of exons
					int temp = 0;
					overlap_Flag = -1;  // clear the flag
				}
				else if(overlap_Flag==0)
				{
					changed_Flag = true;   // change of chromosome
					break;
				}
			}
#pragma endregion
		}
		else  // there is no change of orientation
		{
#pragma region
			if(exon_stop>m_exonStop)
				m_exonStop = exon_stop;  // update the right-end of exon

			if(m_exonStop-m_exonStart>jump_Thresh)    // distance between two adjacent splice sites is larger than the threshold
			{
				p2 = junc_start;     // the starting point of the current junction, i.e., the stopping point of exon
				locate_Flag = false;
				pre_line = line;     // record the current line
				break;
			}
			else
			{
				exon_Stop.push_back(exon_stop);
				exon_Start.push_back(exon_start);
			}

			if(exon_start>m_exonStart)
				m_exonStart = exon_start;
#pragma endregion 
		}
		getline(infile,line);
	}

	/*if(infile.eof())
		cout<<"Reach file end!"<<endl;*/

	exon_Type.clear(); exon_Score.clear();
	if(!exon_Start.empty())
	{	
		// ending_site1 = _max(ending_site1,exon_stop); 
		if((intra_Flag==-1&&overlap_Flag==-1)||intra_Flag==1||overlap_Flag==1)
			exon_Stop.push_back(exon_stop);
		else // intra_Flag=0,2 or overlap_Flag=0,2 , gene locus breaks in the reverse orientation
			exon_Stop.push_back(ending_site1);		
		exon_Predict(exon_Start, exon_Stop, starting_site, ending_site1, exon_Assembly, exon_Type, exon_Score);  // exon assembly
	}

#pragma region cases of gene overlap
	vector<Exon_Node>().swap(rev_exonAssembly);
	vector<float>().swap(rev_exonScore);
	vector<int>().swap(rev_exonType);
	int num1 = reverse_Pair[0].size();
	if(num1>0)   // there are gene-overlaps
	{
		int num2 = reverse_Pair[1].size();
		if(num1<num2)   // stopping point of exon need to be added to the array
		{
			reverse_Pair[0].push_back(reverse_Pair[1][num2-1]+1200);
		}
		int site = reverse_Pair[1][0];
		exon_Predict(reverse_Pair[1], reverse_Pair[0], site, ending_site2, rev_exonAssembly, rev_exonType, rev_exonScore);  // exon assembly
		rev_ExonAssembly.push_back(rev_exonAssembly);
		rev_ExonScore.push_back(rev_exonScore); rev_ExonType.push_back(rev_exonType);
	}
#pragma endregion

#pragma region cases of multiple gene-overlaps
	int intra_Num = rev_intraPairs.size();
	for(int i=0; i<intra_Num; i++)
	{
		vector<VEC_INT>& reverse_Pair1 = rev_intraPairs[i];
		VEC_NODE rev_exonAssembly;
		vector<float> rev_exonScore;
		vector<int> rev_exonType;		
		int num1 = reverse_Pair1[0].size(), num2 = reverse_Pair1[1].size();
		if(num1>0)   // there are gene-overlaps
		{
			if(num1<num2)   // stopping point of exon need to be added to the vector
			{
				reverse_Pair1[0].push_back(reverse_Pair1[1][num2-1]+1200);
			}
			exon_Predict(reverse_Pair1[1], reverse_Pair1[0], reverse_Pair1[1][0], ending_site2, rev_exonAssembly, rev_exonType, rev_exonScore);  // 外显子组装
		}
		rev_ExonAssembly.push_back(rev_exonAssembly);
		rev_ExonScore.push_back(rev_exonScore); rev_ExonType.push_back(rev_exonType);
	}
#pragma endregion

#pragma region clear the arrays and close the files
	vector<int>().swap(exon_Len);
	for(int i=0; i<intra_Num; i++)
	{
		vector<VEC_INT>().swap(rev_intraPairs[i]);	
	}
	rev_intraPairs.clear();
#pragma endregion

	return changed_Flag;
}

// Delimitation of exons
void Gene_Batch::coverage_Assemble(ifstream& infile, Rd_Set& rd_Set)
{
	if(!infile)
		cout<<"Open file error!"<<endl;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Format of BAM file：1.read name 2.flag for description of alignment 3.name of reference sequence 
	// 4.starting point of left read 5.mapping quality 6.description of mapping quality
	// 7.marker for paired matching or single matching 8.starting point of the other read 9.length of segment
	// 10.read sequence 11.read quality 12.other descriptions of the mapping
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma region Declaration of variables
	int x1 = 0, x2 = 0, z1 = 0, z2 = 0, temp_Len = 0; // end points of the paired-end reads; lenght of fragment
	int read_len = rd_Len;  // read length

	stringstream ss;
	ss<<read_len; 
	string len_R = ss.str()+"M"; // string indicating exact total match of the read
	string line, word, words[6], frag_Name;  
	int serial = 0;     // read serial 
	string tPar = "N";
	string::size_type pos, pos0;
	vector<int>::iterator itr;

	int junc_Len = 0, cover_Len = 0, p1 = 0, p2 = 0;
	int x_1 = 0, x_2 = 0, z_1 = 0, z_2 = 0;
	int pre_Left = 0, pre_Right = 0;
	bool anno_flag = false; 
	int line_cnt = 0, line_cnt1 = 0, line_cnt2 = 0;
	bool junc_Left = false, junc_Right = false, flag_Right = false; 
	vector<int> right_Boundary;      // right boundaries of fragments
	int acc_cnt1 = 0;
	int thresh1 = 20, thresh2 = 20, len_Thresh = 5;
	int num = exon_Assembly.size();  // number of exons

	int seg_Mark = 0;
	int exon_start = exon_Assembly[0].bound1, exon_stop = exon_Assembly[0].bound2;
	int pre_start = 0, pre_stop = 0;
	vector<int> seg_Start, seg_Stop, exon_Start, exon_Stop;   // starting points and stopping points of segments
	vector<bool> modified_Flag(2*num-1);   // flag for modification of the segment

	for(int j=0; j<2*num-1; j++)
		modified_Flag[j] = false;

	int max1 = 50, max2 = 50, min1 = 100, min2 = 100;
	float jump_Rate1 = 0.25, jump_Rate2 = 0.25;
	site_1 = exon_Assembly[0].bound1-1;
	site_2 = exon_Assembly[num-1].bound2+1;
	bool locate_Flag = false;
	bool ctr_Flag = false;
	seg_Mark = 0;
	vector<Rd>& rd_Vec = rd_Set.rd_Vec;

#pragma endregion

#pragma region record the segments
	for(seg_Mark=0; seg_Mark<2*num-1; seg_Mark++)
	{
		int index = 0;
		if(seg_Mark%2==0) // segments for exons
		{
			index = (int)seg_Mark/2;
			exon_start = exon_Assembly[index].bound1-2;
			exon_stop = exon_Assembly[index].bound2+2;	
			exon_Start.push_back(exon_start+2);
			exon_Stop.push_back(exon_stop-2);
		}
		else  // segments for introns
		{
			index = (int)(seg_Mark+1)/2;
			exon_start = exon_Assembly[index-1].bound2+3;
			exon_stop = exon_Assembly[index].bound1-3;
		}
		seg_Start.push_back(exon_start);
		seg_Stop.push_back(exon_stop);
	}

	for(int k=0; k<num-1; k++)
	{
		vector<int> seg1, inter_Seg1;
		rd_Left.push_back(seg1);
		rd_Left.push_back(inter_Seg1);
		vector<int> seg2, inter_Seg2;
		rd_Right.push_back(seg2);
		rd_Right.push_back(inter_Seg2);
	}
	vector<int> seg1, seg2;
	rd_Left.push_back(seg1);
	rd_Right.push_back(seg2);
#pragma endregion 

	// Identify the locations of the reads
	vector<int> single_rBoundary;
	read_Locate(seg_Start, seg_Stop, right_Boundary, single_rBoundary, exon_Start, exon_Stop, infile, rd_Set);

#pragma region Cluster the right boundaries of the fragments
	sort(right_Boundary.begin(),right_Boundary.end());  // sort the right boundaries of the fragments
	int num_1 = line_cnt1;   // number of the left reads
	int num_2 = right_Boundary.size();  // number of the right boundaries
	int right_boundary = num_2>0?(right_Boundary[num_2-1]+1):0;
	int num_3 = single_rBoundary.size();
	if(num_3>0)
		right_boundary = _max(single_rBoundary[num_3-1],right_boundary);

	if(exon_Assembly[num-1].len<500)  // the exon length is approximate
		exon_Score[num-1] = 0.8;

	exon_Stop[exon_Stop.size()-1] = right_boundary;
	seg_Stop[seg_Stop.size()-1] = right_boundary;

	seg_Mark = 0;
	exon_start = seg_Start[seg_Mark];
	exon_stop = seg_Stop[seg_Mark];
	int Num = 2*num-1;

	for(int k=0; k<num_2; k++)
	{
		x1 = right_Boundary[k];
		if(x1>=exon_start&&x1<=exon_stop)  // the starting point of the exon falls in the given exon
		{
			rd_Right[seg_Mark].push_back(x1);
		}
		else
		{
			seg_Mark++;
			locate_Flag = false;
			while(seg_Mark<2*num-1&&!locate_Flag)
			{
				exon_start = seg_Start[seg_Mark];
				exon_stop = seg_Stop[seg_Mark];
				if(x1>=exon_start&&x1<=exon_stop)  // the starting point of the exon falls in the given exon
				{
					locate_Flag = true;
					rd_Right[seg_Mark].push_back(x1);
					break;
				}
				seg_Mark++;
			}
		}
	}
#pragma endregion

#pragma region clear the arrays and close the files
	vector<int>().swap(right_Boundary);
	vector<int>().swap(exon_Start);
	vector<int>().swap(exon_Stop);
	vector<int>().swap(seg_Start);
	vector<int>().swap(seg_Stop);
	vector<bool>().swap(modified_Flag);
#pragma endregion
}

// Examine the current results of exon delimitation and assemble exons into genes 
bool Gene_Batch::assemble_SubGraph(vector<vector<int>>& seg_Modify)
{
	int exon_Num = exon_Assembly.size();  // number of exons
	vector<vector<PAIR_INT>> start_New;   // possible starting points of exons
	vector<vector<PAIR_INT>> stop_New;    // possible stopping points of exons

	int Num = 2*exon_Num-1;           // number of segments
	vector<float> jump_thresh1(Num);  // threshold for distance between adjacent splice sites
	vector<float> jump_thresh2(Num);  // threshold for distance between adjacent splice sites
	int rd_num1 = 0, rd_num2 = 0, rd_num = 0, exon_Len = 0, seg_Len = 0;

#pragma region Parameters 
	float seglen_threshHigh = 1000;   // threshold of exon length for relaxing control
	float seglen_threshLow = 600;     // threshold of exon length for restricted control
	float seglen_Rate = 1.8;   // threshold for identifying exon as boundary exon
	float jump_rate = 0.18;    // ratio of splice site distance to exon length
	float score_thresh = 0.8; 
	float rate = 5.0;
	float jump_Thresh = 0;
	// when there are a small number of exons, the control is more restricted
	if(exon_Num<7){
		seglen_threshHigh = 1200;
		seglen_Rate = 2;
		jump_rate = 0.2;
	}
#pragma endregion

#pragma region Adjust the scores of exons according to segment lengths
	float aver_segLen = 200;  // average length of high-confidence exons average length of high-confidence exons
	int valid_Cnt = 0;        
	for(int i=1; i<exon_Num-1; i++)
	{
		if(exon_Score[i]>0.75&&exon_Assembly[i].len<4000)
		{
			aver_segLen += exon_Assembly[i].len;
			valid_Cnt++;
		}
	}
	aver_segLen = valid_Cnt>0?aver_segLen/valid_Cnt:200; 

	float seglen_thresh = _min(seglen_threshHigh,aver_segLen*seglen_Rate);  // threshold of segment length
	seglen_thresh = _max(seglen_threshLow, seglen_thresh);

	for(int i=1; i<exon_Num-1; i++)
	{
		// exon length is beyond the threshold while both of its ends can be identified by junctions
		if(exon_Assembly[i].len>seglen_thresh&&exon_Type[i]==3)
		{
			exon_Score[i] = score_thresh*seglen_thresh/exon_Assembly[i].len;
#pragma region identify whether the exon is intern exon
			bool inner_Flag = false;
			int start = exon_Assembly[i].bound1;
			int stop = exon_Assembly[i].bound2;
			int junction_Num = junction.size(); 
			for(int j=0; j<junction_Num; j++)
			{
				if(stop<junction[j].second&&start>junction[j].first)
				{
					inner_Flag = true;
					break;
				}
				else if(start<junction[j].first)
					break;
			}
#pragma endregion
			if(inner_Flag==false)
				exon_Type[i] = 0;
			else
				exon_Type[i] = 4;
		}
	}
#pragma endregion

#pragma region Calculate the thresholds
	float min1 = 120, min2 = 0, min3 = 0, max1 = 1200, t1 = 0, t2 = 0;
	for(int j=0; j<Num-1; j+=2)
	{
		rd_num1 = rd_Left[j].size();
		rd_num2 = rd_Right[j].size();
		rd_num = rd_num1 + rd_num2;
		int serial = j/2;
		exon_Len = exon_Assembly[serial].len;
		min2 = exon_Len*jump_rate;
		min3 = rate*exon_Len*1.0/(rd_num1+1);
		t1 = _max(min1,min2);
		t2 = _max(t1,min3);                // lower threshold
		jump_thresh1[j] = _min(max1,t2);   // upper threshold
		min3 = rate*exon_Len*1.0/(rd_num2+1);
		t2 = _max(t1,min3);
		jump_thresh2[j] = _min(max1,t2);
	}
	// Relax threshold for exons at the boundaries of genes
	int j = Num-1;
	rd_num1 = rd_Left[j].size();
	rd_num2 = rd_Right[j].size();
	rd_num = rd_num1 + rd_num2;
	int serial = j/2;
	exon_Len = exon_Assembly[serial].len;
	jump_rate = 0.05;
	min2 = exon_Len*jump_rate;
	min3 = rate*exon_Len*1.0/(rd_num1+1);
	t1 = _max(min1,min2);
	t2 = _max(t1,min3);                // lower threshold
	jump_thresh1[j] = _min(max1,t2);   // upper threshold
	min3 = rate*exon_Len*1.0/(rd_num2+1);
	t2 = _max(t1,min3);
	jump_thresh2[j] = _min(max1,t2);
#pragma endregion

	vector<int> Rank;   // record the serial change of exons
	for(int i=0; i<exon_Num; i++)
		Rank.push_back(i);

	int add_Num = intron_Coverage(Rank);  // examine read distribution across the introns
	/*if(add_Num>0)
		cout<<"Intron added"<<endl;*/

	exon_Num = exon_Assembly.size();  // number of current identified exons
	Interval interval;
	for(int i=0; i<exon_Num; i++)
	{
		int ori_Serial = Rank[i];
		if(ori_Serial>=0&&exon_Score[ori_Serial]<score_thresh)
		{
			vector<int> pair1, pair2, cnt1, cnt2, gap_1, gap_2;
			int j = 2*ori_Serial; 
			seg_Len = exon_Assembly[i].len;

			int cnt_thresh = _min(10,0.05*seg_Len);
			if(seg_Len<50)   // not considered when read length is too small
				continue;

#pragma region Search possible gene boundaries
			vector<int>& temp1 = rd_Left[j];   // record the segment the left read belongs to
			vector<int>& temp2 = rd_Right[j];  // record the segment the right read belongs to
			int pre_Left = 0, pre_Right = 0, temp = 0;

			rd_num1 = temp1.size();   // number of left reads in the current segment
			rd_num2 = temp2.size();   // number of right reads in the current segment

			rd_num = rd_num1+rd_num2;
			int exon_start = exon_Assembly[i].bound1;
			int exon_stop = exon_Assembly[i].bound2;
			float jump_Thresh = 0, boundary_Thresh = 0;

#pragma region discontinuties for the left points of reads		
			if(rd_num1>0)
			{
				jump_Thresh = jump_thresh1[j];
				boundary_Thresh = 50+fragLen_Mean;
				locate_Jump(temp1, rd_num1, j, exon_start, exon_stop, jump_Thresh, 
					boundary_Thresh, pair1, cnt1, gap_1);
			}
#pragma endregion 

#pragma region discontinuties for the right points of reads	
			if(rd_num2>0)
			{
				jump_Thresh = jump_thresh2[j];
				boundary_Thresh = 50;
				locate_Jump(temp2, rd_num2, j, exon_start, exon_stop, jump_Thresh, 
					boundary_Thresh, pair2, cnt2, gap_2);
			}
#pragma endregion

#pragma endregion Search possible gene boundaries
			if(pair1.size()>0||pair2.size()>0)
			{
				int change_type = exon_Type[i];
				interval.ch_type.push_back(change_type);
				interval.site_New1.push_back(pair1); 
				interval.site_New2.push_back(pair2); 
				interval.cnt_New1.push_back(cnt1);
				interval.cnt_New2.push_back(cnt2);
				interval.gap1.push_back(gap_1);
				interval.gap2.push_back(gap_2);
				interval.serial.push_back(i); 
			}
		}
	}

	interval_Locate(interval,seg_Modify);

#pragma region Clear the arrays
	int n1 = rd_Left.size();
	for(int j=0; j<n1; j++)
		vector<int>().swap(rd_Left[j]);
	int n2 = rd_Right.size();
	for(int j=0; j<n2; j++)
		vector<int>().swap(rd_Right[j]);
	vector<vector<int>>().swap(rd_Left);
	vector<vector<int>>().swap(rd_Right);
#pragma endregion

	return true;
}

// Modify the predicted gene structure in the current position
// Return: the number of assembled genes in the cluster
int Gene_Batch::gene_Modify(vector<Gene>& gene_Set, vector<vector<int>>& seg_Modify)
{
	int num = seg_Modify.size();  // number of the modified exons
	int serial = 0, i = 0, p1 = 0, p2 = 0, p3 = 0, p4 = 0;
	int exon_Num = exon_Assembly.size();  // number of exons

	vector<Exon_Node>::iterator exon_itr = exon_Assembly.begin();  // iterator of the exon vector
	vector<PAIR_INT> mark;   // marks for the exons
	int gene_Cnt = 1;        // number of genes
	vector<int> add_Serial;  // new serial for added gene

	vector<int> Rank(exon_Num);  // serials of exons
	for(int j=0; j<exon_Num; j++)
		Rank[j] = j;

#pragma region modify the exons
	// deal with the breaking points due to segmentation of the exons
	vector<int>().swap(starting_site_vec);  
	vector<int>().swap(ending_site_vec);    
	int end_Serial = exon_Assembly.size()-1;
	vector<int> exon_Type_1 = exon_Type;    // exon type 

	for(int k=0; k<num; k++)
	{
		vector<int>& points = seg_Modify[k];  // breakpoints in the segments
		serial = points[0];
		int i = Rank[serial];   // ranking in the current exon combination 
		int point_num = (points.size()-1)/2;  // number of breaking points
		p1 = points[1]; 
		p2 = points[2];
		if(point_num>1)
		{
			p3 = points[3];
			p4 = points[4];
			if(p1>p3)  // rank the breaking points
			{
				int temp1 = p1, temp2 = p2;
				p1 = p3; p2 = p4;
				p3 = temp1; p4 = temp2;
			}
		}

		if(exon_Type[serial]==0)  // exon type: exon is too long; it may be boundary exon
		{
#pragma region cases of long exons on the boundaries of genes
			int temp = exon_Assembly[i].bound2;
			if(serial>0&&serial<exon_Num-1)    // internal exo is to be breaked and gene novel genes
			{
				exon_Assembly[i].bound2 = p1;  // modify the boundaries of the previous exon
				exon_Assembly[i].len = exon_Assembly[i].bound2 - exon_Assembly[i].bound1 + 1;
				Exon_Node cur_Exon;
				cur_Exon.bound1 = p2;
				cur_Exon.bound2 = temp;
				cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;      // identfiy a new exon
				exon_Assembly.insert(exon_Assembly.begin()+i+1,cur_Exon);  // insert new predicted exon into the exon assembly vector
				exon_Type_1.insert(exon_Type_1.begin()+i+1,0);   // new exon at the boundary of gene
				gene_Cnt++;
				for(int j=serial+1; j<exon_Num; j++)  // change of exon index
					Rank[j]++;
				add_Serial.push_back(serial);         // serial for the added exon
				starting_site_vec.push_back(p2);
				ending_site_vec.push_back(p1);
			}
			else if(serial==exon_Num-1)  // for the last exon of the gene
			{
				exon_Assembly[i].bound2 = p1;
				starting_site = p2;      // current transcription start site
				ending_site = p1;        // current transcription stop site
				exon_Assembly[i].len = exon_Assembly[i].bound2 - exon_Assembly[i].bound1 + 1;
				exon_Type_1[i] = 0;
			}
#pragma endregion
		}
		else if(exon_Type[serial]==1)  // only the right end of exon is identified; the left end is uncertain
		{
#pragma region cases of exons with uncertain left ends
			int temp = exon_Assembly[i].bound2;
			if(temp-p2>50)  // length satisfying the condition
			{
				exon_Assembly[i].bound1 = p2;				
				exon_Assembly[i].len = exon_Assembly[i].bound2 - exon_Assembly[i].bound1 + 1;
				exon_Type_1[i] = 1;   // exon type
				// cout<<"Left boundary modified"<<endl;
			}
#pragma endregion
		}
		else if(exon_Type[serial]==2)  // only the left end of exon is identified; the right end is uncertain
		{
#pragma region cases of exons with uncertain right ends
			if(abs(p1-exon_Assembly[i].bound1)>50)  // junction satisfying the constraint for distance between splice sites
			{
				exon_Assembly[i].bound2 = p1;
				exon_Assembly[i].len = p1 - exon_Assembly[i].bound1 + 1;
				exon_Type_1[i] = 2;   // exon type
				// cout<<"Right boundary modified"<<endl;
			}
#pragma endregion
		}
		else if(exon_Type[serial]==4)   // long internal exons
		{
#pragma region cases of long internal exons
			int temp1 = exon_Assembly[i].bound1, temp2 = exon_Assembly[i].bound2;
			int add_Start = p2;
			int k = 0;
			if(p1-temp1>50)   // length satisfying constraints
			{
				exon_Assembly[i].bound2 = p1;
				exon_Assembly[i].len = p1 - exon_Assembly[i].bound1 + 1;
				exon_Type_1[i] = 4;  // exons broken into different parts
			}
			if(point_num>1)   // there are multiple possible breaking points
			{
				Exon_Node cur_Exon;
				cur_Exon.bound1 = p2;
				cur_Exon.bound2 = p3;
				cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
				if(cur_Exon.len>50)  // exon length satisfying constraints
				{
					add_Start = p4;
				}
			}
			if(temp2-add_Start>50)
			{
				Exon_Node cur_Exon;
				cur_Exon.bound1 = add_Start;
				cur_Exon.bound2 = temp2;
				cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;				
				exon_Assembly.insert(exon_Assembly.begin()+i+1+k,cur_Exon);
				exon_Type_1.insert(exon_Type_1.begin()+i+1+k,4);  // exons broken into different parts
				k++;
			}
			for(int j=serial+1; j<exon_Num; j++)  // change of exon index
				Rank[j] += k;
#pragma endregion
		}
	}
#pragma endregion

#pragma region Modify gene set
	int add_Num = add_Serial.size();       // number of added genes 
	int k = 0, k1 = 0, k2 = 0, start_Serial = 0;
	float thresh = fragLen_Mean*1.5;
	non_Link.clear();
	if(add_Num>=1)   // re-delimitation of the exons
	{
		add_Serial.push_back(end_Serial);  // add the last exon
		add_Num++;
		k1 = 0;
		for(int j=0; j<add_Num; j++)
		{
			Gene new_Gene;
			end_Serial = add_Serial[j];
			k2 = Rank[end_Serial];
			for(int m = k1; m<=k2; m++)
			{
				new_Gene.subexon.push_back(exon_Assembly[m]);
				new_Gene.exon_Type.push_back(exon_Type_1[m]);
			}
			new_Gene.subexon_Num = new_Gene.subexon.size();
			new_Gene.start = new_Gene.subexon[0].bound1;
			new_Gene.stop = new_Gene.subexon[new_Gene.subexon_Num-1].bound2;

			for(int l1 = 0; l1<new_Gene.subexon_Num; l1++)
			{
				if(new_Gene.exon_Type[l1]==1&&l1>0&&(new_Gene.subexon[l1].bound1-new_Gene.subexon[l1-1].bound2)>thresh)
					new_Gene.non_Link.push_back(make_pair(l1-1,l1));
				else if(new_Gene.exon_Type[l1]==2&&l1<new_Gene.subexon_Num-1&&(new_Gene.subexon[l1+1].bound1-new_Gene.subexon[l1].bound2)>thresh)
					new_Gene.non_Link.push_back(make_pair(l1,l1+1));			
			}
			gene_Set.push_back(new_Gene);
#pragma region Allocate the reads into different read vectors

#pragma endregion
			k1 = k2+1;
		}
	}
	else
	{
		int exon_num = exon_Assembly.size();  // number of exons
		for(int l1 = 0; l1<exon_num; l1++)
		{
			if(exon_Type_1[l1]==1&&l1>0&&(exon_Assembly[l1].bound1-exon_Assembly[l1-1].bound2)>thresh)
				non_Link.push_back(make_pair(l1-1,l1));
			else if(exon_Type_1[l1]==2&&l1<exon_num-1&&(exon_Assembly[l1+1].bound1-exon_Assembly[l1].bound2)>thresh)
				non_Link.push_back(make_pair(l1,l1+1));			
		}
	}
#pragma endregion

	// predict the starting site of the next gene based on the current predicted gene locus
	int exon_num = exon_Assembly.size();
	starting_site_1 = exon_Assembly[exon_num-1].bound2 + 50;
	exon_Type = exon_Type_1;
	return gene_Cnt;
}

// Simultaneous transcript assembly and abundance estimation for all the reads
void Gene_Batch::batch_AbunEstimation(char* filename1, char* filename2, char* data_path_output)
{
	ifstream infile1, infile2;
	infile1.open(filename1,ios::in);
	if(!infile1){  
		cout<<"Open file error！(Junction)"<<endl;
		return;
	}
	infile2.open(filename2, ios::in);
	if(!infile2){
		cout<<"Open file error!(Alignment)"<<endl;
		return; 
	}

	rd_Len = rdLen_estimation(filename2);  // obtain read length
	char filename[200], data_path[200];
	ofstream outfile;
	outfile.open(data_path_output,ios::out);
	if(!outfile)  // the filename is invalid
	{
		sprintf(data_path,data_path_output);
		sprintf(filename,"%s\\exply.txt",data_path_output); // temporary filename
		outfile.open(filename,ios::out);
		if(!outfile)
		{
			cout<<"Output path invalid!"<<endl;
			exit(-1);
		}
		outfile.close();
	}
	else  // the filename is valid
	{
		outfile.close();
	}

	int k = 0, s1 = 0, s2 = 0, start = 0, stop = 0;   // starting site and stopping site of the gene
	int gene_Cnt = 0, gene_vCnt = 1;
	char filename3[200];

	method_Mode = 0;  // set to the annotation-free mode

	sprintf(filename3,"%s\\record.txt",data_path);
	record_file.open(filename3,ios::out);
	
	string line, line1, word, seg[10];
	std::streampos pos = infile1.tellg();   // obtain the current file pointer	

	vector<string> chromoname_Vec;
	string chromo_Name = "chr1";
	int chromo_num = chromoname_Vec.size();  // number of chromomsomes
	chromosome_Num = 0;

	getline(infile1,line);   // read the first line: line of caption
	// cout<<line<<endl;

//////////////////////////////////////////////////////
#pragma region Batch mode 
	getline(infile1,line);   // read the second line
	// cout<<line<<endl;
	pre_line = line;
	stringstream stream(line);
	stream>>chromo_Name;     // the name of the first chromosome
	chromosome_Name = chromo_Name;
	next_chromoName = "";
	cout<<"chromosome name: "<<chromosome_Name<<endl;  // chromosome name
#pragma endregion
//////////////////////////////////////////////////////

	/*while(!infile1.eof())
	{*/
		clock_t start_Time = clock();    // start timing
		vector<int> exon_Start, exon_Stop;
		starting_site = 0;  // initialize the transcription start site

		// initialize fragment length and variance
		fgLen_estimation(filename1,filename2);
		read_Number = 0;
		int valid_Cnt = 0;
		int read_Cnt = 0;
		bool changed_Flag = false;  // flag for change of chromosome

		while(!infile1.eof())
		{
			changed_Flag = false;
			clock_t start_Time = clock();
			exon_Start.clear(); exon_Stop.clear(); exon_Assembly.clear();
			
			changed_Flag = exon_AssembleBatch(infile1,infile2,exon_Start,exon_Stop);

			pos = infile1.tellg();
			std::streampos pos1 = infile2.tellg();   // obtain the current file pointer
			int est_cnt = rev_ExonAssembly.size() + 1; 

#pragma region isoform discovery and abundance estimation
			Rd_Set rd_Set;
			rd_Set.chromosome_Name = chromosome_Name;
			rd_Set.rd_Len = rd_Len;
			for(int k = 0; k<est_cnt; k++)
			{
				if(k>0)
				{
					exon_Assembly = rev_ExonAssembly[k-1];
					exon_Score = rev_ExonScore[k-1];
					exon_Type = rev_ExonType[k-1];
				}

				vector<Gene> gene_Set;
				vector<vector<int>> seg_Modify;
				
				coverage_Assemble(infile2, rd_Set);  // assebly the starting sites and stopping sites of the transcripts
				if(rd_Set.rd_Vec.size()==0)
					continue;

				assemble_SubGraph(seg_Modify);
				int loc_geneNum = gene_Modify(gene_Set, seg_Modify); // modify gene structures in the current position

				if(!infile2)
					cout<<"Error！"<<endl;

				// cout<<"read number in the locus: "<<rd_Set.rd_Vec.size()<<endl;
				read_Cnt += rd_Set.rd_Vec.size();
				
				if(!exon_Assembly.empty())
				{
					for(int k=0; k<loc_geneNum; k++)
					{
						gene_Cnt++;
						Graph_Trans gph;
						if(loc_geneNum>1)
						{
							gph.exon_Vec = gene_Set[k].subexon;
							gph.exon_Type = gene_Set[k].exon_Type;
							gph.non_Link = gene_Set[k].non_Link;
						}
						else
						{
							gph.exon_Vec = exon_Assembly;
							gph.exon_Type = exon_Type;
							gph.non_Link = non_Link;
						}
						
						int exon_Num = gph.exon_Vec.size();  // number of exons
						gph.graph_Initial(exon_Num);         // initialize the graph structure

						int locate_Cnt =  assignment_Pose(rd_Set, gph, infile2, chromo_Name);  // load read alignment
						// exon_PrintBatch(outfile_assemble, gene_Cnt, gph.exon_Vec);

						// estimte mean and variance of fragment lengths
						rd_Set.locate_RdNum = locate_Cnt;
						cout<<"Located reads "<<locate_Cnt<<endl;
						rd_Set.mean = fragLen_Mean;
						rd_Set.sigma = fragLen_Var;
						if(locate_Cnt>100||(gene_Cnt==1&&locate_Cnt>0))
						{
							rd_Set.fragLen_Estimation(read_Number,fragLen_Mean,fragLen_Var);
							fragLen_Mean = rd_Set.mean;
							fragLen_Var = rd_Set.sigma;
						}					

#pragma region isoform discovery and abundance estimation
						if(locate_Cnt>=3)
						{
							stringstream ss;
							ss<<gene_vCnt;
							rd_Set.gene_Name = "MaxInfo." + ss.str();
							rd_Set.orientation = '.';
							cout<<gph.exon_Vec[0].bound1<<" "<<gph.exon_Vec[exon_Num-1].bound2<<endl;
							int s1 = rd_Set.locate_start, s2 = rd_Set.locate_stop;
							cout<<rd_Set.rd_Vec[s1].Left1<<" "<<rd_Set.rd_Vec[s2-1].Right2<<endl;
							// gph.print_Graph(); // output the graph structure
							bool ctr_Flag = abun_EstiSingle(rd_Set, gph, data_path, 0);
							if(ctr_Flag==true)
								gene_vCnt++;
							/*else
							cout<<"Too many dataDims"<<endl;*/
							gph.clear_JuncRank();
						}
#pragma endregion
						// rd_Set.~Rd_Set();
						gph.clear();
					}
				}				
			}
#pragma endregion

			record_file<<"fragment length: "<<rd_Set.gene_Name<<"\t"<<fragLen_Mean<<endl;

			rd_Set.~Rd_Set();
			vector<VEC_NODE>().swap(rev_ExonAssembly);
			vector<vector<int>>().swap(rev_ExonType);
			vector<vector<float>>().swap(rev_ExonScore);

			if(changed_Flag==true&&next_chromoName!="")  // there is change of chromosome
			{
				chromosome_Name = next_chromoName;
				cout<<"Changed: "<<chromosome_Name<<endl;
				next_chromoName = "";
				// break;  // single-chromosome mode
			}
			clock_t finish_Time = clock();  // stop timing
			double cal_Time = finish_Time - start_Time;
			// printf("Running time: %.4f ms\n",cal_Time);
		}

		clock_t finish_Time = clock();  // stop timing
		double cal_Time = finish_Time - start_Time;
		printf("Running time: %.4f ms\n",cal_Time);

		record_file<<"Total read number: "<<read_Cnt<<endl;
		if(infile1.eof())
		{
			// cout<<"reach file end"<<endl;
			// break;
		}
		read_file.close();

		// convert the temporary file to standard file
		int total_Number = read_Cnt;
		char filename_1[200];
		sprintf(filename,"%s\\exply.txt",data_path_output);
		sprintf(filename_1,"%s\\results.gtf",data_path_output);
		expression_conver(filename, filename_1, total_Number);

	// } // end while
	
	infile1.close();
	infile2.close();
	// outfile_assemble.close();
	record_file.close();
	// outfile_record.close();
}

// Assembly between two genes
int Gene_Batch::interval_Assembly(ifstream& infile, vector<PAIR_INT>& junction, vector<int>& orient, vector<Gene>& Gene_Set, int& serial)
{
	int num = junction.size();  // number of junctions
	bool ctr_Flag = false, changed_Flag = false;
	int bound1 = Gene_Set[serial].stop;    // right boundary of the previous gene
	int bound2 = Gene_Set[serial+1].start; // left boundary of the next gene
	int exon_start = bound2-10, exon_stop = bound1+10, junc_start = 0, junc_stop = 0, junc_num = 0, temp = 0, m_exonStop = 0, m_exonStart = bound2;
	int line_Cnt = 0, intra_Flag = 1;
	int temp_endingSite = bound2 - 10;
	string line,seg[12],orientation;
	vector<int> exon_Start, exon_Stop;

	if(bound2<bound1)  // there is gene overlap
		return -1;

	while(exon_start<bound2&&exon_stop>bound1)
	{
		if((line_Cnt==0||intra_Flag==1))
		{
			line = pre_line;
		}

		intra_Flag = -1;
		stringstream stream(line);
		for(int i=0; i<12; i++)
			stream>>seg[i];   // read the 12 fields

		if(seg[0]!=chromosome_Name)
		{
			next_chromoName = seg[0];  // change of chromosome name
			pre_line = line;
			record_file<<"Chromosome changed: "<<chromosome_Name<<endl;
			changed_Flag = true;
			break;
		}

		sscanf(seg[1].c_str(),"%d",&junc_start);  // starting site of junction
		sscanf(seg[2].c_str(),"%d",&junc_stop);   // stopping site of junction
		sscanf(seg[4].c_str(),"%d",&junc_num);    // number of reads spanning junctions

#pragma region exon delimitation according to the annotations
		string::size_type pre = 0, pos;
		char tPar = ',';
		if((pos=seg[10].find(tPar,pre))!=string::npos)  // 在字符串中查找","字符
		{
			string str = seg[10].substr(pre,pos-pre);
			sscanf(str.c_str(),"%d",&temp);
			exon_stop = junc_start + temp;    // junction 起点
		}
		pre = 0;
		if((pos=seg[11].find(tPar,pre))!=string::npos)  // search the splitting character ‘,’
		{
			int len = seg[11].length();
			string str = seg[11].substr(pos+1,len-pos-1);
			sscanf(str.c_str(),"%d",&temp);   // exon length
			exon_start = junc_start + temp + 1;
		}
#pragma endregion

		if(exon_stop>bound2)   // beyond the boundary
		{
			pre_line = line;
			changed_Flag = true;
			break;
		}

		if(exon_stop<bound1)   // beyond the boundary
		{
			getline(infile,line);
			pre_line = line;
			changed_Flag = true;
			return 2;
		}
		//////////////////////////////////////////////////////////////////////////
		// error control
		float spr_Thresh_1 = 50000, spr_Thresh_2 = 20000, spr_Thresh_Scale = 10000;
		int jump_Distance = abs(exon_start-exon_stop);
		if(jump_Distance>spr_Thresh_1)  // too large jumping distance indicates possible error
		{
			sus_Junc.push_back(make_pair(exon_stop,exon_start));
			getline(infile,line);
			line_Cnt++;
			pre_line = line;
			continue;
		}
		else if(jump_Distance>spr_Thresh_2)  // too large jumping distance indicates possible error
		{
			if(junc_num<3)
			{
				sus_Junc.push_back(make_pair(exon_stop,exon_start));			
				getline(infile,line);
				line_Cnt++;
				pre_line = line;
				continue;
			}			
		}

		if(line_Cnt==0)
		{
			orientation = seg[5];  // read the orientation
		}
		line_Cnt++;		

		// record the junction
		junction.push_back(make_pair(exon_stop,exon_start));

		if(seg[5]!=orientation)    // there is change of orientation
		{
			changed_Flag = false;
			pre_line = line;       // record the current line
			temp_endingSite = _min(exon_start+2000,exon_stop-50);
			break;
		}
		else  // there isn't change of orientation
		{
#pragma region
			if(exon_stop>m_exonStop)
				m_exonStop = exon_stop;  // updata right boundary of the exon

			if(m_exonStop-m_exonStart>spr_Thresh_Scale)    // expressed segment is longer than expectation
			{
				changed_Flag = false;
				pre_line = line;     // record the current line
				int t1 = (int)((m_exonStart+m_exonStop)/2.0);
				temp_endingSite = _max(m_exonStart+2000,t1);
				break;
			}
			else
			{
				exon_Stop.push_back(exon_stop);
				exon_Start.push_back(exon_start);
			}

			if(exon_start>m_exonStart)
				m_exonStart = exon_start;
#pragma endregion
		}

		getline(infile,line);
	}

	if(infile.eof())
		cout<<"A reach file end!"<<endl;

	exon_Type.clear(); exon_Score.clear();
	char mark = '"';
	starting_site = bound1+50;  // left boundary of the region
	ending_site1 = bound2-50;   // right boundary of the region
	if(!exon_Start.empty())
	{		
		exon_Stop.push_back(temp_endingSite);		
		exon_Predict(exon_Start, exon_Stop, starting_site, ending_site1, exon_Assembly, exon_Type, exon_Score);  // 外显子组装

		Gene cur_Gene;
		cur_Gene.subexon = exon_Assembly;
		int exon_num = exon_Assembly.size();
		int length = Gene_Set[serial].gene_Name.length();
		string gene_Name = Gene_Set[serial].gene_Name.substr(0,length-2)+"_1";
		cur_Gene.gene_Name = gene_Name + mark + ';';
		cur_Gene.subexon_Num = exon_num;
		cur_Gene.exon_Type = exon_Type;
		cur_Gene.orientation = orientation[0];
		cur_Gene.start = exon_Assembly[0].bound1;
		cur_Gene.stop = exon_Assembly[exon_num-1].bound2;
		serial++;
		Gene_Set.insert(Gene_Set.begin()+serial,cur_Gene);		
	}

	if(changed_Flag==false)
		return 0;
	else
		return 1;

}

// Gene assembly based on gene annotations
bool Gene_Batch::assembly_Annotation(ifstream& infile, vector<Gene>& Gene_Set)
{
#pragma region Declaration of Variables
	vector<int> exon_Start, exon_Stop, junc_Start, junc_Stop, orientation;
	vector<vector<int>> junc_Num;
	vector<PAIR_INT> junction;
	vector<Exon_Node> exon_Vec, exon_Assembly;
	int starting_site = 0, ending_site = 0;
#pragma endregion

	bool ctr_Flag = false, loc_Flag = false;
	int i = 0, serial = 0, gene_Num = Gene_Set.size(), gene_Cnt = gene_Num;
	string gene_Name, line;
	getline(infile,line);  // read the line of caption
	getline(infile,line);
	pre_line = line;

	while(serial<gene_Cnt)
	{
		gene_Name = Gene_Set[serial].gene_Name;  // gene name or gene number
		cout<<"No."<<serial+1<<" "<<gene_Name<<endl;
		exon_Start.clear();
		exon_Stop.clear();

		exon_Vec = Gene_Set[serial].subexon;    // vector of exons
		int exon_Num = exon_Vec.size();         // number of exons
		int len = exon_Vec[exon_Num-1].bound2 - exon_Vec[0].bound1 + 1;  // gene length
		loc_Flag = locate_Junction(infile, Gene_Set[serial], exon_Start, exon_Stop,junc_Num);		
		starting_site = Gene_Set[serial].start; // left boundary of gene
		ending_site = Gene_Set[serial].stop;    // right boundary of gene
		exon_Assembly.clear();
		if(loc_Flag==true)
		{
			exon_Predict(exon_Start, exon_Stop, starting_site, ending_site, exon_Assembly, exon_Type, exon_Score);  // exon assembly
			exon_CompareAnno(exon_Assembly,exon_Vec,Gene_Set[i],junc_Num);  // comparison with results based on gene annotations
		}
		ctr_Flag = false;
		int pre_Serial = serial;
		int ctr_Mark = 0;
		while(ctr_Mark<1&&serial<gene_Cnt)
		{
			ctr_Mark = interval_Assembly(infile, junction, orientation, Gene_Set,serial);
			if(ctr_Mark==-1)
			{
				pre_Serial++;
				serial++;  // jump to the next gene
			}
		}

		gene_Cnt = gene_Num + serial - pre_Serial;
		if(ctr_Mark==2)  // outside the gene region
			serial--;

		serial++;
	}

	return ctr_Flag;
}



