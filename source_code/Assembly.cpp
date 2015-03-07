
#include "Batch_Mode.h"
#include "Read_Assignment.h"
#include "Entropy.h"
#include "ISA.h"

///////////////////////////////////////////////////////////////////////////////
// Deal with the sequencing data for isoform identification and abundance estimation 
// initialize the graph
void Gene_Batch::ini_Graph(Graph_Trans& gph, string Gene_Id, int serial, ifstream& infile)
{
	string line, word, seg[20], word_Value;
	bool flag_Ctrl = false;
	int i = 0;
	vector<int> start_Site;  // start point
	vector<int> stop_Site;   // end point
	int start, stop;   // start point and end point of transciption
	int label_idx;     // location of "GENE_ID" field
	Gene cur_Gene;     // current gene

#pragma region find the location of "GENE_ID" field
	std::streampos pos1 = infile.tellg();   // get the position of current file pointer
	getline(infile,line);
	stringstream stream(line);
	stream>>seg[0];  // tax_id
	while(seg[i]!="gene_name"&&seg[i]!="gene_id"&&seg[i]!="Gene_Id"&&seg[i]!="Gene_id")
	{
		i++;
		stream>>seg[i];
	}

	if(seg[i]=="gene_name"||seg[i]=="gene_id"||seg[i]=="Gene_Id"||seg[i]=="Gene_id")  // find "gene_id" field
		label_idx = i+1;   // "gene_id" field
	infile.seekg(pos1);    // return to last line
#pragma endregion

	while(!flag_Ctrl&&getline(infile,line))
	{
#pragma region load "GENE_ID" filed
		stringstream stream(line);
		for(int i=0; i<15; i++)
			stream>>seg[i];
		word_Value = seg[label_idx];
		word = seg[label_idx-1];
		if(word!="gene_name"&&word!="gene_id"&&word!="Gene_Id"&&word!="Gene_id")
		{
			int i = 0;
			while(seg[i]!="gene_name"&&seg[i]!="gene_id"&&seg[i]!="Gene_Id"&&seg[i]!="Gene_id")
				i++;			
			int label_idx1 = i+1;
			word_Value = seg[label_idx1];
		}

		cur_Gene.gene_Name = word_Value; // gene name
		if(word_Value==Gene_Id)
			flag_Ctrl = true;		
#pragma endregion
	}

	if(flag_Ctrl==true)  // locate gene
	{
		while(word_Value==Gene_Id)      // within the range of current gene
		{
			if(seg[2]=="transcript")    // isoform field
			{
				sscanf(seg[3].c_str(),"%d",&start);
				sscanf(seg[4].c_str(),"%d",&stop);
				start_Site.push_back(start);  // start point of transciption
				stop_Site.push_back(stop);    // end point of transciption
			}
			else if(seg[2]=="exon")    // exon assumed
			{
				sscanf(seg[3].c_str(),"%d",&start);
				sscanf(seg[4].c_str(),"%d",&stop);
				exon_Insertion(cur_Gene.subexon,start,stop);  // insert new exon into current exon list
			}

			pos1 = infile.tellg(); // get the position of current file pointer
			getline(infile,line);  // read a line
			if(line=="")           // reach end of file
				break;
			stringstream stream(line);			
			for(int i = 0; i<15; i++)
				stream>>seg[i];
			word_Value = seg[label_idx];
			word = seg[label_idx-1];
			if(word!="gene_name"&&word!="gene_id"&&word!="Gene_Id"&&word!="Gene_id")
			{
				int i = 0;
				while(seg[i]!="gene_name"&&seg[i]!="gene_id"&&seg[i]!="Gene_Id"&&seg[i]!="Gene_id")
					i++;			
				int label_idx1 = i+1;
				word_Value = seg[label_idx1];
			}		
		}
		infile.seekg(pos1);   // move upward for one line

		gph.exon_Vec = cur_Gene.subexon;
		int exon_Num = gph.exon_Vec.size();   // number of exons
		gph.graph_Initial(exon_Num);
	}
}

// write gene structure according to annotation
void Gene_Batch::ini_GraphSpr(vector<Gene>& gene, char* filename1, string chromo_Name)
{
	ifstream myfile;
	myfile.open(filename1,ios::in);
	if(!myfile)
	{
		cout<<"Open file error"<<endl;
		exit(-1);
	}

	string line, word, word_Value;
	string seg[20];
	int seg_Num = 17;  // number of segments
	int line_Num = 0;
	int gene_Serial = -1, pre_Serial = -1;   // gene index
	int chr_idx = 1, type_idx = 2, label_idx = 9, ori_idx = 6, trs_idx = 11;
	bool ctr_Flag = false;
	string Gene_Id, Gene_IdPre;  // gene id

#pragma region find the location of "GENE_ID" field
	std::streampos pos1 = myfile.tellg();   // get the position of current file pointer
	getline(myfile,line);
	stringstream stream(line);
	stream>>seg[0];  // tax_id
	int i = 0;
	string temp_Str;
	while(seg[i]!="gene_name"&&seg[i]!="gene_id"&&seg[i]!="Gene_Id"&&seg[i]!="Gene_id")
	{
		i++;
		stream>>seg[i];
		cout<<seg[i]<<endl;
	}
	if(seg[i]=="gene_name"||seg[i]=="gene_id"||seg[i]=="Gene_Id"||seg[i]=="Gene_id")  // find "gene_id" field
		label_idx = i+1;   // "gene_id" field
	myfile.seekg(pos1);    // return to last line
#pragma endregion

	bool locate_Flag = false;
	while(std::getline(myfile,line))
	{
		stringstream stream(line);
		for(int i = 0; i<15; i++)
			stream>>seg[i];
		if(seg[0].substr(0,2)==chromo_Name)
		{  
			locate_Flag = true;
#pragma region load "GENE_ID" field
			word = seg[label_idx-1];
			if(word!="gene_name"&&word!="gene_id"&&word!="Gene_Id"&&word!="Gene_id")
			{
				int i = 0;
				while(seg[i]!="gene_id"&&seg[i]!="Gene_Id"&&seg[i]!="Gene_id")
					i++;
				int label_idx1 = i+1;
				Gene_Id = seg[label_idx1];
			}
			else
				Gene_Id = seg[label_idx];
#pragma endregion
#pragma region gene field
			if(Gene_Id!="")  // new gene
			{
				string gene_Name = Gene_Id;   // Locus: gene name
				int found = isExisting(gene,gene_Name);  // check whether the gene exists
				if(found>-1)  // if the gene exists, no gene is added
				{
					gene_Serial = found;  // locate current gene
				}
				else
				{		
					Gene cur_Gene;
					gene_Serial = pre_Serial;
					cur_Gene.gene_Name = gene_Name;
					cur_Gene.orientation = seg[ori_idx][0];
					if(gene_Serial>=0)  // if data exists
					{
						int num1 = gene[gene_Serial].subexon.size();  // number of subexons
						if(num1>0)
						{
							gene[gene_Serial].start = gene[gene_Serial].subexon[0].bound1;     // start point of gene
							gene[gene_Serial].stop = gene[gene_Serial].subexon[num1-1].bound2; // start point of gene
							gene[gene_Serial].subexon_Num = gene[gene_Serial].subexon.size();  // subexon number of gene
						}
					}
					gene_Serial++;
					pre_Serial = gene_Serial;
					gene.push_back(cur_Gene);  // record current gene
					cout<<"Gene "<<gene_Serial+1<<" : "<<gene_Name<<endl;
				}					
			}
#pragma endregion gene field
#pragma region isoform field
			if(seg[type_idx]=="transcript")  // isoform
			{
				Trans cur_Trans;  // load isoform
				sscanf(seg[3].c_str(),"%d",&(cur_Trans.start));  // start point
				sscanf(seg[4].c_str(),"%d",&(cur_Trans.stop));   // end point
				cur_Trans.trsName = seg[trs_idx];   // isoform name
				cur_Trans.orientation = seg[ori_idx][0]; // isoform direction
				cur_Trans.length = 0;  // isoform length
				gene[gene_Serial].trans.push_back(cur_Trans);
			}
#pragma endregion isoform field
#pragma region exon field
			else if(seg[type_idx]=="exon")        // exon
			{
				Exon_Node subexon;
				string trs_Name = seg[trs_idx];   // exon name
				int s1, s2;
				sscanf(seg[3].c_str(),"%d",&s1);  // start point
				sscanf(seg[4].c_str(),"%d",&s2);  // end point
				int found = isExisting_Trans(gene[gene_Serial].trans,trs_Name);  // check whether the gene exists
				if(found>-1)
				{
					int start = gene[gene_Serial].trans[found].start;  // start point of isoform
					int len = s2 - s1 + 1;  // exon length
					gene[gene_Serial].trans[found].length += len;
					gene[gene_Serial].trans[found].exon_len.push_back(len);
					gene[gene_Serial].trans[found].exon_start.push_back(s1-start);
				}

#pragma region load exon field and intron field
				Gene& cur_Gene = gene[gene_Serial];  // get current gene

				int subexon_Num = cur_Gene.subexon.size();  // number of current subexon

				int j1 = 0, j2 = 0;
				while(j1<subexon_Num&&cur_Gene.subexon[j1].bound2<s1) // check the left bound of the new exon
					j1++;

				if(j1==subexon_Num) // the current subexon and former subexons don't overlap
				{
					Exon_Node subexon;
					subexon.bound1 = s1; subexon.bound2 = s2; subexon.len = s2 - s1 + 1;
					cur_Gene.subexon.push_back(subexon);   // add new subexon
				}
				else  // the current subexon and former subexons overlap
				{
					vector<Exon_Node>::iterator iter;
					while(j2<subexon_Num&&cur_Gene.subexon[j2].bound2<s2)  // check the right bound of the new exon
						j2++;
					Exon_Node cur_Exon;
					if(j2==subexon_Num)  // exceed the boundary of the last exon
					{
#pragma region exons overlap
						if(s2>cur_Gene.subexon[j2-1].bound2+2)
						{
							cur_Exon.bound1 = cur_Gene.subexon[j2-1].bound2+1;
							cur_Exon.bound2 = s2;
							cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;  // subexon length
							cur_Gene.subexon.push_back(cur_Exon);  // insert subexon
						}

						iter = cur_Gene.subexon.begin() + j2 - 1;
						if(cur_Gene.subexon[j2-1].bound1<s1-1&&cur_Gene.subexon[j2-1].bound2>s1+1)  // two exons overlap
						{
							int pre_bound1 = cur_Gene.subexon[j2-1].bound1;  // right boundary of subexon
							cur_Gene.subexon[j2-1].bound1 = s1;  // update the boundary of subexon
							cur_Gene.subexon[j2-1].len = cur_Gene.subexon[j2-1].bound2 - s1 + 1;  // update the length of subexon
							cur_Exon.bound1 = pre_bound1; cur_Exon.bound2 = s1 - 1;  // devided exon
							cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
							cur_Gene.subexon.insert(iter,cur_Exon);
						}
#pragma endregion
#pragma region one exon contains the other
						else
						{
							int inter_Num = j2 - j1 - 1 ;
							int l1 = cur_Gene.subexon[j1].bound1, r1 = cur_Gene.subexon[j1].bound2;
#pragma region first end exon
							iter = cur_Gene.subexon.begin() + j1;	
							int insert_itr1 = j1 + 1;

							if(s1<r1-1&&s1>l1+1) // devide the first exon
							{
								cur_Gene.subexon[j1].bound1 = s1;   // update the boundary of subexon
								cur_Gene.subexon[j1].len = cur_Gene.subexon[j1].bound2 - cur_Gene.subexon[j1].bound1 + 1;  // update the length of subexon			
								cur_Exon.bound1 = l1; cur_Exon.bound2 = s1-1; // devided exon
								cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
								cur_Gene.subexon.insert(iter,cur_Exon);
								insert_itr1++;
							}
							else if(s1<l1-1)  // add new initial subexon
							{
								cur_Exon.bound1 = s1; 
								cur_Exon.bound2 = l1-1; // devided exon
								cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
								cur_Gene.subexon.insert(iter,cur_Exon);
								insert_itr1++;
							}

#pragma endregion first end exon
							iter = cur_Gene.subexon.begin() + insert_itr1;  // move the pointer forward for 2 steps
							for(int k1 = 0; k1<inter_Num; k1++)  // insert the current exon to the neighboring exons
							{
								int pre_bound2 = cur_Gene.subexon[insert_itr1-1].bound2;  // right boundary of the exon
								cur_Exon.bound1 = pre_bound2 + 1;
								cur_Exon.bound2 = cur_Gene.subexon[insert_itr1].bound1-1;
								cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;															
								if(cur_Exon.len>1)
								{
									cur_Gene.subexon.insert(iter,cur_Exon);
									insert_itr1 += 2;   //  move the pointer forward for 2 steps
								}
								else //if the exon is not inserted
									insert_itr1++;   //  move the pointer forward for 1 step
								// pre_bound2 = cur_Gene.subexon[insert_itr1-1].bound2;  // right boundary of the exon
								iter = cur_Gene.subexon.begin() + insert_itr1;
								// pos++;
							}
						}
#pragma endregion 
					}
					else
					{
#pragma region
						int l1 = cur_Gene.subexon[j1].bound1, r1 = cur_Gene.subexon[j1].bound2;  // left boundary and right boundary of the first exon
						int l2 = cur_Gene.subexon[j2].bound1, r2 = cur_Gene.subexon[j2].bound2;  // left boundary and right boundary of the last exon
						vector<Exon_Node>::iterator iter_bak;
						int inter_Num = j2 - j1;
						int pre_bound2 = cur_Gene.subexon[j1].bound2;  // right boundary of the exon
						Exon_Node cur_Exon;
						int insert_itr = j1 + 1; // insert intron index
						int pos = j2;  // record last exon index
						if(s1>l1+1&&s1<r1-1)  // devide first exon
						{
							cur_Gene.subexon[j1].bound1 = s1;   // update the boundary of the subexon
							cur_Gene.subexon[j1].len = cur_Gene.subexon[j1].bound2 - cur_Gene.subexon[j1].bound1 + 1;  // update the length of the subexon
							iter = cur_Gene.subexon.begin() + j1;					
							cur_Exon.bound1 = l1; cur_Exon.bound2 = s1-1; // devided exons
							cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
							cur_Gene.subexon.insert(iter,cur_Exon);
							insert_itr++;
							pos++;
							j1++;
						}
						else if(s1<l1-1)    // insert new subexons
						{	
							cur_Exon.bound1 = s1; 
							if(s2<l2-1)  // insert nre subexons into the intron
							{
								cur_Exon.bound2 = s2; // devided exons
							}
							else
							{
								cur_Exon.bound2 = l1-1; // devided exons
							}
							cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
							iter = cur_Gene.subexon.begin() + j1;  // insert new exons into the first exon
							cur_Gene.subexon.insert(iter,cur_Exon);
							insert_itr++;
							pos++;
						}
						if(s2>l2+1)   // devide the last exon
						{ 
#pragma region insert subexon		
							iter = cur_Gene.subexon.begin() + insert_itr;
							for(int k1 = 0; k1<inter_Num; k1++)  // insert the current exon to the neighboring exons
							{
								cur_Exon.bound1 = pre_bound2 + 1;
								cur_Exon.bound2 = cur_Gene.subexon[insert_itr].bound1-1;
								cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;															
								// iter++;
								if(cur_Exon.len>1)
								{
									cur_Gene.subexon.insert(iter,cur_Exon);
									pos++;
									insert_itr += 2;   // move the pointer forward for 2 steps
								}
								else//if the exon is not inserted
									insert_itr++;   // move the pointer forward for 1 step
								pre_bound2 = cur_Gene.subexon[insert_itr-1].bound2;  // right boundary of the exon
								iter = cur_Gene.subexon.begin() + insert_itr;
								// pos++;
							}

							if(s2<r2-1)  // devide the last exon into two exons
							{
								iter = cur_Gene.subexon.begin() + pos;
								cur_Exon.bound1 = cur_Gene.subexon[pos].bound1; cur_Exon.bound2 = s2; // devided exons
								cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;	
								cur_Gene.subexon.insert(iter,cur_Exon);  // devided exons
								pos++;
								cur_Gene.subexon[pos].bound1 = s2+1;  // update the boundary of the subexon
								cur_Gene.subexon[pos].len = cur_Gene.subexon[pos].bound2 - cur_Gene.subexon[pos].bound1 + 1;  // update the length of the subexon
							}
#pragma endregion insert subexon	
						}
						else if(l1!=l2)
						{
#pragma region insert subexon	
							iter = cur_Gene.subexon.begin() + pos;
							pos--;  // move the last exon backward for 1 step
							for(int k1 = 0; k1<inter_Num-1; k1++)
							{
								cur_Exon.bound1 = pre_bound2 + 1;
								cur_Exon.bound2 = cur_Gene.subexon[insert_itr].bound1-1;
								cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;							
								// iter++;
								if(cur_Exon.len>1)
								{
									cur_Gene.subexon.insert(iter,cur_Exon);
									pos++;
									insert_itr += 2;
								}
								else
									insert_itr++;   // move the pointer forward for 1 step
								pre_bound2 = cur_Gene.subexon[insert_itr-1].bound2;  // right boundary of the exon
								iter = cur_Gene.subexon.begin() + insert_itr;
							}
							int r2_pre = cur_Gene.subexon[pos].bound2;
							if(s2> r2_pre + 2)
							{
								cur_Exon.bound1 = r2_pre + 1; cur_Exon.bound2 = s2; // devided exons
								cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
								cur_Gene.subexon.insert(iter,cur_Exon);  // add new subexons
							}						
#pragma endregion insert subexon	
						}
#pragma endregion
					}
				}
#pragma endregion non=isoform field
			}
			else
			{
			}
#pragma endregion exon field

		}
		else if(locate_Flag==true)
		{
			break;
		}
		line_Num++;
	}

	if(!gene.empty())
	{
		int num1 = gene[gene_Serial].subexon.size();
		if(num1>0)
		{
			gene[gene_Serial].start = gene[gene_Serial].subexon[0].bound1;  // start point of the gene
			gene[gene_Serial].stop = gene[gene_Serial].subexon[num1-1].bound2; // end point of the gene
			gene[gene_Serial].subexon_Num = num1;  // subexon number of the gene
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
// check the relation between the interval according to the ends of read and the break point
bool Gene_Batch::junc_Linkage(Graph_Trans& gph, vector<int>& points_1, vector<int>& readIdx)
{
	int num = points_1.size();  // number of ends
	vector<int> points = points_1;
	for(int i=0; i<num; i++)
	{
		if(points_1[i]<0)
			points[i] = -points_1[i]/2-1;
	}
	int j1 = -1, j2 = -1, in_ser1 = 0, in_ser2 = 0;
	if(num==2)   // no hops
	{
		j1 = points[0]; j2 = points[1];
		in_ser1 = points_1[0]; in_ser2 = points_1[1];
	}

	if(in_ser1<0)
		readIdx.push_back(in_ser1);
	else
		readIdx.push_back(j1);
	for(int k1=j1+1; k1<j2; k1++)
	{
		readIdx.push_back(k1);  // find the exon which contains the left end read
		if(abs(gph.exon_Vec[k1].bound1-gph.exon_Vec[k1-1].bound2)>5)
			gph.junc_Link1[k1]++;
		if(abs(gph.exon_Vec[k1].bound2-gph.exon_Vec[k1+1].bound1)>5)
			gph.junc_Link2[k1]++;
	}
	if(in_ser2<0)
		readIdx.push_back(in_ser2);
	else
		readIdx.push_back(j2);  //  find the exon which contains the left end read
	if(j2!=j1)
	{
		if(abs(gph.exon_Vec[j2].bound1-gph.exon_Vec[j2-1].bound1)>5)
			gph.junc_Link1[j2]++;
		if(in_ser2<0)
			gph.junc_Link2[j2]++;
	}

	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////
// modify the relation between the graph according to the cover rate of the read

bool Gene_Batch::graph_Modify(Rd_Set& rd_Set, Graph_Trans& gph)
{
	int exon_Num = gph.exon_Vec.size(); // number of exons
	int intron_Num = exon_Num - 1;  // number of introns
	bool inser_Flag = false;  // check whether the read is in the intron
	int rd_Num = rd_Set.rd_Vec.size();  // number of read

	vector<int> inser_sum(intron_Num);  // number of segments in the intron
	vector<int> inser_cnt(intron_Num);  // number of segments in the intron

	vector<int> new_serial(exon_Num);
	inser_sum[0] = gph.intron_Add[0].size();
	for(int k=1; k<intron_Num; k++)
	{
		inser_cnt[k] = gph.intron_Add[k].size();
		inser_sum[k] = inser_sum[k-1] + inser_cnt[k];
	}

	new_serial[0] = 0;
	for(int i=1; i<exon_Num; i++)
		new_serial[i] = i+inser_sum[i-1];

	for(int i=0; i<exon_Num-1; i++)
	{
		if(gph.link_Mtx[i+1][i+2]>0)
		{

		}
	}	

	return true;
}

//////////////////////////////////////////////////////////////////////////////
// deal with the assignment of sequencing data
int Gene_Batch::assignment_Pose(/*char* filename1, char* filename2, */Rd_Set& rd_Set, Graph_Trans& gph, ifstream& myfile, string chromo_Name)
{
	// bam file format£º1.read name 2.flags describing the result of align 3.reference sequence name
	// 4.start point of left end read 5.mapping quality 6.detail description of mapping quality
	// 7.pair end alignment"=",single end alignment"*"  8.start point of another read 9.sequence length 
	// 10.read sequence 11.read quality 12. description of mapping
	stringstream ss;
	ss<<rd_Len;  // transform the length of read into string
	string len_R = ss.str()+"M"; // matching string
	string line, word, words[10], frag_Name, chromo_name;  // string

#pragma region paramenter declaration
	int t_serial = 0, serial1 = 0;    // read index  
	string tPar = "N";  // flag for isoform
	string tPar1 = "I";  // flag for insertion
	string tPar2 = "D";  // flag for deletion
	string::size_type pos, pos0;
	vector<string> unresolved_read;   // unlocated read
	vector<int> unresolved_serial;    // unlocated read index
	vector<int> right_readIdx;
	vector<int>::iterator itr;
	vector<string>::iterator str_itr;
	vector<Rd>& rd_Vec = rd_Set.rd_Vec;	
	int x1 = 0, x2 = 0, z1 = 0, z2 = 0, temp_Len = 0; // start point of the left read and the right read£¬length of fragment
	int read_len = rd_Len;  // length of read

	int junc_Len = 0, cover_Len = 0, p1 = 0, p2 = 0;
	short j1 = -1, j2 = -1, j3 = -1, j4 = -1;  // interval of the for ends of the two-end read
	int x_1 = 0, x_2 = 0, z_1 = 0, z_2 = 0, je = 0, rlen1 = 0, rlen2 = 0, j_3 = 0, j_4 = 0, j_1 = 0, j_2 = 0;
	bool anno_flag = false;   // flag of the boundary of the two-end read according to the matching information of a single read
	int line_cnt = 0;
	bool flag_Ctrl = true;
	bool locate_Flag = false;
	int exon_Num = gph.exon_Vec.size();   // number of exons
	int e_Num = exon_Num;
	int start_Site = gph.exon_Vec[0].bound1;  // the start point of the transcription
	int stop_Site = gph.exon_Vec[exon_Num-1].bound2;  // the end point of the transcription
	std::streampos pos1 = myfile.tellg();     // get the position of the current file pointer
	int rd_Cnt = 0;  // number of reads within the gene
	int rd_Num = rd_Set.rd_Vec.size();  // number of reads within the segment
	rd_Set.locate_start = -1;
	rd_Set.locate_stop = rd_Num;
	int invalid_Cnt = 0, speci_Cnt = 0;
#pragma endregion paramenter declaration
	ofstream outfile;
	char filename[200];
	
	for(int k=0; k<rd_Num; k++)  // load data from the valid gene
	{
		j1 = -1, j2 = -1, j3 = -1, j4 = -1;

		Rd& cur_Rd = rd_Set.rd_Vec[k];  // current read
		frag_Name = cur_Rd.name;
		anno_flag = false;

		temp_Len = cur_Rd.insert_Len;
		x1 = cur_Rd.Left1;
		x2 = cur_Rd.Left2;
		z1 = cur_Rd.Right1;
		z2 = cur_Rd.Right2;

		int s1 = 0, s2 = 0;
		if(z1==0){
			s1 = x1; s2 = x2;
		}
		else{
			s1 = _min(x1,z1); s2 = _max(_max(x1,z1),z2);
		}
		int in_ser1 = 0, in_ser2 = 0, in_ser3 = 0, in_ser4 = 0;		

		if(s1>=start_Site&&s2<stop_Site)   // the start point and the end point of the fragment should be within the range of the gene
		{
			if(rd_Cnt==0)
				rd_Set.locate_start = k;
			// cout<<rd_Cnt<<" "<<x1<<" "<<x2<<" "<<z1<<" "<<z2<<endl;

			if(rd_Cnt>=3)
				locate_Flag = true;

			flag_Ctrl=true;
			j1 = 0;
			cur_Rd.left_readIdx.clear();
			cur_Rd.right_readIdx.clear();
			cur_Rd.l1 = -1; cur_Rd.l2 = -1; cur_Rd.r1 = -1; cur_Rd.r2 = -1;
			while(j1<e_Num&&x1>=gph.exon_Vec[j1].bound1) j1++; j1--; //  find the exon which contains the left boundary of the read

			if(j1<0||x1>gph.exon_Vec[j1].bound2)
			{
				invalid_Cnt++;
				continue;  // current read mismatches within the exons
			}
			cur_Rd.l1 = j1;  // end1 judge1		

			int num1 = cur_Rd.left_pos.size();   // left end read, number of break point
			int num2 = cur_Rd.right_pos.size();  // right end read, number of break point

			// check the interval of the left end read
			if(num1==0)   // the string matches
			{
#pragma region split points are not in the left end read
				j2 = j1;   // the left end reads are in the same exon, neighboring exons are not taken into account
				while(j2<e_Num&&x2>=gph.exon_Vec[j2].bound1) j2++;  j2--; // find the exon which contains the right boundary of the read
				if(x2<gph.exon_Vec[j2].bound2)
				{
					cur_Rd.l2 = j2;  // end2 judge1
					rd_Cnt++;        // at least single end read exists in the isoform, increasing
					if(j1!=j2)
					{
						for(int k1 = j1; k1<j2; k1++)
						{
							if(z1>0)  // two-end read
								gph.link_Mtx[k1+1][k1+2]++;  // each two neighboring exons are connected between j1 and j2
							else      // single-end read
								gph.link_Mtx[k1+1][k1+2] += 0.75;
						}
						if(j2>j1+1)
						{
							for(int k1 = j1; k1<=j2; k1++)
								cur_Rd.left_readIdx.push_back(k1);  //  find the exon which contains the left boundary of the read
						}				
					}
				}				
#pragma endregion
			}
			else  // left end read, hops exist
			{
#pragma region split points are in the left end read
				int junc_Len = 0, cover_Len = 0;
				bool ctr_Flag = false;
				x_1 = cur_Rd.left_pos[0];  // the original location of the other end of the left end read
				int x_1a = cur_Rd.left_pos[1];  // the original location of the other end of the right end read

				j_1 = j1;
				while(j_1<e_Num&&x_1>=gph.exon_Vec[j_1].bound1) j_1++; j_1--;   // find the exon which contains the right boundary of the read
				for(int k1 = j1; k1<=j_1; k1++)
					cur_Rd.left_readIdx.push_back(k1);

				for(int l=1; l<num1-2; l+=2)
				{
#pragma region  check the middle exon
					p1 = cur_Rd.left_pos[l];
					p2 = cur_Rd.left_pos[l+1];
					// while(j_1<e_Num&&abs(gph.exon_Vec[j_1].bound1-p1)>2)
					// j_1++;  // locate middle exons
					while(j_1<e_Num&&(gph.exon_Vec[j_1].bound1-p1>2||gph.exon_Vec[j_1].bound2-p1<-2))
						j_1++;  // locate middle exons

					cover_Len = p2 - p1 + 1;
					if(j_1<e_Num)
					{
						while(cover_Len>-2)  // may contain neighboring exons
						{
							cur_Rd.left_readIdx.push_back(j_1);
							cover_Len -= gph.exon_Vec[j_1].len;		
							j_1++;
							if(j_1>=e_Num||abs(gph.exon_Vec[j_1].bound2-gph.exon_Vec[j_1-1].bound1)>2)
								break;
						}
						j_1 = _max(0,j_1-1);
					}
#pragma endregion
				}
				p1 = cur_Rd.left_pos[num1-1];
				j_2 = j_1;  // find the exon which contains the right end of the left end read
				while(j_2<exon_Num&&p1>=gph.exon_Vec[j_2].bound2) j_2++; // find the exon which contains the right boundary of the read
				j2 = j_2;   // find the exon which contains the right end of the left end read
				while(j2<exon_Num&&x2>=gph.exon_Vec[j2].bound1-1) j2++; j2--; // find the exon which contains the right boundary of the read

				if(x2<=gph.exon_Vec[j2].bound2){
					cur_Rd.l2 = j2;   // end2 judge2
					rd_Cnt++;
					if(j_2!=j_1)
					{
						for(int k1=j_2; k1<j2; k1++)
							cur_Rd.left_readIdx.push_back(k1);
					}
					cur_Rd.left_readIdx.push_back(j2);
				}

				// find the exon which contains the right end of the left end read
				int idnum = cur_Rd.left_readIdx.size(), index1 = 0, index2 = 0;
				for(int k1 = 0; k1<idnum-1; k1++)
				{
					index1 = cur_Rd.left_readIdx[k1];
					index2 = cur_Rd.left_readIdx[k1+1];
					// cout<<index1+1<<" "<<index2+1<<" "<<endl;
					if(index1!=index2)
					{
						if(z1>0)  // two-end read
							gph.link_Mtx[index1+1][index2+1]++;  // each two neighboring exons are connected between j1 and j2
						else      // single-end read
							gph.link_Mtx[index1+1][index2+1] += 0.75;  // the confidence of single-end read is low
					}
				}		
#pragma endregion 
			}

			if(z1>0)  // two-end read
			{
#pragma region check the interval of the right end read
				j3 = j1;
				while(j3<e_Num&&z1>=gph.exon_Vec[j3].bound1) j3++; j3--; // find the exon which contains the left boundary of the right end read
				if(z1<=gph.exon_Vec[j3].bound2)
				{		
					cur_Rd.r1 = j3;   // end3 judge1
					if(num2==0)
					{
#pragma region right end read split point does not exist
						j4 = j3;  // the right end reads are in the same exon, neighboring exons are not taken into account
						while(j4<e_Num&&z2>=gph.exon_Vec[j4].bound1) j4++; j4--; //  find the exon which contains the right boundary of the read
						if(z2<=gph.exon_Vec[j4].bound2)
						{
							cur_Rd.r2 = j4;   // end4 judge1
							if(j3!=j4)
							{
								int part_Len = gph.exon_Vec[j3].bound2 - z1 + z2 - gph.exon_Vec[j4].bound1 + 2;
								if(part_Len==read_len)  // junction only relates two exons
								{
									gph.link_Mtx[j3+1][j4+1]++;  // j3 and j4 are connected
									anno_flag = true;
									///// calculate the probability of the read
								}
							}
						}
#pragma endregion
					}
					else // right end read, break point exists
					{
#pragma region right end read split point exists
						z_1 = cur_Rd.right_pos[0];  // the original location of the first inner end of the right end read
						int z_1a = cur_Rd.right_pos[1];  // the original location of the other end of the left end read

						j_3 = j3;
						while(j_3<e_Num&&z_1>=gph.exon_Vec[j_3].bound1) j_3++; j_3--;   // find the exon which contains the right boundary of the read
						for(int k1 = j3; k1<=j_3; k1++)
							cur_Rd.right_readIdx.push_back(k1);
						for(int l=1; l<num2-2; l+=2)
						{
#pragma region  check the interval of the right end read
							p1 = cur_Rd.right_pos[l];
							p2 = cur_Rd.right_pos[l+1];
							//while(j_3<e_Num&&abs(gph.exon_Vec[j_3].bound1-p1)>2)
							//	j_3++;  // locate middle exons

							while(j_3<e_Num&&(gph.exon_Vec[j_3].bound1-p1>2||gph.exon_Vec[j_3].bound2-p1<-2))
								j_3++;  // locate middle exons

							cover_Len = p2 - p1 + 1;
							if(j_3<e_Num)
							{						
								while(cover_Len>2)  //  may contain neighboring exons
								{		
									cur_Rd.right_readIdx.push_back(j_3);
									cover_Len -= gph.exon_Vec[j_3].len;
									j_3++;
									if(j_3>=e_Num||abs(gph.exon_Vec[j_3].bound2-gph.exon_Vec[j_3-1].bound1)>2)
										break;							
								}
								j_3 = _max(0,j_3-1);
							}
#pragma endregion
						}
						p1 = cur_Rd.right_pos[num2-1];
						j_4 = j_3;  // find the exon which contains the inner end of the right end read
						while(j_4<exon_Num&&p1>=gph.exon_Vec[j_4].bound2) j_4++; // find the exon which contains the right boundary of the read
						j4 = j_4;   // find the exon which contains the right end of the right end read
						while(j4<exon_Num&&z2>=gph.exon_Vec[j4].bound1-1) j4++; j4--; // find the exon which contains the right boundary of the read
						if(z2<=gph.exon_Vec[j4].bound2){
							cur_Rd.r2 = j4;   // end4 judge2
							/*if(j_4!=j4)
							{*/
								for(int k1 = j_4; k1<=j4; k1++)
									cur_Rd.right_readIdx.push_back(k1);
							//}
						}

						// find the exon which contains the end point of the right end read
						int idnum = cur_Rd.right_readIdx.size(), index1 = 0, index2 = 0;
						//cout<<"right_readIdx:"<<endl;
						for(int k1 = 0; k1<idnum-1; k1++)
						{
							int serial = cur_Rd.right_readIdx[k1];
							index1 = rd_Vec[serial].right_readIdx[k1];
							index2 = rd_Vec[serial].right_readIdx[k1+1];
							//cout<<index1+1<<" "<<index2+1<<" "<<endl;
							if(index1!=index2)
								gph.link_Mtx[index1+1][index2+1]++;  // each two neighboring exons are connected between j1 and j2
						}
#pragma endregion right end read split point exists
					}
				}

				// additional connection
				if(cur_Rd.l2>=0&&j3==j2+1)   // connection exists
				{
					gph.link_Mtx[j2+1][j3+1]++;
				}
#pragma endregion check the interval of the right end read
			}
		}
		else if(s1>=stop_Site)  // exceed the valid range of the gene
		{
			flag_Ctrl = false;
			rd_Set.locate_stop = k;
			break;
		}
	}
#pragma endregion two-end matching of string

#pragma region  connection of graph
	// modify the information of exon
	int id1 = rd_Set.locate_start;
	/*cout<<"Invalid cnt: "<<invalid_Cnt<<endl;
	cout<<"Speci cnt: "<<speci_Cnt<<endl;*/
	if(locate_Flag==true&&id1>=0&&id1<rd_Vec.size())
	{	
		// cout<<"mark: "<<endl;
		int s1 = rd_Vec[id1].Left1;  // first left end read
		if(s1<gph.exon_Vec[0].bound2){
			gph.exon_Vec[0].bound1 = s1;
			gph.exon_Vec[0].len = exon_Assembly[0].bound2 - exon_Assembly[0].bound1+1;
		}

	}
#pragma endregion connection of graph

	return rd_Cnt;	
}

// determine the location of the read
int Gene_Batch::read_Locate(vector<int>& seg_Start, vector<int>& seg_Stop, vector<int>& right_Boundary, vector<int>& single_rBoundary,
	vector<int>& exon_Start, vector<int>& exon_Stop, ifstream& infile, Rd_Set& rd_Set)
{
#pragma region paramenter declaration
	int x1 = 0, x2 = 0, z1 = 0, z2 = 0, temp_Len = 0; // start point of the left read and the right read£¬length of the fragment
	int read_len = rd_Len;  // length of the read
	stringstream ss;
	ss<<read_len;  // transform the length of read into string
	string len_R = ss.str()+"M"; // matching string
	string line, word, words[6], frag_Name, chromo_Name;  // string
	int serial = 0, serial1 = -1;     // read index
	string::size_type pos, pos0, pos1;
	string tPar = "N";   // flag of isoform
	string tPar1 = "I";  // flag of insertion
	string tPar2 = "D";  // flag of deletion
	vector<string> unresolved_read;   // unlocated read
	vector<int> unresolved_serial;    // unlocated read index
	int junc_Len = 0, cover_Len = 0, p1 = 0, p2 = 0;
	char j1 = -1, j2 = -1, j3 = -1, j4 = -1;   // the interval of the four ends of the two-end reads
	int x_1 = 0, x_2 = 0, z_1 = 0, z_2 = 0, j_3 = 0, j_4 = 0, j_1 = 0, j_2 = 0;
	int pre_Left = 0, pre_Right = 0, seg_Mark = 0;
	
	int line_cnt = 0, line_cnt1 = 0, line_cnt2 = 0;
	int num = exon_Assembly.size();
	bool junc_Left = false, junc_Right = false, flag_Right = false;  // flag of hops of the left boundary and right boundary of the isoform
	int exon_start = exon_Assembly[0].bound1;  // count from the first interval
	int exon_stop = exon_Assembly[0].bound2;
	vector<Rd>& rd_Vec = rd_Set.rd_Vec;

	bool locate_Flag = false;
	bool ctr_Flag = false;
	bool anno_flag = false;  // get the location of the boundary of the two-end read from single read alignment
	seg_Mark = 0;
	vector<string>::iterator str_itr;   
	vector<int>::iterator itr;
	int bound_Slack = _min(exon_Assembly[num-1].bound2,exon_Assembly[num-1].bound1 + 30); // slack boundary of the last exon
	vector<Rd> rd_bak;
#pragma endregion paramenter declaration

#pragma region load the read
	int bak_Num = rd_Bak.size();
	for(int i=0; i<bak_Num; i++)
	{
		Rd cur_Rd = rd_Bak[i];
		if(cur_Rd.Left1>site_1&&cur_Rd.Right2<site_2)
			rd_Vec.push_back(cur_Rd);
	}
	vector<Rd>().swap(rd_Bak);
#pragma endregion

	while(/*read_Serial<300&&*/std::getline(infile,line))  // test
	{
		// cout<<line_cnt<<endl;
		// cout<<line<<endl;
		line_cnt++;
		Rd cur_Rd;  // exon
		anno_flag = false;
		stringstream stream(line);
		stream>>word; frag_Name = word;
		stream>>word; stream>>word; // load the first 3 strings 
		chromo_Name = word;
		
		if(chromo_Name==next_chromoName)  // turn to the next chormosome
			break;
		if(chromo_Name!=chromosome_Name)  // not current chromosome
			continue;

		for(int i = 0; i<6; i++)
			stream>>words[i];  // 4 start point of left end read 5 matching level 6 matching level 7 flag of pair-end alignment 8 start point of right end read 9 inserted size
		sscanf(words[5].c_str(),"%d",&temp_Len);  // inserted size
		sscanf(words[0].c_str(),"%d",&x1);  // start point of the left end read
		sscanf(words[4].c_str(),"%d",&z1);  // start point of the right end read
		// cout<<seg_Mark<<"\t"<<x1<<"\t"<<z1<<endl;

		int s1 = _min(x1,z1), s2= _max(x1,z1);
		if(words[3]!="="){
			s1 = x1;
			s2 = x1;
		}
		else{
			z2 = 0;
			if(temp_Len>0)
				z2 = x1 + temp_Len - 1;  // end point of the right end read
			s1 = _min(x1,z1), s2 = _max(_max(x1,z1),z2);
		}

		if(words[3]!="=")  // single-end matching
		{
			if(words[3]=="*"&&x1>site_1&&x1<site_2)
				single_rBoundary.push_back(x1+rd_Len);  // defaults for no hops
		}
		else if(s1>site_1&&s1<site_2)  // left end in the first interval, right end in the second interval
		{
			Rd cur_Rd;  // read
			if(temp_Len>=0)   // forward description of the alignment result
			{
#pragma region forward description of the alignment result
#pragma region check the interval of the left end point
				if(line_cnt1==0&&s1<exon_Assembly[0].bound2)
				{
					exon_Assembly[0].bound1 = s1;  // modify the left boundary of the start exon
					exon_Assembly[0].len = exon_Assembly[0].bound2 - exon_Assembly[0].bound1 + 1;
				}
				line_cnt1++;
				// cout<<line_cnt1<<endl;
				/*if(line_cnt1>20000)
					break;*/
				// cout<<line_cnt1<<" "<<x1<<" "<<x2<<" "<<z1<<" "<<z2<<endl;
				if(seg_Mark<2*num-1&&x1>=exon_start&&x1<=exon_stop)  // within the interval
				{
					rd_Left[seg_Mark].push_back(x1);
				}				
				else  // turn to the intron interval or the next interval
				{
					seg_Mark++;
					if(seg_Mark>=2*num-1){
						// cout<<x1<<"\t"<<z1<<" :Error!"<<endl;
						continue;
					}
					locate_Flag = false;
					while(seg_Mark<2*num-1&&!locate_Flag)
					{
						exon_start = seg_Start[seg_Mark];
						exon_stop = seg_Stop[seg_Mark];
						if(x1>=exon_start&&x1<=exon_stop)  // within the interval
						{
							locate_Flag = true;
							rd_Left[seg_Mark].push_back(x1);
							break;
						}
						seg_Mark++;
					}
				}
#pragma endregion check the interval of the left end point

#pragma region check the location of the split point of the left end read
				if(words[2]==len_R)   // matching string
				{
					j2 = j1;   // the right end reads are in the same exon, neighboring exons are not taken into account
					x2 = x1 + rd_Len - 1;   // end point of the left end read
				}
				else
				{
					int junc_Len = 0, cover_Len = 0;
					string::size_type pre = 0;    // record the location of last "N"

#pragma region left end and cover area of left end read
					p1 = x1;
					while((pos=words[2].find(tPar,pre))!=string::npos||
						(pos=words[2].find(tPar1,pre))!=string::npos||
						(pos=words[2].find(tPar2,pre))!=string::npos)  // find "N" or "I" in the isoform
					{
						pos0 = words[2].find("M",pre); // locate the first separator
						if(pos0!=string::npos)
						{ 
							string str = words[2].substr(pos0+1,pos-pos0-1);
							sscanf(str.c_str(),"%d",&junc_Len);   // distance span by the read
							if(words[2].substr(pos,1)==tPar1)     // insertion exists
								junc_Len = 0;
							str = words[2].substr(pre,pos0-pre);
							sscanf(str.c_str(),"%d",&cover_Len);  // length covered by the read 
							p2 = p1 + cover_Len - 1;
							cur_Rd.left_pos.push_back(p2);   // record the break point
							p1 = p2 + junc_Len + 1;
							cur_Rd.left_pos.push_back(p1);   // record the break point
							pre = pos + 1;   // next start point for searching
						}
					}
#pragma endregion left end and cover area of left end read

					// right end of the left end read
					pos0 = words[2].find("M",pre); // loacte the first separator
					string str;
					if(pos0!=string::npos)
					{
						str = words[2].substr(pre,pos0-pre);
						sscanf(str.c_str(),"%d",&cover_Len); // distance span by the read
						x_2 = p1 + cover_Len - 1;  // the origin position of the other end if the right end of the read
						
					}
					x2 = x_2;
				}
#pragma endregion
				if(z2<site_2&&z2>site_1)
					right_Boundary.push_back(z2);  // record the right end of the fragment
				
				// add read vector
				string::size_type pos0 = frag_Name.find(".",0); // locate the first separator
				frag_Name = frag_Name.substr(pos0+1,frag_Name.length()-pos0-1);
				cur_Rd.name = frag_Name;
				cur_Rd.l1 = -1; cur_Rd.l2 = -1; cur_Rd.r1 = -1; cur_Rd.r2 = -1;
				cur_Rd.Left1 = x1; cur_Rd.Left2 = x2;
				cur_Rd.Right1 = z1; cur_Rd.Right2 = z2;
				cur_Rd.insert_Len = temp_Len;  // length

				if(z2>x1)
				{
					if(z2<site_2)
					{
						rd_Set.rd_Vec.push_back(cur_Rd);  // add read to the list
						anno_flag = true;
						serial1++;
						if(x1>bound_Slack&&z2<site_2)  // read may belong to overlapped areas
							rd_bak.push_back(cur_Rd);
					}
					else
						rd_bak.push_back(cur_Rd);					
				}
				else
					rd_Set.singlerd_Vec.push_back(cur_Rd);  // add single-end read to the list

			    int part_Len = z2 - z1 + 1;   // the distance between the two ends of the right end read
				if(abs(part_Len-read_len)>2&&anno_flag==true)  // the alighment of the right end read should be furthur determined
				{
					unresolved_read.push_back(frag_Name);  // record fragment name
					unresolved_serial.push_back(serial1);  // record fragment index
					anno_flag = false;
				}
#pragma endregion
			}
			else  // inverse description of the alignment result
			{
#pragma region inverse description of the alignment result
				int num = unresolved_serial.size();  // number of unlocated right end read
				for(int k1=0; k1<num; k1++)
				{
					
					if(frag_Name==unresolved_read[k1])  // fragment names are the same
					{
						int serial = unresolved_serial[k1];
						x_1 = rd_Vec[serial].Left1;  x2 = rd_Vec[serial].Left2;
						z_1 = rd_Vec[serial].Right1; z2 = rd_Vec[serial].Right2;
						if(x1==z_1&&z1==x_1&&rd_Vec[serial].right_readIdx.size()==0)  // the same read
						{
							x1 = x_1; z1 = z_1;

#pragma region check the junction according to the matching string
							if(words[2]!=len_R)   // hops exist
							{
								junc_Len = 0, cover_Len = 0;
								string::size_type pre = 0; // record the last loacation of "N"
																	

								int p1 = z1, p2 = 0; 
								bool ctr_Flag = false;
								while((pos=words[2].find(tPar,pre))!=string::npos||
									(pos=words[2].find(tPar1,pre))!=string::npos||
									(pos=words[2].find(tPar2,pre))!=string::npos)  // find "N" in the isoform
								{
									pos0 = words[2].find("M",pre); // loacte the first seperator
									if(pos0!=string::npos)
									{ 
										string str = words[2].substr(pos0+1,pos-pos0-1);
										sscanf(str.c_str(),"%d",&junc_Len);   // distance span by the read
										if(words[2].substr(pos,1)==tPar1)     // insertion exists
											junc_Len = 0;

										str = words[2].substr(pre,pos0-pre);
										sscanf(str.c_str(),"%d",&cover_Len);  // length covered by the read 

										p2 = z1 + cover_Len - 1;  // the original location of the other end of the right end of the read
										cur_Rd.right_pos.push_back(p2);
										p1 = p2 + junc_Len + 1;  // the exon may contain the left boundary of the read
										cur_Rd.right_pos.push_back(p1);
										pre = pos + 1;   // next start point for searching 
									}
								}
							}
#pragma endregion check the junction according to the matching string

							str_itr = unresolved_read.begin() + k1;
							unresolved_read.erase(str_itr);  // delete read from the undetermined set
							itr = unresolved_serial.begin() + k1;
							unresolved_serial.erase(itr);  // delete read from the undetermined set
							break;
						}
					}
				}
#pragma endregion 
			}
		}
		else if(s1>=site_2)  // the left end and the right end exceed the limit
		{
			// cout<<x1<<" "<<x2<<" "<<z1<<" "<<z2<<" beyond the range"<<endl;
			break;
		}
		else
		{
			/*cout<<"next range"<<endl;
			infile.seekg(pos1);
			return false;*/
		}
	}

	rd_Bak = rd_bak;  // save the read in the temporary area
	vector<Rd>().swap(rd_bak);
	return rd_Set.rd_Vec.size();  // return the number of the read
}

// insert new exon
void Gene_Batch::exon_Insertion(vector<Exon_Node>& subexon, int start, int stop)
{
	int s1 = start, s2 = stop;
	int d1 = 0, d2 = 0;

	int subexon_Num = subexon.size();  // number of current subexon

	int j1 = 0, j2 = 0;
	while(j1<subexon_Num&&subexon[j1].bound2<s1) // check the left boundary of the new subexon
		j1++;
	if(j1==subexon_Num) // current subexon doesn't overlap the former subexons
	{
		Exon_Node cur_Exon;
		cur_Exon.bound1 = s1; cur_Exon.bound2 = s2; cur_Exon.len = s2 - s1 + 1;
		subexon.push_back(cur_Exon);   // insert new exon
	}
	else
	{
		vector<Exon_Node>::iterator iter;
		while(j2<subexon_Num&&subexon[j2].bound2<s2)  // check the right boundary of the new subexon
			j2++;
		Exon_Node cur_Exon;
		if(j2==subexon_Num)  // exceed the boundary of the last exon
		{
#pragma region exons overlap
			if(subexon[j2-1].bound1<s1)  // two exons overlap
			{
				int pre_bound2 = subexon[j2-1].bound2;  // right boundary of the subexon
				subexon[j2-1].bound2 = s1;  // update the boundary of the exon
				subexon[j2-1].len = s1 - subexon[j2-1].bound1 + 1;  // update the length of the exon
				cur_Exon.bound1 = s1 + 1; cur_Exon.bound2 = pre_bound2;  // splitted exon
				cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
				subexon.push_back(cur_Exon);  // splitted exon
				cur_Exon.bound1 = pre_bound2 + 1; cur_Exon.bound2 = s2;  // splitted exon
				cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
				if(cur_Exon.len>2)
					subexon.push_back(cur_Exon);  // splitted exon
			}
#pragma endregion
#pragma region one exon contains the other
			else
			{
				iter = subexon.begin() + j2 - 2;    // insert subexon
				if(s1<subexon[j2-1].bound1-2)
				{
					cur_Exon.bound1 = s1;
					cur_Exon.bound2 = subexon[j2-1].bound1-1;
					cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;  // the length of subexon
					subexon.insert(iter,cur_Exon);  // insert subexon
				}

				if(s2>subexon[j2-1].bound2+2)
				{
					cur_Exon.bound1 = subexon[j2-1].bound2+1;
					cur_Exon.bound2 = s2;
					cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;  // the length of subexon
					subexon.push_back(cur_Exon);  // insert subexon
				}
			}
#pragma endregion 
		}
		else
		{
			int l1 = subexon[j1].bound1, r1 = subexon[j1].bound2;  // left boundary and right boundary of the first end exon
			int l2 = subexon[j2].bound1, r2 = subexon[j2].bound2;  // left boundary and right boundary of the last end exon
			vector<Exon_Node>::iterator iter_bak;
			int inter_Num = j2 - j1;
			int pre_bound2 = subexon[j1].bound2;  // right boundary of the exon
			Exon_Node cur_Exon;
			int insert_itr = j1 + 1; // insert intron index
			int pos = j2;  // index of the last end exon
			if(s1>l1+1&&s1<r1-1)  // first end exon should be splitted
			{
				subexon[j1].bound1 = s1;   // update the boundary of the subexon
				subexon[j1].len = subexon[j1].bound2 - subexon[j1].bound1 + 1;  // update the length of the subexon
				iter = subexon.begin() + j1;					
				cur_Exon.bound1 = l1; cur_Exon.bound2 = s1-1; // splitted exon
				cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
				subexon.insert(iter,cur_Exon);
				insert_itr++;
				pos++;
				j1++;
			}
			else if(s1<l1-1)    // insert new subexon
			{	
				cur_Exon.bound1 = s1; 
				if(s2<l2-1)  // insert new subexons into the intron
				{
					cur_Exon.bound2 = s2; // splitted exon
				}
				else
				{
					cur_Exon.bound2 = l1-1; // splitted exon
				}
				cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
				iter = subexon.begin() + j1;  // insert new subexons into the first end exon
				subexon.insert(iter,cur_Exon);
				insert_itr++;
				pos++;
			}
			if(s2>l2+1)   // last end exon should be splitted
			{
#pragma region insert subexon		
				iter = subexon.begin() + insert_itr;
				for(int k1 = 0; k1<inter_Num; k1++)  // insert new subexons into the neighboring exons
				{
					cur_Exon.bound1 = pre_bound2 + 1;
					cur_Exon.bound2 = subexon[insert_itr].bound1-1;
					cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;															
					// iter++;
					if(cur_Exon.len>1)
					{
						subexon.insert(iter,cur_Exon);
						pos++;
						insert_itr += 2;   // move the pointer forward for 2 steps
					}
					else
						insert_itr++;   // move the pointer forward for 1 step
					pre_bound2 = subexon[insert_itr-1].bound2;  // right boundary of the exon
					iter = subexon.begin() + insert_itr;
					// pos++;
				}

				if(s2<r2-1)  // split the last end exon into 2 exons
				{
					iter = subexon.begin() + pos;
					cur_Exon.bound1 = subexon[pos].bound1; cur_Exon.bound2 = s2; // splitted exon
					cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;	
					subexon.insert(iter,cur_Exon);  // splitted exon
					pos++;
					subexon[pos].bound1 = s2+1;  // update the boundary of the subexon
					subexon[pos].len = subexon[pos].bound2 - subexon[pos].bound1 + 1;  // update the length of the subexon
				}
#pragma endregion insert subexon	
			}
			else if(l1!=l2)
			{
#pragma region insert subexon	
				iter = subexon.begin() + pos;
				pos--;  // move the last end exon backward for 1 step
				for(int k1 = 0; k1<inter_Num-1; k1++)
				{
					cur_Exon.bound1 = pre_bound2 + 1;
					cur_Exon.bound2 = subexon[insert_itr].bound1-1;
					cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;							
					// iter++;
					if(cur_Exon.len>1)
					{
						subexon.insert(iter,cur_Exon);
						pos++;
						insert_itr += 2;
					}
					else
						insert_itr++;   // move the pointer forward for 1 step
					pre_bound2 = subexon[insert_itr-1].bound2;  // right boundary of the exon
					iter = subexon.begin() + insert_itr;
				}
				int r2_pre = subexon[pos].bound2;
				if(s2> r2_pre + 2)
				{
					cur_Exon.bound1 = r2_pre + 1; cur_Exon.bound2 = s2; // splitted exon
					cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
					subexon.insert(iter,cur_Exon);  // insert new subexons
				}						
#pragma endregion insert subexon	
			}
		}
	}
}

// exon assembly
// input 1£º read data file
// input 2£º junction file
// input 3£º exon array
void Gene_Batch::exon_Assemble(char* filename1, char* filename)
{
	ifstream infile;
	infile.open(filename,ios::in);   // open file
	if(!infile)
		cout<<"Open file error!"<<endl;

	vector<int> exon_Start;   // start point of the exon
	vector<int> exon_Stop;    // end point of the exon
	vector<int> exon_Len;     // length of the exon

	string line, word, seg[15];
	string::size_type pos, pre;
	int exon_start = 0, exon_stop = 1e09, exon_len = 0, temp=0;
	char tPar = ',';
	exon_Start.push_back(0);

	getline(infile,line);     // load first line
	while(std::getline(infile,line))
	{
		stringstream stream(line);
		for(int i=0; i<12; i++)
			stream>>seg[i];   // load 12 field

#pragma region loacte the boundary of the exon
		pre = 0;
		sscanf(seg[1].c_str(),"%d",&exon_start);  // start point of the junction
		sscanf(seg[2].c_str(),"%d",&exon_stop);  // end point of the junction

		if((pos=seg[10].find(tPar,pre))!=string::npos)  // find "," in the string
		{
			string str = seg[10].substr(pre,pos-pre);
			sscanf(str.c_str(),"%d",&temp);   
			exon_Stop.push_back(exon_start+temp);
		}		

		pre = 0;
		if((pos=seg[11].find(tPar,pre))!=string::npos)  // find "," in the string
		{
			int len = seg[11].length();
			string str = seg[11].substr(pos+1,len-pos-1);
			sscanf(str.c_str(),"%d",&temp);   // length of the exon
			exon_Start.push_back(exon_start+temp+1);
		}
#pragma endregion loacte the boundary of the exon
	}

	int num1 = exon_Start.size();  // number of start point of the exon
	int num2 = exon_Stop.size();  // number of end point of the exon

	if(num1>0&&num2>0)
	{
		exon_Start[0] = _max(0,exon_Stop[0] - 5000);
		exon_stop = exon_Start[num1-1] + 10000;
		exon_Stop.push_back(exon_stop);
		num2 = exon_Stop.size();   // number of end point of the exon

		sort(exon_Start.begin(),exon_Start.end());  // sort the start point of the exon
		sort(exon_Stop.begin(),exon_Stop.end());    // sort the end point of the exon

		// handle the duplicate exons
		vector<int>::iterator iter;
		for(int j = 0; j<num1-1; j++)
		{
			if(abs(exon_Start[j]-exon_Start[j+1])<3)   // neighbor of the start point of the exon
			{
				iter = exon_Start.begin()+j+1;
				exon_Start.erase(iter);
				j--;
				num1--;
			}
		}

		for(int j = 0; j<num2-1; j++)
		{
			if(abs(exon_Stop[j]-exon_Stop[j+1])<3)   // neighbor of the start point of the exon
			{
				iter = exon_Stop.begin()+j+1;
				exon_Stop.erase(iter);
				j--;
				num2--;
			}
		}
#pragma region exon assembly
		// loacte the boundary of the exon in the annotation
		Exon_Node cur_Node;
		num1 = exon_Start.size();  // number of start point of the exon
		num2 = exon_Stop.size();   // number of end point of the exon

		int i = 0, k = 0;
		while(i<num2)
		{
			exon_stop = exon_Stop[i];
			int k1 = k;		
			while(k<num1&&exon_Start[k]<exon_stop)
				k++;
			for(int l=k1; l<k-1; l++)
			{
				Exon_Node cur_Node;
				cur_Node.bound1 = exon_Start[l];
				cur_Node.bound2 = exon_Start[l+1]-1;
				cur_Node.len = cur_Node.bound2-cur_Node.bound1+1;
				exon_Assembly.push_back(cur_Node);	
			}
			Exon_Node cur_Node;
			cur_Node.bound1 = exon_Start[k-1];
			cur_Node.bound2 = exon_stop;
			cur_Node.len = cur_Node.bound2-cur_Node.bound1+1;
			exon_Assembly.push_back(cur_Node);

			int i1 = i;
			while(i<num2-1&&exon_Start[k]>exon_Stop[i+1])   // connected right end of the exon
				i++;
			for(int l=i1; l<i; l++)
			{
				Exon_Node cur_Node;
				cur_Node.bound1 = exon_Stop[l]+1;
				cur_Node.bound2 = exon_Stop[l+1];
				cur_Node.len = cur_Node.bound2-cur_Node.bound1+1;
				exon_Assembly.push_back(cur_Node);  // insert new exons
			}
			i++;
		}
#pragma endregion
	}
	else
	{
		exon_Stop.push_back(10000);
		Exon_Node cur_Node;
		cur_Node.bound1 = exon_Start[0]; cur_Node.bound2 = exon_Stop[exon_Stop.size()-1]; 
		cur_Node.len = cur_Node.bound2 - cur_Node.bound1 + 1;
		exon_Assembly.push_back(cur_Node);
	}

	ifstream infile1;
	infile1.open(filename1,ios::in);  // load file
	if(!infile1)
		cout<<"Open file error!"<<endl;
	bounding_Assemble(infile1, exon_Start, exon_Stop); // find the boundary of the isoform

#pragma region clear array and close file
	vector<int>().swap(exon_Start);
	vector<int>().swap(exon_Stop);
	vector<int>().swap(exon_Len);
	infile.close();
	infile1.close();
#pragma endregion

}

/////////////////////////////////////////////////////////////////////////////////////////
// start point and end point of transcript assembly
int Gene_Batch::locus_Assemble(ifstream& infile, Rd_Set& rd_Set, Graph_Trans& gph, vector<Exon_Node>& exon_Assembly)
{
	if(!infile)
		cout<<"Open file error!"<<endl;

	// bam file format£º1.read name 2.flags describing the result of align 3.reference sequence name
	// 4.start point of left end read 5.mapping quality 6.detail description of mapping quality
	// 7.pair end alignment"=",single end alignment"*"  8.start point of another read 9.sequence length 
	// 10.read sequence 11.read quality 12. description of mapping

#pragma region parameter declaration
	int x1 = 0, x2 = 0, z1 = 0, z2 = 0, temp_Len = 0; // start point of the left read and the right read£¬length of fragment
	int read_len = rd_Len;  // length of read

	stringstream ss;
	ss<<read_len;  // transform the length of read into string
	string gene_Name = rd_Set.gene_Name;  // gene name
	string len_R = ss.str()+"M"; // matching string
	string line, word, words[6], frag_Name;  // string
	int serial = 0;      // read index
	string tPar = "N";   // flag of isoform
	string tPar1 = "I";  // flag of insertion
	string tPar2 = "D";  // flag of deletion
	string::size_type pos, pos0;
	vector<int>::iterator itr;
	vector<string>::iterator str_itr;
	vector<string> unresolved_read;   // unlocated read
	vector<int> unresolved_serial;    // unlocated read index
	vector<int> right_readIdx;	
	vector<Rd>& rd_Vec = rd_Set.rd_Vec;
	short j1 = -1, j2 = -1, j3 = -1, j4 = -1;  // the interval of the four ends of the two-end reads
	int junc_Len = 0, cover_Len = 0, p1 = 0, p2 = 0;
	int x_1 = 0, x_2 = 0, z_1 = 0, z_2 = 0, quality = 0;
	int je = 0, rlen1 = 0, rlen2 = 0, j_3 = 0, j_4 = 0, j_1 = 0, j_2 = 0;
	int line_Thresh = 10000, quality_Thresh = 1, pre_quality = -1;
	string pre_fragName = "";
	bool single_Flag = true;
	
	int pre_Left = 0, pre_Right = 0;
	bool anno_flag = false;  // get the location of the boundary of the two-end read from single read alignment
	int line_cnt = 0, line_cnt1 = 0, line_cnt2 = 0, serial1 = 0;
	bool junc_Left = false, junc_Right = false, flag_Right = false;  // flag of hops of the left and right boundary of isoform
	vector<int> right_Boundary;  // right boundary of isoform

	int acc_cnt1 = 0;
	int thresh1 = 20, thresh2 = 20, len_Thresh = 5;

	vector<vector<int>> inser_Junc;   // possible hop point in the intron
	int num = exon_Assembly.size(), e_Num = num;  // number of exon
	
	int seg_Mark = 0;
	int exon_start = exon_Assembly[0].bound1, exon_stop = exon_Assembly[0].bound2;  // count from the first interval
	int pre_start = 0, pre_stop = 0;
	vector<int> seg_Start, seg_Stop, exon_Start, exon_Stop; // start point and end point of the interval
	vector<bool> modified_Flag(2*num-1);   // flag for whether the interval is modified
	for(int j=0; j<2*num-1; j++)
		modified_Flag[j] = false;
	int site1 = exon_Assembly[0].bound1-1, site2 = exon_Assembly[exon_Assembly.size()-1].bound2+1;
	
	int max1 = 50, max2 = 50, min1 = 100, min2 = 100;
	float jump_Rate1 = 0.25, jump_Rate2 = 0.25;

	for(seg_Mark=0; seg_Mark<2*num-1; seg_Mark++)
	{
		int index = 0;
		if(seg_Mark%2==0) // interval of the exon
		{
			index = (int)seg_Mark/2;
			exon_start = exon_Assembly[index].bound1-2;
			exon_stop = exon_Assembly[index].bound2+2;
			exon_Start.push_back(exon_start+2);
			exon_Stop.push_back(exon_stop-2);
		}
		else
		{
			index = (int)(seg_Mark+1)/2;
			exon_start = exon_Assembly[index-1].bound2+3;
			exon_stop = exon_Assembly[index].bound1-3;
		}
		seg_Start.push_back(exon_start);
		seg_Stop.push_back(exon_stop);
	}

	bool locate_Flag = false;
	bool ctr_Flag = false;
	bool flag_Ctrl = false;
	int flag = 0;
	seg_Mark = 0;

	int exon_Num = gph.exon_Vec.size();

	char path[200];
	char filename[200];
	sprintf(path,"%s","L:\\data");
	
	ofstream outfile;  // if the read data is in the gene£¬output to file
#pragma endregion

	std::streampos pos1 = infile.tellg();  // get the position of current file pointer
	string chr_Name;   // chromosome name
	while(std::getline(infile,line))  // test
	{
		line_cnt++;
		anno_flag = false;
		stringstream stream(line);
		stream>>word; frag_Name = word;
		stream>>word; stream>>word;    // load the first 3 strings 
		chr_Name = word;               // chromosome name
		if(chr_Name!=chromosome_Name)  // chromosome name is changed
		{
			next_chromoName = chr_Name;
			return -1;
		}
		for(int i = 0; i<6; i++)
			stream>>words[i];  // 4 start point of left end read 5 matching level 6 matching level 7 flag of pair-end alignment 8 start point of right end read 9 inserted size
		sscanf(words[5].c_str(),"%d",&temp_Len);  // inserted size
		sscanf(words[0].c_str(),"%d",&x1);  // start point of the left end read
		sscanf(words[4].c_str(),"%d",&z1);  // start point of the right end read
		sscanf(words[1].c_str(),"%d",&quality);
		z2 = 0;
		if(temp_Len>0)
			z2 = x1 + temp_Len - 1;  // end point of the right end read
		
		int s1 = _min(x1,z1), s2 = _max(_max(x1,z1),z2);
		int in_ser1 = 0, in_ser2 = 0, in_ser3 = 0, in_ser4 = 0;

		if(s1>site1&&s2<site2)     // read is within the interval of gene
		{
#pragma region check the quality of read
			if(line_cnt1>line_Thresh)
			{
				quality_Thresh = 2;   // increase threshold
				single_Flag = false;
			}
			if(quality<quality_Thresh)  // ignore the read with low quality  
			{
				if(line_cnt1>line_Thresh)
					continue;
				else if(quality==pre_quality&&frag_Name==pre_fragName)  // the first one is saved
					continue;
			}

			pre_fragName = frag_Name;
			pre_quality = quality;

#pragma endregion check the quality of read

			string::size_type pos0 = frag_Name.find(".",0); // locate the first separator
			frag_Name = frag_Name.substr(pos0+1,frag_Name.length()-pos0-1);
			flag_Ctrl = true;
			flag = 1;

			if(temp_Len>=0)        //  forward description of the alignment result
			{
				if(line_cnt1==0)
				{
#pragma region initialize the interval
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
#pragma endregion initialize the interval

					gph.graph_Initial_1(exon_Num);  // initial the graph
				}
				line_cnt1++;

#pragma region left end read calculate the location of the ends 
				j1 = 0;
				while(j1<exon_Num&&x1>=gph.exon_Vec[j1].bound1) j1++; j1--; // find the exon which contains the left boundary

				if(j1<0||x1>gph.exon_Vec[j1].bound2)
					continue;  // current read doesn't match the exon

				// output the location
				// cout<<frag_Name<<" "<<s1<<" "<<s2<<" "<<words[2]<<endl;
				Rd cur_Rd;  // read
				cur_Rd.l1 = j1; cur_Rd.l2 = -1; cur_Rd.r1 = -1; cur_Rd.r2 = -1;

				if(x1>gph.exon_Vec[j1].bound2)   // read is in the intron
				{
					seg_Mark = 2*j1 + 1;   // intron index
					rd_Left[seg_Mark].push_back(x1);
					in_ser1 = -2*(j1+1);
				}
				else   // read is in the exon
				{
					seg_Mark = 2*j1;   // exon index
					rd_Left[seg_Mark].push_back(x1);
				}

				if(words[2]==len_R)   // matching string
				{
#pragma region matching string
					j2 = j1;   // the left end reads are in the same exon, neighboring exons are not taken into account
					x2 = x1 + rd_Len - 1;   // end point of the left end read
					while(j2<exon_Num&&x2>=gph.exon_Vec[j2].bound1) j2++;  j2--; // find the exon which contains the right boundary
					

					if(x2>gph.exon_Vec[j2].bound2)   // read is in the intron
					{
						seg_Mark = 2*j2 + 1;   // intron index
						rd_Left[seg_Mark].push_back(x2);
						in_ser2 = -2*(j2+1);

					}
					else
					{
						seg_Mark = 2*j2;  // read is in the exon
						cur_Rd.l2 = j2;   // end2 judge1
						rd_Left[seg_Mark].push_back(x2);
						
						if(x1<=gph.exon_Vec[j1].bound2)
						{
							for(int k1 = j1; k1<j2; k1++)
								gph.link_Mtx[k1+1][k1+2]++;  // each two neighboring exons are connected between j1 and j2

							/*if(j2>j1+1)
							{*/
							for(int k1 = j1; k1<=j2; k1++)
								cur_Rd.left_readIdx.push_back(k1);  // find the exon which contains the left end read
							//}
						}
					}
#pragma endregion
				}
				else    // splicing junction exists
				{
					int junc_Len = 0, cover_Len = 0;
					string::size_type pre = 0;    // find the location of last "N"
					bool ctr_Flag = false;
					j2 = j1;

#pragma region location of ends of left end read

					junc_Len = 0;
					cover_Len = 0;
					pre = 0;    // find the location of last "N"

#pragma region left end and covered area of left end read
					p1 = x1;
					while((pos=words[2].find(tPar,pre))!=string::npos||
						(pos=words[2].find(tPar1,pre))!=string::npos||
						(pos=words[2].find(tPar2,pre))!=string::npos)  // find "N" and "I" in the isoform
					{
						pos0 = words[2].find("M",pre); // locate the first separator
						if(pos0!=string::npos)
						{ 
							string str = words[2].substr(pos0+1,pos-pos0-1);
							sscanf(str.c_str(),"%d",&junc_Len);   // distance span by the read
							if(words[2].substr(pos,1)==tPar1)     // insertion exists
								junc_Len = 0;
							str = words[2].substr(pre,pos0-pre);
							sscanf(str.c_str(),"%d",&cover_Len);  // length covered by the read 
							p2 = p1 + cover_Len - 1;
							cur_Rd.left_pos.push_back(p2);   // record the break point
							p1 = p2 + junc_Len + 1;
							cur_Rd.left_pos.push_back(p1);   // record the break point
							pre = pos + 1;   // next start point for searching 
						}
					}
#pragma endregion 

					// right end of the left end read
					pos0 = words[2].find("M",pre); // loacte the first separator
					string str;
					if(pos0!=string::npos)
					{
						str = words[2].substr(pre,pos0-pre);
						sscanf(str.c_str(),"%d",&cover_Len); // distance span by the read
						x_2 = p1 + cover_Len - 1;  // the origin position of the other end if the right end of the read
						}
					x2 = x_2;

					int part_Len = z2 - z1 + 1;   // the distance between the two ends of the right end read
					if(abs(part_Len-read_len)>2&&anno_flag==true)  // the alighment of the right end read should be furthur determined
					{
						unresolved_read.push_back(frag_Name);  // record fragment name
						unresolved_serial.push_back(serial1);  // record fragment index
						anno_flag = false;
					}
#pragma endregion

#pragma region check the intervals of the ends of left end read
					ctr_Flag = false;
					x_1 = cur_Rd.left_pos[0];  // the original location of the other end of the left end of the read
					j_1 = j1;
					while(j_1<e_Num&&x_1>=gph.exon_Vec[j_1].bound1) j_1++; j_1--;   // find the exon which contains the right boundary
					for(int k1 = j1; k1<=j_1; k1++)
						cur_Rd.left_readIdx.push_back(k1);

					int num1 = cur_Rd.left_pos.size();
					for(int l=1; l<num1-2; l+=2)
					{
#pragma region  check the middle exon
						p1 = cur_Rd.left_pos[l];
						p2 = cur_Rd.left_pos[l+1];

						while(j_1<e_Num&&(gph.exon_Vec[j_1].bound1-p1>2||gph.exon_Vec[j_1].bound2-p1<-2))
							j_1++;  // locate the middle exons

						cover_Len = p2 - p1 + 1;
						if(j_1<e_Num)
						{
							while(cover_Len>-2)  // may contain neighboring exons
							{
								cur_Rd.left_readIdx.push_back(j_1);
								cover_Len -= gph.exon_Vec[j_1].len;		
								j_1++;
								if(j_1>=e_Num||abs(gph.exon_Vec[j_1].bound2-gph.exon_Vec[j_1-1].bound1)>2)
									break;
							}
							j_1 = _max(0,j_1-1);
						}
#pragma endregion
					}
					p1 = cur_Rd.left_pos[num1-1];
					j_2 = j_1;  // find the exon which contains the right end of the left end read
					while(j_2<exon_Num&&p1>=gph.exon_Vec[j_2].bound2) j_2++; // find the exon which contains the right boundary of the read
					j2 = j_2;   // find the exon which contains the right end of the left end read
					while(j2<exon_Num&&x2>=gph.exon_Vec[j2].bound1-1) j2++; j2--; // find the exon which contains the right boundary of the read

					if(x2<=gph.exon_Vec[j2].bound2){
						cur_Rd.l2 = j2;   // end2 judge2
						if(j_2!=j_1){
							for(int k1=j_2; k1<j2; k1++)
								cur_Rd.left_readIdx.push_back(k1);
						}
						cur_Rd.left_readIdx.push_back(j2);
					}

					// find the exon which contains the right end of the left end read
					int idnum = cur_Rd.left_readIdx.size(), index1 = 0, index2 = 0;
					for(int k1 = 0; k1<idnum-1; k1++)
					{
						index1 = cur_Rd.left_readIdx[k1];
						index2 = cur_Rd.left_readIdx[k1+1];
						// cout<<index1+1<<" "<<index2+1<<" "<<endl;
						if(index1!=index2)
						{
							if(z1>0)  // two-end read
								gph.link_Mtx[index1+1][index2+1]++;  //  each two neighboring exons are connected between j1 and j2
							else      // single-end read
								gph.link_Mtx[index1+1][index2+1] += 0.75;  // the confidence of single-end read is low
						}
					}		
#pragma endregion 
					
				}
#pragma endregion left end read calculate the location of the ends 

				anno_flag = true;
				if(words[3]=="=")   // determine the right boundary of the isoform according to the two-end matching
				{
#pragma region two-end matching string
					j3 = j1;
					while(j3<exon_Num&&z1>=gph.exon_Vec[j3].bound1) j3++; j3--; // find the exon which contains the left end of the right end read

					if(z1>gph.exon_Vec[j3].bound2)  // read is in the intron
					{
						seg_Mark = 2*j3+1;
						rd_Right[seg_Mark].push_back(z1);
						in_ser3 = -2*(j3+1);
					}
					else  // read is in the exon
					{		
						seg_Mark = 2*j3;
						rd_Right[seg_Mark].push_back(z1);

						cur_Rd.r1 = j3;   // end3 judge1
						j4 = j3;  // the left end reads are in the same exon, neighboring exons are not taken into account
						while(j4<exon_Num&&z2>=gph.exon_Vec[j4].bound1) j4++; j4--; // find the exon which contains the right boundary
						
						if(z2>gph.exon_Vec[j4].bound2)  // read is in the intron
						{
							seg_Mark = 2*j4+1;
							rd_Right[seg_Mark].push_back(z2);
							in_ser4 = -2*(j4+1);
						}
						else  // read is in the exon
						{
							seg_Mark = 2*j4;
							rd_Right[seg_Mark].push_back(z2);

							cur_Rd.r2 = j4;   // end4 judge1
							if(j3==j4)
								cur_Rd.right_readIdx.push_back(j3);
							else
							{
								int part_Len = gph.exon_Vec[j3].bound2 - z1 + z2 - gph.exon_Vec[j4].bound1 + 2;
								if(abs(part_Len-read_len)<3)  // junction only relates two exons
								{
									gph.link_Mtx[j3+1][j4+1]++;  // j3 and j4 are connected
									anno_flag = true;
									cur_Rd.right_readIdx.push_back(j3);
									cur_Rd.right_readIdx.push_back(j4);
								}
								else if(part_Len<read_len-3)  // the alighment of the right end read should be furthur determined
								{
									unresolved_read.push_back(frag_Name);  // record fragment name
									unresolved_serial.push_back(serial);   // record fragment index
									anno_flag = false;
								}
							}
						}
					}

#pragma region connection
					if(x2<=gph.exon_Vec[j2].bound2&&z1<=gph.exon_Vec[j3].bound2) // the end of the read is in the exon
					{
						if(j3==j2+1)   // two exons are adjacent
						{
							gph.link_Mtx[j2+1][j3+1]++;
							// cout<<x1<<" "<<x2<<" "<<z1<<" "<<z2<<":"<<j2+1<<" "<<j3+1<<" Neighbor linkage"<<endl;
						}
					}
#pragma endregion connection				

#pragma endregion two-end matching string
				}

#pragma region add read

				cur_Rd.name = frag_Name;
				// cur_Rd.l1 = j1; cur_Rd.l2 = j2; cur_Rd.r1 = j3; cur_Rd.r2 = j4;
				cur_Rd.Left1 = x1; cur_Rd.Right2 = z2;
				cur_Rd.Left2 = x2; cur_Rd.Right1 = z1;
				cur_Rd.insert_Len = temp_Len;  // length

				if(z2>x1)  // two-end read
				{

						rd_Set.rd_Vec.push_back(cur_Rd);  // add read to the list
						anno_flag = true;
						serial++;
				
				}
				else if(words[3]=="*"&&single_Flag==true)
				{
					cur_Rd.insert_Len = fragLen_Mean;
					rd_Set.rd_Vec.push_back(cur_Rd);   // add read to the list
					serial++;

				}
				
#pragma endregion
			}
			else  // inverse description of the alignment result
			{
#pragma region calculate the right end read and location of its ends
				if(words[3]=="=")   // two-end matching
				{
					int num = unresolved_serial.size();  // number of unlocated right end read
					for(int k1=0; k1<num; k1++)          // find the matching left end read
					{

						if(frag_Name==unresolved_read[k1])   // fragment names are the same
						{
							int serial1 = unresolved_serial[k1];
							x_1 = rd_Vec[serial1].Left1;  x2 = rd_Vec[serial1].Left2;
							z_1 = rd_Vec[serial1].Right1; z2 = rd_Vec[serial1].Right2;

							if(x1==z_1&&z1==x_1&&rd_Vec[serial1].right_readIdx.size()==0)  // the same read
							{
								x1 = x_1; z1 = z_1;
								j1 = rd_Vec[serial1].l1;  j2 = rd_Vec[serial1].l2;  j3 = rd_Vec[serial1].r1;  j4 = rd_Vec[serial1].r2;
								Rd& cur_Rd = rd_Vec[serial1];
								cur_Rd.right_readIdx.clear();
#pragma region check junction according to the matching string and calculate the location
								if(words[2]!=len_R)   // hops exist
								{
									junc_Len = 0, cover_Len = 0;
									string::size_type pre = 0; // record the last loacation of "N"

									int p1 = z1, p2 = 0; 
									bool ctr_Flag = false;
									while((pos=words[2].find(tPar,pre))!=string::npos||
										(pos=words[2].find(tPar1,pre))!=string::npos||
										(pos=words[2].find(tPar2,pre))!=string::npos)  // find "N" in the isoform
									{
										pos0 = words[2].find("M",pre); // loacte the first seperator
										if(pos0!=string::npos)
										{ 
											string str = words[2].substr(pos0+1,pos-pos0-1);
											sscanf(str.c_str(),"%d",&junc_Len);   // distance span by the read
											if(words[2].substr(pos,1)==tPar1)     // insertion exists
												junc_Len = 0;

											str = words[2].substr(pre,pos0-pre);
											sscanf(str.c_str(),"%d",&cover_Len);  // length covered by the read 

											p2 = p1 + cover_Len - 1;  // the original location of the other end of the right end of the read
											cur_Rd.right_pos.push_back(p2);
											p1 = p2 + junc_Len + 1;  // the exon may contain the left boundary of the read
											cur_Rd.right_pos.push_back(p1);
											pre = pos + 1;   // next start point for searching   
										}
									}
								}
#pragma endregion check junction according to the matching string and calculate the location

#pragma region calculate the location of the read
								int num2 = cur_Rd.right_pos.size();  // number of split points of the right end read

#pragma region check the interval of the right end of the read
								if(z1<=gph.exon_Vec[j3].bound2)
								{
									cur_Rd.r1 = j3;   // end3 judge1
									if(num2>0)  // right end read, break point exists
									{
#pragma region right end read split points exists
										z_1 = cur_Rd.right_pos[0];  // the original location of the first inner end of the right end read
										j_3 = j3;
										
										while(j_3<e_Num&&z_1>=gph.exon_Vec[j_3].bound1) j_3++; j_3--;   // find the exon which contains the right boundary of the read
										for(int k1 = j3; k1<=j_3; k1++)
											cur_Rd.right_readIdx.push_back(k1);
										for(int l=1; l<num2-2; l+=2)
										{
#pragma region  check middle exon
											p1 = cur_Rd.right_pos[l];
											p2 = cur_Rd.right_pos[l+1];
											//while(j_3<e_Num&&abs(gph.exon_Vec[j_3].bound1-p1)>2)
											//	j_3++;  //  locate middle exons

											while(j_3<e_Num&&(gph.exon_Vec[j_3].bound1-p1>2||gph.exon_Vec[j_3].bound2-p1<-2))
												j_3++;  //  locate middle exons

											cover_Len = p2 - p1 + 1;
											if(j_3<e_Num)
											{						
												while(cover_Len>2)  // may contain neighboring exons
												{		
													cur_Rd.right_readIdx.push_back(j_3);
													cover_Len -= gph.exon_Vec[j_3].len;
													j_3++;
													if(j_3>=e_Num||abs(gph.exon_Vec[j_3].bound2-gph.exon_Vec[j_3-1].bound1)>2)
														break;							
												}
												j_3 = _max(0,j_3-1);
											}
#pragma endregion
										}
										p1 = cur_Rd.right_pos[num2-1];
										j_4 = j_3;  // find the exon which contains the inner end of the right end read
										while(j_4<exon_Num&&p1>=gph.exon_Vec[j_4].bound2) j_4++; // find the exon which contains the right boundary of the read
										j4 = j_4;   // find the exon which contains the right end of the right end read
										while(j4<exon_Num&&z2>=gph.exon_Vec[j4].bound1-1) j4++; j4--; // find the exon which contains the right boundary of the read
										if(z2<=gph.exon_Vec[j4].bound2){
											cur_Rd.r2 = j4;   // end4 judge2
											for(int k1 = j_4; k1<=j4; k1++)
												cur_Rd.right_readIdx.push_back(k1);
										}

										// find the exon which contains the end point of the right end read
										int idnum = cur_Rd.right_readIdx.size(), index1 = 0, index2 = 0;
										//cout<<"right_readIdx:"<<endl;
										for(int k1 = 0; k1<idnum-1; k1++)
										{
											int serial = cur_Rd.right_readIdx[k1];
											index1 = rd_Vec[serial1].right_readIdx[k1];
											index2 = rd_Vec[serial1].right_readIdx[k1+1];
											//cout<<index1+1<<" "<<index2+1<<" "<<endl;
											if(index1!=index2)
												gph.link_Mtx[index1+1][index2+1]++;  // each two neighboring exons are connected between j1 and j2
										}
#pragma endregion right end read split points exists
									}
									else  // right end read, break point doesn't exist
									{
										if(j3>=0&&j4>=0){
											for(int k1=j3; k1<j4; k1++)
											{
												cur_Rd.right_readIdx.push_back(k1);
												gph.link_Mtx[k1+1][k1+2]++;
											}
											cur_Rd.right_readIdx.push_back(j4);
										}
									}
								}

								// additional connection
								if(cur_Rd.l2>=0&&cur_Rd.r1>=0&&j3==j2+1)   // connection exists
								{
									gph.link_Mtx[j2+1][j3+1]++;
								}
#pragma endregion check the interval of the right end read

#pragma endregion  calculate the location of the read
							}
						}

						str_itr = unresolved_read.begin() + k1;
						unresolved_read.erase(str_itr);  // delete read from the undetermined set
						itr = unresolved_serial.begin() + k1;
						unresolved_serial.erase(itr);  // delete read from the undetermined set
						break;
					}
				}
#pragma endregion calculate the location of the ends of the read
			}
		}
		else if(s1>=site2&&s2>=site2)  // the left end and the right end exceed the limit
		{
			// break;
			cout<<"beyond the range"<<endl;
			infile.seekg(pos1);
			break;
			// return flag_Ctrl;
		}
		else
		{
			/*infile.seekg(pos1);
			return false;*/
		}
		pos1 = infile.tellg();
 	}
	// cout<<"Assignment loaded"<<endl;

	rd_Set.locate_start = 0;
	rd_Set.locate_stop = rd_Vec.size();


#pragma region clear array and close file
	vector<int>().swap(right_Boundary);
	vector<int>().swap(exon_Start);
	vector<int>().swap(exon_Stop);
	vector<int>().swap(seg_Start);
	vector<int>().swap(seg_Stop);
	int n1 = rd_Left.size();
	for(int j=0; j<n1; j++)
		vector<int>().swap(rd_Left[j]);
	int n2 = rd_Right.size();
	for(int j=0; j<n2; j++)
		vector<int>().swap(rd_Right[j]);
	vector<vector<int>>().swap(rd_Left);
	vector<vector<int>>().swap(rd_Right);
	vector<bool>().swap(modified_Flag);
#pragma endregion

	return flag;
}

// start point and the end point of the isoform
void Gene_Batch::bounding_Assemble(ifstream& infile, vector<int>& exon_Start1, vector<int>& exon_Stop1)
{
	//ifstream infile;
	//infile.open(filename,ios::in);  // load file
	if(!infile)
		cout<<"Open file error!"<<endl;

	// bam file format£º1.read name 2.flags describing the result of align 3.reference sequence name
	// 4.start point of left end read 5.mapping quality 6.detail description of mapping quality
	// 7.pair end alignment"=",single end alignment"*"  8.start point of another read 9.sequence length 
	// 10.read sequence 11.read quality 12. description of mapping

	int x1 = 0, x2 = 0, z1 = 0, z2 = 0, temp_Len = 0; // start point of the left read and the right read£¬length of fragment
	int read_len = rd_Len;  // length of read 

	stringstream ss;
	ss<<read_len;  // transform the length of read into string
	string len_R = ss.str()+"M"; //matching string
	string line, word, words[6], frag_Name;  // string
	int serial = 0;     // read index
	string tPar = "N";  // flag for isoform
	string::size_type pos, pos0;
	vector<int>::iterator itr;

	int junc_Len = 0, cover_Len = 0, p1 = 0, p2 = 0;
	int x_1 = 0, x_2 = 0, z_1 = 0, z_2 = 0;
	int pre_Left = 0, pre_Right = 0;
	bool anno_flag = false;  // flag of the boundary of the two-end read according to the matching information of a single read
	int line_cnt = 0, line_cnt1 = 0, line_cnt2 = 0;
	bool junc_Left = false, junc_Right = false, flag_Right = false;  // flag the jump in the isoform
	vector<int> right_Boundary;  // right boundary of the isoform
	//vector<int> left_Boundary;   // right boundary of the isoform
	int acc_cnt1 = 0;
	int thresh1 = 20, thresh2 = 20, len_Thresh = 5;
	vector<vector<int>> rd_Left;  // number of reads within the exon
	vector<vector<int>> rd_Right; // number of reads within the exon
	int num = exon_Assembly.size();  // number of exons

	int seg_Mark = 0;
	int exon_start = exon_Assembly[0].bound1, exon_stop = exon_Assembly[0].bound2;  // count from the first interval
	int pre_start = 0, pre_stop = 0;
	vector<int> seg_Start, seg_Stop, exon_Start, exon_Stop; // start point and end point of the interval
	vector<bool> modified_Flag(2*num-1);   // flag for whether the interval is modified
	for(int j=0; j<2*num-1; j++)
		modified_Flag[j] = false;

	int max1 = 50, max2 = 50, min1 = 100, min2 = 100;
	float jump_Rate1 = 0.25, jump_Rate2 = 0.25;
	int site1 = exon_Assembly[0].bound1-1, site2 = exon_Assembly[exon_Assembly.size()-1].bound2+1;

	for(seg_Mark=0; seg_Mark<2*num-1; seg_Mark++)
	{
		int index = 0;
		if(seg_Mark%2==0) // interval of the exon
		{
			index = (int)seg_Mark/2;
			exon_start = exon_Assembly[index].bound1-2;
			exon_stop = exon_Assembly[index].bound2+2;	
			exon_Start.push_back(exon_start+2);
			exon_Stop.push_back(exon_stop-2);
		}
		else
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

	bool locate_Flag = false;
	bool ctr_Flag = false;
	seg_Mark=0;

#pragma region load the left end and right end of the two-end matching end
	while(/*read_Serial<300&&*/std::getline(infile,line))  // test
	{
		// cout<<line_cnt<<endl;
		line_cnt++;
		Rd cur_Rd;  // exon
		anno_flag = false;
		stringstream stream(line);
		stream>>word; frag_Name = word;
		stream>>word; stream>>word; // load the first 3 strings
		for(int i = 0; i<6; i++)
			stream>>words[i];  // 4 start point of left end read 5 matching level 6 matching level 7 flag of pair-end alignment 8 start point of right end read 9 inserted size
		sscanf(words[5].c_str(),"%d",&temp_Len);  // inserted size
		sscanf(words[0].c_str(),"%d",&x1);  // start point of the left end read
		sscanf(words[4].c_str(),"%d",&z1);  // start point of the right end read
		// cout<<seg_Mark<<"\t"<<x1<<"\t"<<z1<<endl;

		int s1 = _min(x1,z1), s2= _max(x1,z1);
		if(words[3]!="="){
			s1 = x1;
			s2 = x1;
		}
		else
		{
			z2 = 0;
			if(temp_Len>=0)
				z2 = x1 + temp_Len - 1;  // end point of the right end read
			s1 = _min(x1,z1), s2 = _max(_max(x1,z1),z2);
		}

		if(s1>site1&&s2<site2)
		{
			if(temp_Len>=0)        // forward description of the alignment result
			{
				//z2 = x1 + temp_Len - 1;  //  end point of the right end read
				//x2 = x1 + read_len - 1;  //  end point of the left end read
				if(line_cnt1==0)  
				{
					exon_Start[0] = x1-1;  // modify the left end of the gene
					seg_Start[0] = x1-1;   // modify the left end of the gene
					exon_Assembly[0].bound1 = x1-1;
					exon_Assembly[0].len = exon_Assembly[0].bound2 - exon_Assembly[0].bound1 + 1; // modify the length of the gene		
					pre_Left = x1;
					exon_start = seg_Start[0];
					exon_stop = seg_Stop[0];
				}
				line_cnt1++;
				//left_Boundary.push_back(x1);  // left bound of the isoform

				if(x1>=exon_start&&x1<=exon_stop)  // within the interval
				{
					// t_Seg.push_back(x1);
					rd_Left[seg_Mark].push_back(x1);
				}
				else
				{
					seg_Mark++;
					if(seg_Mark>=2*num-1){
						cout<<x1<<"\t"<<z1<<" :Error!"<<endl;
						continue;
					}
					locate_Flag = false;
					while(seg_Mark<2*num-1&&!locate_Flag)
					{
						exon_start = seg_Start[seg_Mark];
						exon_stop = seg_Stop[seg_Mark];
						if(x1>=exon_start&&x1<=exon_stop)  // within the interval
						{
							locate_Flag = true;
							rd_Left[seg_Mark].push_back(x1);
							break;
						}
						seg_Mark++;
					}
				}
			}
			else  // inverse description of the alignment result
			{
				if(words[3]=="=")   // determine the right boundary of the isoform according to the two-end matching
				{  
					int temp = z1;
					z1 = x1;    // left point of the right end read
					x1 = temp;  // left point of the left end read
					
					if(words[2]==len_R)   // matching string
					{
						z2 = z1+read_len-1;
					} 
					else
					{ 
						string::size_type pre;   // find the location of the last "N"
						pre = 0;
						int p1 = 0, p2 = 0, l1 = 0, l2 = 0;
						bool ctr_Flag = false;
						junc_Len = 0;
						cover_Len = 0;
#pragma region calculate the distance span of the read
						while((pos=words[2].find(tPar,pre))!=string::npos)  // find "N" in the isoform
						{
							pos0 = words[2].find("M",pre); // locate the first separator
							if(pos0!=string::npos)
							{ 
								string str = words[2].substr(pos0+1,pos-pos0-1);
								sscanf(str.c_str(),"%d",&l1); // distance span by the read
								str = words[2].substr(pre,pos0-pre);
								sscanf(str.c_str(),"%d",&l2); // length covered by the read 
								junc_Len += l1;
								cover_Len += l2;
								pre = pos + 1;
							}
						}
						if((pos0 = words[2].find("M",pre))!=string::npos)  // locate the first separator
						{
							string str1 = words[2].substr(pre,pos0-pre);
							sscanf(str1.c_str(),"%d",&l2); // length covered by the read 
							cover_Len += l2;
						}
#pragma endregion 
						z2 = z1+junc_Len+read_len-1;
					}
				}
				right_Boundary.push_back(z2);   // right boundry of the isoform
			}
		}
		else if(s1>=site2&&s2>=site2)  // left end and right end both exceed the boundary
		{
			cout<<"beyond the range"<<endl;
			break;
		}
		else
		{
			/*cout<<"next range"<<endl;
			infile.seekg(pos1);
			return false;*/
		}
	}

#pragma endregion 

#pragma region right end of the isoform
	sort(right_Boundary.begin(),right_Boundary.end());  // sort the right boundry of the isoform
	int num_1 = line_cnt1;   // number of the left end read
	int num_2 = right_Boundary.size();  // number of the right end read
	int right_boundary = right_Boundary[num_2-1]+1;
	exon_Assembly[num-1].bound2 = right_boundary;  // right end of the last exon
	exon_Assembly[num-1].len = exon_Assembly[num-1].bound2-exon_Assembly[num-1].bound1+1;

	exon_Stop[exon_Stop.size()-1] = right_boundary;
	seg_Stop[seg_Stop.size()-1] = right_boundary;

	seg_Mark = 0;
	exon_start = seg_Start[seg_Mark];
	exon_stop = seg_Stop[seg_Mark];

	for(int k=0; k<num_2; k++)
	{
		x1 = right_Boundary[k];
		if(x1>=exon_start&&x1<=exon_stop)  // within the interval
		{
			// t_Seg.push_back(x1);
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
				if(x1>=exon_start&&x1<=exon_stop)  // within the interval
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

#pragma region left boundary
	int seg_Num = 2*num-1;
	int interval = 20;     // length of the interval
	float inter_len = interval*1.0;
	int start = 0, stop = 0;
	vector<int> site_New1, site_New1a;  // identify possbile new start points
	vector<int> site_New2; // identify possbile new boundaries of the exons
	vector<Exon_Node> sus_Intron;      // possible introns
	vector<Exon_Node>::iterator exon_itr;  // vector of exons
	vector<int> seg_Exons;  // number of exons in each interval

	int insert_idx = 0;   // location for inserting exons
	int Num = 2*num-1;    // number of intervals
	int* cnt_seg = new int[Num];  // number of exons in the former exon area
	for(int j=0; j<Num; j+=2)
		cnt_seg[j] = 1;   // number of exons in the former exon area
	for(int j=1; j<Num; j+=2)
		cnt_seg[j] = 0;   // number of introns in the former exon area

#pragma region find the possible exons in the current exon
	for(int k = 0; k<seg_Num; k +=2)  // current exons
	{
		site_New1.clear(); site_New1a.clear(); site_New2.clear();		
		insert_idx += cnt_seg[k];  // number of exons

		vector<int>& temp = rd_Left[k];
		vector<int>& temp1 = rd_Right[k];
		int rd_num = temp.size();     // number of left end read in the current interval
		int rd_num1 = temp1.size();   // number of right end read in the current interval
		exon_start = seg_Start[k];    // start point and end point of the exon
		exon_stop = seg_Stop[k];
		int seg_len  = exon_stop-exon_start+1;
		int cnt_thresh = _min(10,0.05*seg_len);
		if(rd_num<cnt_thresh&&rd_num1<cnt_thresh)  // delete current exon if the number of read is small
		{
			exon_Assembly.erase(exon_Assembly.begin()+k/2);
			cnt_seg[k/2] = 0;
			insert_idx--;  // modify the exon index
			continue;
		}
		if(temp.size()<10||seg_len<10)  // ignore the read if the number of read is small or the length of the read is short
			continue;
		int inter_Num = _max(1,(int)(seg_len*1.0/inter_len+0.5));  // number of intervals
		if(abs(temp[0]-exon_start)<5)
			pre_Left = temp[0];   // start point of the exon
		else if(k>0)
			pre_Left = exon_start;
		else
		{
		}
		// float inter_len = seg_len*1.0/inter_Num;  // length of intervals
		junc_Left = false;  // flag of hops
		int* seg_Abun = new int[inter_Num];  // number of read within the interval
		//vector<int> seg_Abun(inter_Num);
		for(int j=0; j<inter_Num; j++)
			seg_Abun[j] = 0;

#pragma region identify the possible boundary of exon and start point according to the discontinuity of abundance
		int inter_Idx = 0;
		start = exon_start+inter_len*inter_Idx-1;    // start point of the interval
		stop = exon_start+inter_len*(inter_Idx+1)+1; // end point of the interval
		float jump_Thresh = _min(min1,seg_len*jump_Rate1);  // threshold of the hop
		if(seg_len<200)
			jump_Thresh = 50;
		for(int k1=0; k1<rd_num; k1++)
		{
			x1 = temp[k1];
			//cout<<k1<<"\t"<<x1<<endl;
			if(abs(x1-pre_Left)>jump_Thresh)
			{
				junc_Left = true;
				site_New1.push_back(x1);  // possible new start point
				site_New1a.push_back(pre_Left);
				//cout<<pre_Left<<"\t"<<x1<<endl;
			}
			pre_Left = x1;
			if(x1<=stop&&x1>=start&&inter_Idx<inter_Num)
				seg_Abun[inter_Idx]++;
			else
			{
				while(x1>stop&&inter_Idx<inter_Num)
				{
					inter_Idx++;
					if(inter_Idx==inter_Num)  // avoid exceeding the boundary
						break;
					start += inter_len;   // start point of the interval
					stop += inter_len;    // end point of the interval
				}
				if(x1<=stop&&x1>=start&&inter_Idx<inter_Num)
					seg_Abun[inter_Idx]++;
			}
		}

		int inter_Num1 = 0;
		for(int j=0; j<inter_Num; j++)
		{
			if(seg_Abun[j]>0)
				inter_Num1++;
		}

		float aver = rd_num*1.0/inter_Num1;   // average over the interval
		float thresh = aver*0.5;
		if(inter_Num>1)  // number of the interval
		{
			float* diff = new float[inter_Num-1]; // number of the interval		
			for(int j=0; j<inter_Num-2; j++)
			{
				diff[j] = abs(seg_Abun[j+1]-seg_Abun[j]);
				if(seg_Abun[j]>2&&diff[j]*1.0>aver/*diff[j]*1.0/(seg_Abun[j]+1e-09)>0.5*/)  // frequency of the read significantly increases
					site_New2.push_back(exon_start+inter_len*(j+1));  // possible new start point of the exon
			}
			delete[] diff;
		}
#pragma endregion identify the possible boundary of exon and start point according to the discontinuity of abundance

		int n1 = site_New1.size();  // start point of the newly located exon
		int n2 = site_New2.size();  // start point of the newly located exon

#pragma region check boundary of the exon according to the breakpoint

		for(int i=0; i<n1; i++)  
		{
			int x = site_New1[i];
			int pre_x = site_New1a[i];
			int j1 = 0, j = 0, mark = 0;
			while(j1<inter_Num&&(exon_start+inter_len*(j1+1))<pre_x)
				j1++;
			j = j1;
			while(j<inter_Num&&(exon_start+inter_len*(j+1))<x)
				j++;

			int bound2 = exon_Stop[k/2];  // initial boundary
			int bound1 = exon_Start[k/2]; // initial boundary
			// cout<<bound1<<"\t"<<bound2<<endl;
			float thresh_1 = _max(thresh,1);
			// cout<<j<<"\t"<<inter_Num<<"\t"<<seg_Abun[j]<<"\t"<<seg_Abun[j+1]<<"\t"<<thresh<<endl;
			float acc = -1;
			if(j<inter_Num-1)
				acc = seg_Abun[j]+seg_Abun[j+1];
			if((acc>thresh&&(abs(pre_x-exon_start)<5||seg_Abun[j1]>thresh))||(j>=inter_Num-1/*&&seg_Abun[j]>thresh_1*/)) // the threshold of the last interval is decreased for noise
			{
				modified_Flag[k] = true;
				//cout<<k+1<<" Modify"<<endl;
				Exon_Node cur_Exon;      // added exon
				cur_Exon.bound1 = x;   
				cur_Exon.bound2 = bound2;
				cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
				int exon_idx = k/2;
				Exon_Node cur_Exon1;
				cur_Exon1.bound1 = pre_x+1;
				cur_Exon1.bound2 = x-1;
				cur_Exon1.len = cur_Exon1.bound2 - cur_Exon1.bound1 + 1;
				sus_Intron.push_back(cur_Exon1);   // deleted exon
				int z1 = 0, z2 = 0;
				bool del_Flag = check_Intron(pre_x+1,x-1,z1,z2,rd_Right[k]);   // check whether the exon should be deleted
				int pre_idx = insert_idx-1;   // index od the previous exon
				if(del_Flag&&cur_Exon.len>10)
				{
					if(k/2<exon_Assembly.size())
					{
						exon_Assembly.insert(exon_Assembly.begin()+insert_idx,cur_Exon);  // add new exon
						insert_idx++;
					}
					else
						exon_Assembly.push_back(cur_Exon);   // add new exon
					cnt_seg[k]++;					
				}
				else if(!del_Flag)
				{
					cur_Exon1.bound1 = abs(z1-pre_x)>5?z1:pre_x+1;
					cur_Exon1.bound2 = abs(z2-x)>5?z2:x-1;
					cur_Exon1.len = cur_Exon1.bound2 - cur_Exon1.bound1 + 1;

					if(k/2<exon_Assembly.size()-1)
					{
						//insert_idx++;   // location for inserting the new exon
						if(cur_Exon1.len>len_Thresh){
							exon_Assembly.insert(exon_Assembly.begin()+insert_idx,cur_Exon1);  // insert new exon
							insert_idx++; cnt_seg[k]++;
						}
						if(cur_Exon.len>len_Thresh){
							exon_Assembly.insert(exon_Assembly.begin()+insert_idx,cur_Exon);  // insert new exon
							insert_idx++; cnt_seg[k]++;
						}
					}
					else
					{
						if(cur_Exon1.len>len_Thresh){
							exon_Assembly.push_back(cur_Exon1);   // insert new exon
							cnt_seg[k]++;
						}
						if(cur_Exon.len>len_Thresh){
							exon_Assembly.push_back(cur_Exon);
							cnt_seg[k]++;
						}
					}
				}

				// find the nearest junction
				j = 0;
				int n1 = exon_Stop.size();  // junction detected
				// exon_Assembly[pre_idx].bound1 = bound1;
				while(j<n1&&abs(exon_Stop[j]-pre_x)>10)
					j++;
				if(j<n1)  // find the nearest junction
					exon_Assembly[pre_idx].bound2 = exon_Stop[j];  // ends of the exon identified by the junction
				else
					exon_Assembly[pre_idx].bound2 = pre_x;      // relocate the end of the exon
				exon_Assembly[pre_idx].len = exon_Assembly[pre_idx].bound2 - exon_Assembly[pre_idx].bound1 + 1;
				if(exon_Assembly[pre_idx].len<=5)  // if the length of the exon is too short, delete the exon
				{
					exon_Assembly.erase(exon_Assembly.begin()+pre_idx);
					insert_idx--;
					cnt_seg[k]--;  // decrease the count
				}
			}
		}		
#pragma endregion


		delete[] seg_Abun;
		// vector<int>().swap(seg_Abun);
	}
#pragma endregion find the exons in the current exon

#pragma region find the exons in the current intron

	for(int k=1; k<seg_Num; k+=2)
	{
		int intron_idx = (k-1)/2;  // intron index
		int site1 = seg_Start[k]-2, site2 = seg_Stop[k]+2;  // compensate the misalienment
		int len = site2-site1+1;

		vector<int>& temp1 = rd_Left[k];
		vector<int>& temp2 = rd_Right[k];
		int n1 = temp1.size(), n2 = temp2.size();
		int s1 = 0;
		for(int l=0; l<=k; l++)
			s1 += cnt_seg[l];

		int x1 = 1e09, x2 = 1e09, y1 = 0, y2 = 0;
		float thresh3 = _max(5,0.03*len);   // threshold of the intron
		bool flag1 = true, flag2 = true;
		if(n1>0){
			x1 = temp1[0]; y1 = temp1[n1-1];
		}
		if(n2>0){
			x2 = temp2[0]; y2 = temp2[n2-1];
		}
		if((n1+n2)>thresh3&&(flag1||flag2)) // exons in the former introns
		{
			int start = _min(x1,x2);
			int stop = _max(y1,y2);
			Exon_Node cur_Exon;  // add exons

			if(abs(site1-start)<5)
				cur_Exon.bound1 = site1;
			else
				cur_Exon.bound1 = start;
			if(abs(site2-stop)<5)
				cur_Exon.bound2 = site2;
			else
				cur_Exon.bound2 = stop;

			cur_Exon.len = cur_Exon.bound2-cur_Exon.bound1+1;
			exon_Assembly.insert(exon_Assembly.begin()+s1,cur_Exon);
			cnt_seg[k]++;  // increase the number of the exon in the interval 
		}
	}
#pragma endregion find the exons in the current intron

#pragma endregion

#pragma region right boundary
	// sort(right_Boundary.begin(),right_Boundary.end());  // sort the right boundry of the isoform
	int num2 = exon_Stop.size();  // number of the right end of the exons
	int rd_num1 = right_Boundary.size();  // detected number of the right end of the isoform
	num = exon_Assembly.size();   // number of exons

	pre_Right = right_Boundary[0];
	int thresh = 20;
	int acc_cnt = 0;
	vector<int> new_exonStop, new_exonStart;  // start point and end point of the new exon
	int exon_Num = exon_Assembly.size();  // number of exons
	site_New1.clear(); site_New1a.clear(); site_New2.clear();   // clear the array
	interval = 20;
	inter_len = interval*1.0;
	vector<int> acc_vcnt;
	acc_cnt1 = 0;

	for(int k=0; k<seg_Num; k+=2)  // current exons
	{
		site_New1.clear(); site_New1a.clear(); site_New2.clear(); 
		acc_vcnt.clear();

		if(/*modified_Flag[k]==false*/cnt_seg[k]==1)   // exons only split once
		{
			vector<int>& temp = rd_Right[k];
			exon_start = seg_Start[k]; // start point and end point of the exon
			exon_stop = seg_Stop[k];
			int seg_len  = exon_stop-exon_start+1;
			int cnt_thresh = _min(10,0.1*seg_len);
			if(temp.size()<cnt_thresh||seg_len<10)  // ignore the read if the number of read is small or the length of the read is short
				continue;

			temp.push_back(exon_stop);  // end point of the added exon
			int rd_num = temp.size();   // number of read within the current interval
			//int inter_Num = (int)num*1.0/interval;  // number of interval 
			//float inter_len = seg_len*1.0/inter_Num;  // length of interval 
			int inter_Num = _max(1,(int)(seg_len/inter_len+0.5));
			junc_Right = false;  // flag of hops
			int* seg_Abun = new int[inter_Num];  // number of read within the current interval
			//vector<int> seg_Abun(inter_Num);
			for(int j=0; j<inter_Num; j++)
				seg_Abun[j] = 0;

#pragma region  identify the possible boundary of exon and start point according to the discontinuity of abundance
			int pre = temp[0];
			int inter_Idx = 0;
			start = exon_start+inter_len*inter_Idx-1;    // start point of the interval
			stop = exon_start+inter_len*(inter_Idx+1)+1; // end point of the interval
			float jump_Thresh = _min(min2,seg_len*jump_Rate2);  // threshold of the hops
			if(seg_len<200)
				jump_Thresh = 50;
			pre_Right = temp[0];

			for(int k1=0; k1<rd_num; k1++)
			{
				x1 = temp[k1];
				if(abs(x1-pre_Right)>jump_Thresh)
				{
					acc_vcnt.push_back(acc_cnt1);
					acc_cnt1 = 0;
					junc_Right = true;
					site_New1.push_back(pre_Right);  // possible new end point
					site_New1a.push_back(x1);
				}
				else
					acc_cnt1++;
				pre_Right = x1;
				if(x1<=stop&&x1>=start&&inter_Idx<inter_Num)
					seg_Abun[inter_Idx]++;
				else
				{
					while(x1>stop&&inter_Idx<inter_Num)
					{
						inter_Idx++;
						if(inter_Idx==inter_Num)  // avoid exceeding the boundary
							break;
						start += inter_len;   // start point of the interval
						stop += inter_len;    // end point of the interval
					}
					if(x1<=stop&&x1>=start&&inter_Idx<inter_Num)
						seg_Abun[inter_Idx]++;
				}
			}

			int inter_Num1 = 0;
			for(int j=0; j<inter_Num; j++)
			{
				if(seg_Abun[j]>0)
					inter_Num1++;
			}

			float aver = rd_num*1.0/inter_Num1;    // average over the interval
			float thresh = aver*0.2;
			if(inter_Num>1)
			{
				float* diff = new float[inter_Num-1];  // number of the interval	
				for(int j=0; j<inter_Num-2; j++)
				{
					diff[j] = abs(seg_Abun[j+1]-seg_Abun[j]);
					if(seg_Abun[j]>2&&diff[j]>aver)    // frequency of the read significantly increases
						site_New2.push_back(exon_start+inter_len*(j+1));  // possible new start point of the exon
				}
				delete[] diff;
			}
#pragma endregion 

			int n1 = site_New1.size();  //  start point of the newly located exon
			int n2 = site_New2.size();  //  start point of the newly located exon

#pragma region check boundary of the exon according to the breakpoint
			insert_idx = 0;
			for(int l=0; l<=k; l++)
				insert_idx += cnt_seg[l];  // sum over the number of the exon of each interval

			for(int i=0; i<n1; i++)  
			{
				int pre_x = site_New1[i];
				int x = site_New1a[i];
				int j = 0;
				while(j<inter_Num&&(exon_start+inter_len*(j+1))<pre_x)  // find the segment which contains the break point
					j++;
				int bound2 = exon_Stop[k/2];   // initial boundary
				int bound1 = exon_Start[k/2];  // initial boundary
				if((j<inter_Num-1&&acc_vcnt[i]>thresh)||j==0) // the threshold of the last interval is decreased for noise
				{
					Exon_Node cur_Exon;      // added exon
					cur_Exon.bound1 = x;   
					cur_Exon.bound2 = bound2;
					cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
					int exon_idx = k/2;
					Exon_Node cur_Exon1;
					cur_Exon1.bound1 = pre_x;
					cur_Exon1.bound2 = x-1;
					cur_Exon1.len = cur_Exon1.bound2 - cur_Exon1.bound1 + 1;
					sus_Intron.push_back(cur_Exon1);   // deleted exon
					int z1 = -1, z2 = -1;
					bool del_Flag = check_Intron(pre_x+1,x-1,z1,z2,rd_Left[k]);  // check whether the exon should be deleted
					int pre_idx = insert_idx-1;   // index od the previous exon
					if(del_Flag&&cur_Exon.len>10)
					{
						if(k/2<exon_Assembly.size()-1)
						{
							exon_Assembly.insert(exon_Assembly.begin()+insert_idx,cur_Exon);  //  add new exon
						}
						else
							exon_Assembly.push_back(cur_Exon);   //  add new exon
						cnt_seg[k]++;    // the number of the exons increases
					}
					else if(!del_Flag)
					{
						cur_Exon1.bound1 = abs(z1-pre_x)>5?z1:pre_x+1;
						cur_Exon1.bound2 = abs(z2-x)>5?z2:x-1;
						cur_Exon1.len = cur_Exon1.bound2 - cur_Exon1.bound1 + 1;

						if(k/2<exon_Assembly.size()-1)
						{
							//insert_idx++;   // location for inserting the exon
							if(cur_Exon1.len>len_Thresh){
								exon_Assembly.insert(exon_Assembly.begin()+insert_idx,cur_Exon1);  // insert the exon
								insert_idx++; cnt_seg[k]++;
							}
							if(cur_Exon.len>len_Thresh){
								exon_Assembly.insert(exon_Assembly.begin()+insert_idx,cur_Exon);  // insert the exon
								insert_idx++; cnt_seg[k]++;
							}
						}
						else
						{
							if(cur_Exon1.len>len_Thresh){
								exon_Assembly.push_back(cur_Exon1);   // insert the exon
								cnt_seg[k]++;
							}
							if(cur_Exon.len>len_Thresh){
								exon_Assembly.push_back(cur_Exon);
								cnt_seg[k]++;
							}
						}
					}					

					// find the nearest junction
					j = 0;
					int n1 = exon_Stop.size();  // junction detected
					exon_Assembly[pre_idx].bound1 = bound1;
					while(j<n1&&abs(exon_Stop[j]-pre_x)>10)
						j++;
					if(j<n1)  // find the nearest junction
						exon_Assembly[pre_idx].bound2 = exon_Stop[j];  // ends of the exon identified by the junction
					else
						exon_Assembly[pre_idx].bound2 = pre_x;      // relocate the end of the exon
					exon_Assembly[pre_idx].len = exon_Assembly[pre_idx].bound2 - exon_Assembly[pre_idx].bound1 + 1;
					if(exon_Assembly[pre_idx].len<=5)  // if the length of the exon is too short, delete the exon
					{
						exon_Assembly.erase(exon_Assembly.begin()+pre_idx);
						insert_idx--;
						cnt_seg[k]--;  // decrease the count
					}
				}
			}		
#pragma endregion

			delete[] seg_Abun;
			// vector<int>().swap(seg_Abun);
		}
	}
#pragma endregion 

#pragma region clear array and close file
	//vector<int>().swap(left_Boundary);
	vector<int>().swap(right_Boundary);
	vector<int>().swap(exon_Start);
	vector<int>().swap(exon_Stop);
	vector<int>().swap(seg_Start);
	vector<int>().swap(seg_Stop);
	vector<Exon_Node>().swap(sus_Intron);
	int n1 = rd_Left.size();
	for(int j=0; j<n1; j++)
		vector<int>().swap(rd_Left[j]);
	int n2 = rd_Right.size();
	for(int j=0; j<n2; j++)
		vector<int>().swap(rd_Right[j]);
	vector<vector<int>>().swap(rd_Left);
	vector<vector<int>>().swap(rd_Right);
	vector<bool>().swap(modified_Flag);
	//infile.close();
#pragma endregion
}

// calculate the hop score
int Gene_Batch::cal_Rate(vector<int>& site, vector<int>& cnt, vector<int>& gap, int num, vector<float>& rate, int& m_idx1, int& m_idx2)
{
	m_idx1 = -1; m_idx2 = -1;
	int m_gap = 0;
	float m_rate = 0;
	for(int k=0; k<num; k++)
	{
		rate[k] = (cnt[k]+cnt[k+1]);  // number of read in the neighboring interval
		if(gap[k]>m_gap)
		{
			m_gap = gap[k];  // maximum hop
			m_idx2 = k;
		}
	}
	for(int k=0; k<num; k++)  // calculate the rate of the cover rate and the hop
	{
		rate[k] = rate[k]*exp(0.5*(gap[k]/m_gap-1));
		if(rate[k]>m_rate)
		{
			m_rate = rate[k];
			m_idx1 = k;
		}
	}
	return m_gap;
}

// calculate the hop score
float Gene_Batch::cal_Rate_1(int serial, vector<int>& site, vector<int>& cnt, vector<int>& gap, int num, vector<float>& rate, int& m_idx1, int& m_idx2)
{
	m_idx2 = -1;
	int m_gap = -1;
	float m_rate = -1;  // maximum hop
	float thresh = 10;
	
	for(int k=0; k<num; k++)
	{
		rate[k] = (cnt[k]+cnt[k+1]);  // number of read in the neighboring interval
		if(gap[k]>m_gap)
		{
			m_gap = gap[k];  // maximum hop
			m_idx1 = k;
		}
	}
	
	if(num>0)
	{
		if(abs(site[0]-exon_Assembly[serial].bound1)<thresh&&cnt[0]<5)  // hop point of the left boundary
			rate[0] = _max(rate[0],2*cnt[1]);
		if(abs(site[num-1]-exon_Assembly[serial].bound2)<thresh&&cnt[num]<5)  // hop point of the right boundary
			rate[num-1] = _max(rate[num-1],2*cnt[num-1]);
	}
	for(int k=0; k<num; k++)
	{
		rate[k] = rate[k]*exp(0.5*(gap[k]/m_gap-1));
		if(rate[k]>m_rate)
		{
			m_rate = rate[k];
			m_idx2 = k;
		}
	}
	return m_rate;
}

// find the possible hop and locate the hop
// input£ºdetected hop point in each interval
// output£ºpossible location of hop
// match_Identify() check whether the left end read and right end read match
int Gene_Batch::interval_Locate(Interval& seg, vector<vector<int>>& seg_Modify)
{
	vector<PAIR_INT> pair_Vec1, pair_Vec2;
	int num = seg.serial.size();  // number of exons which have possible break point
	int serial = 0, change_type = 0;
	int exon_Num = exon_Assembly.size();  // number of exons 
	
	for(int i=0; i<num; i++)
	{
		serial = seg.serial[i];   // interval index
		change_type = seg.ch_type[i];
		int p1 = 0, p2 = 0, p3 = 0, p4 = 0;       // new start point and end point

		vector<int>& temp1 = seg.site_New1[i];  // break point indentified by the left end read
		vector<int>& temp2 = seg.site_New2[i];  // break point indentified by the right end read

		int m_idx1 = 0, m_idx2 = 0, m_idx3 = 0, m_idx4 = 0, m_gap1 = 0, m_gap2 = 0;
		float m_rate1 = 0, m_rate2 = 0;
		int s1 = exon_Assembly[serial].bound1, s2 = exon_Assembly[serial].bound2;  // start point and end point of the gene

		float thresh = exon_Assembly[i].len*0.75;   // threshold of the length of exons
		float thresh_s = 50;   // threshold of the length of exons

		int n1 = temp1.size()/2, n2 = temp2.size()/2;   // number of hop points
		vector<float> rate1(n1), rate2(n2);

		vector<int> match1(n1), match2(n2);  // matching hop points
		/*vector<float> fide1(n1), fide2(n2);
		for(int j=0; j<n1; j++)
		{ fide1[j] = -1; match1[j] = -1;}
		for(int j=0; j<n2; j++)
		{ fide2[j] = -1; match2[j] = -1;}*/

		int add_Cnt = 0;   // number of added exons
		bool flag1 = false;
		vector<int> sel_IdxL, sel_IdxR;  // index of selected hops
		pair_Vec1.clear(); pair_Vec2.clear();
		if(n1>0)
		{
#pragma region find possible split point
			for(int j = 0; j<n1; j++)
				pair_Vec1.push_back(make_pair(j,seg.gap1[i][j]));
			sort(pair_Vec1.begin(),pair_Vec1.end(),cmp_1);  // sort by gap
			m_rate1 = cal_Rate_1(serial,temp1,seg.cnt_New1[i],seg.gap1[i],n1,rate1,m_idx1,m_idx2);
			m_idx1 = pair_Vec1[0].first;  // index of maximum hop
			m_gap1 = pair_Vec1[0].second; // length of maximum hop
			int cnd_idx = -1;
			if(n1>1)
				cnd_idx = pair_Vec1[1].first;
			
			// sel_IdxL[0] = m_idx1;  // index of the first hop	point
			sel_IdxL.push_back(m_idx1);
			// select between the maximum hop and the maximum cover rate
			if(m_idx1==m_idx2)   // the maximum hop and the maximum cover rate are the same
			{
				add_Cnt++;
				if(n1>1)
					sel_IdxL.push_back(cnd_idx); // add the second maximum hop
			}
			else
			{
				sel_IdxL.push_back(m_idx2);
				if(n1>1&&m_idx2!=cnd_idx)
					sel_IdxL.push_back(cnd_idx);
			}
#pragma endregion	
		}
		if(n2>0)   // identify by the right end
		{
#pragma region find possible split point
			for(int j=0; j<n2; j++)
				pair_Vec2.push_back(make_pair(j,seg.gap2[i][j]));
			sort(pair_Vec2.begin(),pair_Vec2.end(),cmp_1);  // sort by gap
			m_rate2 = cal_Rate_1(serial,temp2,seg.cnt_New2[i],seg.gap2[i],n2,rate2,m_idx3,m_idx4);
			m_idx3 = pair_Vec2[0].first;  // index of maximum hop
			m_gap2 = pair_Vec2[0].second; // length of maximum hop
			int cnd_idx = -1;
			if(n2>1)
				cnd_idx = pair_Vec2[1].first;

			sel_IdxR.push_back(m_idx3);
			// select between the maximum hop and the maximum cover rate
			if(m_idx3==m_idx4)   // the maximum hop and the maximum cover rate are the same
			{
				if(n2>1)
					sel_IdxR.push_back(cnd_idx);
			}
			else
			{
				sel_IdxR.push_back(m_idx4);
				if(n2>1&&m_idx4!=cnd_idx)
					sel_IdxR.push_back(cnd_idx);
			}
#pragma endregion find possible split point
		}

#pragma region check whether the possible hop points identified by the left end read and right end read match
		int m1 = sel_IdxL.size(), m2 = sel_IdxR.size();  // most possible hop point in the left end and the right end
		vector<int> match(m1);
		vector<int> valid_Cnt;
		match.clear();
		for(int j=0; j<m1; j++)
		{
			bool match_Flag = false;
			for(int l=0; l<m2; l++)
			{
				match_Flag = match_Identify(sel_IdxL[j],sel_IdxR[l],exon_Type[serial],temp1,temp2);
				if(match_Flag==true)
				{
					match.push_back(sel_IdxR[l]);
					valid_Cnt.push_back(j);
					break;
				}
			}
			if(match_Flag==false)
				match.push_back(-1);
		}
#pragma endregion

#pragma region find the possible boundary of the exon according to the hop point
		vector<int> points_L, points_R, points;
		int idx = -1, idx_2 = -1, idx_s = -1;
		if(valid_Cnt.size()>0)  // matching of the hop point in the left end and the right end 
		{	
#pragma region matching of the hop points of the left end read and right end read
			if(valid_Cnt.size()==1)
			{
				idx = sel_IdxL[valid_Cnt[0]];   // find the first location
				idx_2 = match[valid_Cnt[0]];    // index in the right end
				points_Locate(idx,idx_2,serial,s1,temp1,temp2,points);
			}
			else
			{
				bool ctr_Flag = true;
				int valid_Num = valid_Cnt.size();
				if(serial==exon_Num-1)
				{
#pragma region  last end exon
					ctr_Flag = false;
					vector<PAIR_INT> idx_vec;
					int idx = -1, idx_2 = -1;				
					for(int l=0; l<valid_Num; l++)
						idx_vec.push_back(make_pair(valid_Cnt[l],sel_IdxL[valid_Cnt[l]]));  // index in the left end
					sort(idx_vec.begin(),idx_vec.end(),cmp_1);  // sort by the index
					for(int l=valid_Num-1; l>=0; l--)  // select from the nearest
					{
						int t_idx = idx_vec[l].second;
						if(abs(exon_Assembly[serial].bound1-temp1[2*t_idx])>100&&rate1[t_idx]>m_rate1*0.5)  // the distance to the boundary should be far and the cover rate should be high
						{
							idx = t_idx;
							idx_2 = match[idx_vec[l].first];
							break;
						}
					}
					if(idx>=0)
						points_Locate(idx,idx_2,serial,s1,temp1,temp2,points);
					else    // normal mode
					{
						ctr_Flag = true;
					}
#pragma endregion 
				}
				if(ctr_Flag==true)    // normal mode
				{	
#pragma region last end exon
					idx = -1; idx_2 = -1;
					int l = 0;
					for(l=0; l<valid_Num; l++)  // find the first hop point which satisfy the cover rate
					{
						idx = sel_IdxL[valid_Cnt[l]];
						idx_2 = match[valid_Cnt[l]];
						if(rate1[idx]>m_rate1*0.6)
							break;
					}
					if(l>=valid_Num)  // select by the hop value
					{
						idx = sel_IdxL[valid_Cnt[0]];   // the first location
						idx_2 = match[valid_Cnt[0]];    // index of the right end
					}
					points_Locate(idx,idx_2,serial,s1,temp1,temp2,points);
					if(change_type!=0)    // only one end of the exon is identified
					{
						idx_s = -1;
						for(l=0; l<valid_Num; l++)  // find the first hop point which satisfy the cover rate
						{
							idx_s = sel_IdxL[valid_Cnt[l]];
							idx_2 = match[valid_Cnt[l]];
							int temp = _min(1000,seg.gap1[i][idx]*0.8);
							if((idx_s!=idx)&&(rate1[idx_s]>m_rate1*0.5)&&(seg.gap1[i][idx_s]>temp))
								break;
						}
						if(l<valid_Num)
						{
							vector<int> points_1;
							points_Locate(idx_s,idx_2,serial,s1,temp1,temp2,points_1);
							points.push_back(points_1[1]); points.push_back(points_1[2]);
						}
					}
#pragma endregion
				}
			}
			seg_Modify.push_back(points);
#pragma endregion
		}
		else   // hop points of the left end and the right end mismatch
		{
#pragma region hop points identified by the left and right end read mismatch
			int idx = -1, chosen = -1;
			float gap_thresh = 300, cnt_thresh = 10;  // select the hop with high confidence
			bool mark = false;
			if(n1==0)
			{
				if(serial==exon_Num-1||(m_gap2>gap_thresh&&(seg.cnt_New2[i][m_idx3]+seg.cnt_New2[i][m_idx3+1]>cnt_thresh)))
				{
					idx = m_idx3; chosen = 1;
				}
			}
			else 
			{
				if(serial==exon_Num-1||(m_gap1>gap_thresh&&(seg.cnt_New1[i][m_idx1]+seg.cnt_New1[i][m_idx1+1]>cnt_thresh)))
				{
					idx = m_idx1; chosen = 0;
					if(n2>0)
					{
						if(m_gap2>m_gap1&&seg.cnt_New2[i][m_idx3]>seg.cnt_New1[i][m_idx1])
						{
							idx = m_idx3; chosen = 1;
						}
					}
				}				
			}
			if(idx>-1)
			{
				if(chosen==0)
				{
					p1 = temp1[2*idx]; p2 = temp1[2*idx+1];
				}
				else
				{
					p1 = temp2[2*idx]; p2 = temp2[2*idx+1];
				}
				points.push_back(serial); points.push_back(p1); points.push_back(p2);
				seg_Modify.push_back(points);
			}
#pragma endregion
		}

		if(serial==exon_Num-1)
		{
			starting_site = points[2];  // select the boundary of the next gene
		}
#pragma endregion
	}

	return seg_Modify.size();
}

// check whether the two hop points match
bool Gene_Batch::match_Identify(int idx_1, int idx_2, int change_type, vector<int>& temp1, vector<int>& temp2)
{
	float match_thresh1 = 100, match_thresh2 = 120;
	float match_thresh3 = 0.7*fragLen_Mean, match_thresh4 = 0.8*fragLen_Mean;
	int l_1 = temp1[2*idx_1], l_2 = temp1[2*idx_1+1];    // hop point
	int r_1 = temp2[2*idx_2], r_2 = temp2[2*idx_2+1], min_idx = 0;
	bool ctr_Flag = false;
	float dis_L = r_1 - l_1, dis_R = r_2 - l_2;
	float dis_1 = abs(dis_L-0.7*fragLen_Mean), dis_2 = abs(dis_R-0.7*fragLen_Mean);
	float dis1 = _min(abs(fragLen_Mean-dis_L),abs(dis_L));
	float dis2 = _min(abs(fragLen_Mean-dis_R),abs(dis_R));

	//  one end is strict, and the other slack
	if((dis_1<match_thresh3&&dis_2<match_thresh4)||(dis_2<match_thresh3&&dis_1<match_thresh4))
		return true;
	else if((dis1<match_thresh1&&dis2<match_thresh2)||(dis2<match_thresh1&&dis1<match_thresh2)) // two hop points match
		return true;	

	if(change_type==0||change_type==2)   // accurately locate 5' end
	{
		if(dis_1<match_thresh3||dis1<match_thresh1)
			return true;
	}
	else if(change_type==1)
	{
		if(dis_2<match_thresh3||dis2<match_thresh1)
			return true;
	}

	return ctr_Flag;
}

// locate split point
void Gene_Batch::points_Locate(int idx_1, int idx_2, int serial, int s1, vector<int>& temp1, vector<int>& temp2, vector<int>& points)
{
	int left1 = temp1[2*idx_1], left2 = temp1[2*idx_1+1];
	int right1 = temp2[2*idx_2], right2 = temp2[2*idx_2+1];
	int p1 = (left1+fragLen_Mean)*0.4+right1*0.6;   // get the left hop point by weight
	int temp = right2 - fragLen_Mean;
	int p2 = 0;
	if(temp>s1)
		p2 = temp*0.4 + left2*0.6;
	else
		p2 = left2;
	points.push_back(serial); points.push_back(p1); points.push_back(p2);
}

// find the break point
// input£ºpoints£¬number of points£¬interval index£¬interval boundary£¬threshold of hops£¬threshold of the distance to the boundary
// output£ºdetected hop points£¬number of points in the neighboring points of the hop point, distance of hop points
bool Gene_Batch::locate_Jump(vector<int>& points, int rd_Num, int idx, int bound1, int bound2, float jump_Thresh, float boundary_Thresh,
	                         vector<int>& pair_Jump, vector<int>& cnt, vector<int>&gap)
{
	int pre = 0, temp = 0;
	if(abs(points[0]-bound1)<5||idx==0)
		pre = points[0];   // start point of the exon
	else
		pre = bound1;

	bool add_Flag = false;
	int serial = idx/2;
	if(abs(points[rd_Num-1]-exon_Assembly[serial].bound2)>boundary_Thresh){
		points.push_back(exon_Assembly[serial].bound2);   // add the right end of the interval
		rd_Num++;
		add_Flag = true;
	}

	bool junc = false;   // flag of whether the hops exist
	int num = 0, x = 0;
	for(int k=0; k<rd_Num; k++)
	{
		x = points[k];
		if((temp=x-pre)>jump_Thresh)
		{
			junc = true;						
			pair_Jump.push_back(pre);
			pair_Jump.push_back(x);  // possible new start point
			cnt.push_back(num);
			gap.push_back(temp);
			num = 0;
		}
		else
			num++;
		pre = x;
	}
	cnt.push_back(num);  // number of read in the last interval
	if(add_Flag)
		points.pop_back();  // delete the last element
	
	return true;
}

// transcipt assembly
bool Gene_Batch::exon_Predict(vector<int>& exon_Start, vector<int>& exon_Stop, int start_site, int end_site, 
	                          vector<Exon_Node>& exon_Vec, vector<int>& exon_type, vector<float>& exon_score)
{
	int num1 = exon_Start.size();  // number of start point of the exon
	int num2 = exon_Stop.size();   // number of end point of the exon

	int exon_stop = 0, exon_start = 0;
	sort(exon_Start.begin(),exon_Start.end());  // sort the start point of the exon
	sort(exon_Stop.begin(),exon_Stop.end());    //  sort the end point of the exon

#pragma region duplicate boundries of exon
	vector<int>::iterator iter;
	for(int j = 0; j<num1-1; j++)
	{
		if(abs(exon_Start[j]-exon_Start[j+1])<3)   // neighbor of the start point of the exon
		{
			iter = exon_Start.begin()+j+1;
			exon_Start.erase(iter);
			j--;
			num1--;
		}
	}

	for(int j = 0; j<num2-1; j++)
	{
		if(abs(exon_Stop[j]-exon_Stop[j+1])<3)   // neighbor of the start point of the exon
		{
			iter = exon_Stop.begin()+j+1;
			exon_Stop.erase(iter);
			j--;
			num2--;
		}
	}
#pragma endregion

	if(exon_Stop[num2-1]<exon_Start[num1-1])   // miss the boundary of the exon
		exon_Stop.push_back(exon_Start[num1-1]+2000);

#pragma region transcript assembly
	// loacte the boundary of the exon in the annotation
	Exon_Node cur_Node;	

	// clear the array of exon score
	exon_score.clear();
	exon_type.clear();
	exon_Vec.clear();

	if(exon_Start[0]>exon_Stop[0])
		exon_Start.insert(exon_Start.begin(),exon_Stop[0]-1200);

	num1 = exon_Start.size();  // number of start point of the exon
	num2 = exon_Stop.size();   // number of end point of the exon

	int i = 0, k = 0, exon_idx = 0;
	while(i<num2)
	{
		exon_stop = exon_Stop[i];
		int k1 = k;		
		while(k<num1&&exon_Start[k]<exon_stop)
			k++;
		for(int l=k1; l<k-1; l++)
		{
			Exon_Node cur_Node;
			cur_Node.bound1 = exon_Start[l];
			cur_Node.bound2 = exon_Start[l+1]-1;
			cur_Node.len = cur_Node.bound2-cur_Node.bound1+1;
			exon_Vec.push_back(cur_Node);
			exon_idx++;
		}
		Exon_Node cur_Node;
		cur_Node.bound1 = exon_Start[k-1];
		cur_Node.bound2 = exon_stop;
		cur_Node.len = cur_Node.bound2-cur_Node.bound1+1;
		exon_Vec.push_back(cur_Node);

		int n = k-k1;			
		for(int l=0; l<n-1; l++)
		{
			exon_score.push_back(0.5);  // the left end is determined, while the right end not
			exon_type.push_back(2);
		}
		exon_score.push_back(1);
		exon_type.push_back(3);

		int i1 = i;
		while(i<num2-1&&exon_Start[k]>exon_Stop[i+1])   // right end of the exon is continuous
			i++;
		for(int l=i1; l<i; l++)
		{
			Exon_Node cur_Node;
			cur_Node.bound1 = exon_Stop[l]+1;
			cur_Node.bound2 = exon_Stop[l+1];
			cur_Node.len = cur_Node.bound2-cur_Node.bound1+1;
			exon_Vec.push_back(cur_Node);  // add new exons
		}

		int m = i-i1;
		for(int l=0; l<m; l++)
		{
			exon_score.push_back(0.5);  // the right end is determined, while the left end not
			exon_type.push_back(1);
		}

		i++;
	}
#pragma endregion
	int exon_Num = exon_Vec.size();
	// score of the last exon
	exon_score[exon_Num-1] = 0.5;
	exon_type[exon_Num-1] = 0;

	if(exon_Num>0)
	{
		if(start_site-exon_Vec[0].bound2>fragLen_Mean*1.5)
		{
			/*vector<Exon_Node>::iterator exon_iter = exon_Vec.begin();
			exon_Vec.erase(exon_iter);
			exon_Num--;*/
		}
		else
		{
			int temp1 = exon_Vec[exon_Num-1].bound1;
			int temp2 = _max(start_site,exon_Vec[0].bound2-1200);
			if(temp2<exon_Vec[0].bound2)
				exon_Vec[0].bound1 = temp2;   // temporary start point of the gene
			else
				exon_Vec[0].bound1 = _max(0,exon_Vec[0].bound2-1200);
			exon_Vec[0].len = exon_Vec[0].bound2 - exon_Vec[0].bound1 + 1;
			if(exon_Num>0)
				exon_Vec[exon_Num-1].bound2 = _min(temp1+3000,exon_Vec[exon_Num-1].bound2)/*(int)(0.5*(temp+p2)*/;  // temporary end point of the gene
			exon_Vec[exon_Num-1].len = exon_Vec[exon_Num-1].bound2 - exon_Vec[exon_Num-1].bound1 + 1;
		}
	}	

	return true;
}

// check whether the genes overlap
int Gene_Batch::gene_Overlap(ifstream& infile, string orientation, int& pre_meStart, int& pre_meStop, int site1, int site2, 
	                          vector<vector<int>>& reverse_Pair, vector<int>& exon_Start, vector<int>& exon_Stop)
{
#pragma region paramenter declaration
	string word, line, seg[12];
	float jump_Thresh = 3000;    // forward threshold of the distance of jump
    float jump_Thresh_1 = 2500;  // inverse threshold of the distance of jump
 	int junc_start = 0, junc_stop = 0, exon_start = 0, exon_stop = 0, pre = 0, pos = 0, temp = 0;
	int s1 = pre_meStart;    // start point of the maximum gene
	int link_Flag = 2;
	int m_exonStart2 = site2, m_exonStop2 = site1;  // start point and end point of the exon inserted in the gene
	int& m_exonStart1 = pre_meStart;
	int& m_exonStop1 = pre_meStop;  // start point and end point of the former exon 
		
	char tPar = ',';
	if(stack_exonBoundary.size()==0)
	{
		vector<int> start, stop;
		stack_exonBoundary.push_back(stop);
		stack_exonBoundary.push_back(start);
	}
	vector<int>& exon_start1 = stack_exonBoundary[1];
	vector<int>& exon_stop1 = stack_exonBoundary[0];  // start point and end point of the exon 
	
	ending_site1 = _max(pre_meStart+200,site1);  // end point of transcipt
	ending_site2 = site2+200;   // end point of transcipt
	
	int num1 = 1, num2 = 2;  // number of hops in the inverse direction
	int ori_cnt = 0, rev_cnt = 0;
	bool pre_state = false;  // direction of the last gene

#pragma endregion
	
	int reverse_num = 0, reverse_end = 0;
	bool restart_Flag = true;
	if(intra_activeFlag==true||reverse_Pair[0].size()>0)   // overlap exists
	{
		// intra_activeFlag = true;
		reverse_num = reverse_Pair[1].size();
		reverse_end = reverse_Pair[1][reverse_num-1];  // end point
		if(site1-reverse_end<jump_Thresh)
		{
			intra_activeFlag = true;
			reverse_Pair[0].pop_back();  // delete last element
			reverse_Pair[0].push_back(site1);
			reverse_Pair[1].push_back(site2);
			reverse_num++;
			reverse_end = exon_start;
			restart_Flag = false;
		}
		else
		{
			intra_activeFlag = false;   // last gene finished
		}
	}

	if(restart_Flag==true)
	{
		exon_stop1.push_back(site1);
		int temp = _max(pre_meStart,site1-1000);
		exon_start1.push_back(temp);   // add first end
		exon_start1.push_back(site2);
	}

	while(std::getline(infile,line))
	{
		stringstream stream(line);
		// cout<<line<<endl;
		for(int i=0; i<12; i++)
			stream>>seg[i];   // load 12 fields

		if(seg[0]!=chromosome_Name)
		{
			// check the next chromosome
			pre_line = line;           // current line
			next_chromoName = seg[0];  // modify chromosome name
			record_file<<"Chromosome name changed: "<<chromosome_Name<<endl;
			ending_site1 = m_exonStart1 + 200;
			return 0;
		}

		sscanf(seg[1].c_str(),"%d",&junc_start);  // start point of junction
		sscanf(seg[2].c_str(),"%d",&junc_stop);   // end point of junction
		
#pragma region locate the boundary of the exon in the annotation
		pre = 0;
		if((pos=seg[10].find(tPar,pre))!=string::npos)  // find "," in the string
		{
			string str = seg[10].substr(pre,pos-pre);
			sscanf(str.c_str(),"%d",&temp);
			exon_stop = junc_start + temp;      // start point of junction
			// exon_Stop.push_back(junc_start+temp);
		}

		pre = 0;
		if((pos=seg[11].find(tPar,pre))!=string::npos)  // find "," in the string
		{
			int len = seg[11].length();
			string str = seg[11].substr(pos+1,len-pos-1);
			sscanf(str.c_str(),"%d",&temp);   // length of the exon
			exon_start = junc_start + temp + 1;
			// exon_Start.push_back(junc_start+temp+1);
		}
#pragma endregion locate the boundary of the exon in the annotation

		//////////////////////////////////////////////////////////////////////////
		// error control
		if(abs(exon_start-exon_stop)>100000)  // hop distance is too large, error may exist
		{
			sus_Junc.push_back(make_pair(exon_stop,exon_start));
			record_file<<"mark: "<<line<<endl;
			continue;
		}

		if(seg[5]!=orientation)  // inverse direction
		{
			if(pre_state==false)
				rev_cnt++;   // inverse direction continuously exists
			else
				rev_cnt = 1;
			pre_state = false;

			if(rev_cnt>=3)
			{
				jump_Thresh = 1500;
				jump_Thresh_1 = 800;
			}

			if(exon_stop>m_exonStop2)
				m_exonStop2 = exon_stop;  // update the right boundary of exon

			if(exon_stop-m_exonStart1>jump_Thresh_1)
			{
				if(intra_activeFlag==true)   // overlap exists
				{
					if(exon_stop-reverse_end<jump_Thresh)
					{
						reverse_Pair[0].push_back(exon_stop);
						reverse_Pair[1].push_back(exon_start);
						reverse_num++;
						reverse_end = exon_start;
					}
					else
					{
						ending_site2 = reverse_end+1200;
						if(reverse_num>0)
						{
							ending_site2 = _min(ending_site2,exon_stop);
							reverse_Pair[0].push_back(ending_site2);
						}
	
						ending_site1 = _min(m_exonStart1+1200,ending_site2-200);
						break;
					}
				}
				else
				{
					break;    // overlap doesn't exist
				}
			}
			else
			{
				exon_start1.push_back(exon_start);
				exon_stop1.push_back(exon_stop);
			}

			if(exon_stop>m_exonStart2)
				m_exonStart2 = exon_stop;  // update the right boundary of exon
		}
		else    // positive direction
		{
			pre_state = true;

			if(exon_stop>m_exonStop1)
				m_exonStop1 = exon_stop;  // update the right boundary of exon

			if(m_exonStop1-m_exonStart1>jump_Thresh)    // distance span is too long
			{

				ending_site2 = m_exonStop1;
				break;
			}
			else
			{
				pre_meStart = _max(exon_start,m_exonStart1);  // update the maximum start point of the exon
				pre_meStop = m_exonStop1;     // update the end point of maximum exon
				// exon with inverse direction
				int num1 = exon_stop1.size();
				int num2 = exon_start1.size();
				exon_Start.push_back(exon_start);
				exon_Stop.push_back(exon_stop);
				
				if(num1>0&&num2>0)
				{
					if(reverse_Pair[1].size()==0)
					{
						/*int temp = _max(s1,exon_stop1[0]-1000);
						reverse_Pair[1].push_back(temp);*/
						start_sites.push_back(exon_start1[0]);
					}
					for(int k=0; k<num1; k++)
					{
						reverse_Pair[0].push_back(exon_stop1[k]);
					}
					for(int k=0; k<num2; k++)
					{
						reverse_Pair[1].push_back(exon_start1[k]);
					}
					ending_site2 = _min(exon_start1[num2-1]+1000,exon_stop);
					reverse_Pair[0].push_back(ending_site2);
				}
				link_Flag = 1;  // return true
			    ending_site1 = m_exonStop1;
				stack_exonBoundary[0].clear();
				stack_exonBoundary[1].clear();
				break;            
			}

			if(exon_start>m_exonStart1)
				m_exonStart1 = exon_start;
		}
	} 

	pre_line = line;  // current line
	record_file<<"gene_Overlap: "<<pre_line<<endl;
	return link_Flag;
}

// check whether the genes overlap
// function£ºidentify the ends of exons until the direction is changed or the distance span is too long
// identify the ends of exons until the direction is changed or the distance span is too long
int Gene_Batch::intra_Overlap(ifstream& infile, string orientation, int& pre_meStart, int& pre_meStop, int site1, int site2, 
	                           vector<vector<int>>& reverse_Pair, vector<int>& exon_Start, vector<int>& exon_Stop)
{
	int link_Flag = 2;
	string word, line, seg[12];
	float jump_Thresh = 3000;    // threshold of the hop distance in the positive dirction
	float jump_Thresh_1 = 2500;  // threshold of the hop distance in the inverse dirction
	int junc_start = 0, junc_stop = 0, exon_start = 0, exon_stop = 0, pre = 0, pos = 0, temp = 0;
	int p1 = 0, p2 = 0;
	int s1 = pre_meStart;        // start point of the maximum gene in the positive direction
	char tPar = ',';	
	int m_exonStart2 = site2, m_exonStop2 = site1;  // start point and end point of the exons in the inserted gene
	int& m_exonStart1 = pre_meStart;
	int& m_exonStop1 = pre_meStop;   // start point and end point of the gene

	ending_site1 = pre_meStart+200;  // start point of transciption
	ending_site2 = site2+100;        // end point of transciption
	int ending_site2_bak = site2+5000;  // candidates of end point of transciption

	int num1 = reverse_Pair[0].size(), num2 = reverse_Pair[1].size();  // number of hop points in the inverse direction
	int ori_cnt = 1, rev_cnt = 1;
	bool pre_state = false;    // direction of the last gene
	double spr_Thresh = 200000;
	
	while(std::getline(infile,line))
	{
		stringstream stream(line);
		// cout<<line<<endl;
		for(int i=0; i<12; i++)
			stream>>seg[i];   // load 12 field

		if(seg[0]!=chromosome_Name)
		{
			pre_line = line;
			next_chromoName = seg[0];  // modify chromosome name
			record_file<<"Chromosome name changed: "<<chromosome_Name<<endl;
			ending_site1 = m_exonStart1+200;
			return 0;
		}

		sscanf(seg[1].c_str(),"%d",&junc_start);  // start point of junction
		sscanf(seg[2].c_str(),"%d",&junc_stop);   // end point of junction

#pragma region locate the boundary of the exon in the annotation
		pre = 0;
		if((pos=seg[10].find(tPar,pre))!=string::npos)  // find "," in the string
		{
			string str = seg[10].substr(pre,pos-pre);
			sscanf(str.c_str(),"%d",&temp);
			exon_stop = junc_start + temp;      //start point of junction
			// exon_Stop.push_back(junc_start+temp);
		}

		pre = 0;
		if((pos=seg[11].find(tPar,pre))!=string::npos)  // find "," in the string
		{
			int len = seg[11].length();
			string str = seg[11].substr(pos+1,len-pos-1);
			sscanf(str.c_str(),"%d",&temp);   // length of the exon
			exon_start = junc_start + temp + 1;
			// exon_Start.push_back(junc_start+temp+1);
		}
#pragma endregion locate the boundary of the exon in the annotation

		//////////////////////////////////////////////////////////////////////////
		// error control
		if(abs(exon_start-exon_stop)>spr_Thresh)  // hop distance is too large, error may exist
		{
			sus_Junc.push_back(make_pair(exon_stop,exon_start));
			record_file<<"mark: "<<line<<endl;
			continue;
		}

		if(seg[5]!=orientation)   // inverse direction
		{
			if(pre_state==false)
				rev_cnt++;   //  inverse direction continuously exists
			else
				rev_cnt = 1;
			pre_state = false;

			if(exon_stop>m_exonStop2)
				m_exonStop2 = exon_stop;  // update the right boundary of exon

			if(ori_cnt>=3)
			{
				jump_Thresh_1 = 800; // adjust threshold according to the probability of split
				jump_Thresh = 1500; 
			}

			if(m_exonStop2-m_exonStart2>jump_Thresh)  // distance of hop is too long in the inverse direction
			{
				pre_line = line; // current line
				if(num1>0&&num1<num2)
				{
					int temp = _min(ending_site2_bak,reverse_Pair[1][num2-1]+1200);
					temp = _max(ending_site2,temp);    // avoid the end of the exons from exceeding the limit
					reverse_Pair[0].push_back(temp);   // add the end of the exons 
					ending_site2 = temp;
					//ending_site2 = _min(m_exonStop2,reverse_Pair[1][num2-1]+1000);
				}
				record_file<<"intra_Overlap: "<<pre_line<<endl;
				intra_activeFlag = false;  // set the flag to be inactive
				rev_breakPoints.push_back(num2);
				return 1;   // retur nthe flag
			}
			else
			{

				reverse_Pair[0].push_back(exon_stop);
				reverse_Pair[1].push_back(exon_start);
				num1++;
				num2++;
			}
			if(exon_start>m_exonStart2)
				m_exonStart2 = exon_start;
		}
		else    // positive direction
		{
			if(pre_state==true)
				ori_cnt++;   // positive direction continuously exists
			else
				ori_cnt = 1;
			pre_state = true;

			if(exon_stop>m_exonStop1)
				m_exonStop1 = exon_stop;  // update the right boundary of exon

			if(exon_stop>m_exonStart2&&ori_cnt==1)  // direction is changed
				ending_site2_bak = exon_stop;
			
			if(rev_cnt>=3)
			{
				jump_Thresh_1 = 800; //  adjust threshold according to the probability of split
				jump_Thresh = 1500; 
			}

			if(m_exonStop1-m_exonStart2>jump_Thresh_1)    // distance span is too long
			{

				p2 = junc_start;     // end point of the exon
				pre_line = line;
				if(m_exonStop1-m_exonStart1>jump_Thresh)
				{
					link_Flag = 2;
					ending_site1 = m_exonStop1;  // temporary end point of transciption
				}
				else
					link_Flag = 1;

				intra_activeFlag = false;  // set the flag to be inactive
				rev_breakPoints.push_back(reverse_Pair[0].size());
				break;
			}
			else
			{
				exon_Stop.push_back(exon_stop);
				exon_Start.push_back(exon_start);
			}
			if(exon_start>m_exonStart1)
				m_exonStart1 = exon_start;
		}
	}

	record_file<<"intra_Overlap: "<<pre_line<<endl;

	ending_site1 = _min(m_exonStop1,pre_meStart+500);
	ending_site2 = m_exonStop2;
	return link_Flag;
}

// locate the possible start site and end site of transcription
void Gene_Batch::locate_transcriptionSites(Gene& gene, vector<PAIR_INT>& sites)
{
	int trs_Num = gene.trans.size();  // number of isoform
	sites.clear();
	for(int i=0; i<trs_Num; i++)
	{
		int exon_Num = gene.trans[i].exon_start.size();  // number of exons
		int s1 = gene.trans[i].exon_start[0];   // start point of transcription
		int s2 = gene.trans[i].exon_stop[exon_Num-1];  // end point of transcription
		sites.push_back(make_pair(s1,s2));
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
// modify the interval index of read according to the data
bool Gene_Batch::serial_Modify(Rd_Set& rd_Set, Graph_Trans& gph)
{	
	int exon_Num = gph.exon_Vec.size(); // number of exons
	// int intron_Num = gph.intron_Add.size();  // number of introns
	int intron_Num = exon_Num - 1;  // number of introns
	bool inser_Flag = false;  // check whether the read is in the intron
	int rd_Num = rd_Set.rd_Vec.size();  // number of read

	vector<int> inser_sum(intron_Num);  // number of segments in the intron
	vector<int> inser_cnt(intron_Num);  // number of segments in the intron

	inser_sum[0] = gph.intron_Add[0].size();
	for(int k=1; k<intron_Num; k++)
	{
		inser_cnt[k] = gph.intron_Add[k].size();
		inser_sum[k] = inser_sum[k-1] + inser_cnt[k-1];
	}

	int inser_serial = 0, serial = 0;
	int pos = 0;
	for(int i=0; i<rd_Num; i++)
	{
		Rd& cur_Rd = rd_Set.rd_Vec[i];
		int n1 = cur_Rd.left_readIdx.size();
		int n2 = cur_Rd.right_readIdx.size();
		vector<int> cnt;
		vector<int>& readIdx = cur_Rd.left_readIdx;
		cnt.push_back(n1); cnt.push_back(n2);

		for(int l=0; l<2; l++)
		{
			if(l==1)
				readIdx = cur_Rd.right_readIdx;
			for(int k=0; k<cnt[l]; k++)
			{
#pragma region modify the index
				if(readIdx[k]<0)
				{
					inser_Flag = true;
					inser_serial = -readIdx[k]/2-1;
					pos = cur_Rd.bi_pos[l][k];  // locatin of the split point
					if(inser_cnt[inser_serial]>0)
					{
						vector<Exon_Node>& node_vec = intron_Add[inser_serial];
						int num = inser_cnt[inser_serial];
						for(int j=0; j<num; j++)
						{
							if(pos<node_vec[j].bound2&&pos>node_vec[j].bound1)
								cur_Rd.left_readIdx[k] = inser_serial+inser_sum[j]-(num-j-1);
						}
					}
				}
				else
				{
					serial = cur_Rd.left_readIdx[k];
					if(serial>0)
						cur_Rd.left_readIdx[k] = serial+inser_sum[serial-1];
				}
#pragma endregion modify the index
			}
		}
	}

	return true;
}

// check whether the gene is valid
bool Gene_Batch::check_Valid(Rd_Set& rd_Set, Graph_Trans& gph)
{
	bool ctr_Flag = true;
	int rd_Num = rd_Set.rd_Vec.size();  // number of read
	float thresh_number = 10;
	float thresh_rate = 0.2;
	int exon_cover = 0;
	int seg_Num = rd_Left.size();
	if(rd_Num<3)
		return false;

	return ctr_Flag;
}

// check the retention of the intron
bool Gene_Batch::intron_Retention(Rd_Set& rd_Set, Graph_Trans& gph, vector<vector<int>>& intron_Junc1,
	                              vector<vector<int>>& intron_Junc2, vector<vector<int>>& rd_Left, vector<vector<int>>& rd_Right)
{
	int seg_Num = rd_Left.size();  // number of intervals
	int exon_Num = (seg_Num+1)/2;    // number of exons
	int intron_Num = (seg_Num-1)/2;  // number of introns
	vector<int> intron_Cnt(intron_Num);
	for(int i=0; i<intron_Num; i++)
		intron_Cnt[i] = 0;

	int threshold = 10;
	int add_serial = 0;
	for(int i=1; i<seg_Num; i=i+2)
	{
		add_serial=0;
		int in_serial = (i-1)/2;
		vector<Exon_Node>& intron_Vec = intron_Add[in_serial];
		vector<int>& serial = serial_Add[in_serial];

		vector<int> seg1 = rd_Left[i];
		vector<int> seg2 = rd_Right[i];
		bool ctr_Flag = false;
		bool inser_Flag = false;
		int num1 = seg1.size(), num2 = seg2.size();  // number of read in the interval
		int s1 = 0, s2 = -1, z1 = 0, z2 = -1;
		int site1 = gph.exon_Vec[in_serial].bound2;
		int site2 = gph.exon_Vec[in_serial+1].bound1;
		int cnt1 = 0, cnt2 = 0;
#pragma region calculate the cover rate of the read in the intron
		if(num1+num2>threshold)   
		{
			for(int k=0; k<seg1.size(); k++)
			{
				if(seg1[k]>site1+5&&seg1[k]<site2-5)
					cnt1++;
			}
			for(int k=0; k<seg2.size(); k++)
			{
				if(seg2[k]>site1+5&&seg2[k]<site2-5)
					cnt2++;
			}
			if(cnt1+cnt2>threshold)
			{
				inser_Flag = true;
				if(num2>0)
				{
					sort(seg2.begin(),seg2.end());
					z1 = seg2[0]; z2 = seg2[num2-1];
				}
				s1 = num1>0?_min(seg1[0],z1):z1;
				s2 = num1>0?_max(seg1[num1-1],z2):z2;
			}
		}
#pragma endregion calculate the cover rate of the read in the intron
		int seg_Len = gph.exon_Vec[i].bound1-gph.exon_Vec[i-1].bound2-1;	

		vector<int>& in_Junc1 = intron_Junc1[in_serial];
		vector<int>& in_Junc2 = intron_Junc2[in_serial];
		int n1 = in_Junc1.size();  // number of juctions in the intron
		int n2 = in_Junc2.size();  // number of juctions in the intron
		int p5 = n1>0?in_Junc1[n1-1]:(gph.exon_Vec[in_serial].bound2+1);
		int p3 = n2>0?in_Junc2[0]:(gph.exon_Vec[in_serial+1].bound1-1);

		cnt1 = 0;
		cnt2 = 0;
#pragma region calculate the cover rate of the middle read
		for(int k=0; k<num1; k++)
		{
			if(seg1[k]>=p5&&seg1[k]<=p3)
				cnt1++;
		}
		for(int k=0; k<num2; k++)
		{
			if(seg2[k]>=p5&&seg2[k]<=p3)
				cnt2++;
		}
		if(cnt1+cnt2>threshold)
			ctr_Flag = true;
#pragma endregion calculate the cover rate of the middle read

#pragma region add exons
		int inter1 = 0, inter2 = 0;
		if(n1<2&&n2<2&&inser_Flag==true)
		{
#pragma region intron retention
			Exon_Node node;
			float thresh1 = seg_Len*0.2; 
			inter1 = s1-site1;
			inter2 = site2-s2;
			if(inter1>thresh1&&inter1>inter2)
			{
				node.bound1 = s1;
				node.bound2 = site2-1;
				cout<<"Intron Retention 1"<<endl;
			}
			else if(inter2>thresh1&&inter2>inter1)
			{
				node.bound1 = site1+1;
				node.bound2 = s2;
				cout<<"Intron Retention 2"<<endl;
			}
			else if(inter1<thresh1&&inter2<thresh1)
			{
				node.bound1 = site1+1;
				node.bound2 = site2-1;
				cout<<"Intron Retention 0"<<endl;
			}
			node.len = node.bound2-node.bound1+1;
			intron_Vec.push_back(node);
			serial.push_back(in_serial);
			add_serial++;
#pragma endregion
		}
		else
		{
#pragma region alternative 5',3'
			bool add_Right = false;
			if(n1>=2)  // possible 5' end in the intron
			{
				if(p5<p3)
				{
					Exon_Node node;
					node.bound1 = site1+1;  // left boundary of the new exon
					node.bound2 = p5;  // right boundary of the new exon
					node.len = node.bound2-node.bound1+1;
					if(node.len>10){
						intron_Vec.push_back(node);
						serial.push_back(in_serial);
						add_serial++;
					}
					if(ctr_Flag==true)
					{
						Exon_Node node;
						node.bound1 = p5+1;
						node.bound2 = p3-1;
						node.len = node.bound2-node.bound1+1; // length of the exon
						intron_Vec.push_back(node);
						serial.push_back(in_serial);
						add_serial++;
					}
					if(n2>=2)
						add_Right = true;
				}
				else if(n2>=2) // p5>p3
				{
					Exon_Node node;
					node.bound1 = p3;
					node.bound2 = p5;
					node.len = node.bound2-node.bound1+1;
					if(node.len>20){
						intron_Vec.push_back(node);
						serial.push_back(in_serial);
						add_serial++;
					}
					if(node.len>20){
						intron_Vec.push_back(node);
						serial.push_back(in_serial);
						add_serial++;
					}
				}
			}
			else if(n2>=2||add_Right==true)
			{
				Exon_Node node;
				node.bound1 = p3;
				node.bound2 = site2-1;
				node.len = node.bound2-node.bound1+1;
				intron_Vec.push_back(node);
				serial.push_back(in_serial);
				add_serial++;
			}			
#pragma endregion alternative 5',3'
		}
		intron_Cnt[in_serial] = add_serial;  // nodes inserted into the intron
#pragma endregion 
	}
	return true;
}

// check wheter the exon should be deleted
bool Gene_Batch::check_Intron(int x1, int x2, int& z1, int& z2, vector<int>& seg)
{
	int idx = 0, x = 0;
	int num = seg.size();  // number of reads in the interval
	while(idx<num&&(x=seg[idx])<x1+2)
		idx++;
	int k1 = idx;
	while(idx<num&&(x=seg[idx])<x2-2)
		idx++;
	int count = idx-k1;

	if(count<5)   // transcript rate is low
		return true; //intron
	else
	{
		z2 = seg[idx-1];  // start point of the right end read in the interval
		z1 = seg[k1];     // end point of the right end read in the interval
		return false;
	}
}

// Comparison with results based on gene annotations: modify the results based on gene annotations
void Gene_Batch::exon_CompareAnno(vector<Exon_Node>& exon_Vec, vector<Exon_Node>& exon_Vec_Anno, Gene& cur_Gene, vector<vector<int>>& junc_Num)
{
	int num1 = exon_Vec.size();       // number of exons
	int num2 = exon_Vec_Anno.size(); 
	int num3 = junc_Num.size();       // number of junctions
	float thresh_1 = 20, thresh_2 = 0;
	float sum = 0, aver = 0;
	for(int i=0; i<num3; i++)
		sum += junc_Num[i][2];

	aver = sum*1.0/num3;
	thresh_2 = _max(thresh_1,aver*0.5);
	int exon_start = 0, exon_stop = 0;
	int proceed1 = 0, proceed2 = 0, ser1 = 0, ser2 = 0, ser3 = 0, k1 = 0, k2 = 0;

	vector<int> mark1, mark2;
	for(int k=0; k<num1; k++)
	{
		mark1.push_back(0);
		mark2.push_back(0);
	}

#pragma region Check the validity of the left boundary and right boundary of the exon
	for(int k=0; k<num1-1; k++)
	{
		exon_stop = exon_Vec[k].bound2;  // right boundary of the exon
		bool ctr_Flag = false;
		for(k2 = ser2; k2<num2; k2++)
		{
			int temp = exon_Vec_Anno[k2].bound2;  // right boundary of the exon
			if(abs(exon_stop-temp)<5)     // the right boundaries coincide approximately
			{
				ctr_Flag = true;
				mark2[k] = 1;
				break;
			}
		}
		if(ctr_Flag==true)
			ser2 = k2;
		else   // check the confidence of the exon
		{
			int valid = 0;
			for(int j=0; j<num3; j++)
			{
				if(abs(exon_stop-junc_Num[j][0])<3)    // the junctions coincide
				{
					valid += junc_Num[j][2];   // count the valid split reads
				}
			}
			if(valid>thresh_2)  // the junction is valid
			{
				mark2[k] = 2;
			}
		}
	}

	for(int i=1; i<num1; i++)
	{
		exon_start = exon_Vec[i].bound1;  // left boundary of the exon
		bool ctr_Flag = false;
		for(int k1 = ser1; k1<num2; k1++)
		{
			int temp = exon_Vec_Anno[k1].bound1;  // left boundary of the exon
			if(abs(exon_start-temp)<10)     // the left boundaries coincide approximately
			{
				ctr_Flag = true;
				mark1[i] = 1;
				break;
			}
		}
		if(ctr_Flag==true)
			ser1 = k1;
		else   // check the confidence of the exon
		{
			int valid = 0;
			for(int j=0; j<num3; j++)
			{
				if(abs(exon_start-junc_Num[j][1])<3)    // the junctions coincide
				{
					valid += junc_Num[j][2];   // count the valid split reads
				}
			}
			if(valid>thresh_2)  // the junction is valid
			{
				mark1[i] = 2;
			}
		}
	}
#pragma endregion

#pragma region  choose the exon to be added

	int exon_cnt = exon_Vec_Anno.size(); // number of exons
	vector<int> pos, serial;
	for(int i=0; i<num1; i++)
	{
		if(mark1[i]==2||mark2[i]==2)  // either the left boundary or the right boundary meets the requirement
		{
			while(ser3<exon_cnt&&exon_Vec_Anno[ser3].bound2>exon_Vec[i].bound2)
				ser3++;
			if(ser3<exon_cnt)
			{
				if(exon_Vec_Anno[ser3].bound2>exon_Vec[i].bound1)
					exon_Vec[i].bound1 = exon_Vec_Anno[ser3].bound2+1;  // modify the boundary of the exon
				pos.push_back(i);
				serial.push_back(ser3);
			}
		}
	}

	int num = pos.size();  // number of added exons
	for(int k=0; k<num; k++)
	{
		int ser = serial[k];  // the serial of the exon
		int insert_ser = pos[k];
		exon_Vec_Anno.insert(exon_Vec_Anno.begin()+ser+k,exon_Vec[insert_ser]);
	}

	cur_Gene.subexon_Num = cur_Gene.subexon.size();  // update the number of exons
#pragma endregion

}
