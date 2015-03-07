
//////////////////////////////////////////////////////////////////////////
//// TRANSCRIPT ASSEMBLY AND ABUNDANCE ESTIMATION BASED ON RNA-SEQ DATA

#include "Entropy.h"
#include "Batch_Mode.h"
#include "ISA.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Deal with the sequencing data containing intron-retention
//// Examine the length of read
int Gene_Batch::rdLen_estimation(char* filename)
{
	ifstream infile;
	infile.open(filename,ios::in);
	if(!infile)
	{
		cout<<"Open file error!"<<endl;
		exit(-1);
	}

	string line, word, seg[12];
	int read_Len = 75;
	string::size_type pos0 = 0;
	string substring;
	while(getline(infile,line))
	{
		stringstream stream(line);
		for(int k=0; k<6; k++)
			stream>>seg[k];
		pos0 = seg[5].find('M',0);   // locate the first splitting character
		if(pos0<seg[5].length()-1)
			continue;
		else
		{
			substring = seg[5].substr(0,pos0);
			sscanf(substring.c_str(),"%d",&read_Len);
			break;
		}
	}

	rd_Len = read_Len;
	return read_Len;
}

// Estimate the mean and variance of fragment length
void Gene_Batch::fgLen_estimation(char* filename1, char* filename2)
{
	ifstream infile1, infile2;
	infile1.open(filename1,ios::in);
	if(!infile1){  
		cout<<"Open file error£¡(Junction)"<<endl;
		return;
	}
	infile2.open(filename2, ios::in);
	if(!infile2){
		cout<<"Open file error!(Alignment)"<<endl;
		return;
	}

	char filename[200];
	int k = 0, s1 = 0, s2 = 0, start = 0, stop = 0;   // starting site and stopping site of gene
	int gene_Cnt = 0;
	char filename3[200];

	string line, line1, word, seg[10];
	std::streampos pos = infile1.tellg();   // obtain the current file pointer
	string chromo_Name = "chr1";
	chromosome_Num = 0;

	//////////////////////////////////////////////////////
#pragma region batch mode 
	getline(infile1,line);   // read the first line: line of caption
	getline(infile1,line);   // read the second line
	cout<<line<<endl;
	pre_line = line;
	stringstream stream(line);
	stream>>chromo_Name;     // name of the first chromosome
	chromosome_Name = chromo_Name;
	next_chromoName = "";
	cout<<"chromosome name: "<<chromosome_Name<<endl;  // chromosome name
#pragma endregion
	//////////////////////////////////////////////////////

	vector<int> exon_Start, exon_Stop;
	starting_site = 0;  // initialize the transcription start site

	// initialize the mean and variance of the fragment length
	fragLen_Mean = DEFAULT_MEAN;
	fragLen_Var = DEFAULT_VAR;
	int valid_Cnt = 0, read_Cnt = 0;
	bool changed_Flag = false;
	vector<int> read_Length;
	while(!infile1.eof())
	{
		clock_t start_Time = clock();
		exon_Start.clear(); exon_Stop.clear(); exon_Assembly.clear();
		changed_Flag = exon_AssembleBatch(infile1,infile2,exon_Start,exon_Stop);

		pos = infile1.tellg();
		std::streampos pos1 = infile2.tellg();   // obtain the current file pointer
		int est_cnt = rev_ExonAssembly.size() + 1;

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

			coverage_Assemble(infile2, rd_Set);  // assemble the starting site and stopping site of the transcripts
			if(rd_Set.rd_Vec.size()==0)
				continue;

			assemble_SubGraph(seg_Modify);
			int loc_geneNum = gene_Modify(gene_Set, seg_Modify); // modify genes in the current position

			if(!infile2)
				cout<<"Error£¡"<<endl;

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

					stringstream ss;
					ss<<gene_Cnt;
					rd_Set.gene_Name = ss.str();
					int locate_Cnt =  assignment_Pose(rd_Set, gph, infile2, chromo_Name);  // load read alignment

					// estimate the mean and variance of the fragment length
					rd_Set.locate_RdNum = locate_Cnt;
					// cout<<"Located reads "<<locate_Cnt<<endl;
					rd_Set.mean = fragLen_Mean;
					rd_Set.sigma = fragLen_Var;
					if(locate_Cnt>100||(gene_Cnt==1&&locate_Cnt>0))
					{
						valid_Cnt += rd_Set.fragLen_Estimation_Initial(read_Length);
					}					
				}
			}
		}

		rd_Set.~Rd_Set();
		vector<VEC_NODE>().swap(rev_ExonAssembly);
		vector<vector<int>>().swap(rev_ExonType);
		vector<vector<float>>().swap(rev_ExonScore);

		if(valid_Cnt>200) 
			break;
	}

	double sum1 = 0, mean1 = 0, var1;
	for(int i=0; i<valid_Cnt; i++)
		sum1 += read_Length[i];
	mean1 = sum1*1.0/valid_Cnt;  // initial estimation of the mean
	sum1 = 0;
	for(int i=0; i<valid_Cnt; i++)
		sum1 += (read_Length[i]-mean1)*(read_Length[i]-mean1);
	var1 = sqrt(sum1*1.0/(valid_Cnt-1));

	if(var1<20)
		var1 = 20;
	else if(var1>50)
		var1 = 50;

	fragLen_Mean = mean1;
	fragLen_Var = var1;
}

bool Gene_Batch::assignment_Pose_1(Gene& cur_Gene, Rd_Set& rd_Set, Graph_Trans& gph, ifstream& myfile, string chromo_Name)
{
	if(!myfile){
		cout << "Unable to open myfile";
		exit(1); // terminate with error
	}

	int x1 = 0, x2 = 0, z1 = 0, z2 = 0, temp_Len = 0; // ends of left reads and right reads, and fragment length
	int read_len = rd_Len;  // read length

	stringstream ss;
	ss<<rd_Len;
	string len_R = ss.str()+"M";  // string indicating exact match of the whole read 
	string line, word, words[10], frag_Name, chromo_name;

#pragma region Declaration of Variables
	int serial = 0;       // read serial
	string tPar = "N";    // marker for mismatch
	string tPar1 = "I";   // marker for insertion 
	string tPar2 = "D";   // marker for deletion
	string::size_type pos, pos0;
	vector<string> unresolved_read;   // reads not located 
	vector<int> unresolved_serial;    // serials of reads not located
	vector<int> right_readIdx;
	vector<int>::iterator itr;
	vector<string>::iterator str_itr;  // the indices of the segments in which the four ends of the two reads are located
	vector<Rd>& rd_Vec = rd_Set.rd_Vec;
	char j1 = -1, j2 = -1, j3 = -1, j4 = -1;  // the indices of the segments to which the four ends of the two reads belong

	int junc_Len = 0, cover_Len = 0, p1 = 0, p2 = 0;
	int x_1 = 0, x_2 = 0, z_1 = 0, z_2 = 0, je = 0, rlen1 = 0, rlen2 = 0, j_3 = 0, j_4 = 0, j_1 = 0, j_2 = 0;
	bool anno_flag = false;
	int line_cnt = 0;
	bool flag_Ctrl = true;
	bool locate_Flag = false;
	int exon_Num = gph.exon_Vec.size();       // number of exons
	int e_Num = exon_Num;
	int start_Site = gph.exon_Vec[0].bound1;  // transcription start site of the gene
	int stop_Site = gph.exon_Vec[exon_Num-1].bound2;  // transcription stop site of the gene
	std::streampos pos1 = myfile.tellg();     // obtain the present file pointer
	int rd_Cnt = 0;      // counting the reads falling into the genetic segments

#pragma endregion Declaration of Variables
	while(flag_Ctrl==true&&std::getline(myfile,line))  // load data in the valid genetic range
	{
		line_cnt++;
		Rd cur_Rd;
		j1 = -1, j2 = -1, j3 = -1, j4 = -1;

		anno_flag = false;
		stringstream stream(line);
		stream>>word; frag_Name = word;
		stream>>word; stream>>word; // read the first three fields
		chromo_name = word;
		for(int i = 0; i<10; i++)
			stream>>words[i];  // meaning of the fields: 4 starting site of left read 5 matching quality 6 matching quality 7 marker for paired match 8 starting site of right read 9 inserted size
		sscanf(words[5].c_str(),"%d",&temp_Len);  // inserted size
		sscanf(words[0].c_str(),"%d",&x1);  // starting site of left read
		sscanf(words[4].c_str(),"%d",&z1);  // starting site of right read
		z2 = 0;
		if(temp_Len>0)
			z2 = x1 + temp_Len - 1;  // stopping site of right read
		else   // single-end read
			z1 = x1;

		int s1 = _min(x1,z1), s2 = _max(_max(x1,z1),z2);
		int in_ser1 = 0, in_ser2 = 0, in_ser3 = 0, in_ser4 = 0;

#pragma region both of the paired reads can be matched
		if(s1>=start_Site&&s2<stop_Site&&chromo_name==chromo_Name)   // the starting site and stopping site of the fragment are in the range of the gene
		{
			rd_Cnt++;
			if(rd_Cnt>10)
				locate_Flag = true;
			if(temp_Len>=0)        // alignment result is shown in left-right order
			{
				x2 = x1 + rd_Len - 1;    // stopping site of left read stopping site of left read st
				if(x1>=start_Site&&z2<=stop_Site)
				{
					flag_Ctrl=true;
#pragma region identify the exons the ends of both reads fall in
					j1 = 0;
					while(j1<e_Num&&x1>=gph.exon_Vec[j1].bound1) j1++; j1--; // exon in which the left end of left read falls

					if(j1<0)
						continue;  // the read isn't matched into the identified exon regions

#pragma region identify junctions according to read alignment
					if(words[2]==len_R)   // exact match of the whole read
					{
						j2 = j1;   // both ends of the left read fall into the same exon
						x2 = x1 + rd_Len - 1;   // stopping site of the left read
						while(j2<e_Num&&x2>=gph.exon_Vec[j2].bound1) j2++;  j2--; // exon in which the right end of the left read falls

						if(j1!=j2)
						{
							for(int k1 = j1; k1<j2; k1++)
								gph.link_Mtx[k1+1][k1+2]++;  // each two adjacent exons between j1 and j2 are linked
							if(j2>j1+1)
							{
								for(int k1 = j1; k1<=j2; k1++)
									cur_Rd.left_readIdx.push_back(k1);  // exons in which the left read falls
							}				
						}
					}
					else
					{
						int junc_Len = 0, cover_Len = 0;
						string::size_type pre = 0;    // record the last position of "N"

#pragma region identify the segments in which the left read falls according to the identified junctions
						bool ctr_Flag = false;
						j2 = j1;
#pragma region left end of the left read and the middle matched part of the read
						while((pos=words[2].find(tPar,pre))!=string::npos||
							(pos=words[2].find(tPar1,pre))!=string::npos||
							(pos=words[2].find(tPar2,pre))!=string::npos)  // search character 'N' and 'I' in the alignment description
						{
							pos0 = words[2].find("M",pre); // locate the first splitting character
							if(pos0!=string::npos)
							{ 
								string str = words[2].substr(pos0+1,pos-pos0-1);
								sscanf(str.c_str(),"%d",&junc_Len);   // the distance spanned by the read
								if(words[2].substr(pos,1)==tPar1)     // there is insertion
									junc_Len = 0;

								str = words[2].substr(pre,pos0-pre);
								sscanf(str.c_str(),"%d",&cover_Len);  // the length covered by the read
								j_1 = j1;

								if(!ctr_Flag)
								{
#pragma region identification of the left end 		
									x_1 = x1 + cover_Len - 1;  // the location of the internal end of the left read
									j_1 = j1;
									while(j_1<e_Num&&x_1>=gph.exon_Vec[j_1].bound1) j_1++; j_1--;   // exon in which the internal end of the left read falls

									for(int k1 = j1; k1<=j_1; k1++)
										cur_Rd.left_readIdx.push_back(k1);

#pragma endregion identification of the left end
									p1 = x1 + cover_Len + junc_Len;  // the left end of the exon which is probably covered by the read
									ctr_Flag = true;
								}
								else
								{
#pragma region identfication of the middle exon covered by the read
									p2 = p1 + cover_Len - 1;
									while(j_1<e_Num&&abs(gph.exon_Vec[j_1].bound1-p1)>2)
										j_1++;  // locate the middle exon covered by the read
									while(abs(cover_Len)>2)   // adjacent exons are probably contained
									{
										cur_Rd.left_readIdx.push_back(j_1);
										cover_Len -= gph.exon_Vec[j_1].len;
										j_1++;
										if(j_1>=e_Num-1||abs(gph.exon_Vec[j_1].bound2-gph.exon_Vec[j_1+1].bound1)>2)
											break;
									}
#pragma endregion identfication of the middle exon covered by the read
									p1 = p2 + junc_Len + 1; 
								}
								pre = pos + 1;   // search after this position next time
							}
						}
#pragma endregion

						// right end(internal end) of the left read
						pos0 = words[2].find("M",pre);    // locate the first splitting character
						string str;
						if(pos0!=string::npos)
						{
							str = words[2].substr(pre,pos0-pre);
							sscanf(str.c_str(),"%d",&cover_Len); // the distance spanned by the read
							x_2 = p1 + cover_Len - 1;   // the position of the right end of the second part of the left read
							j_2 = j_1;  // exon in which the left end of the left read falls
							while(j_2<exon_Num&&p1>=gph.exon_Vec[j_2].bound2) j_2++; 
							j2 = j_2;   // exon in which the right end of the left read falls
							while(j2<exon_Num&&x_2>=gph.exon_Vec[j2].bound1-1) j2++; j2--;
							if(j_2!=j_1)
							{
								for(int k1 = j_2; k1<j2; k1++)
									cur_Rd.left_readIdx.push_back(k1);
							}
							cur_Rd.left_readIdx.push_back(j2);
						}
						x2 = x_2;
#pragma endregion identify the segments in which the left read falls according to the identified junctions 
						// search the exon the right end of the left read falls
						int idnum = cur_Rd.left_readIdx.size(), index1 = 0, index2 = 0;
						for(int k1 = 0; k1<idnum-1; k1++)
						{
							index1 = cur_Rd.left_readIdx[k1];
							index2 = cur_Rd.left_readIdx[k1+1];
							gph.link_Mtx[index1+1][index2+1]++;  // each two adjacent exons between j1 and j2 are linked
						}		
					}
#pragma endregion identify junctions according to read alignment

					if(words[3]=="=")   // paired alignment, hence another alignment result can be referred to
					{
						j3 = j1;
						while(j3<e_Num&&z1>=gph.exon_Vec[j3].bound1) j3++; j3--;  // exon in which the left end of the right read falls
						j4 = j3;        // both ends of the left read fall into the same exon
						while(j4<e_Num&&z2>=gph.exon_Vec[j4].bound1) j4++; j4--;  // exon in which the right end of the right read falls
						if(j3!=j4)
						{
							int part_Len = gph.exon_Vec[j3].bound2 - z1 + z2 - gph.exon_Vec[j4].bound1 + 2;
							if(part_Len==read_len)  // junction is only related to two exons
							{
								gph.link_Mtx[j3+1][j4+1]++;  // j3 is linked with j4
								anno_flag = true;
							}
							else if(part_Len<read_len)  // the exact alignment of the right read needs further examination
							{
								unresolved_read.push_back(frag_Name);  // record the name of the fragment
								unresolved_serial.push_back(serial);   // record the serial of the fragment
								anno_flag = false;
							}
						}
					}

					if(j3==j2+1)   // the two exons are adjacent
					{
						if(j1==j2)  // junction there is no junction within the left read
						{
							gph.link_Mtx[j2+1][j3+1] += rd_Len*1.0/gph.exon_Vec[j2].len;  // j1 and j2 are linked
						}
						else
						{
							gph.link_Mtx[j2+1][j3+1]++;
						}
					}

					cur_Rd.name = frag_Name;
					cur_Rd.l1 = j1; cur_Rd.l2 = j2; cur_Rd.r1 = j3; cur_Rd.r2 = j4;
					cur_Rd.Left1 = x1; cur_Rd.Left2 = x2;
					cur_Rd.Right1 = z1; cur_Rd.Right2 = z2;
					cur_Rd.insert_Len = temp_Len;  // fragment length

					if(z2>x1)
						rd_Set.rd_Vec.push_back(cur_Rd);  // add read to the vector
					else
						rd_Set.singlerd_Vec.push_back(cur_Rd);  // add single-end read to the vector

					serial++;
#pragma endregion 
				}
			}
			else if(temp_Len<0)   // alignment is shown in the reverse (right-left) way
			{
				if(z1>=start_Site&&x1<stop_Site-rd_Len)
				{
#pragma region alignment is shown in the reverse (right-left) way
					if(words[3]=="=")   // paired alignment
					{  
						int num = unresolved_serial.size();  // number of right reads not located
						for(int k1=0; k1<num; k1++)
						{
							int basic_line = frag_Name.length()-2;
							if(frag_Name.substr(0,basic_line)==unresolved_read[k1].substr(0,basic_line))   // the fragment names are the same
							{
								int serial = unresolved_serial[k1];
								x_1 = rd_Vec[serial].Left1;  x2 = rd_Vec[serial].Left2;
								z_1 = rd_Vec[serial].Right1; z2 = rd_Vec[serial].Right2;
								if(x1==z_1&&z1==x_1&&rd_Vec[serial].right_readIdx.size()==0)  // the same pair of reads
								{
									x1 = x_1; z1 = z_1;
									j1 = rd_Vec[serial].l1;  j2 = rd_Vec[serial].l2;  j3 = rd_Vec[serial].r1;  j4 = rd_Vec[serial].r2;
#pragma region identify junction according to the results of read alignment
									if(words[2]==len_R)   // exact match of the whole read
									{
										j4 = j3;   // both ends of the right read fall in the same exon
										while(j4<e_Num&&z2>=gph.exon_Vec[j4].bound1) j4++;  j4--; // the exon in which the right end of the right read falls
										if(j3!=j4)
										{
											for(int k1 = j3; k1<j4; k1++)
												gph.link_Mtx[k1+1][k1+2]++;  // each two adjacent exons between j3 and j4 are linked
											if(j4-j3>1)
											{
												for(int k1=j3; k1<=j4; k1++)
													rd_Vec[serial].right_readIdx.push_back(k1);
											}
										}
									} 
									else
									{ 
										junc_Len = 0, cover_Len = 0;
										rlen2 = z2 - gph.exon_Vec[j4].bound1 + 1;
										rlen1 = gph.exon_Vec[j3].bound2 - z1 + 1;
										string::size_type pre;    // record the last position of character 'N'
										pre = 0;
										right_readIdx.clear();
#pragma region identify the segments in which the right read falls according to the identified junctions
										int p1 = 0, p2 = 0; 
										bool ctr_Flag = false;
										j_3 = j3;
#pragma region left end of the right read and the middle matched part of the read
										while((pos=words[2].find(tPar,pre))!=string::npos||
											(pos=words[2].find(tPar1,pre))!=string::npos||
											(pos=words[2].find(tPar2,pre))!=string::npos)  // search the character 'N'
										{
											pos0 = words[2].find("M",pre); // loate the first splitting character 'M'
											if(pos0!=string::npos)
											{ 
												string str = words[2].substr(pos0+1,pos-pos0-1);
												sscanf(str.c_str(),"%d",&junc_Len);   // distance spanned by the read
												if(words[2].substr(pos,1)==tPar1)     // there is insertion
													junc_Len = 0;

												str = words[2].substr(pre,pos0-pre);
												sscanf(str.c_str(),"%d",&cover_Len);  // part of the length covered by the read

												if(!ctr_Flag)
												{
#pragma region identfiy the left end	
													z_1 = z1 + cover_Len - 1;  // the position of the internal end of the right read
													j_3 = j3;
													while(j_3<e_Num&&z_1>=gph.exon_Vec[j_3].bound1) j_3++; j_3--;   // the exon in which the right end of the right read falls

													for(int k1 = j3; k1<=j_3; k1++)
														rd_Vec[serial].right_readIdx.push_back(k1);
#pragma endregion identfiy the left end	
													p1 = z1 + cover_Len + junc_Len;  // the left boundary of the exon possibly covered by the read
													ctr_Flag = true;
												}
												else
												{
#pragma region identify the middle exons
													int j_pre = j_3;
													p2 = p1 + cover_Len - 1;
													while(j_3<e_Num&&abs(gph.exon_Vec[j_3].bound1-p1)>2)
														j_3++;  // locate the middle exon covered by the read

													if(j_3!=j_pre)
														right_readIdx.push_back(j_3);
													cover_Len -= gph.exon_Vec[j_3].len;
													while(cover_Len>2)  // adjacent exons are probably contained
													{												
														j_3++;
														if(j_3>=e_Num-1||abs(gph.exon_Vec[j_3].bound2-gph.exon_Vec[j_3+1].bound1)>2)
															break;
														cover_Len -= gph.exon_Vec[j_3].len;
														right_readIdx.push_back(j_3);																								
													}
#pragma endregion identify the middle exons
													p1 = p2 + junc_Len;
												}
												pre = pos + 1;   // search after this position next time  
											}
										}
#pragma endregion

										// right end of the right read
										pos0 = words[2].find("M",pre);   // locate the first splitting character 'M'					
										string str;
										if(pos0!=string::npos)
										{
											str = words[2].substr(pre,pos0-pre);
											sscanf(str.c_str(),"%d",&cover_Len); // distance spanned by the read
											z_2 = z2 - cover_Len + 1;  // the position of the outer end of the right read
											j_4 = j4;
											while(j_4>=0&&z_2<gph.exon_Vec[j_4].bound1-1) j_4--; // the exon the right end of the right read falls in
											if(j_4!=j_3)
												rd_Vec[serial].right_readIdx.push_back(j_4);

											for(int k1 = j_4+1; k1<=j4; k1++)
												rd_Vec[serial].right_readIdx.push_back(k1);
										}

#pragma endregion
										// search the exon in which the right end of the left read falls in 
										int idnum = rd_Vec[serial].right_readIdx.size(), index1 = 0, index2 = 0;
										for(int k1 = 0; k1<idnum-1; k1++)
										{
											index1 = rd_Vec[serial].right_readIdx[k1];
											index2 = rd_Vec[serial].right_readIdx[k1+1];
											gph.link_Mtx[index1+1][index2+1]++;  //  each two adjacent exons between j1 and j2 are linked
										}										
									}
#pragma endregion 

									str_itr = unresolved_read.begin() + k1;
									unresolved_read.erase(str_itr);  // delete read from the uncertainty set
									itr = unresolved_serial.begin() + k1;
									unresolved_serial.erase(itr);    // delete read from the uncertainty set
									break;
								}
							}
						}
					}
#pragma endregion
				}
			}
		}
		else if(words[3]=="="&&s1>=stop_Site)  // beyond the valid range of the gene
		{
			flag_Ctrl = false;
			// back to the previous line
			myfile.seekg(pos1);
			break;
		}
		else
		{
		}
		pos1 = myfile.tellg();   // obtain the current file pointer

#pragma endregion
	}

#pragma region  modify the linkage in the graphs
	// modify the information of the exons
	if(rd_Vec.size()>0){
		int s1 = rd_Vec[0].Left1;  // left end of the first read
		exon_Assembly[0].bound1 = s1;
		exon_Assembly[0].len = exon_Assembly[0].bound2 - exon_Assembly[0].bound1+1;
		gph.exon_Vec[0] = exon_Assembly[0];
		starting_site = exon_Assembly[exon_Num-1].bound2+2;

		if(locate_Flag==true)
		{	
			int num = gph.exon_Vec.size();  // number of exons
			for(int i = 0; i<num-1; i++)
			{
				if(abs(gph.exon_Vec[i+1].bound1-gph.exon_Vec[i].bound2)<2)  // two exons are adjacent
				{
					if(gph.link_Mtx[i+1][i+2]==0)
						gph.link_Mtx[i+1][i+2] = 1;  // the two nodes are linked
				}
			} 
		}
	}
#pragma endregion

	return locate_Flag;	
}

// Examine the current exon assembly and conduct refinement
int Gene_Batch::intron_Coverage(vector<int>& rev_Rank)
{
	int exon_Num = exon_Assembly.size();  // number of exons
	int ori_exonNum = exon_Num;
	vector<vector<PAIR_INT>> start_New;   // possible starting points of exons
	vector<vector<PAIR_INT>> stop_New;    // possible stopping points of exons

	int Num = 2*exon_Num-1;           // number of segments
	vector<float> jump_thresh1(Num);  // threshold for read distribution discontinuity in the segments
	vector<float> jump_thresh2(Num);
	int rd_num1 = 0, rd_num2 = 0, rd_num = 0, exon_Len = 0, seg_Len = 0;

	Interval interval;
	vector<vector<int>> seg_Modify;
	float cnt_Thresh = 0, jump_Thresh = 0, boundary_Thresh = 0;
	for(int j=1; j<Num-1; j+=2)
	{
		int i = (j-1)/2;   // serial of the segment
		vector<int>& temp1 = rd_Left[j];   // record the segments of the left points
		vector<int>& temp2 = rd_Right[j];  // record the segments of the right points
		rd_num1 = temp1.size();   // number of left reads in the current segment
		rd_num2 = temp2.size();   // number of right reads in the current segment

		rd_num = rd_num1+rd_num2;
		float temp = (rd_Left[j-1].size()+rd_Left[j-1].size()+rd_Right[j+1].size()+rd_Right[j+1].size())*0.2;
		cnt_Thresh = _max(temp,100);

		if(rd_num>cnt_Thresh)
		{
#pragma region Search possible breaking points
			int exon_start = exon_Assembly[i].bound2+1;   // starting point of exon
			int exon_stop = exon_Assembly[i+1].bound1-1;  // stopping point of exon
			vector<int> pair1, pair2, cnt1, cnt2, gap_1, gap_2;
			seg_Len = exon_Assembly[i+1].bound1 - exon_Assembly[i].bound2;
			float jump_Thresh = _max(100,seg_Len*0.25);
			jump_Thresh = _min(1000,jump_Thresh);

#pragma region discontinuity of distribution of left reads		
			if(rd_num1>0)
			{
				jump_Thresh = jump_thresh1[j];
				boundary_Thresh = fragLen_Mean;
				locate_Jump(temp1, rd_num1, i, exon_start, exon_stop, jump_Thresh, 
					boundary_Thresh, pair1, cnt1, gap_1);
			}
#pragma endregion 

#pragma region discontinuity of distribution of right reads
			if(rd_num2>0)
			{
				jump_Thresh = jump_thresh2[j];
				boundary_Thresh = 50;
				locate_Jump(temp2, rd_num2, i, exon_start, exon_stop, jump_Thresh, 
					boundary_Thresh, pair2, cnt2, gap_2);
			}
#pragma endregion

#pragma endregion search possible breaking points
			if(pair1.size()>0||pair2.size()>0)
			{
				Interval interval;
				int change_type = exon_Type[i];
				interval.ch_type.push_back(change_type);
				interval.site_New1.push_back(pair1); 
				interval.site_New2.push_back(pair2); 
				interval.cnt_New1.push_back(cnt1);
				interval.cnt_New2.push_back(cnt2);
				interval.gap1.push_back(gap_1);
				interval.gap2.push_back(gap_2);
				interval.serial.push_back(i);    // record segment serial
				interval_Locate_1(interval,seg_Modify);
			}
			else
			{
				vector<int> points;
				points.push_back(i);
				seg_Modify.push_back(points);
			}
		}
	}

#pragma region Modify exons
	int num = seg_Modify.size();
	int p1 = 0, p2 = 0;
	vector<int> Rank(exon_Num);  // exon serial
	for(int j=0; j<exon_Num; j++)
		Rank[j] = j;
	for(int k=0; k<num; k++)
	{
		vector<int>& points = seg_Modify[k];  // breaking points in the segments
		int serial = points[0];
		int i = Rank[serial];
		int point_num = (points.size()-1)/2;  // number of breaking points
		if(point_num>0)
		{
			p1 = points[1]; 
			p2 = points[2];
			float dis1 = p1-exon_Assembly[i].bound2;
			float dis2 = exon_Assembly[i+1].bound1 - p2;
			if(dis2>dis1)
			{
				Exon_Node cur_Exon;
				cur_Exon.bound1 = p2;
				cur_Exon.bound2 = exon_Assembly[i+1].bound1-1;
				cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;				
				exon_Assembly.insert(exon_Assembly.begin()+i+1,cur_Exon);
				exon_Type.insert(exon_Type.begin()+i+1,5);  // the back end is linked with the exon
				for(int j=serial+1; j<exon_Num; j++)        // change of exon index
					Rank[j]++;
			}
			else
			{
				Exon_Node cur_Exon;
				cur_Exon.bound1 = exon_Assembly[i].bound2+1;
				cur_Exon.bound2 = p1;
				cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;				
				exon_Assembly.insert(exon_Assembly.begin()+i+1,cur_Exon);
				exon_Type.insert(exon_Type.begin()+i+1,6);  // the front end is linked with the exon
				for(int j=serial+1; j<exon_Num; j++)        // change of exon index
					Rank[j]++;
			}
		}
		else  // insert complete intron
		{
			Exon_Node cur_Exon;
			cur_Exon.bound1 = exon_Assembly[i].bound2+1;
			cur_Exon.bound2 = exon_Assembly[i+1].bound1-1;
			cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;				
			exon_Assembly.insert(exon_Assembly.begin()+i+1,cur_Exon);
			exon_Type.insert(exon_Type.begin()+i+1,7);  // insert complete exon
			for(int j=serial+1; j<exon_Num; j++)        // change of exon index
				Rank[j]++;
		}
	}
#pragma endregion

	vector<vector<int>>().swap(seg_Modify);
	int exon_num = exon_Assembly.size();
	if(exon_num>ori_exonNum)
	{
		rev_Rank.clear();
		for(int l=0; l<exon_num; l++)
			rev_Rank.push_back(-1);
		for(int l=0; l<ori_exonNum; l++)
		{
			int ser = Rank[l];
			rev_Rank[ser] = l;
		}
	}

	return (exon_num-ori_exonNum);

}

// Search possible jumping points for cases of introns
// Locate possible jumping points from detected discontinuities
// Input: detected points of discontinuities in each segment
// Output: possible jumping points
// match_Identify() check the consistency of the jumping positions identified by left reads and right reads
int Gene_Batch::interval_Locate_1(Interval& seg, vector<vector<int>>& seg_Modify)
{
	vector<PAIR_INT> pair_Vec1, pair_Vec2;
	int num = seg.serial.size();          // number of exons which possibly contain breaking points
	int serial = 0, change_type = 0;
	int exon_Num = exon_Assembly.size();  // number of exons

	for(int i=0; i<num; i++)
	{
		serial = seg.serial[i];           // segment serial
		change_type = seg.ch_type[i];
		int p1 = 0, p2 = 0, p3 = 0, p4 = 0;     // starting points and stopping points

		vector<int>& temp1 = seg.site_New1[i];  // breaking points identified by left reads
		vector<int>& temp2 = seg.site_New2[i];  // breaking points identified by right reads

		int m_idx1 = 0, m_idx2 = 0, m_idx3 = 0, m_idx4 = 0, m_gap1 = 0, m_gap2 = 0;
		float m_rate1 = 0, m_rate2 = 0;
		int s1 = exon_Assembly[serial].bound2, s2 = exon_Assembly[serial+1].bound1;  // starting site and stopping site of a gene

		float thresh = exon_Assembly[i].len*0.75;      // threshold of exon length
		float thresh_s = 50;

		int n1 = temp1.size()/2, n2 = temp2.size()/2;  // number of jumping points
		vector<float> rate1(n1), rate2(n2);

		vector<int> match1(n1), match2(n2);  // matched jumping points

		int add_Cnt = 0;    // number of added exons
		bool flag1 = false;
		vector<int> sel_IdxL, sel_IdxR;      // index of selected jumping point
		pair_Vec1.clear(); pair_Vec2.clear();
		if(n1>0)
		{
#pragma region search possible breaking points
			for(int j = 0; j<n1; j++)
				pair_Vec1.push_back(make_pair(j,seg.gap1[i][j]));
			sort(pair_Vec1.begin(),pair_Vec1.end(),cmp_1);  // rank according to gaps
			m_rate1 = cal_Rate_1(serial,temp1,seg.cnt_New1[i],seg.gap1[i],n1,rate1,m_idx1,m_idx2);
			m_idx1 = pair_Vec1[0].first;      // index of the largest gap
			m_gap1 = pair_Vec1[0].second;     // value of the largest gap
			int cnd_idx = -1;
			if(n1>1)
				cnd_idx = pair_Vec1[1].first;

			sel_IdxL.push_back(m_idx1);
			// chose between the largest gap and the largest covering percentage
			if(m_idx1==m_idx2)   // identification results from the largest gap and the largest covering percentage coincide
			{
				add_Cnt++;
				if(n1>1)
					sel_IdxL.push_back(cnd_idx); // add the index of the second largest gap
			}
			else
			{
				sel_IdxL.push_back(m_idx2);
				if(n1>1&&m_idx2!=cnd_idx)
					sel_IdxL.push_back(cnd_idx);
			}
#pragma endregion	
		}
		if(n2>0)   // identification based on the right ends
		{
#pragma region search possible breaking points
			for(int j=0; j<n2; j++)
				pair_Vec2.push_back(make_pair(j,seg.gap2[i][j]));
			sort(pair_Vec2.begin(),pair_Vec2.end(),cmp_1);   // rank according to the values of gaps
			m_rate2 = cal_Rate_1(serial,temp2,seg.cnt_New2[i],seg.gap2[i],n2,rate2,m_idx3,m_idx4);
			m_idx3 = pair_Vec2[0].first;   // index of the largest gap
			m_gap2 = pair_Vec2[0].second;  // value of the largest gap
			int cnd_idx = -1;
			if(n2>1)
				cnd_idx = pair_Vec2[1].first;

			sel_IdxR.push_back(m_idx3);
			// chose between the largest gap and the largest covering percentage
			if(m_idx3==m_idx4)   // identification results from the largest gap and the largest covering percentage coincide
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
#pragma endregion
		}

#pragma region check the consistency of the jumping positions identified by left reads and right reads
		int m1 = sel_IdxL.size(), m2 = sel_IdxR.size();  // jumping points identified by left reads and right reads
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

#pragma region Estimate possible exon boundaries according to identified jumping points
		vector<int> points_L, points_R, points;
		int idx = -1, idx_2 = -1, idx_s = -1;
		if(valid_Cnt.size()>0)  // there are consistent jumping points identified by both left reads and right reads
		{	
#pragma region jumping points identified by left reads and right reads coincide
			if(valid_Cnt.size()==1)
			{
				idx = sel_IdxL[valid_Cnt[0]];
				idx_2 = match[valid_Cnt[0]];    // responding index in the vector for right reads
				points_Locate(idx,idx_2,serial,s1,temp1,temp2,points);
			}
			else
			{
				bool ctr_Flag = true;
				int valid_Num = valid_Cnt.size();
#pragma region cases of internal exons
				idx = -1; idx_2 = -1;
				int l = 0;
				for(l=0; l<valid_Num; l++)  // search the first jumping point satisfying the requirement of read coverage
				{
					idx = sel_IdxL[valid_Cnt[l]];
					idx_2 = match[valid_Cnt[l]];
					if(rate1[idx]>m_rate1*0.6)
						break;
				}
				if(l>=valid_Num)  // choose only according to gaps
				{
					idx = sel_IdxL[valid_Cnt[0]];   // the first selected position 
					idx_2 = match[valid_Cnt[0]];    // responding index in the vector for the right reads
				}
				points_Locate(idx,idx_2,serial,s1,temp1,temp2,points);

#pragma endregion
			}
			seg_Modify.push_back(points);
#pragma endregion
		}
		else  // jumping points identified by left reads and right reads don't match
		{
#pragma region 
			int idx = -1, chosen = -1;
			float gap_thresh = 300, cnt_thresh = 10;   // choose high-confidence gap
			bool mark = false;
			if(n1==0)
			{
				idx = m_idx3; chosen = 1;
			}
			else 
			{
				idx = m_idx1; chosen = 0;	
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
#pragma endregion
	}

	return seg_Modify.size();
}

// Locate the junctions in the gene
bool Gene_Batch::locate_Junction(ifstream& infile, Gene& cur_Gene, vector<int>& exon_Start, vector<int>& exon_Stop, vector<vector<int>>& junc_Num)
{
	stringstream stream(pre_line);
	string word, seg[12], line;
	for(int i=0; i<12; i++)
		stream>>seg[i];   // read the 12 fields

	int num = cur_Gene.subexon.size();  // number of subexons
	int bound1 = cur_Gene.subexon[0].bound1, bound2 = cur_Gene.subexon[num-1].bound2;  // gene boundary
	int junc_start = 0, junc_stop = 0, junc_num = 0, exon_start = 0, exon_stop = 0, temp = 0;
	bool ctr_Flag = false, loc_Flag = false;
	string::size_type pre, pos;
	char tPar = ',';

	exon_Start.push_back(bound1);  // add starting site
	while(ctr_Flag==false)
	{
		sscanf(seg[1].c_str(),"%d",&junc_start);  // starting site of junction
		sscanf(seg[2].c_str(),"%d",&junc_stop);   // stopping site of junction
		sscanf(seg[4].c_str(),"%d",&junc_num);    // number of reads spanning the junctions
		
#pragma region delimit exon boundaries based on the annotations
		pre = 0;
		if((pos=seg[10].find(tPar,pre))!=string::npos)  // search the splitting character ','
		{
			string str = seg[10].substr(pre,pos-pre);
			sscanf(str.c_str(),"%d",&temp);
			exon_stop = junc_start + temp;    // starting site of junction, or stopping site of exon
		}

		pre = 0;
		if((pos=seg[11].find(tPar,pre))!=string::npos)  // search the splitting character ','
		{
			int len = seg[11].length();
			string str = seg[11].substr(pos+1,len-pos-1);
			sscanf(str.c_str(),"%d",&temp);   // length of exon
			exon_start = junc_start + temp + 1;
		}
#pragma endregion

		if(exon_stop>bound1-10&&exon_start<bound2+10)   // the predicted exon is inside the gene
		{
			loc_Flag = true;
			exon_Start.push_back(exon_start);
			exon_Stop.push_back(exon_stop);
			vector<int> temp;
			temp.push_back(exon_stop); temp.push_back(exon_start); temp.push_back(junc_num);
			
			junc_Num.push_back(temp);
			if(infile.eof())
				break;
			getline(infile,line);
			loc_Flag = true;			
			stringstream stream(line);
			for(int i=0; i<12; i++)
				stream>>seg[i];   // read 12 fields
			pre_line = line;
		}
		else if(exon_stop>=bound2+10||exon_start>bound2)  // beyond the range of the gene
		{
			break;
		}	
		else if(exon_start<=bound1-10||exon_stop<=bound1-10)  // outside the gene region
		{
			if(infile.eof())
				break;
			getline(infile,line);
			stringstream stream(line);
			for(int i=0; i<12; i++)
				stream>>seg[i];   // read 12 fields
			pre_line = line;
		}
	}

	return loc_Flag;
}

// Locate junctions between the genes
void Gene_Batch::locate_Junction_Interval(ifstream& infile, vector<Gene>& Gene_Set, int serial, vector<PAIR_INT> junction, vector<int>& orientation)
{
	stringstream stream(pre_line);
	string word, seg[12],line;
	for(int i=0; i<12; i++)
		stream>>seg[i];   // read the 12 fields

	int bound1 = Gene_Set[serial].stop;     // right boundary of the previous gene
	int bound2 = Gene_Set[serial+1].start;  // left boundary of the next gene
	int exon_start = 0, exon_stop = 0, junc_start = 0, junc_stop = 0, junc_num = 0, ori = 0, temp = 0;
	bool ctr_Flag = false;

	orientation.clear();
	while(ctr_Flag==false)
	{
		sscanf(seg[1].c_str(),"%d",&junc_start);  // starting site of junction
		sscanf(seg[2].c_str(),"%d",&junc_stop);   // stopping site of junction
		sscanf(seg[4].c_str(),"%d",&junc_num);    // number of reads spanning junctions
		
#pragma region delimit exon boundaries based on the annotations
		string::size_type pre = 0, pos;
		char tPar = ',';
		if((pos=seg[10].find(tPar,pre))!=string::npos)  // search the splitting character ','
		{
			string str = seg[10].substr(pre,pos-pre);
			sscanf(str.c_str(),"%d",&temp);
			exon_stop = junc_start + temp;    // starting site of junction, or stopping site of exon
		}

		pre = 0;
		if((pos=seg[11].find(tPar,pre))!=string::npos)  // search the splitting character ','
		{
			int len = seg[11].length();
			string str = seg[11].substr(pos+1,len-pos-1);
			sscanf(str.c_str(),"%d",&temp);   // exon length
			exon_start = junc_start + temp + 1;
		}
#pragma endregion

		if(exon_stop>bound1&&exon_start<bound2)   // the predicted exon is inside the gene
		{
			junction.push_back(make_pair(exon_stop,exon_start));
			if(seg[5]=="+")
				orientation.push_back(1);  // forward orientation
			else
				orientation.push_back(0);  // reverse orientation

			getline(infile,line);
			stringstream stream(pre_line);
			for(int i=0; i<12; i++)
				stream>>seg[i];   // read the 12 fields
			pre_line = line;
		}
		else if(exon_start>=bound2)  // beyond the range of the gene
		{
			break;
		}		
	}
}



