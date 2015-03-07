
#include "Entropy.h"

///////////////////////////////////////////////////////////////////////////////
// Load the results of read alignment
void Rd_Set::load_Assignment(char* filename, Graph_Trans &gph, int exon_Num)
{
	// initialize the graph
	gph.graph_Initial(exon_Num);
	// gph.print_Graph();

	ifstream myfile;
	myfile.open(filename,ios::in);

	if(!myfile){
		cout << "Unable to open myfile";
		exit(1); // terminate with error
	}

	int x1 = 0, x2 = 0, z1 = 0, z2 = 0, temp_Len = 0; // starting sites of the left read and right read; fragment length
	int read_len = rd_Len;  // read length

	stringstream ss;
	ss<<rd_Len;                  
	string len_R = ss.str()+"M"; // string indicating exact match of the whole read 
	string line, word, words[6], frag_Name;
	int serial = 0;     // read serial
	string tPar = "N";  // marker for mismatch
	string::size_type pos, pos0;
	vector<string> unresolved_read;   // reads not located
	vector<int> unresolved_serial;    // serials of reads not located
	vector<int> right_readIdx;
	vector<int>::iterator itr;
	vector<string>::iterator str_itr;

	char j1 = -1, j2 = -1, j3 = -1, j4 = -1;  // the indices of the segments to which the four ends of the two reads belong
	int e_Num = exon_Num;     // number of exons
	int junc_Len = 0, cover_Len = 0, p1 = 0, p2 = 0;
	int x_1 = 0, x_2 = 0, z_1 = 0, z_2 = 0, je = 0, rlen1 = 0, rlen2 = 0, j_3 = 0, j_4 = 0, j_1 = 0, j_2 = 0;
	bool anno_flag = false;
	int line_cnt = 0;

	while(std::getline(myfile,line))
	{
		line_cnt++;
		Rd cur_Rd; 
		j1 = -1, j2 = -1, j3 = -1, j4 = -1;
		anno_flag = false;

		stringstream stream(line);
		stream>>word; frag_Name = word;
		stream>>word; stream>>word; // read the first three fields
		for(int i = 0; i<6; i++)
			stream>>words[i];  // meaning of the fields: 4 starting site of left read 5 matching quality 6 matching quality 7 marker for paired match 8 starting site of right read 9 inserted size
		sscanf(words[5].c_str(),"%d",&temp_Len);  // inserted size
		sscanf(words[0].c_str(),"%d",&x1);  // starting site of left read
		sscanf(words[4].c_str(),"%d",&z1);  // starting site of right read
		if(temp_Len>=0)        
		{
			z2 = x1 + temp_Len - 1;  // stopping site of right read
			x2 = x1 + rd_Len - 1;    // stopping site of left read	

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
				string::size_type pre;    // record the last position of "N"
				pre = 0;
#pragma region identify the segments in which the left read falls according to the identified junctions
				bool ctr_Flag = false;
				j2 = j1;
#pragma region left end of the left read and the middle matched part of the read
				while((pos=words[2].find(tPar,pre))!=string::npos)  // search character 'N'
				{
					pos0 = words[2].find("M",pre); // locate the first splitting character
					if(pos0!=string::npos)
					{ 
						string str = words[2].substr(pos0+1,pos-pos0-1);
						sscanf(str.c_str(),"%d",&junc_Len);   // the distance spanned by the read
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

							p1 = x1 + cover_Len + junc_Len;  // the left end of the exon which is probably covered by the read
							ctr_Flag = true;
#pragma endregion identification of the left end
						}
						else
						{
#pragma region identfication of the middle exon covered by the read
							p2 = p1 + cover_Len - 1;
							while(j_1<e_Num&&abs(gph.exon_Vec[j_1].bound1-p1)>2)
								j_1++;  // locate the middle exon covered by the read
							while(abs(cover_Len)>2)  // adjacent exons are probably contained
							{
								cur_Rd.left_readIdx.push_back(j_1);
								cover_Len -= gph.exon_Vec[j_1].len;
								j_1++;
								if(j_1>=e_Num-1||abs(gph.exon_Vec[j_1].bound2-gph.exon_Vec[j_1+1].bound1)>2)
									break;
							}
#pragma endregion  identfication of the middle exon covered by the read
							p1 = p2 + junc_Len + 1; 
						}
						pre = pos + 1;   // search after this position next time 
					}
				}
#pragma endregion

				// right end(internal end) of the left read
				pos0 = words[2].find("M",pre); // locate the first splitting character
				string str;
				if(pos0!=string::npos)
				{
					str = words[2].substr(pre,pos0-pre);
					sscanf(str.c_str(),"%d",&cover_Len); // the distance spanned by the read
					x_2 = p1 + cover_Len - 1;  // the position of the right end of the second part of the left read
					// cout<<cover_Len<<endl;
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

#pragma region paired match
			if(words[3]=="=")   //  paired alignment, hence another alignment result can be referred to
			{
				j3 = j2;
				while(j3<e_Num&&z1>=gph.exon_Vec[j3].bound1) j3++; j3--; // exon in which the left end of the right read falls
				j4 = j3;  // both ends of the left read fall into the same exon
				while(j4<e_Num&&z2>=gph.exon_Vec[j4].bound1) j4++; j4--;
				if(j3!=j4)
				{
					int part_Len = gph.exon_Vec[j3].bound2 - z1 + z2 - gph.exon_Vec[j4].bound1 + 2; 
					if(part_Len==read_len)  // junction is only related to two exons
					{
						gph.link_Mtx[j3+1][j4+1]++;  // j3 is linked with j4
						anno_flag = true;
					}
					else if(part_Len<read_len)   // the exact alignment of the right read needs further examination
					{
						unresolved_read.push_back(frag_Name);  // record the name of the fragment
						unresolved_serial.push_back(serial);   // record the serial of the fragment
						anno_flag = false;
					}
				}
			}
#pragma endregion

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
			cur_Rd.Left1 = x1; cur_Rd.Right2 = z2;
			cur_Rd.Left2 = x2; cur_Rd.Right1 = z1;

			rd_Vec.push_back(cur_Rd);

			serial++;
#pragma endregion 
		}
		else  // alignment is shown in the reverse (right-left) way
		{
			if(words[3]=="=")   // paired alignment
			{  
				int num = unresolved_serial.size();     // number of right reads not located
				for(int k1=0; k1<num; k1++)
				{
					if(frag_Name==unresolved_read[k1])  // the fragment names are the same
					{
						int serial = unresolved_serial[k1];
						x_1 = rd_Vec[serial].Left1;  x2 = rd_Vec[serial].Left2;
						z_1 = rd_Vec[serial].Right1; z2 = rd_Vec[serial].Right2;

						if(x1==z_1&&z1==x_1&&rd_Vec[serial].right_readIdx.size()==0)  //  the same pair of reads
						{
							x1 = x_1; z1 = z_1;
							j1 = rd_Vec[serial].l1;  j2 = rd_Vec[serial].l2;  j3 = rd_Vec[serial].r1;  j4 = rd_Vec[serial].r2;
#pragma region identify junction according to the results of read alignment
							if(words[2]==len_R)   // exact match of the whole read
							{
								j4 = j3;   // both ends of the right read fall in the same exon
								while(j4<e_Num&&z2>=gph.exon_Vec[j4].bound1) j4++;  j4--;
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
								while((pos=words[2].find(tPar,pre))!=string::npos)  // search the character 'N'
								{
									pos0 = words[2].find("M",pre); // search the character 'M'
									if(pos0!=string::npos)
									{ 
										string str = words[2].substr(pos0+1,pos-pos0-1);
										sscanf(str.c_str(),"%d",&junc_Len);   // distance spanned by the read
										str = words[2].substr(pre,pos0-pre);
										sscanf(str.c_str(),"%d",&cover_Len);  // length covered by the read

										if(!ctr_Flag)
										{
#pragma region identfiy the left end		
											z_1 = z1 + cover_Len - 1;  // the position of the internal end of the right read
											j_3 = j3;
											while(j_3<e_Num&&z_1>=gph.exon_Vec[j_3].bound1) j_3++; j_3--;

											for(int k1 = j3; k1<=j_3; k1++)
												rd_Vec[serial].right_readIdx.push_back(k1);
#pragma endregion identfiy the left end
											p1 = z1 + cover_Len + junc_Len; // the left boundary of the exon possibly covered by the read
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
								//cout<<"right_readIdx:"<<endl;
								for(int k1 = 0; k1<idnum-1; k1++)
								{
									index1 = rd_Vec[serial].right_readIdx[k1];
									index2 = rd_Vec[serial].right_readIdx[k1+1];
									//cout<<index1+1<<" "<<index2+1<<" "<<endl;
									gph.link_Mtx[index1+1][index2+1]++;  // each two adjacent exons between j1 and j2 are linked
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
		}
	}

#pragma region  modify the linkage in the graphs
	int num = gph.exon_Vec.size();
	for(int i = 0; i<num-1; i++)
	{
		if(abs(gph.exon_Vec[i+1].bound1-gph.exon_Vec[i].bound2)<2)  // two exons are adjacent
		{
			if(gph.link_Mtx[i+1][i+2]==0)
				gph.link_Mtx[i+1][i+2] = 1;  // the two nodes are linked
		}
	} 

#pragma endregion
}

// Calculate the generation probability of the read
double Rd_Set::rd_Entro(Rd &rd, Path &path, vector<Exon_Node> &exon_Vec, int id1, int id2, int id3, int id4)
{
	int frag_Len = 0;  // fragment length
	int sub_Len = 0;
	double pro_Len = 0, entro = 0;

	if(id3>=0)  // paired-end match
	{
		if(id4>id1)  // the left read and right read aren't in the same exon
		{
			for(int k = id1+1; k<id4;k++)
			{
				sub_Len += path.bound[2*k+1]-path.bound[2*k]+1;
			}
			frag_Len = sub_Len + path.bound[2*id1+1]-rd.Left1 + rd.Right2 - path.bound[2*id4] + 2;
		}
		else  // the left read and right read are in the same exon
		{
			frag_Len = rd.Right2 - rd.Left1 + 1;
		}

		// suppose the fragment length follows Gaussian distribution
		int l = frag_Len>Lim1?frag_Len:Lim1;
		int len = l<Lim2?l:Lim2;

		int inter_Idx = int((len-Lim1+1)/10-0.005);  // the interval the fragment length belongs to 
		pro_Len = fragLen_Interval[inter_Idx];       // the probability of the fragment length

		int start_len = path.length - len + 1;
		if(start_len>0)
			entro = pro_Len*1.0/start_len;

		if(start_len<0)
			printf("wrong start length: %f %d %d %d\n",pro_Len,path.length,len,rd.Left1);
	}
	else
	{
		int temp = path.length - rd.Left1 + path.bound[0];  // the possible max length
		int min_Len = temp>Lim1?temp:Lim1; 
		int aver_Len = mean;    // mean of fragment length
		int max_Len = temp<Lim2?temp:aver_Len;
		double sum = 0;
		double thresh = 1e-06;
		for(int l = min_Len; l<max_Len; l++)
		{
			int inter_Idx = int((l-Lim1+1)/10-0.005);  // the interval the fragment length belongs to 
			pro_Len = fragLen_Interval[inter_Idx];     // the probability of the fragment length
			entro += pro_Len/(path.length-l+1);        // entropy of a single read
		}
	}
	return entro;
}

// Normalize p(x|y)
void Rd_Set::rd_Normalize()
{	
	double min1 = 1e-09;
	int p_Num = path.size();	
	for(int i = 0; i<p_Num; i++)
	{
		double sum = 0;
		for(int j = 0; j<rd_Num; j++)
		{
			sum += rd_Pro[j][i];
		}

		for(int j = 0; j<rd_Num; j++)
		{
			if(sum>1e-05)
				rd_Pro[j][i] /= sum;
			else
				rd_Pro[j][i] = 1.0/rd_Num;
		}
	}

	double sum1 = 0; double sum = 0;
	double eps = 1e-09;
	int id_Num = idx.size();
	for(int i = 0; i<id_Num; i++)
	{
		sum1 = 0;
		sum = 0;
		int k = idx[i];    // index of the valid transcript
		for(int j = 0; j<rd_Num; j++)
		{
			if(rd_Pro[j][k]>0)
			{
				sum1 += rd_Pro[j][k]*log(rd_Pro[j][k]);
				sum += log(rd_Pro[j][k]);
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////
// EM algorithm: calculate the expectation of transcript abundance
float Rd_Set::abun_Estimation(float* abun_Est, float* abun, int interval)
{
	int id_Num = sel_Idx.size();  // selected valid transcript
	float valid_count = 0;        // number of valid transcripts

	double tsum = 0, temp = 0;
	for(int j = 0; j<rd_Num; j += interval)
	{		
		tsum = 0;
		for(int i = 0; i<id_Num; i++)
		{
			int k = sel_Idx[i];   // index of the valid path in the vector of sel_Idx
			int k1 = idx[k];      // index of the valid path
			temp = rd_Pro[j][k1]*abun[k];
			tsum += temp;	
			rd_Pro1[j][k1] = temp;
		}
		if(tsum>0)
		{
			for(int i = 0; i<id_Num; i++)      // normalization
			{
				int k = sel_Idx[i];            // index of the valid path
				int k1 = idx[k];
				rd_Pro1[j][k1] /= tsum;
				abun_Est[k] += rd_Pro1[j][k1]; // accumulate the probability 
			}
			valid_count++;
		}		
	}

	return valid_count;
}

/////////////////////////////////////////////////////////////////////////////////
// EM algorithm: calculate the expectation of transcript abundance
void Rd_Set::abun_Estimation_Re(float* abun_Est, float* abun_Est_Re, float* abun, int interval)
{
	int id_Num = sel_Idx.size();  // selected valid transcript
	double tsum = 0, temp = 0;    

	for(int j = 0; j<rd_Num; j += interval)
	{		
		tsum = 0;
		for(int i = 0; i<id_Num; i++)
		{
			int k = sel_Idx[i];  // index of the valid path in the vector of sel_Idx
			int k1 = idx[k];
			temp = rd_Pro[j][k1]*abun[k];
			tsum += temp;	
			rd_Pro1[j][k1] = temp;
		}
		if(tsum>0)
		{
			for(int i = 0; i<id_Num; i++)  // normalization
			{
				int k = sel_Idx[i];        // index of the valid path
				int k1 = idx[k];
				rd_Pro1[j][k1] /= tsum;    //   
				if(rd_Pro1[j][k1]>0)
				{
					abun_Est[k] += rd_Pro1[j][k1];  // accumulate the probability
					abun_Est_Re[k] += rd_Pro1[j][k1]*log(rd_Pro1[j][k1]);  // accumulate the entropy
				}
			}
		}		
	}
}

////////////////////////////////////////////////////////////////////////////////////
//// Calculate mutual information
// Calculate I(X,Y) = H(X) - H(X|Y)
double Rd_Set::Ixy(float* abun)
{
	double px = 0, pxy = 0;
	int id_Num = sel_Idx.size(); // number of valid paths
	for(int j = 0; j<id_Num; j++)
	{
		int k = sel_Idx[j];      // index of the valid path in the vector of sel_Idx
		int k1 = idx[k];         // index of the valid path
		pxy += path[k1].entro_Abun1*abun[k];
	}
	px = Hx1(abun);

	double I_xy = px-pxy;

	return I_xy;
}

////////////////////////////////////////////////////////////////////////////////////
//// Calculate mutual information
// Calculate I(X,Y) = H(X) - H(X|Y) 
double Rd_Set::Ixy1(float* abun, vector<int>& sel, vector<int>& nsel)
{
	double px = 0, pxy = 0;
	int id_Num = sel.size(); // number of valid paths
	for(int j = 0; j<id_Num; j++)
	{
		int k = sel[j];      // index of the valid path in the vector of sel_Idx
		int k1 = idx[k];     // index of the valid path
		pxy += path[k1].entro_Abun1*abun[k];
	}
	px = Hx1(abun);

	double I_xy = px-pxy;

	return I_xy;
}

// Calculate H(X)
double Rd_Set::Hx1(float* abun)
{
	int id_Num = sel_Idx.size();	
	double sum = 0, temp = 0, sum1 = 0;
	double eps = 1e-09;
	for(int i = 0; i<rd_Num; i++)
	{
		sum1 = 0;
		for(int j = 0; j<id_Num; j++)
		{
			int k = sel_Idx[j];  // index of the valid path in the vector of sel_Idx
			int k1 = idx[k];     // index of the valid path
			sum1 += rd_Pro[i][k1]*abun[k];
		}
		if(sum1>0)
			sum += sum1*log(sum1);
	}
	return -sum;
}

// Calculate H(X)
double Rd_Set::Hx1(float* abun, vector<int>& sel)
{
	int id_Num = sel.size();	
	double sum = 0, temp = 0, sum1 = 0;
	double eps = 1e-09;
	for(int i = 0; i<rd_Num; i++)
	{
		sum1 = 0;
		for(int j = 0; j<id_Num; j++)
		{
			int k = sel[j];
			int k1 = idx[k];
			sum1 += rd_Pro[i][k1]*abun[k];
		}
		if(sum1>0)
			sum += sum1*log(sum1);
	}
	return -sum;
}

// Calculate P(X)
double Rd_Set::Px(float* abun, int interval)
{
	int id_Num = sel_Idx.size();
	double sum = 0, sum1 = 0, temp = 0;
	double eps = 1e-09;
	ofstream outfile;
	outfile.open("D:\\Splicing graph\\px.txt",ios::out);
	for(int i = 0; i<rd_Num; i+=interval)
	{		
		sum = 0;
		for(int j = 0; j<id_Num; j++)
		{
			int k = sel_Idx[j]; 
			int k1 = idx[k];     
			temp = rd_Pro[i][k1]*abun[k];
			if(rd_Pro[i][k1]<0)
				cout<<i<<" "<<k1<<" "<<rd_Pro[i][k1]<<endl;

			sum += temp;
		}
		sum1 += log(sum+eps);
	}

	outfile.close(); 
	return sum1;
}

// Calculate H(Y)
double Rd_Set::Hy(float* abun)
{
	double hy = 0;
	int id_Num = sel_Idx.size();
	for(int i = 0; i<id_Num; i++)
	{
		int k = sel_Idx[i];
		if(abun[k]>1e-06)
			hy -= abun[k]*log(abun[k]);
	}

	return hy;
}

// Calculate I(Y,X)=H(Y)-H(Y|X)
double Rd_Set::Iyx1(float* abun, int interval)
{
	double hy = 0, hyx = 0;
	int id_Num = sel_Idx.size();
	for(int i = 0; i<id_Num; i ++)
	{
		int k = sel_Idx[i];
		if(abun[k]>1e-06)
			hy -= abun[k]*log(abun[k]);
	}

	hyx = Hyx(abun,interval);

	if(abs(hy)>1e-06)
		cxy = 1-hyx/hy;
	else
		cxy = 0;

	return (hy-hyx);
}

// Calculate I(Y,X)=H(Y)-H(Y|X)
double Rd_Set::Iyx2(float* abun)
{
	double hy = 0, hyx = 0;
	int id_Num = sel_Idx.size();
	for(int i = 0; i<id_Num; i++)
	{
		int k = sel_Idx[i];
		hy -= path[k].entro_Abun2*log(path[k].entro_Abun2);
	}

	hyx = Hyx(abun); // calculate conditional entropy

	return (hy-hyx);
}

// Calculate H(Y|X)
double Rd_Set::Hyx(float* abun, int interval)
{
	int id_Num = sel_Idx.size();
	double tsum = 0, temp = 0;
	for(int j = 0; j<rd_Num; j += interval)
	{		
		tsum = 0;
		for(int i = 0; i<id_Num; i++)
		{
			int k = sel_Idx[i];  // index of the valid path
			int k1 = idx[k];     // index of the valid path
			temp = rd_Pro[j][k1]*abun[k];  // joint probability of the read and the transcript 
			tsum += temp;
			rd_Pro1[j][k1] = temp;
		}
		if(abs(tsum)>1e-09)
		{
			for(int i = 0; i<id_Num; i++)  // ¹éÒ»»¯
			{
				int k = sel_Idx[i];  // index of the valid path
				int k1 = idx[k];
				rd_Pro1[j][k1] = rd_Pro1[j][k1]/tsum;
			}
		}
		else
		{
			for(int i = 0; i<id_Num; i++)
			{
				int k = sel_Idx[i];
				int k1 = idx[k];
				rd_Pro1[j][k1] = 1.0/id_Num;  // equal probabilities
			}
		}
	}

	double sum = 0, sum1 = 0;
	for(int i = 0; i<id_Num; i++)
	{
		int k = sel_Idx[i];
		int k1 = idx[k];
		sum1 = 0;
		for(int j = 0; j<rd_Num; j += interval)
		{
			sum1 += rd_Pro1[j][k1];
			if(rd_Pro1[j][k1]>0)
			{
				sum += (rd_Pro[j][k1]*abun[k])*log(rd_Pro1[j][k1]);
			}
		}
		path[k1].entro_Abun2 = sum1/rd_Num;
	}
	return -sum;
}
