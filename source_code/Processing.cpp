
#include "Entropy.h"

///////////////////////////////////////////////////////////////////////////////
// Function for file processing
// Record the annotaion of each gene separately
// Select the genes according to their transcripts and subexons
void Rd_Set::transcript_Select(char* root_path, char* filename)
{
	ifstream infile;
	infile.open(filename,ios::in);
	if(!infile)
		cout<<"Open file error!"<<endl;

	char file1[200];
	sprintf(file1,"%s\\test_gene1.txt",root_path);  // file for saving the selected genes

	ofstream outfile;
	outfile.open(file1,ios::out);
	if(!outfile)
		cout<<"Open file error!"<<endl;

	string line, word;
	string seg[6];

	int trs_Num = 0, exon_Num = 0;  // isoform number and exon number
	string gene_Name;  // gene name

	while(std::getline(infile,line))
	{
		stringstream stream(line);
		if(line.length()>0&&line[0]=='@')   // line of caption
		{
			for(int k=0; k<6; k++)
				stream>>seg[k];

			gene_Name = seg[1];   // gene name
			gene_Name = valid_filename(gene_Name);  // convert to valid filename

			sscanf(seg[4].c_str(),"%d",&exon_Num);  // starting site
			sscanf(seg[5].c_str(),"%d",&trs_Num);   // stopping site

			// select genes with a given number of exons and transcripts
			if(exon_Num<=20&&trs_Num>=2&&trs_Num<=10)
			{
				outfile<<gene_Name<<'\t'<<trs_Num<<'\t'<<exon_Num<<endl;
			}
		}	
	}

	outfile.close();
}

// Load gene annotations
void Rd_Set::load_exonAnnotation(char* filename,char* path)
{
	ifstream myfile;
	myfile.open(filename,ios::in);

	std::streampos pos;

	if(!myfile)
		cout<<"Open file error"<<endl;

	ofstream outfile;
	string line, word, gene_Name;
	int subexon_Num = 0, trs_Num = 0;  // number of subexons and number of transcripts

	char filename1[100];
	int cnt = 0;

	while(std::getline(myfile,line))
	{
		if(line[0]=='@')
		{
			cnt++;			
			stringstream stream(line);
			stream>>word; // gene serial
			stream>>gene_Name;
			cout<<cnt<<" "<<gene_Name<<endl;

			gene_Name = valid_filename(gene_Name);   // convert to valid filename

			sprintf(filename1,"%s\\%s.txt",path,gene_Name.c_str());   // file for saving the exons
			outfile.open(filename1,ios::out);

			stream>>word; stream>>word; stream>>word;  // convert to valid filename
			sscanf(word.c_str(),"%d",&subexon_Num);    // number of subexons	
			stream>>word;
			sscanf(word.c_str(),"%d",&trs_Num);        // number of transcripts

			getline(myfile,line);
			pos = myfile.tellg();   // record the file pointer
			getline(myfile,line);
			myfile.seekg(pos);      // back to the start of the line

			if(trs_Num>0&&line[0]!='@')
			{				
				for(int k = 0; k<subexon_Num; k++)
				{
					getline(myfile,line);
					outfile<<line<<endl;
				}	
			}
			outfile.close();
		}
	}
}

// Generate the abundances of transcripts
void Rd_Set::gen_Exply(char* path, char* file1)
{
	ifstream myfile;
	myfile.open(file1,ios::in);  // load the information of genes and transcripts

	ofstream outfile;
	char filename[200];
	string line, word;
	string words[20];
	string gene_Name;
	int subexon_Num = 0, trs_Num = 0;
	vector<int> trs_Len;
	int group_id = 1;
	int gen_Serial = 0;

	srand(time(NULL));

	while(std::getline(myfile,line))
	{
		if(line.length()>0&&line[0]=='@')
		{
			trs_Len.clear();
			stringstream stream(line);
			stream>>word;  // gene serial
			stream>>gene_Name;
			stream>>word;
			sscanf(word.c_str(),"%d",&subexon_Num);  // number of subexons
			stream>>word;
			sscanf(word.c_str(),"%d",&trs_Num);
			float* explv = new float[trs_Num];

			// convert to valid filename
			gene_Name = valid_filename(gene_Name);

			cout<<"Gene "<<gen_Serial+1<<" : "<<gene_Name<<endl;  // gene name			

			sprintf(filename,"%s\\trs_num%d\\%s.txt",path,trs_Num,gene_Name.c_str());
			outfile.open(filename,ios::out);

			float sum = 0;
			for(int k = 0; k<trs_Num; k++)
			{
				explv[k] = 2*UnifoRand()+1.0;
				sum += explv[k] ;
			}
			for(int k = 0; k<trs_Num; k++)
				explv[k] = explv[k]*100/sum;   // normalization of the transcript abundance

			outfile<<"#ID	Length	Dir	Exons	Position	GroupID	NIsoformInGroup	Explv"<<endl;
			for(int i = 0; i<trs_Num; i++)
			{
				stream>>word;
				int temp = 0;
				sscanf(word.c_str(),"%d",&temp);
				trs_Len.push_back(temp);
			}
			int s1 = 0, s2 = 0;
			for(int i = 0; i<trs_Num; i++)
			{
				getline(myfile,line);
				stringstream stream(line);
				for(int k = 0; k<10; k++)
					stream>>words[k]; // fields: 1.transcirpt serial 2.  starting site 3：stopping site 4：transcript name 5：quality 6：orientation 10: number of exons
				sscanf(words[9].c_str(),"%d",&subexon_Num);
				sscanf(words[1].c_str(),"%d",&s1);
				sscanf(words[2].c_str(),"%d",&s2);
				outfile<<words[3]<<'\t'<<trs_Len[i]<<'\t'<<words[5]<<'\t'<<subexon_Num<<'\t'
					<<"3R:"<<s1+1<<"-"<<s2<<'\t'<<group_id<<" "<<trs_Num<<'\t'<<explv[i]<<endl;
			}	

			gen_Serial++;
			delete[] explv;
			outfile.close();

		}
	}

	myfile.close();
	outfile.close();
}

// Output the gene information to files
void Rd_Set::output_transBed(vector<Gene>& gene, char* filename)
{
#pragma region output transcript information 
	// output the exon information to files
	ofstream outfile;
	outfile.open(filename,ios::out);
	ofstream outfile1;

	char file1[200];

	int cnt[25];   // counting by exon number from 1 to 25
	for(int i = 0; i<25; i++)
		cnt[i] = 0;

	int gene_Num = gene.size();  // number of loaded genes
	for(int k = 0; k<gene_Num; k++)
	{
		int exon_Num = gene[k].subexon.size();
		outfile<<"@"<<k+1<<" "<<gene[k].gene_Name<<" "<<exon_Num<<" "
			<<gene[k].trans.size()<<" ";   
		int trans_Num = gene[k].trans.size(); 
		for(int i = 0; i<trans_Num; i++)
			outfile<<gene[k].trans[i].length<<" ";
		outfile<<endl;

#pragma region
		string gene_Name = gene[k].gene_Name;
		string tPar = ":"; 
		string::size_type pos, pos0;
		string::size_type pre;    // record the last position of characer 'N'
		pre = 0;

		// convert to valid filename
		if(gene_Name=="aux"||gene_Name=="com"||gene_Name=="com2")
			gene_Name = gene_Name + "1";
		string gene_Name_out = gene_Name;
		while((pos=gene_Name_out.find(tPar,pre))!=string::npos)  // search the character ':'
		{
			gene_Name = gene_Name_out.replace(pos, tPar.length(), "-");
			pre = pos + 1;   // search after the position next time
			gene_Name_out = gene_Name;
		}
#pragma endregion

		int itemRgb = 0;
		sprintf(file1,"%s\\%s.bed","G:\\My documents\\FragmentStreaming\\D Melanogaster\\chr3R\\seq_gene\\single_genes",gene_Name.c_str());

		for(int i = 0; i<trans_Num; i++)
		{
			int blockSize = gene[k].trans[i].exon_start.size();  // number of exons contained by the transcript

			outfile<<i+1<<'\t'<<gene[k].trans[i].start<<'\t'<<gene[k].trans[i].stop<<'\t'
				<<gene[k].trans[i].trsName<<'\t'<<0<<'\t'<<gene[k].trans[i].orientation<<'\t'
				<<gene[k].trans[i].start<<'\t'<<gene[k].trans[i].stop<<'\t'
				<<itemRgb<<'\t'<<blockSize<<'\t'; // output exons

			if(blockSize>0)
			{
				int len = 0;
				cnt[blockSize]++;
				// output length of the exon
				for(int l = 0; l<blockSize; l++)
				{
					int tLen = gene[k].trans[i].exon_len[l];
					outfile<<tLen<<",";
					len += tLen;
				}
				outfile<<'\t';

				// output starting sites of the subexons
				for(int l = 0; l<blockSize; l++)
				{
					outfile<<gene[k].trans[i].exon_start[l]<<",";
				}

				outfile<<endl;	
			}
			else
			{
				int len = gene[k].trans[i].stop - gene[k].trans[i].start + 1;
				outfile<<len<<","<<"\t"<<0<<endl;
				gene[k].trans[i].length = len;
				cnt[1]++;
			}

		}
	}

	for(int i = 0; i<25; i++)
		outfile<<i<<":"<<cnt[i]<<" ";
	outfile<<endl;

	outfile.close();	

#pragma endregion
}

// Load annotations of transcripts
void Rd_Set::load_GTF(char* file1, char* file2)
{
	ifstream myfile;
	myfile.open(file1,ios::in);
	if(!myfile)
	{
		cout<<"Open file error"<<endl;
		exit(-1);
	}

	ofstream outfile;
	outfile.open(file2,ios::out);

	string line, word;
	int cnt = 0;

	while(cnt<100&&std::getline(myfile,line))
	{
		stringstream stream(line);
		cout<<line<<endl;
		outfile<<line<<endl;
		cnt++;
	}

	myfile.close();
	outfile.close();
}

// Generate file of GTF format according to the annotations
void Rd_Set::gen_transBed(char* file1, char* file2)
{
	ifstream myfile;
	myfile.open(file1,ios::in);
	if(!myfile)
	{
		cout<<"Open file error"<<endl;
		exit(-1);
	}

	string line, word;
	string seg[20];
	int seg_Num = 17;   // number of fields
	int line_Num = 0;

	vector<Gene> gene;  // vector of genes
	int gene_Serial = -1, pre_Serial = -1;   // gene serial

	int chr_idx = 1, type_idx = 11, label_idx = 9, ori_idx = 4, trs_idx = 15;

	bool ctr_Flag = false;

	while(std::getline(myfile,line))
	{
		if(line[0]!='#'&&line.length()>0)
		{
			stringstream stream(line);
			stream>>seg[0];  // tax_id
			stream>>seg[1];  // chr_Name
			if(seg[1].substr(0,2)==chromosome_Name)
			{
				ctr_Flag = true;
				for(int k = 2; k<seg_Num; k++)
				{
					stream>>seg[k];
					// cout<<k<<" "<<seg[k]<<" "<<endl;
				}
#pragma region the field is Gene
				if(seg[type_idx]=="GENE")  // new gene
				{
					string gene_Name = seg[label_idx]; 
					int found = isExisting(gene,gene_Name);  // check whether the gene exists

					if(found>-1)   // the gene exists, not need to add a new gene
					{
						gene_Serial = found;  // locate to the gene
					}
					else
					{		
						Gene cur_Gene;
						sscanf(seg[2].c_str(),"%d",&(cur_Gene.start));
						sscanf(seg[3].c_str(),"%d",&(cur_Gene.stop));
						gene_Serial = pre_Serial;
						cur_Gene.gene_Name = gene_Name;		
						cur_Gene.orientation = seg[ori_idx][0];
						gene_Serial++;
						pre_Serial = gene_Serial;
						gene.push_back(cur_Gene);
						cout<<"Gene "<<gene_Serial+1<<" : "<<gene_Name<<endl;
					}					
				}
#pragma endregion
#pragma region the field is Transcript
				else if(seg[type_idx]=="RNA")
				{
					Trans cur_Trans; 
					sscanf(seg[2].c_str(),"%d",&(cur_Trans.start));
					sscanf(seg[3].c_str(),"%d",&(cur_Trans.stop));
					cur_Trans.trsName = seg[label_idx];
					cur_Trans.orientation = seg[ori_idx][0];
					cur_Trans.length = 0;
					gene[gene_Serial].trans.push_back(cur_Trans);
				}
#pragma endregion
#pragma region the field is Exon
				else if(seg[type_idx]=="UTR"||seg[type_idx]=="CDS")
				{
					Exon_Node subexon;
					string trs_Name = seg[trs_idx];
					int s1, s2;
					sscanf(seg[2].c_str(),"%d",&s1);
					sscanf(seg[3].c_str(),"%d",&s2);
					int found = isExisting_Trans(gene[gene_Serial].trans,trs_Name);
					if(found>-1)
					{
						int start = gene[gene_Serial].trans[found].start;
						int len = s2 - s1 + 1;
						gene[gene_Serial].trans[found].length += len;
						gene[gene_Serial].trans[found].exon_len.push_back(len);
						gene[gene_Serial].trans[found].exon_start.push_back(s1-start);
					}

#pragma region read the fields of Exon and Intron
					Gene& cur_Gene = gene[gene_Serial];

					int subexon_Num = cur_Gene.subexon.size();

					int j1 = 0, j2 = 0;
					while(j1<subexon_Num&&cur_Gene.subexon[j1].bound2<s1)
						j1++;

					if(j1==subexon_Num) // not overlap with existing exons
					{
						Exon_Node subexon;
						subexon.bound1 = s1; subexon.bound2 = s2; subexon.len = s2 - s1 + 1;
						cur_Gene.subexon.push_back(subexon);	
					}
					else // overlap with existing exons
					{
						vector<Exon_Node>::iterator iter;
						while(j2<subexon_Num&&cur_Gene.subexon[j2].bound2<s2)
							j2++;
						Exon_Node cur_Exon;
						if(j2==subexon_Num)
						{
#pragma region exons overlap
							if(s2>cur_Gene.subexon[j2-1].bound2+2)
							{
								cur_Exon.bound1 = cur_Gene.subexon[j2-1].bound2+1;
								cur_Exon.bound2 = s2;
								cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
								cur_Gene.subexon.push_back(cur_Exon);
							}

							iter = cur_Gene.subexon.begin() + j2 - 1;
							if(cur_Gene.subexon[j2-1].bound1<s1-1&&cur_Gene.subexon[j2-1].bound2>s1+1)  // two exons overlap
							{
								int pre_bound1 = cur_Gene.subexon[j2-1].bound1;
								cur_Gene.subexon[j2-1].bound1 = s1;
								cur_Gene.subexon[j2-1].len = cur_Gene.subexon[j2-1].bound2 - s1 + 1; 
								cur_Exon.bound1 = pre_bound1; cur_Exon.bound2 = s1 - 1; 
								cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
								cur_Gene.subexon.insert(iter,cur_Exon);
							}
#pragma endregion
#pragma region
							else
							{
								int inter_Num = j2 - j1 - 1 ;
								int l1 = cur_Gene.subexon[j1].bound1, r1 = cur_Gene.subexon[j1].bound2;
#pragma region deal with the first exon
								iter = cur_Gene.subexon.begin() + j1;	
								int insert_itr1 = j1 + 1;

								if(s1<r1-1&&s1>l1+1) // the first exon need to be split
								{
									cur_Gene.subexon[j1].bound1 = s1; 
									cur_Gene.subexon[j1].len = cur_Gene.subexon[j1].bound2 - cur_Gene.subexon[j1].bound1 + 1;	
									cur_Exon.bound1 = l1; cur_Exon.bound2 = s1-1;
									cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
									cur_Gene.subexon.insert(iter,cur_Exon);
									insert_itr1++;
								}
								else if(s1<l1-1)     // insert new starting subexon
								{
									cur_Exon.bound1 = s1; 
									cur_Exon.bound2 = l1-1;
									cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
									cur_Gene.subexon.insert(iter,cur_Exon);
									insert_itr1++;
								}

#pragma endregion
								iter = cur_Gene.subexon.begin() + insert_itr1; 
								for(int k1 = 0; k1<inter_Num; k1++)
								{
									int pre_bound2 = cur_Gene.subexon[insert_itr1-1].bound2;
									cur_Exon.bound1 = pre_bound2 + 1;
									cur_Exon.bound2 = cur_Gene.subexon[insert_itr1].bound1-1;
									cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;															
									if(cur_Exon.len>1)
									{
										cur_Gene.subexon.insert(iter,cur_Exon);
										insert_itr1 += 2;
									}
									else
										insert_itr1++;
									// pre_bound2 = cur_Gene.subexon[insert_itr1-1].bound2; 
									iter = cur_Gene.subexon.begin() + insert_itr1;
									// pos++;
								}
							}
#pragma endregion 
						}
						else
						{
							int l1 = cur_Gene.subexon[j1].bound1, r1 = cur_Gene.subexon[j1].bound2;  // 首端外显子左、右边界
							int l2 = cur_Gene.subexon[j2].bound1, r2 = cur_Gene.subexon[j2].bound2;  // 末端外显子左、右边界
							vector<Exon_Node>::iterator iter_bak;
							int inter_Num = j2 - j1;
							int pre_bound2 = cur_Gene.subexon[j1].bound2;  // 外显子右边界
							Exon_Node cur_Exon;
							int insert_itr = j1 + 1; // 插入的内含子的索引
							int pos = j2;            // 记录末端外显子的索引	
							if(s1>l1+1&&s1<r1-1)     // 首端外显子需要被分开
							{
								cur_Gene.subexon[j1].bound1 = s1;             // 更新亚外显子边界
								cur_Gene.subexon[j1].len = cur_Gene.subexon[j1].bound2 - cur_Gene.subexon[j1].bound1 + 1;  // 更新亚外显子长度
								iter = cur_Gene.subexon.begin() + j1;					
								cur_Exon.bound1 = l1; cur_Exon.bound2 = s1-1; // 分裂出的外显子
								cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
								cur_Gene.subexon.insert(iter,cur_Exon);
								insert_itr++;
								pos++;
								j1++;
							}
							else if(s1<l1-1)    // 
							{	
								cur_Exon.bound1 = s1; 
								if(s2<l2-1)     // 在原内含子区域插入新的亚外显子
								{
									cur_Exon.bound2 = s2;
								}
								else
								{
									cur_Exon.bound2 = l1-1;
								}
								cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
								iter = cur_Gene.subexon.begin() + j1;
								cur_Gene.subexon.insert(iter,cur_Exon);
								insert_itr++;
								pos++;
							}
							if(s2>l2+1)  // the last exon need to be split
							{ 
#pragma region insert new subexon	
								iter = cur_Gene.subexon.begin() + insert_itr;
								for(int k1 = 0; k1<inter_Num; k1++)
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
										insert_itr++; 
									pre_bound2 = cur_Gene.subexon[insert_itr-1].bound2;  // right boundary of the exon
									iter = cur_Gene.subexon.begin() + insert_itr;
									// pos++;
								}

								if(s2<r2-1)  // the last exon is split into two exons
								{
									iter = cur_Gene.subexon.begin() + pos;
									cur_Exon.bound1 = cur_Gene.subexon[pos].bound1; cur_Exon.bound2 = s2;  // split exon
									cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;	
									cur_Gene.subexon.insert(iter,cur_Exon);  // split exon
									pos++;
									cur_Gene.subexon[pos].bound1 = s2+1;     // update the boundary of the new subexon
									cur_Gene.subexon[pos].len = cur_Gene.subexon[pos].bound2 - cur_Gene.subexon[pos].bound1 + 1;  // update the length of the subexon
								}
#pragma endregion
							}
							else if(l1!=l2)
							{
#pragma region insert new subexon
								iter = cur_Gene.subexon.begin() + pos;
								pos--;				
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
										insert_itr++;
									pre_bound2 = cur_Gene.subexon[insert_itr-1].bound2;
									iter = cur_Gene.subexon.begin() + insert_itr;
								}
								int r2_pre = cur_Gene.subexon[pos].bound2;
								if(s2> r2_pre + 2)
								{
									cur_Exon.bound1 = r2_pre + 1; cur_Exon.bound2 = s2; // split exon
									cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
									cur_Gene.subexon.insert(iter,cur_Exon);             // insert new subexon
								}						
#pragma endregion
							}
						}
					}
#pragma endregion
				}
				else
				{
				}
#pragma endregion

				line_Num++;
			}
			else if(ctr_Flag)
				break;
		}
	}

#pragma region Output the information of exons of genes
	ofstream outfile;
	outfile.open(file2,ios::out);
	int gene_Num = gene.size();
	int subexon_cnt[50];
	for(int i = 0; i<50; i++)
		subexon_cnt[i] = 0;
	int trs_cnt[50];
	for(int i = 0; i<50; i++)
		trs_cnt[i] = 0;

	for(int k = 0; k<gene_Num; k++)
	{
		int exon_Num1 = gene[k].subexon.size();
		int trs_Num = gene[k].trans.size();
		int exon_Num = exon_Num1>0?exon_Num1:1;
		outfile<<"@"<<k+1<<" "<<gene[k].gene_Name<<" "<<gene[k].start<<" "<<gene[k].stop<<" "<<exon_Num<<" "
			<<gene[k].trans.size()<<endl;

		subexon_cnt[exon_Num]++;
		trs_cnt[gene[k].trans.size()]++;

		for(int i = 0; i<trs_Num; i++)
		{
			outfile<<gene[k].trans[i].length<<" ";
		}
		outfile<<endl;

		for(int i = 0; i<exon_Num1; i++)
		{
			outfile<<i+1<<" "<<gene[k].subexon[i].bound1<<" "<<gene[k].subexon[i].bound2<<" "<<gene[k].subexon[i].len<<endl; // 输出外显子信息
		}
	}

	for(int i = 0; i<50; i++)
		outfile<<i<<":"<<subexon_cnt[i]<<" ";
	outfile<<endl;
	for(int i = 0; i<50; i++)
		outfile<<i<<":"<<trs_cnt[i]<<" ";
	outfile<<endl;

	outfile.close();

#pragma endregion
}

// Obtain information of exons and transcripts according to the annotations
void Rd_Set::load_Annotation(char* path, char* filename, char* ext)
{
	ifstream myfile;
	char file1[200];
	sprintf(file1,"%s\\%s.%s",path,filename,ext);

	myfile.open(file1,ios::in);
	if(!myfile)
	{
		cout<<"Open file error"<<endl;
		exit(-1);
	}

	int line_Number = 0, vline_Number = 0; 
	string word, line;
	int seg_Num = 6;   // number of fields
	vector<Gene> gene; // vector of genes
	Gene cur_Gene1;
	int gene_Serial = -1, pre_Serial = -1;   // gene serial
	Trans cur_Trans;
	char split = '\t';

	string seg[6];     // fields：1.start 2.stop 3.Accession 4.Oritation 5.Align quality 6.Description

	bool ctr_Flag = true;

	std::streampos pos1 = myfile.tellg();
	std::streampos pos2;

	while(std::getline(myfile,line))  // to reach the non-caption lines
	{			
		pos2 = myfile.tellg();
		cout<<line<<endl;
		if(line_Number>0&&line.length()>0&&line[0]!='#')
			break;	
		line_Number++;
		pos1 = myfile.tellg();
	}

	myfile.seekg(pos1);

#pragma region Load information of the genes
	while(!myfile.eof())
	{
		for(int k = 0; k<seg_Num; k++)
			myfile>>seg[k];

#pragma region read the field of Transcript
		if(seg[2]!="exon"&&seg[2]!="intron")
		{
			Trans cur_Trans;
			cur_Trans.trsName = seg[2];
			string gene_Name = seg[3];
			int found = isExisting(gene,gene_Name);  // check whether the gene already exists

			if(found>-1)      // the gene exists
			{
				gene_Serial = found;
			}
			else
			{
				gene_Serial = pre_Serial;
				cur_Gene1.gene_Name = gene_Name;				
				gene_Serial++;
				pre_Serial = gene_Serial;
				gene.push_back(cur_Gene1);
			}
			sscanf(seg[0].c_str(),"%d",&(cur_Trans.start));
			sscanf(seg[1].c_str(),"%d",&(cur_Trans.stop));

			gene[gene_Serial].trans.push_back(cur_Trans);

		}
#pragma endregion
#pragma region deal with fields of Exon and Intron
		else if(seg[2]=="intron")
		{
			int s1 = 0, s2 = 0;
		}
		else 
		{
			int s1 = 0, s2 = 0;
			sscanf(seg[0].c_str(),"%d",&s1);  // left boundary of exon
			sscanf(seg[1].c_str(),"%d",&s2);  // right boundary of exon
			Gene& cur_Gene = gene[gene_Serial];

			int subexon_Num = cur_Gene.subexon.size();

			int j1 = 0, j2 = 0;
			while(j1<subexon_Num&&cur_Gene.subexon[j1].bound2<s1)
				j1++;
			if(j1==subexon_Num)
			{
				Exon_Node subexon;
				subexon.bound1 = s1; subexon.bound2 = s2; subexon.len = s2 - s1 + 1;
				cur_Gene.subexon.push_back(subexon);   // 加入新的亚外显子				
			}
			else
			{
				vector<Exon_Node>::iterator iter;
				while(j2<subexon_Num&&cur_Gene.subexon[j2].bound2<s2)
					j2++;
				Exon_Node cur_Exon;
				if(j2==subexon_Num)
				{
#pragma region exons overlap
					if(cur_Gene.subexon[j2-1].bound1<s1)  // the twon exons overlap
					{
						int pre_bound2 = cur_Gene.subexon[j2-1].bound2;
						cur_Gene.subexon[j2-1].bound2 = s1;
						cur_Gene.subexon[j2-1].len = s1 - cur_Gene.subexon[j2-1].bound1 + 1;
						cur_Exon.bound1 = s1 + 1; cur_Exon.bound2 = pre_bound2;
						cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
						cur_Gene.subexon.push_back(cur_Exon);
						cur_Exon.bound1 = pre_bound2 + 1; cur_Exon.bound2 = s2;
						cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
						if(cur_Exon.len>2)
							cur_Gene.subexon.push_back(cur_Exon);
					}
#pragma endregion
#pragma region
					else
					{
						iter = cur_Gene.subexon.begin() + j2 - 2;
						if(s1<cur_Gene.subexon[j2-1].bound1-2)
						{
							cur_Exon.bound1 = s1;
							cur_Exon.bound2 = cur_Gene.subexon[j2-1].bound1-1;
							cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
							cur_Gene.subexon.insert(iter,cur_Exon);
						}

						if(s2>cur_Gene.subexon[j2-1].bound2+2)
						{
							cur_Exon.bound1 = cur_Gene.subexon[j2-1].bound2+1;
							cur_Exon.bound2 = s2;
							cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
							cur_Gene.subexon.push_back(cur_Exon);
						}
					}
#pragma endregion 
				}
				else
				{
					int l1 = cur_Gene.subexon[j1].bound1, r1 = cur_Gene.subexon[j1].bound2; 
					int l2 = cur_Gene.subexon[j2].bound1, r2 = cur_Gene.subexon[j2].bound2; 
					vector<Exon_Node>::iterator iter_bak;
					int inter_Num = j2 - j1;
					int pre_bound2 = cur_Gene.subexon[j1].bound2;
					Exon_Node cur_Exon;
					int insert_itr = j1 + 1;
					int pos = j2; 
					if(s1>l1+1&&s1<r1-1)   // the first exon need to be split
					{
						cur_Gene.subexon[j1].bound1 = s1;
						cur_Gene.subexon[j1].len = cur_Gene.subexon[j1].bound2 - cur_Gene.subexon[j1].bound1 + 1;
						iter = cur_Gene.subexon.begin() + j1;					
						cur_Exon.bound1 = l1; cur_Exon.bound2 = s1-1;
						cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
						cur_Gene.subexon.insert(iter,cur_Exon);
						insert_itr++;
						pos++;
						j1++;
					}
					else if(s1<l1-1)    // need to insert new subexon
					{	
						cur_Exon.bound1 = s1; 
						if(s2<l2-1)     // new subexon need to be inserted in the intron region
						{
							cur_Exon.bound2 = s2; 
						}
						else
						{
							cur_Exon.bound2 = l1-1;
						}
						cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
						iter = cur_Gene.subexon.begin() + j1; 
						cur_Gene.subexon.insert(iter,cur_Exon);
						insert_itr++;
						pos++;
					}
					if(s2>l2+1)
					{
#pragma region insert new subexon		
						iter = cur_Gene.subexon.begin() + insert_itr;
						for(int k1 = 0; k1<inter_Num; k1++)
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
								insert_itr++;
							pre_bound2 = cur_Gene.subexon[insert_itr-1].bound2;
							iter = cur_Gene.subexon.begin() + insert_itr;
						}

						if(s2<r2-1)
						{
							iter = cur_Gene.subexon.begin() + pos;
							cur_Exon.bound1 = cur_Gene.subexon[pos].bound1; cur_Exon.bound2 = s2;
							cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;	
							cur_Gene.subexon.insert(iter,cur_Exon);
							pos++;
							cur_Gene.subexon[pos].bound1 = s2+1;
							cur_Gene.subexon[pos].len = cur_Gene.subexon[pos].bound2 - cur_Gene.subexon[pos].bound1 + 1;  // 更新亚外显子长度
						}
#pragma endregion
					}
					else if(l1!=l2)
					{
#pragma region insert new subexon
						iter = cur_Gene.subexon.begin() + pos;
						pos--;				
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
								insert_itr++; 
							pre_bound2 = cur_Gene.subexon[insert_itr-1].bound2;
							iter = cur_Gene.subexon.begin() + insert_itr;
						}
						int r2_pre = cur_Gene.subexon[pos].bound2;
						if(s2> r2_pre + 2)
						{
							cur_Exon.bound1 = r2_pre + 1; cur_Exon.bound2 = s2;
							cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
							cur_Gene.subexon.insert(iter,cur_Exon);
						}						
#pragma endregion
					}
				}
			}
		}
#pragma endregion

		myfile.seekg(pos1);
		getline(myfile,line);
		cout<<line<<endl; 
		pos1 = myfile.tellg();
		line_Number++; 
	} // end while
#pragma endregion

	locate_transSub(gene,path,filename,ext);  // determine the structure of the transcript
}

// Determine the exons contained by each transcript
void Rd_Set::locate_transSub(vector<Gene>& gene, char* path, char* filename, char* ext)
{
	int line_Number = 0, vline_Number = 0;
	string word, line;
	int seg_Num = 6;
	Gene cur_Gene1;
	int gene_Serial = -1;
	Trans cur_Trans;
	int subexon_Num = 0;
	string seg[6];     // fields：1.start 2.stop 3.Accession 4.Oritation 5.Align quality 6.Description
	bool ctr_Flag = true;
	ifstream myfile;

	char file1[200];
	sprintf(file1,"%s\\%s.%s",path,filename,ext);

	myfile.open(file1);
	if(!file1)
	{
		cout<<"Open file error"<<endl;
		exit(-1);
	}

	std::streampos pos1 = myfile.tellg();
	std::streampos pos2;

	while(std::getline(myfile,line))
	{			
		pos2 = myfile.tellg();
		cout<<line<<endl;
		if(line_Number>0&&line.length()>0&&line[0]!='#')
			break;	
		line_Number++;
		pos1 = myfile.tellg();
	}

	myfile.seekg(pos1);

	int gene_Idx = -1;
	int trs_Idx = -1;

#pragma region read gene data
	while(!myfile.eof())
	{
		for(int k = 0; k<seg_Num; k++)
			myfile>>seg[k];

		if(seg[2]!="exon"&&seg[2]!="intron") 
		{
			Trans cur_Trans;           
			cur_Trans.trsName = seg[2];
			string gene_Name = seg[3];
			gene_Idx = isExisting(gene,gene_Name);

			if(gene_Idx>-1)   // the gene already exists
			{
				gene_Serial = gene_Idx;
				cout<<"Gene_Serial: "<<gene_Serial<<endl;
				trs_Idx = indexTrans(gene[gene_Serial],seg[2]);  // index of the transcript
				subexon_Num = gene[gene_Serial].subexon.size();  // number of subexons
			}
		}
		else if(gene_Serial>-1&&trs_Idx>-1&&seg[2]=="exon")
		{
#pragma region read the fields of intron and exon
			int s1 = 0, s2 = 0;
			sscanf(seg[0].c_str(),"%d",&s1);
			sscanf(seg[1].c_str(),"%d",&s2);

			int i = 0;
			bool cds_Flag = false;
			if(seg[4]=="CDS")
				cds_Flag = true;

			int subexon_Num = gene[gene_Serial].subexon.size();

			while(i<subexon_Num&&abs(gene[gene_Serial].subexon[i].bound1-s1)>2)
				i++;

			for(int k = i; k<subexon_Num; k++)
			{
				if(gene[gene_Serial].subexon[k].bound2>s2+2)
					break;
				gene[gene_Serial].trans[trs_Idx].subexon_Serial.push_back(k);
				gene[gene_Serial].trans[trs_Idx].cds.push_back(cds_Flag);
			}

#pragma endregion
		}

		myfile.seekg(pos1); 
		getline(myfile,line);
		cout<<line<<endl;   
		pos1 = myfile.tellg();
		line_Number++; 

	} // end while
#pragma endregion

}

// Downsampling of sequencing data
void Rd_Set::subSampling(char* filename1, char* filename2, int subFactor)
{
	ifstream infile;
	infile.open(filename1,ios::in);
	if(!infile)
		cout<<"Open file error!"<<endl;
	ofstream outfile;
	outfile.open(filename2,ios::out);
	if(!outfile)
		cout<<"Open file error!"<<endl;

	string line, word;
	int line_Cnt = 0;
	while(std::getline(infile,line))
	{
		stringstream stream(line);
		if(line_Cnt%5==0)
			outfile<<line<<endl;
		line_Cnt++;
	}
	infile.close();
	outfile.close();
}

// Merge the sequencing data
void Rd_Set::merge(char* root_path, char* filename1, vector<string>& gene_Name)
{
	int gene_Num = gene_Name.size();
	char filename[200];
	string gene_name,line;
	ofstream outfile;
	outfile.open(filename1,ios::out);   // output to the file
	if(!outfile) cout<<"Open file error!"<<endl;

	for(int k = 0; k<gene_Num; k++)
	{
		gene_name = gene_Name[k];
		cout<<gene_name<<endl;
		sprintf(filename,"%s\\trs_num3\\%s\\read_2.fa",root_path,gene_name.c_str());
		ifstream infile;
		infile.open(filename,ios::in);
		if(!infile) cout<<"Open file error!"<<endl;
		while(getline(infile,line))
			outfile<<line<<endl;
		infile.close();
	}

	outfile.close();
}

// Separted original sequencing data into the form of paired-end reads
void Rd_Set::separate_pairedData(char* filename1, char* filename2, char* filename3)
{
	ifstream infile;
	infile.open(filename1,ios::in);      // file of original sequencing data
	if(!infile)
		cout<<"Open file error!"<<endl;
	ofstream outfile1, outfile2;
	outfile1.open(filename2,ios::out);   // output file directory for the separated data of left reads
	if(!outfile1)
		cout<<"Open file error!"<<endl;
	outfile2.open(filename3,ios::out);   // output file directory for the separted data of right reads
	if(!outfile2)
		cout<<"Open file error!"<<endl;

	int line_cnt = 0;
	int rd_len = 76;
	string line1,word,seg[3],line[4];
	while(!infile.eof())
	{
		for(int k=0;k<4;k++)
		{
			getline(infile,line[k]);
		}
		stringstream stream(line[0]);
		for(int i=0; i<3; i++)
			stream>>seg[i];
		string temp1 = seg[0]+".1";
		string temp2 = seg[0]+".2";
		temp1 = temp1+" "+seg[1]+" "+seg[2];
		temp2 = temp2+" "+seg[1]+" "+seg[2];
		outfile1<<temp1<<endl;
		outfile2<<temp2<<endl;

		temp1 = line[1]+line[2]+line[3];
		line1 = temp1.substr(0,rd_len);
		outfile1<<line1<<endl;
		line1 = temp1.substr(rd_len,rd_len);
		outfile2<<line1<<endl;
		line_cnt++;
		if(line_cnt%5000==0)
			cout<<line_cnt<<endl;
	}

	infile.close();
	outfile1.close();
	outfile2.close();
}

// Output the gene annotations to files to provide information of exons
void Rd_Set::output_GeneAnnotation(char* filename, vector<Gene>& gene)
{
	ofstream outfile;
	outfile.open(filename,ios::out);
	if(!outfile)
		cout<<"Open file(output file) error!"<<endl;
	int gene_Num = gene.size();

	for(int k = 0; k<gene_Num; k++)
	{
		int exon_Num1 = gene[k].subexon.size();  // number of exons
		int trs_Num = gene[k].trans.size();      // number of transcripts
		int exon_Num = exon_Num1;
		outfile<<"@"<<k+1<<" "<<gene[k].gene_Name<<" "<<gene[k].start<<" "<<gene[k].stop<<" "<<exon_Num<<" "
			<<trs_Num<<endl;   // output gene serial, gene name, starting site, number of subexons and number of transcripts

		for(int i = 0; i<exon_Num; i++)
		{
			gene[k].subexon[i].len = gene[k].subexon[i].bound2 - gene[k].subexon[i].bound1 + 1; // length of subexon
			outfile<<i+1<<" "<<gene[k].subexon[i].bound1<<" "<<gene[k].subexon[i].bound2<<" "<<gene[k].subexon[i].len<<endl; // output information of subexons
		}
	}
	outfile.close();
}

// Merge the data of paired-end reads into one file
void Rd_Set::Merge_pairedData(char* filename1, char* filename2, char* filename3)
{
	ifstream infile1, infile2;
	ofstream outfile;
	infile1.open(filename1,ios::in);
	if(!infile1)
		cout<<"Open file (file1) error!"<<endl;
	infile2.open(filename2,ios::in);
	if(!infile2)
		cout<<"Open file (file2) error!"<<endl;
	outfile.open(filename3,ios::out);
	if(!outfile)
		cout<<"Open file (outfile) error!"<<endl;

	string line1, line2, line;
	std::streampos pos;
	char mark1 = '>', mark2 = '@';

	int num = 1;
	char temp;
	while(getline(infile1,line1))
	{
		stringstream stream(line1);
		outfile<<line1<<endl;
		pos = infile1.tellg();
		while(getline(infile1,line1)&&line1[0]!=mark1&&line1[0]!=mark2)
		{
			outfile<<line1<<endl;
			pos = infile1.tellg();
			num++;
		}
		for(int k=0; k<num; k++)
		{
			getline(infile2,line2);
			outfile<<line2<<endl;
		}
		if(!infile1.eof()) 
		{
			line = line1;
			break;
		}
	}

	outfile<<line1<<endl;
	for(int k=1; k<num; k++){
		getline(infile1,line1);
		outfile<<line1<<endl;
	}
	for(int k=0; k<num; k++){
		getline(infile2,line2);
		outfile<<line2<<endl;
	}
	while(!infile1.eof())
	{
		for(int k=0; k<num; k++){
			getline(infile1,line1);
			outfile<<line1<<endl;
		}
		for(int k=0; k<num; k++){
			getline(infile2,line2);
			outfile<<line2<<endl;
		}
	}

	infile1.close();
	infile2.close();
	outfile.close();
}

// Extract the read data in the specific gene
void Rd_Set::extract_pairedData(char* filename1, char* filename2, int site1, int site2, string chromo_Name)
{
	ifstream infile;
	infile.open(filename1,ios::in);     // file of read data
	if(!infile)
		cout<<"Open file error!"<<endl;

	ofstream outfile;
	outfile.open(filename2,ios::out);   // output file directory
	if(!outfile)
		cout<<"Open file error!"<<endl;

	string line, word, frag_Name, chromosome_Name;
	string words[6];
	int x1 = 0, x2 = 0, z1 = 0, z2 = 0, temp_Len = 0; // starting sites of left read and right read; fragment length
	int line_Cnt1 = 0;
	bool locate_Flag = false;
	int valid_Cnt = 0;

	while(std::getline(infile,line))  // test a given number of genes
	{    
		Rd cur_Rd;  // read
		stringstream stream(line);
		stream>>word; frag_Name = word;
		stream>>word; stream>>word; // read the first three fields
		chromosome_Name = word;     // chromosome name

		if(chromosome_Name!=chromo_Name)  // the specific chromosome not reached
		{
			continue;
		}

		for(int i = 0; i<6; i++)
			stream>>words[i];  // 4 starting site of left read 5 match quality 6 match quality 7 marker for paired match 8 starting site of right read 9 inserted size
		sscanf(words[5].c_str(),"%d",&temp_Len);  // inserted size
		sscanf(words[0].c_str(),"%d",&x1);  // starting site of the left read
		sscanf(words[4].c_str(),"%d",&z1);  // stopping site of the right read
		z2 = 0;
		if(temp_Len>=0)
			z2 = x1 + temp_Len - 1;  // right end of the right read

		int s1 = _min(x1,z1), s2 = _max(_max(x1,z1),z2);		

		if(s1>50000)
		{
			line_Cnt1++;
			if(line_Cnt1%100000==0)
				cout<<s1<<" "<<s2<<endl;
		}

		if(s1==0&&s2>=site1&&s2<=site2)
			locate_Flag = true;
		else if(s2==0&&s1>=site1&&s1<=site2)
			locate_Flag = true;
		else if(s1>=site1&&s2<=site2)
			locate_Flag = true;
		else if((s1>=site2))
		{
			cout<<line<<endl;
			cout<<valid_Cnt<<endl;
			break;
		}

		if(locate_Flag==true)
		{
			valid_Cnt++;
			outfile<<line<<endl;			
		}
	}
	for(int i=0; i<1000; i++)
	{
		getline(infile,line);
		outfile<<line<<endl;
	}
	cout<<"Extracted"<<endl;
	outfile.close();
}

// Load transcript information
void Rd_Set::load_Trans(char* filename)
{
	ifstream myfile;
	myfile.open(filename,ios::in);
	if(!myfile)
	{
		cout<<"can not open the file."<<endl;
		exit(-1);
	}

	string word, line;
	int e_Num = 0, temp = 0;
	double abun = 0;
	while(std::getline(myfile,line))
	{
		Path cur_Path;
		stringstream stream(line);
		stream>>word;  // transcript serial
		stream>>word;  // transcript abundance
		istringstream ss(word);
		ss>>abun;
		cur_Path.entro_Abun1 = abun;
		stream>>word;  // read the transcript length
		sscanf(word.c_str(),"%d",&(cur_Path.length));
		stream>>word;  // read the number of exons
		sscanf(word.c_str(),"%d",&e_Num);
		cur_Path.exon_Number = e_Num;
		cur_Path.exon_Node = new int[e_Num];
		for(int k = 0; k<e_Num; k++)
		{
			stream>>word;
			sscanf(word.c_str(),"%d",&temp);
			cur_Path.exon_Node[k] = temp;
		}
		for(int i = 0; i<cur_Path.exon_Number; i++)
			cout<<cur_Path.exon_Node[i]<<" ";
		cout<<endl;
		path.push_back(cur_Path);
	}
}

// Print transcript information
void Rd_Set::print_Trans()
{
	int p_Num = path.size();
	for(int k = 0; k<p_Num; k++)
	{
		cout<<k<<" "<<path[k].entro_Abun1<<" "<<path[k].length<<" "<<path[k].exon_Number<<" ";
		for(int i = 0; i<path[k].exon_Number; i++)
			cout<<path[k].exon_Node[i]<<" ";
		cout<<endl;
	}
}

// Output the transcript information to the file
void Rd_Set::output_Trans(char* filename)
{
	ofstream outfile;
	outfile.open(filename,ios::app);
	if(!outfile)
	{
		std::cout<<"can not open the file."<<endl;
		exit(-1);
	}
	int num = path.size();

	for(int i = 0; i<num; i++)
	{
		int e_Num = path[i].exon_Number;
		outfile<<i<<" "<<path[i].entro_Abun1<<" "<<path[i].length<<" "<<e_Num<<" ";		
		for(int k = 0; k<e_Num; k++)
			outfile<<path[i].exon_Node[k]<<" ";
		outfile<<endl;
	}

	outfile.close();
}

// Load read information
void Rd_Set::load_Read(char* filename)
{
	ifstream myfile;
	myfile.open(filename,ios::in);
	if(!myfile)
	{
		cout<<"can not open the file."<<endl;
		exit(-1);
	}

	string word, line;
	int temp = 0;
	while(std::getline(myfile,line))
	{
		Rd cur_Rd;
		stringstream stream(line);
		stream>>word;
		stream>>word;  // transcript length
		// load the left boundary and the right boundary
		stream>>word; sscanf(word.c_str(),"%d",&(cur_Rd.Left1));
		stream>>word; stream>>word;
		stream>>word; sscanf(word.c_str(),"%d",&(cur_Rd.Right2));
		// load the indices of exons in which the boundaries of the read fall
		stream>>word; sscanf(word.c_str(),"%d",&temp); cur_Rd.l1 = temp - 1;
		stream>>word; sscanf(word.c_str(),"%d",&temp); cur_Rd.l2 = temp - 1;
		stream>>word; sscanf(word.c_str(),"%d",&temp); cur_Rd.r1 = temp - 1;
		stream>>word; sscanf(word.c_str(),"%d",&temp); cur_Rd.r2 = temp - 1;
		rd_Vec.push_back(cur_Rd);
	}

	rd_Num = rd_Vec.size();  // number of reads
}

// Output the read information to the file
void Rd_Set::output_Read(char* filename)
{
	ofstream outfile;
	outfile.open(filename,ios::out);

	if(!outfile)
	{
		cout<<"can not open the file."<<endl;
		exit(-1);
	}
	int num = rd_Vec.size();

	for(int i = 0; i<num; i++)
	{
		outfile<<rd_Vec[i].name<<" "<<rd_Vec[i].Left1<<" "<<rd_Vec[i].Left2<<" "<<rd_Vec[i].Right1<<" "<<rd_Vec[i].Right2<<" ";
		int n1 = rd_Vec[i].left_readIdx.size();   // the segment in which the boundary of the left read falls
		int n2 = rd_Vec[i].right_readIdx.size();  // the segment in which the boundary of the right read falls
		if(n1>0)
		{
			for(int k = 0; k<n1; k++)
				outfile<<rd_Vec[i].left_readIdx[k]+1<<" ";
			outfile<<"\t";
		}
		else
			outfile<<(int)rd_Vec[i].l1+1<<" "<<(int)rd_Vec[i].l2+1<<"\t";
		if(n2>0)
		{
			for(int k = 0; k<n2; k++)
				outfile<<rd_Vec[i].right_readIdx[k]+1<<" ";
			outfile<<endl;
		}
		else
			outfile<<(int)rd_Vec[i].r1+1<<" "<<(int)rd_Vec[i].r2+1<<endl;
	}

	outfile.close();
}

