
/************************************************************************************
 * Assembly.cpp
 * This file is part of the software program of MaxInfo
 * It can be redistributed and/or modified for uncommercial purposes
 *
 * First created on: 2014-06-27
 * Last modified on: 2015-02-11
 * Author: Yang Yang
 ***********************************************************************************/

//////////////////////////////////////////////////////////////////////////
//// TRANSCRIPT ASSEMBLY AND ABUNDANCE ESTIMATION BASED ON RNA-SEQ DATA

#include "Entropy.h"
#include "Batch_Mode.h"
#include "ISA.h"

// abundance estimation
// Input: file of genes for abundance estimation; file of read alignment
void Gene_Batch::abun_Estimation(char* filename1, char* filename2, char* data_path, string chromo_Name)
{
	char path[200], data_path1[200];
	char filename[200];
	ifstream infile1;
	infile1.open(filename1,ios::in);   // load the genes
	if(!infile1) cout<<"Open file(annotation file) error!"<<endl;

	ifstream infile2;
	infile2.open(filename2,ios::in);   // load the read alignment data for the genes
	if(!infile2) cout<<"Open file(assignment file) error!"<<endl;

	sprintf(data_path1,"%s",data_path);
	int serial = -1;   // gene serial
	int Number = 0;
	string line, word, seg[20]; 
	string gene_Name, gene_Name1; // gene name
	bool flag_Ctrl = false;
	string gene_NamePre = "";
	std::streampos pos1 = infile1.tellg();  // obtain the current file pointer

	// gene serial
	while(std::getline(infile1,line))  // test a certain number of genes
	{
		stringstream stream(line);
		int i = 0;
		stream>>seg[i];
		while(seg[i]!="gene_name"&&seg[i]!="gene_id"&&seg[i]!="Gene_Id"&&seg[i]!="Gene_id")
		{
			i++;
			stream>>seg[i];
		}

		if(seg[i]=="gene_name"||seg[i]=="gene_id"||seg[i]=="Gene_Id"||seg[i]=="Gene_id")  // search the field "gene_id"
		{
			stream>>gene_Name;
		}

		// gene_Name = valid_filename(gene_Name1);   // valid filename
		if(gene_Name!=gene_NamePre)
		{
			Graph_Trans gph;     // directed graph
			Rd_Set rd_Set;
			serial++;
			infile1.seekg(pos1); // file pointer back to the original position
			gene_NamePre = gene_Name;  // gene name

			ini_Graph(gph, gene_Name, serial, infile1);  // initialize the graph structure
			assignment_Pose(rd_Set, gph, infile2, chromo_Name); // deal with the data of read alignment
			rd_Set.gene_Name = gene_Name;
			rd_Set.rd_Len = rd_Len;

			abun_EstiSingle(rd_Set, gph, data_path1);
			gph.clear();
			rd_Set.~Rd_Set();
		}

		pos1 = infile1.tellg();     // obtain the current file pointer
	}
}

// Abundance estimation: assemble exons, and conduct isoform identification and abundance estimation based on the data of read alignment 
void Gene_Batch::abun_EstimationSpr(char* filename1, char* filename2, char* data_path, string chromo_Name)
{
	char data_path1[200];
	sprintf(data_path1,"%s",data_path);
	
#pragma region Declaration of Variables
	char path[200];
	char filename[200];

	ifstream infile2;
	infile2.open(filename2,ios::in);  // load the data of read alignment of the genes for abundance estimation
	if(!infile2) cout<<"Open file(assignment file) error!"<<endl;

	int Number = 0;
	string line, word, seg[20], chromo_name;
	string gene_Name, gene_Name1; // gene name
	bool flag_Ctrl = false;
	int serial = -1;   // gene serial
	std::streampos pos1 = infile2.tellg();  // obtain the current file pointer

	//////////////////////////////////////////////////////////////////////////
	// located to the specific chromosome
	while(getline(infile2,line))
	{
		stringstream stream(line);
		stream>>word;
		stream>>word;  // field related to matching algorithms
		stream>>chromo_name;
		if(chromo_name==chromo_Name)
			break;		
		pos1 = infile2.tellg();
	}

	// cout<<line<<endl;
	infile2.seekg(pos1);
#pragma endregion

	vector<Gene> Gene_Set;   // set of genes
	ini_GraphSpr(Gene_Set, filename1, chromo_Name);  // load gene structure according to annotations

	int gene_Num = Gene_Set.size();  // gene number	
	gene_Num = 20;
	for(int i=0; i<gene_Num; i++)
	{
		gene_Name = Gene_Set[i].gene_Name;  // gene name
		cout<<"No."<<i+1<<" "<<gene_Name<<endl;
		Graph_Trans gph;
		Rd_Set rd_Set;
		gph.exon_Vec = Gene_Set[i].subexon; // exon
		int exon_Num = gph.exon_Vec.size(); // number of exons
		string gene_name = Gene_Set[i].gene_Name;  // gene name or gene number
		gph.graph_Initial(exon_Num);
		rd_Set.gene_Name = gene_name;  // gene name
		rd_Set.rd_Len = rd_Len;        // read length
		// load_assignment deal with read alignment data
		bool locate_Flag = assignment_Pose_1(Gene_Set[i], rd_Set, gph, infile2, chromo_Name);

		if(locate_Flag==true)  // sequencing data located in the gene
		{
			abun_EstiSingle(rd_Set, gph, data_path1);
			gph.clear();
			rd_Set.~Rd_Set();
		}
	}
}

// Isoform identification and abundance estimation in gene loci
bool Gene_Batch::abun_EstiSingle(Rd_Set& rd_Set, Graph_Trans& gph, char* data_path, int method_mode)
{
	char filename[200];
	ofstream outfile_record;  // file recording the running status each time
	/*sprintf(filename,"%s\\run_record.txt",data_path);
	outfile_record.open(filename,ios::app);*/

#pragma region obtain path information from graph structure
	// gph.graph_Modify();   // modify the graph
	// gph.print_Graph();    // output the graph structure
	bool flag = false;
	if(method_mode==0)
		gph.graph_Constraint_Initial();    // initial constraints for the linkage in the graph
	else
		gph.graph_Constraint_Annotated();  // initial constraints for the linkage in the graph

	if(gph.exon_Vec.size()<=20)
	{
		flag = gph.path_DFS1();
	}
	if(flag==false)
	{
		// gph.graph_Constraint();   // enhance the constraints for the first node and the last node
		// gph.print_Graph();        // output the graph structure
		// gph.path_DFS1_constraint();
		gph.path_DFS1_Adaptive();
	}

	// cout<<"Path searched"<<endl;
	gph.junc_Index();
	gph.cost_FirstLastNode();

	// mark the information of the isoforms
	gph.path_Speci(rd_Set.path);
	// rd_Set.print_Trans();

	if(rd_Set.path.size()>5000)
	{
		/*cout<<"Too many paths: "<<rd_Set.path.size()<<endl;
		int exon_Num = gph.exon_Vec.size();
		outfile_record<<"0	"<<rd_Set.path.size()<<" "<<gph.exon_Vec[0].bound1<<" "<<gph.exon_Vec[exon_Num-1].bound2<<endl;
		outfile_record.close();*/
		return false;
	}
#pragma endregion

#pragma region Calculate the entropy of the sequencing data
	int trans_Num = rd_Set.path.size();  // number of transcripts
	rd_Set.initial_RdPro(trans_Num);    // array initialization

	double mean = rd_Set.mean;   // mean of fragment length
	double sigma = rd_Set.sigma; // standard deviation of fragment length
	
	// estimate the probability distribution of fragment length based on discrete intervals
	rd_Set.fragLen_DisInterval(mean, sigma, Lim1, Lim2);

	sprintf(filename,"%s\\%s_a.txt",data_path,rd_Set.gene_Name.c_str());
	rd_Set.path_Entro1(gph,filename);   // calculate the probability and output it to the file
	int rd_num = rd_Set.rd_Vec.size();  // number of valid reads
	int dataDim = rd_Set.idx.size();    // number of valid transcripts

	if(dataDim>5000||dataDim==0)
	{
		/*cout<<"Too many paths: "<<dataDim<<endl;
		int exon_Num = gph.exon_Vec.size();
		outfile_record<<"0	"<<dataDim<<" "<<gph.exon_Vec[0].bound1<<" "<<gph.exon_Vec[exon_Num-1].bound2<<endl;
		outfile_record.close();*/
		return false;
	}
	rd_Set.cost_Linkage(gph);      // calcuate the connectivity cost of the path 
#pragma endregion

#pragma region Isoform identification and abundance estimatino based on gene structures and read alignment of sequencing data

	if(dataDim>0)
	{
		int K = 1;
		float sigma_1 = 5;
		int numSample = 20;
		int iteration = 50;

		int initSampleNum = numSample;

		if(dataDim==1)
		{
			numSample = 1; initSampleNum = 1; iteration = 1;
		}
		for(int t = 0; t<K; t++)
		{
			ISA* pISA = new ISA(dataDim,sigma_1,numSample,iteration);
			float* ini_Sample = new float[initSampleNum*dataDim];
			pISA->initial_Sampling(rd_Set,ini_Sample,initSampleNum);
			pISA->setMode(method_mode);
			pISA->Optimize_IteCtr(ini_Sample,rd_Set,gph,data_path,initSampleNum,iteration);
			pISA->~ISA();

			delete[] ini_Sample;
		}
	}
#pragma endregion

	rd_Set.clear_RdPro(rd_num);
	return true;
}

// Transcript assembly and abundance estimation based on sequencing data and read alignment with annotations
// Input: filename1: file of genome annotations; filename2: file of read alignment; data_path: file directory of output files
void Gene_Batch::abun_EstimationSpr_1(char* filename1, char* filename2, char* data_path)
{
	char data_path1[200];
	sprintf(data_path1,"%s",data_path);

#pragma region Declaration of Variables
	char path[200];
	char filename[200];

	ofstream outfile;
	sprintf(filename,"%s\\exply.txt",data_path1);
	outfile.open(filename,ios::out);
	outfile.close();
	sprintf(filename,"%s\\record-a.txt",data_path1);
	record_file.open(filename,ios::out);

	ifstream infile2;
	infile2.open(filename2,ios::in);  // load read alignment results
	if(!infile2) cout<<"Open file(assignment file) error!"<<endl;

	//sprintf(filename,"%s\\gene_Locate.txt",data_path1);
	//outfile.open(filename,ios::out);  // output the gene location results

	vector<Gene> Gene_Set; // set of genes
	int gene_Num = 0, trsnum_Ctrl = 0, read_Cnt = 0;

	int Number = 0;
	string line, word, seg[20], chromo_Name;	
	string gene_Name, gene_Name1;
	bool flag_Ctrl = false;
	int serial = -1;   // gene serial
	std::streampos pos1 = infile2.tellg();     // obtain the current file pointer

#pragma endregion

	//////////////////////////////////////////////////////////////////////////
	rdLen_estimation(filename2);  // query read length

	//////////////////////////////////////////////////////////////////////////
	// locate specific chromosome
	// initialize the mean and variance of fragment length
	fragLen_Mean = DEFAULT_MEAN;
	fragLen_Var = DEFAULT_VAR;
	read_Number = 0;

	getline(infile2,line);
	stringstream stream(line);
	stream>>word;
	stream>>word;  // field related to matching algorithms
	stream>>chromo_Name;
	chromosome_Name = chromo_Name;
	bool changed_Flag = true;

	while(getline(infile2,line))
	{
		int trsnum_Ctrl = 0;
		if(changed_Flag==true)
		{
			Gene_Set.clear();    // clear the arrays
			load_Annotation(Gene_Set, filename1, chromo_Name);  // load gene structures based on the annotations
			changed_Flag = false;
			gene_Num = Gene_Set.size();  // number of genes
		}
		if(changed_Flag==false)  // gene annotations are loaded
		{
#pragma region Isoform discovery and abundance estimation
			int gene_Num = Gene_Set.size();  // number of genes
			int num = 0, gene_Cnt = 0, i = 0;
			bool ctr_Flag = false;
			for(i=0; i<gene_Num; i++)
			{
				gene_Name = Gene_Set[i].gene_Name;  // gene name or gene number
				cout<<"No."<<i+1<<" "<<gene_Name<<endl;
				Graph_Trans gph;
				Rd_Set rd_Set;
				rd_Set.rd_Len = rd_Len;
				rd_Set.chromosome_Name = chromo_Name;
				gph.exon_Vec = Gene_Set[i].subexon; // vector of exons
				int exon_Num = gph.exon_Vec.size(); // number of exons
				int len = gph.exon_Vec[exon_Num-1].bound2 - gph.exon_Vec[0].bound1 + 1;  // gene length

				cout<<Gene_Set[i].subexon[0].bound1<<" "<<Gene_Set[i].subexon[exon_Num-1].bound2<<" "<<len<<endl;
				gene_Name = gene_Name.substr(1,gene_Name.length()-3);  // extract gene name
				rd_Set.gene_Name = gene_Name;
				rd_Set.orientation = Gene_Set[i].orientation;
				int locate_Flag = locus_Assemble(infile2, rd_Set, gph, Gene_Set[i].subexon); 

				if(locate_Flag==1)  // sequencing data is located in the current gene
				{
					gene_Cnt++;
					cur_Gene = Gene_Set[i];
					//outfile<<gene_Cnt<<"\t"<<i+1<<"\t"<<gene_Name<<"\t"<<rd_Set.rd_Vec.size()<<endl;
					locate_transcriptionSites(Gene_Set[i], transcription_sites);

					// estimate the mean and variance of fragment length
					rd_Set.fragLen_Estimation(read_Number,fragLen_Mean,fragLen_Var);
					fragLen_Mean = rd_Set.mean;
					fragLen_Var = rd_Set.sigma;
					record_file<<rd_Set.rd_Vec.size()<<"\t"<<fragLen_Mean<<endl;
					num++;
					record_file<<num<<" "<<gene_Name<<endl;
					ctr_Flag = check_Valid(rd_Set, gph);

					if(ctr_Flag==true)
					{
						cout<<"read number: "<<rd_Set.rd_Vec.size()<<endl;  // output number of reads	
						read_Cnt += rd_Set.rd_Vec.size();
						gph.cur_Gene = cur_Gene;
						bool est_Flag = abun_EstiSingle(rd_Set, gph, data_path1,2);
						if(est_Flag==true)
							cout<<"True"<<endl;
						else
							cout<<"Too many dataDim"<<endl;	
						gph.clear();
					}
#pragma region clear the arrays
					int n1 = rd_Left.size();
					for(int j=0; j<n1; j++)
						vector<int>().swap(rd_Left[j]);
					int n2 = rd_Right.size();
					for(int j=0; j<n2; j++)
						vector<int>().swap(rd_Right[j]);
					vector<vector<int>>().swap(rd_Left);
					vector<vector<int>>().swap(rd_Right);
					rd_Left.clear();
					rd_Right.clear();
#pragma endregion
				}
				else if(locate_Flag==-1&&next_chromoName!="")  // there is change of chromosome
				{
					chromosome_Name = next_chromoName;
					changed_Flag = true;
					next_chromoName = "";
					rd_Set.~Rd_Set();
					break;
				}
				else    // there is no sequencing data located in the current gene
				{
					infile2.seekg(pos1);  // located to the original position
				}

				rd_Set.~Rd_Set();
			}
			if(changed_Flag==false&&i==gene_Num)
			{
				while(getline(infile2,line))  // test a certain number of genes
				{
					stringstream stream(line);
					stream>>word;
					stream>>word; stream>>word; // read and skip the first three fields
					string chr_Name = word;     // chromosome name
					if(chr_Name!=chromosome_Name)    // there is change of chromosome name
					{
						chromosome_Name = chr_Name;  // located to the next chromosome
						changed_Flag = true;
						next_chromoName = "";
						break;
					}
				}
			}
#pragma endregion
		}
	}

	// convert the temporary file to standard file
	int total_Number = read_Cnt;
	char filename_1[200];
	sprintf(filename,"%s\\exply.txt",data_path);
	sprintf(filename_1,"%s\\results.gtf",data_path);
	expression_conver(filename, filename_1, total_Number);
	infile2.close();
	//outfile.close();
	record_file.close();
}

// Transcript assembly and abundance estimation based on sequencing data and read alignment with annotations
// Input: filename1: file of genome annotations; filename_2: file of junctions; filename2: file of read alignment; data_path: file directory of output files
void Gene_Batch::abun_EstimationSpr_2(char* filename1, char* filename_2, char* filename2, char* data_path)
{
	char data_path1[200];
	sprintf(data_path1,"%s",data_path);
	
#pragma region Declaration of Variables
	char path[200];
	char filename[200];

	sprintf(filename,"%s\\exply.txt",data_path1);
	ofstream outfile;
	outfile.open(filename,ios::out);
	outfile.close();
	sprintf(filename,"%s\\record-a.txt",data_path1);
	record_file.open(filename,ios::out);

	ifstream infile1;
	infile1.open(filename_2,ios::in);  // load the identified junctions of read alignment results
	if(!infile1) cout<<"Open file(assignment file) error!"<<endl;

	ifstream infile2;
	infile2.open(filename2,ios::in);   // load read alignment results
	if(!infile2) cout<<"Open file(assignment file) error!"<<endl;

	//sprintf(filename,"%s\\gene_Locate.txt",data_path1);
	//outfile.open(filename,ios::out);  // output the gene location results

	method_Mode = 1;   // set to the annotation-dependent mode

	int Number = 0, gene_Num = 0, read_Cnt = 0;
	string line, word, seg[20], chromo_name;
	string gene_Name, gene_Name1; // gene name
	bool flag_Ctrl = false;
	int serial = -1;   // gene serial 
	std::streampos pos1 = infile2.tellg();     // obtain the current file pointer

	//////////////////////////////////////////////////////////////////////////
	rdLen_estimation(filename2);  // query read length

	// located to the specific chromosome
	// initialize the mean and variance of fragment length
	fragLen_Mean = DEFAULT_MEAN;
	fragLen_Var = DEFAULT_VAR;
	read_Number = 0;

	getline(infile2,line);
	stringstream stream(line);
	stream>>word;
	stream>>word;  // field related to matching algorithms
	stream>>chromo_name;
	chromosome_Name = chromo_name;
	bool changed_Flag = true;

#pragma endregion

	vector<Gene> Gene_Set;   // set of genes
	while(getline(infile2,line))
	{
		int trsnum_Ctrl = 0;
		if(changed_Flag==true)
		{
			Gene_Set.clear();    // clear the arrays
			load_Annotation(Gene_Set, filename1, chromosome_Name);
			assembly_Annotation(infile1, Gene_Set);
			cout<<"Assembled!"<<endl;
			changed_Flag = false;
			gene_Num = Gene_Set.size();  // number of genes
		}
		if(changed_Flag==false)  // gene annotations are loaded
		{
#pragma region isoform identification and abundance estimation
			int num = 0, gene_Cnt = 0;
			bool ctr_Flag = false;
			int i = 0;
			for(i=0; i<gene_Num; i++)
			{
				gene_Name = Gene_Set[i].gene_Name;  // gene name or gene number
				cout<<"No."<<i+1<<" "<<gene_Name<<endl;
				Graph_Trans gph;
				Rd_Set rd_Set;
				rd_Set.rd_Len = rd_Len;
				rd_Set.chromosome_Name = chromosome_Name;
				gph.exon_Vec = Gene_Set[i].subexon; // vector of exons
				int exon_Num = gph.exon_Vec.size(); // number of exons
				int len = gph.exon_Vec[exon_Num-1].bound2 - gph.exon_Vec[0].bound1 + 1;  // gene length

				record_file<<Gene_Set[i].subexon[0].bound1<<" "<<Gene_Set[i].subexon[exon_Num-1].bound2<<" "<<len<<endl;
				cout<<Gene_Set[i].subexon[0].bound1<<" "<<Gene_Set[i].subexon[exon_Num-1].bound2<<" "<<len<<endl;
				gene_Name = gene_Name.substr(1,gene_Name.length()-3);  // extract gene name
				rd_Set.gene_Name = gene_Name;
				rd_Set.orientation = Gene_Set[i].orientation;
				int locate_Flag = locus_Assemble(infile2, rd_Set, gph, Gene_Set[i].subexon); 

				if(locate_Flag==1)  // sequencing data is located in the current gene
				{
					gene_Cnt++;
					cur_Gene = Gene_Set[i];
					// outfile<<gene_Cnt<<"\t"<<i+1<<"\t"<<gene_Name<<"\t"<<rd_Set.rd_Vec.size()<<endl;
					locate_transcriptionSites(Gene_Set[i], transcription_sites);

					// estimate the mean and variance of fragment length
					rd_Set.fragLen_Estimation(read_Number,fragLen_Mean,fragLen_Var);
					fragLen_Mean = rd_Set.mean;
					fragLen_Var = rd_Set.sigma;
					record_file<<rd_Set.rd_Vec.size()<<"\t"<<fragLen_Mean<<endl;
					num++;
					record_file<<num<<" "<<gene_Name<<endl;
					ctr_Flag = check_Valid(rd_Set, gph);

					if(ctr_Flag==true)
					{
						cout<<"read number: "<<rd_Set.rd_Vec.size()<<endl;  // output read number
						read_Cnt += rd_Set.rd_Vec.size();  
						gph.cur_Gene = cur_Gene;
						bool est_Flag = abun_EstiSingle(rd_Set, gph, data_path1,1);
						if(est_Flag==true)
							cout<<"True"<<endl;
						else
							cout<<"Too many dataDim"<<endl;	
					}			

					if(ctr_Flag==true)
					{
						gph.clear();
					}
#pragma region clear the arrays
					int n1 = rd_Left.size();
					for(int j=0; j<n1; j++)
						vector<int>().swap(rd_Left[j]);
					int n2 = rd_Right.size();
					for(int j=0; j<n2; j++)
						vector<int>().swap(rd_Right[j]);
					vector<vector<int>>().swap(rd_Left);
					vector<vector<int>>().swap(rd_Right);
					rd_Left.clear();
					rd_Right.clear();
#pragma endregion
				}
				else if(locate_Flag==-1)  // there is change of chromosome
				{
					chromosome_Name = next_chromoName;
					changed_Flag = true;
					next_chromoName = "";
					rd_Set.~Rd_Set();
					break;
				}
				else    // there is no sequencing data located in the current gene
				{
					infile2.seekg(pos1);  // located to the original position
				}

				rd_Set.~Rd_Set();
			}
			if(changed_Flag==false&&i==gene_Num)
			{
				while(getline(infile2,line))  // test a certain number of genes
				{
					stringstream stream(line);
					stream>>word;
					stream>>word; stream>>word; // read and skip the first three fields
					string chr_Name = word;     // chromosome name
					if(chr_Name!=chromosome_Name)    // there is change of chromosome name
					{
						chromosome_Name = chr_Name;  // located to the next chromosome
						changed_Flag = true;
						next_chromoName = "";
						break;
					}
				}
			}
#pragma endregion
		}
	}

	// convert the temporary file to standard file
	int total_Number = read_Cnt;
	char filename_1[200];
	sprintf(filename,"%s\\exply.txt",data_path);
	sprintf(filename_1,"%s\\results.gtf",data_path);
	expression_conver(filename, filename_1, total_Number);

	infile2.close();
	outfile.close();
	record_file.close();
}


