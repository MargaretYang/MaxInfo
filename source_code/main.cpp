
/************************************************************************************
 * Assembly.cpp
 * This file is part of the software program of MaxInfo
 * It can be redistributed and/or modified for uncommercial purposes
 *
 * First created on: 2014-06-27
 * Last modified on: 2015-02-11
 * Author: Yang Yang
 ***********************************************************************************/

//////////////////////////////////////////////////////////////////////
//// Main Function
//// Isoform Identification and Abundance Estimation

#include "Read_Assignment.h"
#include "Entropy.h"
#include "Batch_Mode.h"
#include "ISA.h"

int main(int argc, char* argv[])
{
	// read the transcripts
	char filename[200];
	char filename1[200], filename2[200], filename3[200], filename4[200];
	Gene_Batch gene_Batch;
	string chromo_Name = "1";
	vector<string> chr_Name;

	char data_path[200];
	int mode = 0;  // mode indicating whether annotations are used and whether novel isoforms are discovered		
	int opt_Num = argc - 2;  // number of arguments
	
	if(strcmp(argv[1],"-d")==0)  // requiring annotations and discovering novel isoforms
		mode = 1;
	else if(strcmp(argv[1],"-t")==0)
		mode = 2;  // requiring annotations and not discovering novel isoforms
	else if(strcmp(argv[1],"-i")==0)
		mode = 0;  // annotation-independent mode
	else 
		cout<<"Please input the right mode!"<<endl;
	
	if(opt_Num==3&&mode==0)          // without genome annotations
	{
		sprintf(filename1,argv[2]);  // file of junctions
		sprintf(filename2,argv[3]);  // file of read assignment
		sprintf(data_path,argv[4]);  // output path
		gene_Batch.batch_AbunEstimation(filename1,filename2,data_path);
	}
	else if(opt_Num==4&&mode==1)     // with genome annotations: novel transcripts are also to be discovered
	{
		sprintf(filename1,argv[2]);  // file of genome annotations
		sprintf(filename2,argv[3]);  // file of junctions                                                                                                                                                                                                                                                                                                                                                                                                       
		sprintf(filename3,argv[4]);  // file of read assignment
		sprintf(data_path,argv[5]);  // output path
		gene_Batch.abun_EstimationSpr_2(filename1,filename2,filename3,data_path);
	}
	else if(opt_Num==3&&mode==2)     //with genome annotations: not to discover novel transcripts
	{
		sprintf(filename1,argv[2]);  // file of genome annotations
		sprintf(filename2,argv[3]);  // file of read assignment                                                                                                                                                                                                                                                                                                                                                                                                     
		sprintf(data_path,argv[4]);  // output path
		gene_Batch.abun_EstimationSpr_1(filename1,filename2,data_path);
	}
	else
		cout<<"Please input right filenames!"<<argv[1]<<" "<<mode<<" "<<opt_Num<<endl;

	return 0;
}
