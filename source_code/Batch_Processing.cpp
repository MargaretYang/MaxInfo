
/////////////////////////////////////////////////////////////////////////////
////  File Processing Function
/////////////////////////////////////////////////////////////////////////////

#include "Batch_Mode.h"
#include "Read_Assignment.h"
#include "Entropy.h"
#include "ISA.h"

// write gene annotation to exon file write gene annotation to exon file
void Gene_Batch::output_GeneAnnotation(char* filename, vector<Gene>& gene)
{
	ofstream outfile;
	outfile.open(filename,ios::out);
	if(!outfile)
		cout<<"Open file(output file) error!"<<endl;

	int cnt[25];   // count the number of exons from 1 to 25 
	for(int i = 0; i<25; i++)
		cnt[i] = 0;

#pragma region output exon and isoform information
	int gene_Num = gene.size();  // number of input gene 
	for(int k = 0; k<gene_Num; k++)
	{
		int exon_Num = gene[k].subexon.size();  // number of exons 
		int trans_Num = gene[k].trans.size();   // number of isoform 

		outfile<<"@"<<k+1<<" "<<gene[k].gene_Name<<" "<<gene[k].start<<" "<<gene[k].stop<<" "<<exon_Num<<" "
			<<trans_Num<<endl;   // output gene index, name, start point, number of exons, number of isoform

		for(int i = 0; i<exon_Num; i++)
		{
			gene[k].subexon[i].len = gene[k].subexon[i].bound2 - gene[k].subexon[i].bound1 + 1; // length of subexons
			outfile<<i+1<<" "<<gene[k].subexon[i].bound1<<" "<<gene[k].subexon[i].bound2<<" "<<gene[k].subexon[i].len<<endl; // output exons' information
		}
	}

	outfile.close();

#pragma endregion
}

// write gene information to file
void Gene_Batch::output_transBed(char* path1, char* filename, vector<Gene>& gene)
{
#pragma region output isoform information
	// write gene information to file
	ofstream outfile;
	outfile.open(filename,ios::out);
	if(!outfile)
	{
		cout<<"Open file error!"<<endl;
		return;
	}

	int cnt[100];   // count the number of exons from 1 to 25 
	int trs_cnt[100];   // count the number of isoform
	for(int i = 0; i<100; i++)
		cnt[i] = 0;
	for(int i = 0; i<100; i++)
		trs_cnt[i] = 0;

	int gene_Num = gene.size();  // number of input gene
	for(int k = 0; k<gene_Num; k++)
	{
		int exon_Num = gene[k].subexon.size();  // number of exons
		outfile<<"@"<<k+1<<" "<<gene[k].gene_Name<<" "<<exon_Num<<" "
			<<gene[k].trans.size()<<" ";   // output gene index, name, start point, number of exons, number of isoform
		int trans_Num = gene[k].trans.size();   // number of isoform
		if(trans_Num<100)
			trs_cnt[trans_Num]++;
		for(int i = 0; i<trans_Num; i++)
			outfile<<gene[k].trans[i].length<<" ";
		outfile<<endl;


		int itemRgb = 0;
#pragma region output isoform
		for(int i = 0; i<trans_Num; i++)
		{
			int blockSize = gene[k].trans[i].exon_start.size();  // number of exons in the isoform

			outfile<<i+1<<'\t'<<gene[k].trans[i].start-1<<'\t'<<gene[k].trans[i].stop<<'\t'
				<<gene[k].trans[i].trsName<<'\t'<<0<<'\t'<<gene[k].trans[i].orientation<<'\t'
				<<gene[k].trans[i].start<<'\t'<<gene[k].trans[i].stop<<'\t'
				<<itemRgb<<'\t'<<blockSize<<'\t'; // output exons

			if(blockSize>0)
			{
				int len = 0;
				cnt[blockSize]++;  // increase the count
				// output length of subexons
				for(int l = 0; l<blockSize; l++)
				{
					int tLen = gene[k].trans[i].exon_len[l];
					outfile<<tLen<<",";
					len += tLen;
				}
				outfile<<'\t';

				// output the start point of the subexons
				for(int l = 0; l<blockSize; l++)
				{
					outfile<<gene[k].trans[i].exon_start[l]-gene[k].trans[i].start<<",";
				}

				outfile<<endl;	
			}
			else
			{
				int len = gene[k].trans[i].length; // isoform with no exons
				outfile<<len<<","<<"\t"<<0<<endl;
				cnt[1]++;
			}
		}
#pragma endregion

		vector<Exon_Node>().swap(gene[k].subexon);
		vector<Trans>().swap(gene[k].trans);
	}

	ofstream outfile1;
	char filename[200];
	sprintf(filename,"%s\\number.txt",path1);
	outfile1.open(filename,ios::app);
	
	for(int i = 0; i<100; i++)
		outfile1<<i<<":"<<cnt[i]<<" ";
	outfile1<<endl;
	for(int i = 0; i<100; i++)
		outfile1<<i<<":"<<trs_cnt[i]<<" ";
	outfile1<<endl;

	outfile1.close();
	outfile.close();	

#pragma endregion output isoform information
}

// write gene information to files
void Gene_Batch::output_standardBed(char* path1, char* filename, vector<Gene>& gene, string chromo_Name)
{
#pragma region output isoform information
	// write exons to files
	ofstream outfile;
	outfile.open(filename,ios::out);
	if(!outfile)
	{
		cout<<"Open file error!"<<endl;
		return;
	}

	int cnt[100];   // count the number of exons from 1 to 25 
	int trs_cnt[100];   // count the number of isoforms
	for(int i = 0; i<100; i++)
		cnt[i] = 0;
	for(int i = 0; i<100; i++)
		trs_cnt[i] = 0;

	int gene_Num = gene.size();  // number of input gene
	for(int k = 0; k<gene_Num; k++)
	{
		int exon_Num = gene[k].subexon.size();  // number of exons
		//outfile<<"@"<<k+1<<" "<<gene[k].gene_Name<<" "<<exon_Num<<" "
		//	<<gene[k].trans.size()<<" ";   // output gene index, name, start point, number of exons, number of isoform
		int trans_Num = gene[k].trans.size();   // number of isoforms
		if(trans_Num<100)
			trs_cnt[trans_Num]++;
		for(int i = 0; i<trans_Num; i++)
			outfile<<gene[k].trans[i].length<<" ";
		outfile<<endl;


		int itemRgb = 0;
#pragma region output isoform
		for(int i = 0; i<trans_Num; i++)
		{
			int blockSize = gene[k].trans[i].exon_start.size();  // number of the exons in the isoform

			outfile<<chromo_Name<<'\t'<<gene[k].trans[i].start-1<<'\t'<<gene[k].trans[i].stop<<'\t'
				<<gene[k].trans[i].trsName<<'\t'<<0<<'\t'<<gene[k].trans[i].orientation<<'\t'
				<<gene[k].trans[i].start<<'\t'<<gene[k].trans[i].stop<<'\t'
				<<itemRgb<<'\t'<<blockSize<<'\t'; // output extons
			if(blockSize>0)
			{
				int len = 0;
				cnt[blockSize]++;  // increase count
				// output the length of subexons
				for(int l = 0; l<blockSize; l++)
				{
					int tLen = gene[k].trans[i].exon_len[l];
					outfile<<tLen<<",";
					len += tLen;
				}
				outfile<<'\t';

				// output the start point of the exons
				for(int l = 0; l<blockSize; l++)
				{
					outfile<<gene[k].trans[i].exon_start[l]-gene[k].trans[i].start<<",";
				}

				outfile<<endl;	
			}
			else
			{
				int len = gene[k].trans[i].length; // isoform with no exons
				outfile<<len<<","<<"\t"<<0<<endl;
				cnt[1]++;
			}
		}
#pragma endregion

		vector<Exon_Node>().swap(gene[k].subexon);
		vector<Trans>().swap(gene[k].trans);
	}

	outfile.close();	

#pragma endregion output isoform information
}

// output the comparison of isoform identification
void Gene_Batch::output_Comparison(int candidate, string gene_Name, vector<Trans>& trs_SetCtrl, vector<Trans>& trs_SetCase, 
	float precision, float recall, ofstream& outfile)
{
	int trs_num1 = trs_SetCtrl.size();  // isoform number of reference gene
	int trs_num2 = trs_SetCase.size();  // estimate isoform number of current gene

	outfile<<"@"<<candidate<<"\t"<<gene_Name<<"\t"<<trs_num1<<"\t"<<trs_num2<<"\t"<<endl;
	// write comparison result to file as bed format
	for(int i = 0; i<trs_num1; i++)
	{
		int exon_num1 = trs_SetCtrl[i].exon_start.size(); // number of exons
		outfile<<i+1<<"\t"<<trs_SetCtrl[i].start<<"\t"<<trs_SetCtrl[i].stop<<"\t"<<exon_num1<<"\t";	
		for(int j = 0; j<exon_num1; j++)
			outfile<<trs_SetCtrl[i].exon_len[j]<<",";
		outfile<<"\t";
		for(int j = 0; j<exon_num1; j++)
			outfile<<trs_SetCtrl[i].exon_start[j]<<",";
		outfile<<endl;
	}

	for(int i = 0; i<trs_num2; i++)
	{
		int exon_num2 = trs_SetCase[i].exon_start.size(); // number of exons
		outfile<<i+1<<"\t"<<trs_SetCase[i].start<<"\t"<<trs_SetCase[i].stop<<"\t"<<exon_num2<<"\t";
		for(int j = 0; j<exon_num2; j++)
			outfile<<trs_SetCase[i].exon_len[j]<<",";
		outfile<<"\t";
		for(int j = 0; j<exon_num2; j++)
			outfile<<trs_SetCase[i].exon_start[j]<<",";
		outfile<<endl;
	}

	cout<<precision/trs_num2<<"\t"<<recall/trs_num1<<endl;

#pragma endregion write comparison to file

}

// transform the file from refGene format to BED format
void Gene_Batch::refGene_transBed(char* filename1, char* filename2)
{
	ifstream infile;
	infile.open(filename1,ios::in);
	if(!infile)
		cout<<"Open file error!"<<endl;
	ofstream outfile;
	outfile.open(filename2,ios::out);
	if(!outfile)
		cout<<"Open file error!"<<endl;

	string line, word, seg[20];
	string::size_type pos, pre;
	char tPar=',';
	int exon_Len = 0, exon_Start = 0, exon_Stop = 0, exon_Number = 0;
	int start = 0, stop = 0;
	vector<int> exon_len;     // vector the exons' length
	vector<int> exon_start;   // vector the exons' start point
	vector<int> exon_stop;    // vector the exons' end point
	std::streampos pos0;
	int quality = 0;

	// obtaion isoform form annotation
	while(getline(infile,line))
	{
		exon_start.clear(); exon_stop.clear();
		stringstream stream(line);
		// 0: index 1£ºisoform name 2: chromosome name 3: direction 4£ºstart point 5£ºend point 6£ºtranslation start point 7£ºtranslation end point
		// 8£ºnumber of exons 9: start point of exons 10£ºend point of exons
		for(int k = 0; k<11; k++)
			stream>>seg[k];   
		sscanf(seg[4].c_str(),"%d",&start);
		sscanf(seg[5].c_str(),"%d",&stop);
		sscanf(seg[8].c_str(),"%d",&exon_Number);

		cout<<seg[1]<<endl;  // isofrom name

		// extract subexons from annotation
#pragma region extract subexons from annotation
		pre = 0;
		while((pos=seg[9].find(tPar,pre))!=string::npos)  // find ',' from the string
		{
			string str = seg[9].substr(pre,pos-pre);
			sscanf(str.c_str(),"%d",&exon_Start);   // length of exons
			exon_start.push_back(exon_Start);		
			pre = pos + 1;
		}		

		pre = 0;
		while((pos=seg[10].find(tPar,pre))!=string::npos)  // find ',' from the string
		{
			string str = seg[10].substr(pre,pos-pre);
			sscanf(str.c_str(),"%d",&exon_Stop);   // length of exons
			exon_stop.push_back(exon_Stop);		
			pre = pos + 1;
		}
#pragma endregion extract subexons from annotation

#pragma region write isoform information to bed format file
		outfile<<seg[2]<<"\t"<<seg[4]<<"\t"<<seg[5]<<"\t"<<seg[1]<<"\t"<<quality<<"\t"
			<<seg[3]<<"\t"<<seg[6]<<"\t"<<seg[7]<<"\t"<<quality<<"\t"<<exon_Number<<"\t";

		for(int i=0; i<exon_Number; i++)
		{
			outfile<<exon_stop[i]-exon_start[i]<<",";
		}
		outfile<<"\t";
		for(int i=0; i<exon_Number; i++)
		{
			outfile<<exon_start[i]-start<<",";
		}
		outfile<<endl;
#pragma endregion

	}

	infile.close();
	outfile.close();
}

// seperate the file by chromosome
void Gene_Batch::separate_File(char* filename1, char* filename2)
{
	ifstream infile;
	infile.open(filename1,ios::in);
	if(!infile)
		cout<<"Open file error!"<<endl;

	char filename[200];

	string line, word, chrom_Name;
	string::size_type pos, pre = 0;
	char tPar = '|';
	int quality = 0;
	string pre_name = "";

	getline(infile,line);
	stringstream stream(line);
	stream>>chrom_Name;  // load chromosone name
	ofstream outfile;
	if(pre_name=="")
	{
		pre_name = chrom_Name;  // chromosone name
		sprintf(filename,"%s\\%s.bed",filename2,chrom_Name.c_str());
		outfile.open(filename,ios::out);
		if(!outfile)
			cout<<"Open file error!"<<endl;
		outfile<<line<<endl;
	}

	// isoform obtained from the annotation
	while(getline(infile,line))
	{
		stringstream stream(line);
		stream>>chrom_Name;  // load chromosone name
		if(chrom_Name!=pre_name)
		{
			pre_name = chrom_Name;  // chromosone name
			cout<<chrom_Name<<endl;
			if((pos=chrom_Name.find(tPar,pre))!=string::npos)  // find the separator
			{
				string temp_Chr = chrom_Name.substr(pre,pos-pre);
				if(temp_Chr!=pre_name)
				{
					outfile.close();
					sprintf(filename,"%s\\%s.bed",filename2,temp_Chr.c_str());
					outfile.open(filename,ios::app);  // open new file
					pre_name = temp_Chr;
				}
			}
			else
			{
				outfile.close();
				sprintf(filename,"%s\\%s.bed",filename2,chrom_Name.c_str());
				outfile.open(filename,ios::out);  // open new file
			}
			if(!outfile)
				cout<<"Open file error!"<<endl;					
		}
		outfile<<line<<endl;
	}
	outfile.close();
}

// calculate the distribution of the number of exons and isoforms
void Gene_Batch::statistical(char* filename1, char* filename2)
{
	ifstream infile;
	infile.open(filename1,ios::in);
	if(!infile)
		cout<<"Open file error!"<<endl;	
	
	string line, word, seg[100];
	string::size_type pos = 0, pre = 0;
	char tPar = ':';
	int cnt[100], trs_cnt[100];
	for(int i=0; i<100; i++)
	{
		cnt[i] = 0;
		trs_cnt[i] = 0;
	}
	
	int serial = 0, ser = 0, num = 0;
	string serstr, numstr;
	while(getline(infile,line))
	{
		stringstream stream(line);
		for(int k=0; k<100; k++)
		{
			stream>>seg[k];
			if((pos = seg[k].find(tPar,0))!=string::npos)
			{
				serstr = seg[k].substr(0,pos);
				numstr = seg[k].substr(pos+1,seg[k].length()-pos-1);
				sscanf(serstr.c_str(),"%d",&ser);
				sscanf(numstr.c_str(),"%d",&num);
			}
			if(serial%2==0)
				cnt[ser] += num;
			else
				trs_cnt[ser] += num;
		}
		serial++;
	}

	int sum1 = 0, sum2 = 0;
	ofstream outfile;
	outfile.open(filename2,ios::out);
	for(int i=0; i<100; i++)
	{
		outfile<<cnt[i]<<"\t";
		sum1 += cnt[i];  // number of exons
	}
	outfile<<endl;
	for(int i=0; i<100; i++)
	{
		outfile<<trs_cnt[i]<<"\t";
		sum2 += trs_cnt[i];  // number of isoforms
	}
	outfile<<endl;
	outfile<<sum1<<"\t"<<sum2<<endl;

	infile.close();
	outfile.close();
}

// comparison result
void Gene_Batch::estimation_Compare(int candidate, string gene_Name, vector<Trans>& trs_SetCtrl, vector<Trans>& trs_SetCase, float& precision, float& recall)
{
#pragma region compare the exons
	// the strictest compare method
	int trs_num1 = trs_SetCtrl.size();   // annotated isoform number
	int trs_num2 = trs_SetCase.size();   // estimated isoform number

	int* compare_flag1 = new int[trs_num1];    // comparison of the isoform
	for(int i = 0; i<trs_num1; i++)
		compare_flag1[i] = -1;

	int* compare_flag2 = new int[trs_num2];    // comparison of the isoform
	for(int i = 0; i<trs_num2; i++)
		compare_flag2[i] = -1;

	float** score = new float*[trs_num2];
	for(int i = 0; i<trs_num2; i++)
	{
		score[i] = new float[trs_num1];   // the similarity between the annotated isoform and estimated isoform
		for(int j = 0; j<trs_num1; j++)
			score[i][j] = 0;
	}	

	int exon_num1 = 0, exon_num2 = 0, exon_cnt = 0;
	int start1 = 0, start2 = 0, stop1 = 0, stop2 = 0;
	int b1 = 0, b2 = 0;  // junction location
	float len1 = 0, len2 = 0;  // length of exons
	bool coin_flag = false;
	float thresh1 = 30;

#pragma region matching level of the isoform
	for(int k = 0; k<trs_num2; k++)
	{		
		exon_num2 = trs_SetCase[k].exon_start.size();  // estimate the number of exons in the isoform
		exon_cnt = 0;
		start2 = trs_SetCase[k].start;  
		stop2 = trs_SetCase[k].stop;
		for(int l = 0; l<trs_num1; l++)
		{
			exon_num1 = trs_SetCtrl[l].exon_start.size();  // annotated exon number in the isoform
			start1 = trs_SetCtrl[l].start;   
			// stop1 = trs_SetCtrl[l].stop;                // estimated start point and end point may be different from the exons
			
			stop1 = trs_SetCtrl[l].exon_start[exon_num1-1] + trs_SetCtrl[l].exon_len[exon_num1-1] - 1;

			exon_num1 = trs_SetCtrl[l].exon_start.size();  // annotated exon number in the isoform
			if(exon_num1==exon_num2)  // the numbers of the exons are identical
			{
				int s1 = start1-start2, s2 = stop1-stop2;
				if(exon_num1>1/*&&abs(start1-start2)<thresh1&&abs(stop1-stop2)<thresh1*/)  // the numbers of the exons and the start points and end points of the isoform are identical
				{
					// comparison of the first and last exton
					b1 = trs_SetCase[k].exon_stop[0] - trs_SetCtrl[l].exon_stop[0];  // right location difference of the first exon
					b2 = trs_SetCase[k].exon_start[exon_num2-1] - trs_SetCtrl[l].exon_start[exon_num1-1];  // left location difference of the last exon
					len1 = trs_SetCtrl[l].exon_len[0], len2 = trs_SetCtrl[l].exon_len[exon_num1-1];
					float r11 = abs(s1)*1.0/len1, r12 = abs(b1)*1.0/len1;  // error rate of the first exon
					float r22 = abs(s2)*1.0/len2, r21 = abs(b2)*1.0/len2;  // error rate of the last exon
					bool ctr_Flag = true;
					if((len1<500&&r12>0.05)||(len1>500&&abs(b1)>50))
						ctr_Flag = false;
					if((len2<500&&r21>0.05)||(len2>500&&abs(b2)>50))
						ctr_Flag = false;

					if(ctr_Flag)
					{
						score[k][l] = 1 - (r11+r12+r21+r22)*0.5/exon_num1;						
#pragma region comparison of the exons
						for(int j = 1; j<exon_num2-1; j++)  // the boundries of the exon should be indentical
						{
							len1 = trs_SetCtrl[l].exon_len[j];  // length of exons
							len2 = trs_SetCase[k].exon_len[j]; 
							b1 = trs_SetCase[k].exon_start[j] - trs_SetCtrl[l].exon_start[j];  // difference of junction position£ºdifference of the right boundary of the exon
							b2 = trs_SetCase[k].exon_stop[j] - trs_SetCtrl[l].exon_stop[j];  // difference of junction position£ºdifference of the left boundary of the exon

							float rate1 = abs(b1)*1.0/len1;  // error rate of the exon
							float rate2 = abs(b2)*1.0/len1;  // error rate of the exon

							if(rate1<0.05&&rate2<0.05)   // low error rate 
							{
								score[k][l] -= (rate1+rate2)/(2*exon_num1);  // the higher the error rate, the lower the score
							}
							else
							{
								score[k][l] = 0;  // mismatch if the error rate is high
								break;
							}
						}
#pragma endregion comparison of the exons
					}
				}
				else if(exon_num2==1)  // if the isoform consists of 1 exon, it is considered mismatch if the error rate is over 20%
				{
#pragma region isoform with 1 exon
					float len1 = trs_SetCtrl[l].length;   // annotated isoform length
					float len2 = trs_SetCase[k].length;   // estimated isoform length				
					float overlap_rate1 = 0, overlap_rate2 = 0, overlap_len = 0;
					float thresh2 = 0.80;
					if(start2<=stop1&&start1<=stop2)  // two isoforms overlap
					{
						int s1 = _max(start1,start2);
						int s2 = _min(stop1,stop2);
						overlap_len = s2 - s1 + 1;   // overlap length
						overlap_rate1 = overlap_len*1.0/len1;  // overlap length
						overlap_rate2 = overlap_len*1.0/len2;
					}

					if(overlap_rate1>thresh2&&overlap_rate2>thresh2) // high overlap rate
					{
						score[k][l] = overlap_rate1; // matching score
					}

#pragma endregion isoform with 1 exon
				}
			}
		}
	}

#pragma endregion matching level of the isoform

#pragma region find the matching isoform according to the matching level

	for(int k = 0; k<trs_num2; k++)
	{
		float temp_max = 0.8;
		for(int l = 0; l<trs_num1; l++)
		{
			if(score[k][l]>temp_max){
				compare_flag2[k] = l;
				temp_max = score[k][l];
			}
		}
	}

	for(int l = 0; l<trs_num1; l++)
	{
		float temp_max = 0.8;
		for(int k = 0; k<trs_num2; k++)
		{
			if(score[k][l]>temp_max){
				compare_flag1[l] = k;
				temp_max = score[k][l];
			}
		}
	}
#pragma endregion find the matching isoform according to the matching level

	// float precision = 0, recall = 0;
	precision = 0;
	recall = 0;
	float cnt1 = 0, cnt2 = 0;

	for(int i = 0; i<trs_num1; i++){
		if(compare_flag1[i]>-1)
			cnt1 += score[compare_flag1[i]][i];
	}

	for(int i = 0; i<trs_num2; i++){
		if(compare_flag2[i]>-1)
			cnt2 += score[i][compare_flag2[i]];
	}

	// precision = cnt2/trs_num2;
	// recall = cnt1/trs_num1;

	precision = cnt2;
	recall = cnt1;

#pragma region write results to file
	ofstream outfile;
	outfile.open("G:\\My_documents\\FragmentStreaming\\Data_Species\\hg19\\temp_Compare1.txt",ios::app);

	outfile<<"@"<<candidate<<"\t"<<gene_Name<<"\t"<<trs_num1<<"\t"<<trs_num2<<"\t"<<endl;
	// write comparison result to file as bed format
	for(int i = 0; i<trs_num1; i++)
	{
		outfile<<i+1<<"\t"<<trs_SetCtrl[i].start<<"\t"<<trs_SetCtrl[i].stop<<"\t";			
		int exon_num1 = trs_SetCtrl[i].exon_start.size(); // number of exons
		for(int j = 0; j<exon_num1; j++)
			outfile<<trs_SetCtrl[i].exon_len[j]<<",";
		outfile<<"\t";
		for(int j = 0; j<exon_num1; j++)
			outfile<<trs_SetCtrl[i].exon_start[j]<<",";
		//outfile<<"\t"<<trs_SetCtrl[i].length<<"\t"<<trs_SetCtrl[i].abundance<<endl;
		outfile<<endl;
	}

	for(int i = 0; i<trs_num2; i++)
	{
		outfile<<i+1<<"\t"<<trs_SetCase[i].start<<"\t"<<trs_SetCase[i].stop<<"\t";
		int exon_num2 = trs_SetCase[i].exon_start.size(); // number of exons
		for(int j = 0; j<exon_num2; j++)
			outfile<<trs_SetCase[i].exon_len[j]<<",";
		outfile<<"\t";
		for(int j = 0; j<exon_num2; j++)
			outfile<<trs_SetCase[i].exon_start[j]<<",";
		// outfile<<"\t"<<trs_SetCase[i].length<<"\t"<<trs_SetCase[i].abundance<<endl;
		outfile<<endl;
	}
	for(int i = 0; i<trs_num1; i++)
		outfile<<i+1<<":"<<compare_flag1[i]+1<<"\t";  // output from index 1
	outfile<<recall<<endl;

	for(int i = 0; i<trs_num2; i++)
		outfile<<i+1<<":"<<compare_flag2[i]+1<<"\t";  // output from index 1
	outfile<<precision<<endl;

	for(int i = 0; i<trs_num2; i++)
	{
		for(int j = 0; j<trs_num1; j++)
		{
			outfile<<score[i][j]<<"\t";
		}
		outfile<<endl;
	}
#pragma endregion

	delete[] compare_flag1;
	delete[] compare_flag2;
#pragma endregion compare exons
}

// exon comparison
void Gene_Batch::trans_Compare(char* filename1, char* filename2, char* filename3)
{
	ifstream infile1, infile2;
	ofstream outfile;
	infile1.open(filename1,ios::in);
	if(!infile1)
		cout<<"Open file error!"<<endl;
	infile2.open(filename2,ios::in);
	if(!infile2)
		cout<<"Open file error£¡"<<endl;
	outfile.open(filename3,ios::out);
	if(!outfile)
		cout<<"Open file error!"<<endl;

	string::size_type pos0;
	string line, line_pre, word, seg[20], chromo_Name, gene_Name;
	int exon_Num = 0, trs_Num = 0, trs_start = 0, trs_Serial = 0, start = 0, stop = 0;
	int gene_idx = 9;
	vector<Gene> gene;
	vector<int> exon_start, exon_stop, exon_len;
	vector<Trans> trs_SetCase;
	int trsnum_Ctrl = 0, trsnum_Case = 0;
	int gene_Cnt1 = 0, gene_Cnt2 = 0;

	std::streampos pos = infile2.tellg();   // get the position of current file pointer
	getline(infile2,line);  // read a line
	bool flag_ctr = false;
	if(line[0]=='#')  // the first line is note
	{
		// the first line is note
		getline(infile2,line);
		stringstream stream(line);
		stream>>chromo_Name;  // load chromosome name
		trsnum_Ctrl += load_Annotation_1(gene, filename1, chromo_Name);  // load all the gene on the chromosome
		getline(infile2,line);
		stringstream stream1(line);
		stream1>>chromo_Name; // load chromosome name
		infile2.seekg(pos);  // locate return position
		getline(infile2,line);
	}
	else
	{
		stringstream stream(line);
		stream>>chromo_Name;  // load chromosome name
		trsnum_Ctrl += load_Annotation_1(gene, filename1, chromo_Name);  // load all the gene on the chromosome
		infile2.close();
		infile2.open(filename2,ios::in);  // reopen the file
	}

	getline(infile2,line_pre);
	string pre_Name = chromo_Name;  // last chormosome name
	string pre_GeneName;
	int gene_Serial = 0;
	int gene_Num = gene.size();  // number of gene
	bool begin_Flag = false;
	int local_Start = 0, local_Stop = 0;  // local start point and end point
	float precision_cnt = 0, recall_cnt = 0;

	while(flag_ctr==false)  // load isoform estimation
	{
		line = line_pre;
		stringstream stream(line);
		for(int k=0; k<16; k++)
			stream>>seg[k];   // 0: chormosome name 1: function type 2£ºtype 3£ºstart point 4£ºend point 9: gene index 13£ºabundance estimation
		sscanf(seg[3].c_str(),"%d",&start);
		sscanf(seg[4].c_str(),"%d",&stop);

		chromo_Name = seg[0];  // chormosome name
		if(pos0=chromo_Name.find('|',0)!=string::npos)
			chromo_Name = chromo_Name.substr(0,pos0);
		if(chromo_Name!=pre_Name)  // next chormosome
		{
			vector<Gene>().swap(gene);  // calculate the gene on the new segment
			gene_Serial = 0;
			trsnum_Ctrl += load_Annotation_1(gene, filename1, chromo_Name);  // load all the gene on the chromosome
			pre_Name = chromo_Name;
			gene_Num = gene.size();     // number of gene
			cout<<"Chromosome changed: "<<chromo_Name<<endl;
		}


		int len = stop - start + 1;

		if(chromo_Name==pre_Name)
		{
#pragma region load isoform identification result
			if(seg[2]=="transcript")
			{	
				if(seg[1]=="Cufflinks"||seg[1]=="SLIDE")
					gene_Name = seg[9];  // gene name
				else
					gene_Name = seg[6];

				if(begin_Flag==false)
				{
					pre_GeneName = gene_Name;
					begin_Flag = true;
				}
				else if(gene_Name!=pre_GeneName)  // gene transformation
				{
					begin_Flag = false;
				}
#pragma region load isoform
				if(begin_Flag==true)
				{
					Trans cur_Trans;		
					cur_Trans.start = start; cur_Trans.stop = stop;
					trs_start = start;   // start point of isoform
					double temp_Abun = 0, temp_FPKM = 0;

					gene_Name = seg[gene_idx];  // gene name
					if(seg[1]=="Cufflinks"||seg[1]=="SLIDE")
					{
						int temp_Len = seg[15].length();
						string temp = seg[15].substr(1,temp_Len-3);
						sscanf(temp.c_str(),"%lf",&temp_Abun);
						temp_Len = seg[13].length();
						temp = seg[13].substr(1,temp_Len-3);
						sscanf(temp.c_str(),"%lf",&temp_FPKM);
						cur_Trans.FPKM1 = temp_FPKM;
					}
					else
						sscanf(seg[12].c_str(),"%lf",&temp_Abun);
					cur_Trans.abundance = temp_Abun;
					cur_Trans.FPKM = cur_Trans.abundance*1e09/cur_Trans.length;  // tranfer to FPKM value
					if(seg[1]=="new")
						cur_Trans.FPKM1 = cur_Trans.FPKM;

					trs_SetCase.push_back(cur_Trans);	
					// trs_Serial++;
				}
#pragma endregion
			}
			else if(seg[2]=="exon")
			{
				exon_start.push_back(start);   // start point of subexons
				exon_len.push_back(len);       // length of subexons
			}

			if(begin_Flag==true&&!getline(infile2,line_pre)) flag_ctr = true;  // end of file

			if(seg[2]=="transcript"||flag_ctr==true)
			{			
				if(trs_Serial>0||flag_ctr==true)  // add exon information of last isoform
				{
					int trs_Pre = _max(0,trs_Serial-1);
#pragma region combine neighboring exons
					int exon_num = exon_len.size();  // number of exons
					if(trs_Serial==1)
					{
						local_Start = exon_start[0];
						local_Stop = exon_start[exon_num-1] + exon_len[exon_num-1] - 1;
					}
					else
					{
						local_Start = exon_start[0]<local_Start?exon_start[0]:local_Start;
						int temp = exon_start[exon_num-1] + exon_len[exon_num-1] - 1;
						local_Stop = temp>local_Stop?temp:local_Stop;  // local end point
					}
					int k = 0, trs_len = 0;
					while(k<exon_num)
					{
						int j = k;
						while(j<exon_num-1&&(abs(exon_start[j]+exon_len[j]-exon_start[j+1])<2))
							j++;
						len = exon_start[j] - exon_start[k] + exon_len[j];
						trs_SetCase[trs_Pre].exon_start.push_back(exon_start[k]);
						trs_SetCase[trs_Pre].exon_len.push_back(len);
						trs_SetCase[trs_Pre].exon_stop.push_back(exon_start[k]+len-1);
						trs_len += len;
						k = j+1;
					}
#pragma endregion combine neighboring exons
					trs_SetCase[trs_Pre].length = trs_len;  // length of isoform		
					vector<int>().swap(exon_start); 
					vector<int>().swap(exon_stop);
					vector<int>().swap(exon_len);
				}
				trs_Serial++;
			}

			if(begin_Flag==false||flag_ctr==true)   // end of current gene segment
			{
#pragma region find the most similar estimated isoform from the annotated isoform
				trsnum_Case += trs_SetCase.size();  // number of isoform
				int i = gene_Serial;
				int candidate = -1;
				float precision = 0, recall = 0;
				while(i<gene_Num)
				{
					Gene& cur_Gene = gene[i];
					int s1 = cur_Gene.start;
					int s2 = cur_Gene.stop;
					if(s2>local_Start&&s1<local_Stop)  // the gene regions overlap
					{
						float precision_1 = 0, recall_1 = 0;
						gene_Cnt1++;
						estimation_Compare(i, gene[i].gene_Name, gene[i].trans,trs_SetCase, precision_1, recall_1);
						if(precision_1>precision&&recall_1>recall)
						{
							precision = precision_1;
							recall = recall_1;
							cout<<gene[i].gene_Name<<" "<<precision<<" "<<recall<<endl;
							gene_Cnt2++;
							candidate = i;
						}
					}
					else if(s1>local_Start)
					{
						gene_Serial = _max(0,i-5);
						break;
					}
					i++;
				}
#pragma endregion

				if(candidate>=0)  // write comparison result to file
				{
					precision_cnt += precision;  // precision
					recall_cnt += recall; // recall
					output_Comparison(candidate, gene[candidate].gene_Name, gene[candidate].trans, trs_SetCase, precision, recall, outfile);
				}
				vector<Trans>().swap(trs_SetCase);
				trs_Serial = 0;
				local_Start = 0;
				local_Stop = 0;
			}
#pragma endregion
		}
	}

	outfile<<trsnum_Case<<"\t"<<trsnum_Ctrl<<"\t"<<precision_cnt<<"\t"<<recall_cnt<<"\t"<<precision_cnt/trsnum_Case<<"\t"<<recall_cnt/trsnum_Ctrl<<endl;
	outfile<<gene_Cnt1<<"\t"<<gene_Cnt2<<endl;
	infile1.close();
	infile2.close();
	outfile.close();
}

// output the result of exon assembly
void Gene_Batch::exon_Print(char* filename)
{
	ofstream outfile;
	outfile.open(filename,ios::out);
	if(!outfile)
		cout<<"Open file error!"<<endl;

	int exon_Num = exon_Assembly.size();
	for(int i=0; i<exon_Num; i++)
	{
		exon_Assembly[i].len = exon_Assembly[i].bound2-exon_Assembly[i].bound1+1;
		outfile<<i+1<<"\t"<<exon_Assembly[i].bound1<<"\t"<<exon_Assembly[i].bound2<<"\t"<<exon_Assembly[i].len<<endl;
	}
}

// output the result of exon assembly
void Gene_Batch::exon_PrintBatch(ofstream& outfile, int gene_Cnt, vector<Exon_Node>& exon_Vec)
{
	if(!outfile)
	{
		cout<<"Open file error!"<<endl;
		return;
	}

	int exon_Num = exon_Vec.size();
	outfile<<gene_Cnt<<" "<<exon_Num<<endl;
	for(int i=0; i<exon_Num; i++)
	{
		exon_Vec[i].len = exon_Vec[i].bound2-exon_Vec[i].bound1+1;
		outfile<<i+1<<"\t"<<exon_Vec[i].bound1<<"\t"<<exon_Vec[i].bound2<<"\t"<<exon_Vec[i].len<<endl;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
// write the result of gene location to file
void Gene_Batch::gene_LocateOutput(Rd_Set& rd_Set, Graph_Trans& gph, string gene_Name, ofstream& outfile, int serial)
{
	int exon_Num = gph.exon_Vec.size();  // number of exons
	int len = gph.exon_Vec[exon_Num-1].bound2-gph.exon_Vec[exon_Num-1].bound1+1;
	outfile<<serial<<" "<<gene_Name<<" read number: "<<rd_Set.rd_Vec.size()<<endl;  // output the number of read
	outfile<<exon_Num<<" "<<gph.exon_Vec[0].bound1<<" "<<gph.exon_Vec[exon_Num-1].bound2<<" "<<len<<endl;
	for(int i=0; i<exon_Num-1; i++)
	{
		if(gph.link_Mtx[i+1][i+2]>0)
			outfile<<i+1<<" "<<i+2<<" : "<<gph.link_Mtx[i+1][i+2]<<" "<<
			gph.exon_Vec[i].bound2<<" "<<gph.exon_Vec[i+1].bound1<<endl;
	}

	int n1 = gph.junc_Link1.size();  // junction of exons and introns
	int n2 = gph.junc_Link2.size();  // junction of exons and introns
	int num1 = 0, num2 = 0;
	for(int i=0; i<n1; i++)
	{
		if(gph.junc_Link1[i]>0){
			outfile<<i+1<<":"<<gph.junc_Link1[i]<<" ";
			num1++;
		}
	}
	outfile<<num1<<endl;
	for(int i=0; i<n2; i++)
	{
		if(gph.junc_Link2[i]>0){
			outfile<<i+1<<":"<<gph.junc_Link2[i]<<" ";
			num2++;
		}
	}
	outfile<<num2<<endl;

#pragma region count the number of read in the introns
	int seg_Num = rd_Left.size();  // number of segment
	int cnt1 = 0, cnt2 = 0;
	int s1 = 0, s2 = 0, z1 = 0, z2 = 0;
	for(int i=1; i<seg_Num; i+=2)
	{
		s1 = 0; s2 = 0; z1 = 0; z2 = 0;
		cnt1 = rd_Left[i].size();
		cnt2 = rd_Right[i].size();
		if(cnt1>0)
		{
			s1 = rd_Left[i][0];
			s2 = rd_Left[i][cnt1-1];
		}
		if(cnt2>0)
		{
			sort(rd_Right[i].begin(),rd_Right[i].end());
			z1 = rd_Right[i][0];
			z2 = rd_Right[i][cnt2-1];
		}
		if(cnt1>0||cnt2>0)
			outfile<<(i+1)/2<<": "<<s1<<" "<<s2<<" "<<z1<<" "<<z2<<" "<<cnt1<<" "<<cnt2<<endl;
	}
#pragma endregion count the number of read in the introns

}

// load in the gene structures based on the genome annotations
void Gene_Batch::load_Annotation(vector<Gene>& gene, char* filename1, string chromo_Name)
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
	int seg_Num = 17;  // number of fields
	int line_Num = 0;
	int gene_Serial = -1, pre_Serial = -1;   // gene serials
	int chr_idx = 0, type_idx = 2, label_idx = 9, ori_idx = 6, trs_idx = 11;
	bool ctr_Flag = false;
	string Gene_Id, Gene_IdPre;  // index of gene
	int chrName_Len = chromo_Name.length();
	string::size_type pos = 0, pre = 0;
	char tPar = '|';

#pragma region locate the field of GENE_ID
	std::streampos pos1 = myfile.tellg();   // obtain the present file pointer
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
	if(seg[i]=="gene_name"||seg[i]=="gene_id"||seg[i]=="Gene_Id"||seg[i]=="Gene_id")  // search the field of "gene_id"
		label_idx = i+1;
	i = 0;
	stringstream stream1(line);
	stream1>>seg[0];
	while(seg[i]!="transcript_id")
	{
		i++;
		stream1>>seg[i];
		cout<<seg[i]<<endl;
	}
	if(seg[i]=="transcript_id")  // search the field of "transcript_id"
		trs_idx = i+1;
	myfile.close();
	myfile.open(filename1,ios::in);   // back to the first line
#pragma endregion

	bool locate_Flag = false;
	bool iden_Flag = false;
	string str;

	while(std::getline(myfile,line))
	{
		stringstream stream(line);
		for(int i = 0; i<20; i++)
			stream>>seg[i];

		iden_Flag = false;
		pre = 0;
		if((pos=seg[0].find(tPar,pre))!=string::npos)  // search the character '|'
		{
			str = seg[0].substr(pre,pos-pre);
			if(str==chromo_Name)
				iden_Flag = true;
		}
		else if(seg[0]==chromo_Name)
			iden_Flag = true;

		if(iden_Flag==true)
		{  
			locate_Flag = true;
#pragma region read the field of "GENE_ID"
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
#pragma region the field is Gene
			if(Gene_Id!="")   // new gene
			{
				string gene_Name = Gene_Id;   // Locus: gene name
				int found = isExisting(gene,gene_Name);  // examine whether the gene already exists
				if(found>-1)  // the gene already exists: there is no need to add new genes
				{
					gene_Serial = found;  // locate the gene
				}
				else
				{		
					Gene cur_Gene;
					gene_Serial = pre_Serial;
					cur_Gene.gene_Name = gene_Name;
					cur_Gene.orientation = seg[ori_idx][0];
					gene_Serial++;
					pre_Serial = gene_Serial;
					gene.push_back(cur_Gene);  // record the gene
					cout<<"Gene "<<gene_Serial+1<<" : "<<gene_Name<<endl;
				}					
			}
#pragma endregion
#pragma region the field is Transcript
			if(seg[type_idx]=="transcript")
			{
				Trans cur_Trans;  // transcript
				sscanf(seg[3].c_str(),"%d",&(cur_Trans.start));  // starting site of the transcript
				sscanf(seg[4].c_str(),"%d",&(cur_Trans.stop));   // stopping site of the transcript
				if(seg[trs_idx-1]=="transcript_id")
					cur_Trans.trsName  = seg[trs_idx];
				else
				{
					for(int j=trs_idx-3; j<trs_idx+3; j++){
						if(seg[j]=="transcript_id"){
							cur_Trans.trsName  = seg[j+1];  // transcript name
							break;
						}
					}		
				}
				cur_Trans.orientation = seg[ori_idx][0]; // the orientation of the transcript
				cur_Trans.length = 0;  // the length of the transcript
				gene[gene_Serial].trans.push_back(cur_Trans);
			}
#pragma endregion
#pragma region the field is Exon
			else if(seg[type_idx]=="exon")  // exon
			{
				Exon_Node subexon;
				string trs_Name;
				if(seg[trs_idx-1]=="transcript_id")
					trs_Name = seg[trs_idx];
				else
				{
					for(int j=trs_idx-3; j<trs_idx+3; j++){
						if(seg[j]=="transcript_id"){
							trs_Name = seg[j+1];  // transcript name
							break;
						}
					}		
				}
				int s1, s2;
				sscanf(seg[3].c_str(),"%d",&s1);  // starting site
				sscanf(seg[4].c_str(),"%d",&s2);  // stopping site
				int found = isExisting_Trans(gene[gene_Serial].trans,trs_Name);  // examine whether the gene already exists
				if(found>-1)
				{
					int len = s2 - s1 + 1;  // length of the exon
					gene[gene_Serial].trans[found].length += len;
					gene[gene_Serial].trans[found].exon_len.push_back(len);
					gene[gene_Serial].trans[found].exon_stop.push_back(s2);
					gene[gene_Serial].trans[found].exon_start.push_back(s1);
				}
				else  // there is no corresponding transcript
				{
					Trans cur_Trans;   // transcript
					int len = s2 - s1 + 1;
					cur_Trans.trsName = trs_Name;   // transcript name
					cur_Trans.orientation = seg[ori_idx][0]; // orientation of the transcript
					cur_Trans.length = len;    // length of the transcript
					cur_Trans.exon_len.push_back(len);
					cur_Trans.exon_start.push_back(s1);
					cur_Trans.exon_stop.push_back(s2);
					gene[gene_Serial].trans.push_back(cur_Trans);
				}

#pragma region read the fields of Exon and Intron
				Gene& cur_Gene = gene[gene_Serial];  // current gene
				int subexon_Num = cur_Gene.subexon.size();  // number of subexons of the gene
				int j1 = 0, j2 = 0;
				while(j1<subexon_Num&&cur_Gene.subexon[j1].bound2<s1) // identify the left boundary of the new exon
					j1++;

				if(j1==subexon_Num) // not overlap with existing subexons
				{
					Exon_Node subexon;
					subexon.bound1 = s1; subexon.bound2 = s2; subexon.len = s2 - s1 + 1;
					cur_Gene.subexon.push_back(subexon);   // add new subexons
				}
				else  // overlap with existing subexons
				{
					vector<Exon_Node>::iterator iter;
					while(j2<subexon_Num&&cur_Gene.subexon[j2].bound2<s2)  // identify the right boundary of the new exon
						j2++;
					Exon_Node cur_Exon;
					if(j2==subexon_Num)  // beyond the boundary of the last exon
					{
#pragma region cases of exon overlap
						if(s2>cur_Gene.subexon[j2-1].bound2+2)
						{
							cur_Exon.bound1 = cur_Gene.subexon[j2-1].bound2+1;
							cur_Exon.bound2 = s2;
							cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;  // length of subexon
							cur_Gene.subexon.push_back(cur_Exon);  // insert subexon
						}

						iter = cur_Gene.subexon.begin() + j2 - 1;
						if(cur_Gene.subexon[j2-1].bound1<s1-1&&cur_Gene.subexon[j2-1].bound2>s1+1)  // two exons overlap
						{
							int pre_bound1 = cur_Gene.subexon[j2-1].bound1;  // right boundary of exon
							cur_Gene.subexon[j2-1].bound1 = s1;  // update the boundaries of subexons
							cur_Gene.subexon[j2-1].len = cur_Gene.subexon[j2-1].bound2 - s1 + 1;  // update the length of subexon
							cur_Exon.bound1 = pre_bound1; cur_Exon.bound2 = s1 - 1;  // split exon
							cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
							cur_Gene.subexon.insert(iter,cur_Exon);
						}
#pragma endregion
#pragma region cases of one exon containing another
						else
						{
							int inter_Num = j2 - j1 - 1 ;
							int l1 = cur_Gene.subexon[j1].bound1, r1 = cur_Gene.subexon[j1].bound2;
#pragma region deal with the first exon
							iter = cur_Gene.subexon.begin() + j1;	
							int insert_itr1 = j1 + 1;

							if(s1<r1-1&&s1>l1+1) // the first exon need to be split
							{
								cur_Gene.subexon[j1].bound1 = s1;   // update the boundary of the new subexon
								cur_Gene.subexon[j1].len = cur_Gene.subexon[j1].bound2 - cur_Gene.subexon[j1].bound1 + 1;  // update the length of subexon			
								cur_Exon.bound1 = l1; cur_Exon.bound2 = s1-1; // split exon
								cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
								cur_Gene.subexon.insert(iter,cur_Exon);
								insert_itr1++;
							}
							else if(s1<l1-1)  // add new starting subexon
							{
								cur_Exon.bound1 = s1; 
								cur_Exon.bound2 = l1-1; // split exon
								cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
								cur_Gene.subexon.insert(iter,cur_Exon);
								insert_itr1++;
							}

#pragma endregion
							iter = cur_Gene.subexon.begin() + insert_itr1;						
							for(int k1 = 0; k1<inter_Num; k1++)  // insert the part between the previously adjacent exons
							{
								int pre_bound2 = cur_Gene.subexon[insert_itr1-1].bound2;  // right boundary of the exon
								cur_Exon.bound1 = pre_bound2 + 1;
								cur_Exon.bound2 = cur_Gene.subexon[insert_itr1].bound1-1;
								cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;															
								if(cur_Exon.len>1)
								{
									cur_Gene.subexon.insert(iter,cur_Exon);
									insert_itr1 += 2;   // move the pointer by two steps
								}
								else
									insert_itr1++;   // not insert; move the pointer by one step
								iter = cur_Gene.subexon.begin() + insert_itr1;
							}
						}
#pragma endregion 
					}
					else
					{
#pragma region
						int l1 = cur_Gene.subexon[j1].bound1, r1 = cur_Gene.subexon[j1].bound2;  // left and right boundaries of the starting exon
						int l2 = cur_Gene.subexon[j2].bound1, r2 = cur_Gene.subexon[j2].bound2;  // left and right boundaries of the stopping exon
						vector<Exon_Node>::iterator iter_bak;
						int inter_Num = j2 - j1;
						int pre_bound2 = cur_Gene.subexon[j1].bound2;  // right boundary of the exon
						Exon_Node cur_Exon;
						int insert_itr = j1 + 1; // index of the inserted intron
						int pos = j2;   // record the index of the stopping exon
						if(s1>l1+1&&s1<r1-1)     // the first exon need to be split
						{
							cur_Gene.subexon[j1].bound1 = s1;   // update the boundary of the new subexon
							cur_Gene.subexon[j1].len = cur_Gene.subexon[j1].bound2 - cur_Gene.subexon[j1].bound1 + 1;  // update the length of subexon			
							iter = cur_Gene.subexon.begin() + j1;					
							cur_Exon.bound1 = l1; cur_Exon.bound2 = s1-1; // split exon
							cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
							cur_Gene.subexon.insert(iter,cur_Exon);
							insert_itr++;
							pos++;
							j1++;
						}
						else if(s1<l1-1) // need to insert new subexon
						{	
							cur_Exon.bound1 = s1; 
							if(s2<l2-1)  // insert new subexon in the region of intron
							{
								cur_Exon.bound2 = s2; // split exon
							}
							else
							{
								cur_Exon.bound2 = l1-1; // split exon
							}
							cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
							iter = cur_Gene.subexon.begin() + j1;  // insert new exon at the position of the first exon
							cur_Gene.subexon.insert(iter,cur_Exon);
							insert_itr++;
							pos++;
						}
						if(s2>l2+1)   // the last exon need to split
						{ 
#pragma region insert new subexon			
							iter = cur_Gene.subexon.begin() + insert_itr;
							for(int k1 = 0; k1<inter_Num; k1++)  // insert the part between the previously adjacent exons
							{
								cur_Exon.bound1 = pre_bound2 + 1;
								cur_Exon.bound2 = cur_Gene.subexon[insert_itr].bound1-1;
								cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;

								// iter++;
								if(cur_Exon.len>1)
								{
									cur_Gene.subexon.insert(iter,cur_Exon);
									pos++;
									insert_itr += 2;   // move the pointer by two steps
								}
								else
									insert_itr++;   // not insert, mover the pointer by one step
								pre_bound2 = cur_Gene.subexon[insert_itr-1].bound2;  // right boundary of exon
								iter = cur_Gene.subexon.begin() + insert_itr;
							}

							if(s2<r2-1)  // the last exon is split into two exons
							{
								iter = cur_Gene.subexon.begin() + pos;
								cur_Exon.bound1 = cur_Gene.subexon[pos].bound1; cur_Exon.bound2 = s2; // split exon
								cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;	
								cur_Gene.subexon.insert(iter,cur_Exon);  // split exon
								pos++;
								cur_Gene.subexon[pos].bound1 = s2+1;  // update the boundary of subexon
								cur_Gene.subexon[pos].len = cur_Gene.subexon[pos].bound2 - cur_Gene.subexon[pos].bound1 + 1;  // update the length of subexon
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
									insert_itr++;   // move the pointer by one step
								pre_bound2 = cur_Gene.subexon[insert_itr-1].bound2;  // right boundary of the exon
								iter = cur_Gene.subexon.begin() + insert_itr;
							}
							int r2_pre = cur_Gene.subexon[pos].bound2;
							if(s2> r2_pre + 2)
							{
								cur_Exon.bound1 = r2_pre + 1; cur_Exon.bound2 = s2; // split exon
								cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;
								cur_Gene.subexon.insert(iter,cur_Exon);  // add new subexon
							}						
#pragma endregion
						}
#pragma endregion
					}
				}
#pragma endregion
			}
			else
			{
			}
#pragma endregion

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
			gene[gene_Serial].start = gene[gene_Serial].subexon[0].bound1;     // starting site of the gene
			gene[gene_Serial].stop = gene[gene_Serial].subexon[num1-1].bound2; // stopping site of the gene
			gene[gene_Serial].subexon_Num = num1;  // number of subexons of the gene
		}
	}

#pragma region Record the starting sites and stopping sites of the transcripts and genes
	int gene_Num = gene.size();  // number of gene
	for(int i=0; i<gene_Num; i++)
	{
		int site1 = 0, site2 = 0;
		Gene& cur_Gene = gene[i];
		int trs_Num = gene[i].trans.size();
		for(int k=0; k<trs_Num; k++)
		{
			Trans& cur_Trans = gene[i].trans[k];   // current transcript
			int exon_Num = cur_Trans.exon_start.size();  // exon number of transcript
			int s1 = cur_Trans.exon_start[0];
			int s2 = cur_Trans.exon_stop[exon_Num-1];
			if(k==0)
				site1 = s1;
			else if(s1<site1)
				site1 = s1;
			if(s2>site2)
				site2 = s2;
			cur_Trans.start = s1;
			cur_Trans.stop = s2;
		}
		cur_Gene.start = site1;  // left boundary of gene
		cur_Gene.stop = site2;   // rigth boundary of gene
	}
#pragma endregion

}

// Load gene annotations
int Gene_Batch::load_Annotation_1(vector<Gene>& gene, char* filename1, string chromo_Name)
{
	ifstream infile;
	char filename[200];
	sprintf(filename,"%s\\trans_%s.bed",filename1,chromo_Name.c_str());
	infile.open(filename,ios::in);
	if(!infile)
		cout<<"Open file error!"<<endl;

	string line, word, seg[12], trs_Name;
	int gene_Serial = 0, exon_Num = 0, trs_Num = 0, trs_Cnt = 0, exon_Len = 0, exon_Start = 0;
	int bound1 = 0, bound2 = 0;
	int ori_idx = 5;
	char tPar = ',';
	string::size_type pos, pre;
	vector<int> exon_start, exon_len;
	while(getline(infile,line))
	{
		stringstream stream(line);
		if(line.length()>0&&line[0]=='@')   // line of caption
		{
			for(int k=0; k<4; k++)
				stream>>seg[k];

			string gene_Name = seg[1];   // gene name
			sscanf(seg[2].c_str(),"%d",&exon_Num);  // number of transcripts
			sscanf(seg[3].c_str(),"%d",&trs_Num);   // number of exons
			trs_Cnt += trs_Num;

			Gene cur_Gene;
			cur_Gene.gene_Name = gene_Name;			
			for(int i=0; i<trs_Num; i++)
			{	
				vector<int>().swap(exon_start);
				vector<int>().swap(exon_len);
				getline(infile,line); 
				stringstream stream1(line);
				for(int i = 0; i<12; i++)
					stream1>>seg[i];
				int s1, s2;
				sscanf(seg[1].c_str(),"%d",&s1);  // starting site
				sscanf(seg[2].c_str(),"%d",&s2);  // stopping site
				if(i==0){
					bound1 = s1; bound2 = s2;
				}
				else{
					bound1 = s1<bound1?s1:bound1;
					bound2 = s2>bound2?s2:bound2;
				}
				Trans cur_Trans;  // transcript 
				int len = s2 - s1 + 1;
				cur_Trans.trsName = seg[3];   // transcript name
				cur_Trans.orientation = seg[ori_idx][0]; // orientation of transcript
				cur_Trans.length = len;   // length of transcript				
				cur_Trans.start = s1; // starting site
				cur_Trans.stop = s2;  // stopping site
				pre = 0;
				while((pos=seg[10].find(tPar,pre))!=string::npos)  // search the character ','
				{
					string str = seg[10].substr(pre,pos-pre);
					sscanf(str.c_str(),"%d",&exon_Len);   // length of exon
					exon_len.push_back(exon_Len);		
					pre = pos + 1;
				}

				pre = 0;
				while((pos=seg[11].find(tPar,pre))!=string::npos)  // search the character ','
				{
					string str = seg[11].substr(pre,pos-pre);
					sscanf(str.c_str(),"%d",&exon_Start);   // length of exon
					exon_Start++;
					exon_start.push_back(s1+exon_Start);		
					pre = pos + 1;
				}

#pragma region merge adjacent exons
				int exon_num = exon_len.size();
				int k = 0, trs_len = 0;
				while(k<exon_num)
				{
					int j = k;
					while(j<exon_num-1&&(abs(exon_start[j]+exon_len[j]-exon_start[j+1])<2))
						j++;
					int len = exon_start[j] - exon_start[k] + exon_len[j];
					cur_Trans.exon_start.push_back(exon_start[k]);
					cur_Trans.exon_len.push_back(len);
					cur_Trans.exon_stop.push_back(exon_start[k]+len-1);
					trs_len += len;
					k = j+1;
				}
#pragma endregion
				cur_Gene.trans.push_back(cur_Trans);
			}

			cur_Gene.start = bound1;
			cur_Gene.stop = bound2;
			gene.push_back(cur_Gene);  // save the gene
			gene_Serial++;
			cout<<gene_Serial<<endl;
		}	
	}

	return trs_Cnt;
}

// Choose a reasonable number of transcripts
void Gene_Batch::select_Trans(char* path1, char* filename1, vector<string>& chr_Name, int trs1, int trs2, int exon1, int exon2)
{
	int chr_Num = chr_Name.size();  // number of chromosomes
	string::size_type pos = 0, pre = 0;
	char tPar = '\t';

	ofstream outfile;
	outfile.open(filename1,ios::out);
	if(!outfile)
		cout<<"Open file error!"<<endl;

	char filename[200];
	string line, word, seg[20], chromo_Name;
	int exon_Num = 0, trs_Num = 0;
	for(int i=0; i<chr_Num; i++)
	{
		chromo_Name = chr_Name[i];
		sprintf(filename,"%s\\trans_%s.bed",path1,chromo_Name.c_str());
		ifstream infile;
		infile.open(filename,ios::in);
		if(!infile)
			cout<<"Open file error!"<<endl;
		while(getline(infile,line))
		{
			stringstream stream(line); 
			if(line.length()>0&&line[0]=='@')   // line of caption
			{
				for(int k=0; k<6; k++)
					stream>>seg[k];

				string gene_Name = seg[1];   // gene name
				sscanf(seg[2].c_str(),"%d",&exon_Num);  // number of transcripts
				sscanf(seg[3].c_str(),"%d",&trs_Num);   // number of exons

				// choose genes with a reasonable number of exons and transcripts
				if(exon_Num>=exon1&&exon_Num<=exon2&&trs_Num>=trs1&&trs_Num<=trs2)
				{
					cout<<gene_Name<<"\t"<<exon_Num<<"\t"<<trs_Num<<endl;
					for(int i=0; i<trs_Num; i++)
					{
						getline(infile,line);
						pos = line.find(tPar,0);
						string line1 = chromo_Name + line.substr(pos,line.length()-pos);
						outfile<<line1<<endl;
					}
				}
			}	
		}

		infile.close();
	}

	outfile.close();	
}

// Abundance estimation
int Gene_Batch::expression_conver(char* filename1, char* filename2, int total_Number)
{	
	ifstream infile;
	infile.open(filename1,ios::in);
	if(!infile)
		cout<<"Open file error!"<<endl;
	ofstream outfile;
	outfile.open(filename2,ios::out);
	if(!outfile)
		cout<<"Open file error!"<<endl;

	string::size_type pos, pre;
	int base_cnt = 0, coverage = 0;
	int start = 0, stop = 0;
	char filename[200];
	double temp_Abun = 0, temp_FPKM = 0;
	int gene_idx = 9, trans_idx = 11, fpkm_idx = 15, frac_idx = 17, cov_idx = 19;
	char tPar1='"', tPar2 = ';';

	string::size_type pos0;
	string chromo_Name, method_Name, line, line_pre, word, seg[20], seg_1[20];

	bool flag_ctr = false;
	getline(infile,line_pre);
	while(flag_ctr==false)  // load the results of transcript reconstruction and abundance estimation
	{
		line = line_pre;

#pragma region load transcripts and exons
			stringstream stream(line);
			for(int k=0; k<20; k++)
				stream>>seg[k];
			sscanf(seg[3].c_str(),"%d",&start);
			sscanf(seg[4].c_str(),"%d",&stop);
			sscanf(seg[cov_idx].c_str(),"%d",&coverage);

			if(seg[2]=="transcript")
			{
				gene_idx = 9;
				trans_idx = 11;
			}
			else
			{
				gene_idx = 9;
				trans_idx = 11;
			}

			method_Name = seg[1];
			int gene_id = 0, trans_id = 0;
			
			/*int temp_Len = seg[frac_idx].length();
			string temp = seg[frac_idx].substr(1,temp_Len-3);*/
			string temp = seg[frac_idx];
			sscanf(temp.c_str(),"%lf",&temp_Abun);
			/*temp_Len = seg[fpkm_idx].length();
			temp = seg[fpkm_idx].substr(1,temp_Len-3);*/
			temp = seg[fpkm_idx];
			sscanf(temp.c_str(),"%lf",&temp_FPKM);

			bool flag = false;  // locate to the next chromosome
			int exon_Number = 0, exon_Len = 0, start = 0, stop = 0, rd_num = 0;
			vector<string> exon_Line;
			
			while(flag==false&&getline(infile,line))
			{
				stringstream stream(line);
				for(int k=0; k<20; k++)
					stream>>seg_1[k];
				if(seg_1[2]=="exon")
				{
					exon_Number++;
					exon_Line.push_back(line);
					sscanf(seg_1[3].c_str(),"%d",&start);
					sscanf(seg_1[4].c_str(),"%d",&stop);
					sscanf(seg_1[19].c_str(),"%d",&rd_num);  // coverage
					exon_Len += stop-start+1;
				}
				else
				{
					flag = true;
					line_pre = line;
					break;
				}
			}
#pragma endregion

#pragma region write to new files
			// int local_Num = rd_num*temp_Abun;
			int local_Num = coverage*temp_Abun;
			double fpkm = local_Num*1e09/(1.0*exon_Len*total_Number);			

			if(seg[2]=="transcript")
			{
				for(int k=0; k<8; k++)
					outfile<<seg[k]<<"\t";
				for(int k=8; k<11; k=k+2)
					outfile<<seg[k]<<"\t"<<tPar1<<seg[k+1]<<tPar1<<tPar2<<"\t";
				outfile<<seg[14]<<"\t"<<tPar1<<fpkm<<tPar1<<tPar2<<"\t"<<seg[16]<<"\t"<<tPar1<<seg[17]<<tPar1<<tPar2<<endl;
								
				for(int i=0; i<exon_Number; i++)
				{
					line = exon_Line[i];
					stringstream stream(line);
					for(int k=0; k<20; k++)
						stream>>seg_1[k];
					for(int k=0; k<8; k++)
						outfile<<seg_1[k]<<"\t";
					for(int k=8; k<13; k=k+2)
						outfile<<seg_1[k]<<"\t"<<tPar1<<seg_1[k+1]<<tPar1<<tPar2<<"\t";
					outfile<<seg_1[14]<<"\t"<<tPar1<<fpkm<<tPar1<<tPar2<<"\t"<<seg_1[16]<<"\t"<<tPar1<<seg_1[17]<<tPar1<<tPar2<<endl;
				}
			}
#pragma endregion

			if(flag==false) flag_ctr = true;  // reach the end of the file
	}

	infile.close();
	outfile.close();

	return 0;
}
