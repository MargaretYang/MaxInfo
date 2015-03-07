
#include "Entropy.h"

///////////////////////////////////////////////////////////////////////////////
//// Result Analysis
//// Analysis of the isoform identification results
void Rd_Set::quantification_Compare(char* filename1, char* filename2, char* filename3, char* filename4, string gene_Name)
{
#pragma region Declaration of file pointer
	ifstream infile1;
	infile1.open(filename1,ios::in);
	if(!infile1){
		cout<<"Open file error!(file1)"<<endl;
		return;
	}
	ifstream infile2;
	infile2.open(filename3,ios::in);
	if(!infile2){
		cout<<"Open file error!(file2)"<<endl;
		return;
	}

	ofstream outfile;
	outfile.open(filename4,ios::app);
	if(!outfile) cout<<"Open file error!(outfile)"<<endl;

	char abun_filename[200];
	sprintf(abun_filename,"%s\\%s.txt",filename2,gene_Name.c_str());
	ifstream abun_infile;
	abun_infile.open(abun_filename,ios::in);
	if(!abun_infile)
		cout<<"Open file error!(abundance)"<<endl;
#pragma endregion

#pragma region load the information of the transcripts
	string line, line_pre, word, seg1[8];
	getline(abun_infile,line);   // read the line of caption
	vector<float> trs_AbunCtrl;  // transcript abundance
	double sum = 0;
	while(getline(abun_infile,line)&&!line.empty())
	{
		stringstream stream(line);
		for(int i=0; i<8; i++)
			stream>>seg1[i];
		int len = 0;
		double abun = 0;
		sscanf(seg1[7].c_str(),"%lf",&abun);
		trs_AbunCtrl.push_back(abun);
		sum += abun;
	}
	int trs_Num = trs_AbunCtrl.size();
	for(int i=0; i<trs_Num; i++)
		trs_AbunCtrl[i] = trs_AbunCtrl[i]/sum;
	abun_infile.close();
#pragma endregion 

	// search transcript structures from the annotations
	// load data from the BED file
	// Trans cur_Trans;
	vector<Trans> trs_SetCtrl, trs_SetCase;  // annoated transcripts and estimated transcripts
	int start = 0, stop = 0;  // starting site and stopping site of the transcript
	int exon_Number = 0;      // number of exons
	string seg[20];
	string::size_type pos, pre;
	char tPar=',';
	int exon_Len = 0, exon_Start = 0, exon_Stop = 0;  // exon length
	vector<int> exon_len;     // vector of exon lengths
	vector<int> exon_start;   // vector of starting sites of exons
	vector<int> exon_stop;    // vector of stopping sites of exons
	std::streampos pos0;

#pragma region load transcripts from the annotations
	int line_Cnt = 0;
	while(getline(infile1,line))
	{
		if(line=="")   // empty file
			return;
		exon_start.clear(); exon_stop.clear();
		exon_len.clear();

		Trans cur_Trans;
		stringstream stream(line);
		for(int k = 0; k<12; k++)
			stream>>seg[k];   // fields: 1.starting site 2.stopping site 9.exon number 10.exon length 11.starting site of exon
		sscanf(seg[1].c_str(),"%d",&start);
		sscanf(seg[2].c_str(),"%d",&stop);
		sscanf(seg[9].c_str(),"%d",&exon_Number);

		cur_Trans.start = start; cur_Trans.stop = stop;

#pragma region specify exons from the annotations
		pre = 0;
		while((pos=seg[10].find(tPar,pre))!=string::npos)  // search the splitting character ','
		{
			string str = seg[10].substr(pre,pos-pre);
			sscanf(str.c_str(),"%d",&exon_Len);     // exon length
			exon_len.push_back(exon_Len);		
			pre = pos + 1;
		}		

		pre = 0;
		while((pos=seg[11].find(tPar,pre))!=string::npos)  // search the splitting character ','
		{
			string str = seg[11].substr(pre,pos-pre);
			sscanf(str.c_str(),"%d",&exon_Start);   // exon length
			exon_start.push_back(exon_Start);		
			pre = pos + 1;
		}
#pragma endregion

		// merge similar exons
		int exon_num = exon_len.size();  // exon number
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
		cur_Trans.length = trs_len;   // transcript length
		cur_Trans.abundance = trs_AbunCtrl[line_Cnt];  // abundance of the transcript
		trs_SetCtrl.push_back(cur_Trans);   // save transcript in the set
		line_Cnt++;
	}

#pragma endregion

#pragma region reconstructed transcripts
	int trs_Serial = 0;	
	bool flag_ctr = false;
	int trs_start = 0;
	exon_start.clear();  exon_stop.clear();
	exon_len.clear();

	getline(infile2,line_pre);
	if(line_pre.empty())   // empty file
		return;
	while(flag_ctr==false)
	{
		line = line_pre;
		stringstream stream(line);
		for(int k = 0; k<16; k++)
			stream>>seg[k];   // fields: 2. type 3. starting site 4. stopping site 11. abundance estimation
		sscanf(seg[3].c_str(),"%d",&start);
		sscanf(seg[4].c_str(),"%d",&stop);
		start--; stop--;  // Todo: different methods have different ways of recording the starting sites
		int len = stop - start + 1;

		if(seg[2]=="transcript")
		{	
			Trans cur_Trans;			
			cur_Trans.start = start; cur_Trans.stop = stop;
			trs_start = start;   // starting site of the transcript
			double temp_Abun = 0, temp_FPKM = 0;
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
				sscanf(seg[11].c_str(),"%lf",&temp_Abun);
			cur_Trans.abundance = temp_Abun;
			cur_Trans.FPKM = cur_Trans.abundance*1e09/cur_Trans.length;  // convert to FPKM
			if(seg[1]=="new")
				cur_Trans.FPKM1 = cur_Trans.FPKM;
			
			trs_SetCase.push_back(cur_Trans);	
		}
		else
		{
			start -= trs_start;
			stop -= trs_start;
			exon_start.push_back(start);   // starting site of the exon
			exon_len.push_back(len);       // exon length
		}

		pos0 = infile2.tellg();    // obtain the current file pointer
		if(!getline(infile2,line_pre)) flag_ctr = true;  // reach the end of file

		if(seg[2]=="transcript"||flag_ctr==true)
		{			
			if(trs_Serial>0||flag_ctr==true)  // reach the next transcript or reach the end of file
			{
				int trs_Pre = _max(0,trs_Serial-1);
#pragma region merge adjacent exons
				int exon_num = exon_len.size();  // number of exons
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
#pragma endregion
				trs_SetCase[trs_Pre].length = trs_len;  // length of transcript			
				exon_start.clear(); exon_stop.clear();
				exon_len.clear();
			}
			trs_Serial++;
		}
	}	
#pragma endregion

#pragma region comparison of the exons
	// the stricter comparing method
	int trs_num1 = trs_SetCtrl.size();   // number of annotated transcripts
	int trs_num2 = trs_SetCase.size();   // number of reconstructed transcripts
	
	int* compare_flag1 = new int[trs_num1];   // flag for transcript comparison
	for(int i = 0; i<trs_num1; i++)
		compare_flag1[i] = -1;

	int* compare_flag2 = new int[trs_num2];   // flag for transcript comparison
	for(int i = 0; i<trs_num2; i++)
		compare_flag2[i] = -1;

	float** score = new float*[trs_num2];
	for(int i = 0; i<trs_num2; i++)
	{
		score[i] = new float[trs_num1];   // matching score(similarity) of the reconstructed transcript and the annotated transcript
		for(int j = 0; j<trs_num1; j++)
			score[i][j] = 0;
	}	

	int exon_num1 = 0, exon_num2 = 0, exon_cnt = 0;
	int start1 = 0, start2 = 0, stop1 = 0, stop2 = 0;
	int b1 = 0, b2 = 0;        // position of junction
	float len1 = 0, len2 = 0;  // exon length
	bool coin_flag = false;
	float thresh1 = 30;

#pragma region matching confidence of transcripts
	for(int k = 0; k<trs_num2; k++)
	{		
		exon_num2 = trs_SetCase[k].exon_start.size();  // number of exons of the reconstructed transcript
		exon_cnt = 0;
		start2 = trs_SetCase[k].start;  
		stop2 = trs_SetCase[k].stop;
		for(int l = 0; l<trs_num1; l++)
		{
			exon_num1 = trs_SetCtrl[l].exon_start.size();  // number of exons of the annotated transcript
			start1 = trs_SetCtrl[l].start;    
			// annotated transcription start/stop site may be different from those estimated
			stop1 = trs_SetCtrl[l].start + trs_SetCtrl[l].exon_start[exon_num1-1] + trs_SetCtrl[l].exon_len[exon_num1-1] - 1;

			if(exon_num1==exon_num2)  // the numbers of exons are the same
			{
				int s1 = start1-start2, s2 = stop1-stop2;
				if(exon_num1>1&&abs(start1-start2)<thresh1&&abs(stop1-stop2)<thresh1)  // numbers of exons are the same; starting/stopping sites of the annotated transcript and the estimated transcript are the same
				{
					// comparison of the first/last exon
					b1 = start2 + trs_SetCase[k].exon_stop[0] - start1 - trs_SetCtrl[l].exon_stop[0];  // distance between the right boundaries of the first exons
					b2 = start2 + trs_SetCase[k].exon_start[exon_num2-1] - start1 - trs_SetCtrl[l].exon_start[exon_num1-1];  // distance between the left boundaries of the last exons
					
					len1 = trs_SetCtrl[l].exon_len[0], len2 = trs_SetCtrl[l].exon_len[exon_num1-1];
					float r11 = abs(s1)*1.0/len1, r12 = abs(b1)*1.0/len1;  // rate of mismatch of the first exons 
					float r22 = abs(s2)*1.0/len2, r21 = abs(b2)*1.0/len2;  // rate of mismatch of the last exosn
			        bool ctr_Flag = true;
					if((len1<500&&r12>0.05)||(len1>500&&abs(b1)>50))
						ctr_Flag = false;
					if((len2<500&&r21>0.05)||(len2>500&&abs(b2)>50))
						ctr_Flag = false;
					
					if(ctr_Flag)
					{
						score[k][l] = 1 - (r11+r12+r21+r22)*0.5/exon_num1;						
#pragma region comparison of middle exons
						for(int j = 1; j<exon_num2-1; j++)  
						{
							len1 = trs_SetCtrl[l].exon_len[j]; 
							len2 = trs_SetCase[k].exon_len[j]; 
							b1 = start2 + trs_SetCase[k].exon_start[j] - start1 - trs_SetCtrl[l].exon_start[j];  // distance between the right boundaries of the exons
							b2 = start2 + trs_SetCase[k].exon_stop[j] - start1 - trs_SetCtrl[l].exon_stop[j];    // distance between the left boundaries of the exons

							float rate1 = abs(b1)*1.0/len1;  // rate of mismatch of the exons
							float rate2 = abs(b2)*1.0/len1;  // rate of mismatch of the exons

							if(rate1<0.05&&rate2<0.05)   // low mismatch rate
							{
								score[k][l] -= (rate1+rate2)/(2*exon_num1); // match score is reduced due to the mismatch
							}
							else
							{
								score[k][l] = 0;
								break;
							}
						}
#pragma endregion
					}
				}
				else if(exon_num2==1)  // transcript is composed of a single exon; the two single-exon transcripts are regarded as matched if the overlap rate reaches the threshold
				{
#pragma region the transcript is composed of a single exon
					float len1 = trs_SetCtrl[l].length;   // length of annotated transcript
					float len2 = trs_SetCase[k].length;   // length of estimated transcript
					float overlap_rate1 = 0, overlap_rate2 = 0, overlap_len = 0;
					float thresh2 = 0.80;
					if(start2<=stop1&&start1<=stop2)  // two transcripts overlap
					{
						int s1 = _max(start1,start2);
						int s2 = _min(stop1,stop2);
						overlap_len = s2 - s1 + 1;             // length of overlap
						overlap_rate1 = overlap_len*1.0/len1;  // length of overlap
						overlap_rate2 = overlap_len*1.0/len2;
					}

					if(overlap_rate1>thresh2&&overlap_rate2>thresh2) // overlap rate reaches the threshold
					{
						score[k][l] = overlap_rate1; // matching score
					}
#pragma endregion
				}
			}
		}
	}

#pragma endregion

#pragma region identify the matched transcript according to the matching confidence

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
#pragma endregion

	float precision = 0, recall = 0;
	float cnt1 = 0, cnt2 = 0;

	for(int i = 0; i<trs_num1; i++){
		if(compare_flag1[i]>-1)
			cnt1 += score[compare_flag1[i]][i];
	}

	for(int i = 0; i<trs_num2; i++){
		if(compare_flag2[i]>-1)
			cnt2 += score[i][compare_flag2[i]];
	}

	precision = cnt2/trs_num2;
	recall = cnt1/trs_num1;
#pragma endregion

#pragma region Output the comparison results to the files

	outfile<<"@"<<gene_Name<<"\t"<<trs_num1<<"\t"<<trs_num2<<"\t"<<exon_Num<<endl;
	// output the comparison results to the files in the format of BED
	for(int i = 0; i<trs_num1; i++)
	{
		outfile<<i+1<<"\t"<<trs_SetCtrl[i].start<<"\t"<<trs_SetCtrl[i].stop<<"\t";			
		int exon_num1 = trs_SetCtrl[i].exon_start.size();
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
		int exon_num2 = trs_SetCase[i].exon_start.size();
		for(int j = 0; j<exon_num2; j++)
			outfile<<trs_SetCase[i].exon_len[j]<<",";
		outfile<<"\t";
		for(int j = 0; j<exon_num2; j++)
			outfile<<trs_SetCase[i].exon_start[j]<<",";
		outfile<<endl;
	}

	for(int i = 0; i<trs_num1; i++)
		outfile<<i+1<<":"<<compare_flag1[i]+1<<"\t";  // output the comparison results
	outfile<<recall<<endl;

	for(int i = 0; i<trs_num2; i++)
		outfile<<i+1<<":"<<compare_flag2[i]+1<<"\t";  // output the comparison results
	outfile<<precision<<endl;

	for(int i = 0; i<trs_num1; i++)
	{
		trs_SetCtrl[i].FPKM = trs_SetCtrl[i].abundance*1e09/trs_SetCtrl[i].length;
		trs_SetCtrl[i].FPKM1 = trs_SetCtrl[i].FPKM;
		outfile<<i+1<<"\t"<<trs_SetCtrl[i].length<<"\t"<<trs_SetCtrl[i].abundance<<"\t"
		<<setiosflags(ios::fixed)<<setprecision(4)<<trs_SetCtrl[i].FPKM<<"\t"<<trs_SetCtrl[i].FPKM1<<endl;
	}
	for(int i = 0; i<trs_num2; i++)
	{
		trs_SetCase[i].FPKM = trs_SetCase[i].abundance*1e09/trs_SetCase[i].length;
		outfile<<i+1<<"\t"<<trs_SetCase[i].length<<"\t"<<trs_SetCase[i].abundance<<"\t"
		<<setiosflags(ios::fixed)<<setprecision(4)<<trs_SetCase[i].FPKM<<"\t"<<trs_SetCase[i].FPKM1<<endl;
	}
	
	p1 += precision;
	p2 += recall;

#pragma endregion

	outfile.close();
	infile2.close();
	infile1.close();
}

// Evaluate the accuracy of exon assembly
float Rd_Set::exon_assemblyEvaluation(char* filename1, char* filename2, float& pre, float& rec)
{
	ifstream infile1;
	infile1.open(filename1,ios::in);
	if(!infile1)
		cout<<"Open file1 error!"<<endl;
	ifstream infile2;
	infile2.open(filename2,ios::in);
	if(!infile2)
		cout<<"Open file2 error!"<<endl;

	string line, word, seg[3];
	vector<Exon_Node> exons_1;
	vector<Exon_Node> exons_2;
	int exon_start = 0, exon_stop = 0, exon_len = 0;

	while(std::getline(infile1,line))
	{
		stringstream stream(line);
		stream>>word;  // exon serial
		for(int k = 0; k<3; k++)
			stream>>seg[k];
		Exon_Node cur_Exon;
		sscanf(seg[0].c_str(),"%d",&(cur_Exon.bound1));
		sscanf(seg[1].c_str(),"%d",&(cur_Exon.bound2));
		cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;  // exon length
		exons_1.push_back(cur_Exon);
	}
	infile1.close();

	while(std::getline(infile2,line))
	{
		stringstream stream(line);
		stream>>word;  // exon serial
		for(int k = 0; k<3; k++)
			stream>>seg[k];
		Exon_Node cur_Exon;
		sscanf(seg[0].c_str(),"%d",&(cur_Exon.bound1));
		sscanf(seg[1].c_str(),"%d",&(cur_Exon.bound2));
		cur_Exon.len = cur_Exon.bound2 - cur_Exon.bound1 + 1;  // exon length
		exons_2.push_back(cur_Exon);
	}
	infile2.close();

	// Compare the accuracy of exon assembly
	int num1 = exons_1.size();
	int num2 = exons_2.size();
	float* precision = new float[num1];  // precision of exon assembly
	float* recall = new float[num2];     // recall of exon assembly
	int start = 0, stop = 0;
	int k1 = 0;
	float thresh1 = 10, thresh2 = 10;

#pragma region calculate precision
	for(int k = 0; k<num1; k++)
	{
		start = exons_1[k].bound1;
		stop = exons_1[k].bound2;
		int j = k1;
		while(j<num2&&(abs(start-exons_2[j].bound1)>10&&abs(stop-exons_2[j].bound2)>10))
			j++;
		if(j==num2)
		{
			k1 = 0;
			precision[k] = 0;  // matched exon not found
		}
		else
		{
			if((abs(start-exons_2[j].bound1)<10&&abs(stop-exons_2[j].bound2)<10))  // both ends of the exon are matched
				precision[k] = 1-abs(1-exons_1[k].len*1.0/exons_2[j].len);
			else
			{
				if(j<num2-1&&(abs(start-exons_2[j+1].bound1)<10&&abs(stop-exons_2[j+1].bound2)<10)) // test the next exon
				{
					precision[k] = 1-abs(1-exons_1[k].len*1.0/exons_2[j+1].len);
					k1 = j;
				}
				else
				{
					precision[k] = (1-abs(1-exons_1[k].len*1.0/exons_2[j].len))*0.5;
					k1 = _max(0,j-1);
				}

			}
		}
	}
#pragma endregion

#pragma region calculate the recall
	k1 = 0;
	for(int k = 0; k<num2; k++)
	{
		start = exons_2[k].bound1;
		stop = exons_2[k].bound2;
		int j = k1;
		while(j<num1&&(abs(start-exons_1[j].bound1)>10&&abs(stop-exons_1[j].bound2)>10))
			j++;
		if(j==num1)
		{
			k1 = 0;
			recall[k] = 0;  // matched exon not found
		}
		else
		{
			if((abs(start-exons_1[j].bound1)<10&&abs(stop-exons_1[j].bound2)<10))  // both ends of the exon are matched
				recall[k] = 1-abs(1-exons_1[j].len*1.0/exons_2[k].len);
			else
			{
				if(j<num1-1&&(abs(start-exons_1[j+1].bound1)<10&&abs(stop-exons_1[j+1].bound2)<10)) // test the next exon
				{
					recall[k] = 1-abs(1-exons_1[j+1].len*1.0/exons_2[k].len);
					k1 = j;
				}
				else
				{
					recall[k] = (1-abs(1-exons_1[j].len*1.0/exons_2[k].len))*0.5;
					k1 = _max(0,j-1);
				}

			}
		}
	}
#pragma endregion

	for(int i = 0; i<num1; i++)
		pre += precision[i];
	for(int i = 0; i<num2; i++)
		rec += recall[i];

	pre = pre*1.0/num1;
	rec = rec*1.0/num2;
	cout<<pre<<"\t"<<rec<<endl;

	float score = pre*0.5+rec*0.5;

	return score;
}

// Estimate the mean and variance of fragment length based on the sequencing data
bool Rd_Set::fragLen_Estimation(int& read_Number, float len_Mean, float len_Var)
{
	double len1 = 0;
	int valid_Cnt = 0;
	bool changed_Flag = false;
	// rd_Num = locate_stop - locate_start;  // number of valid reads
	rd_Num = rd_Vec.size();                  // number of valid reads
	locate_RdNum = rd_Num;
	mean = len_Mean;
	sigma = len_Var;
	float ori_sq = len_Var*len_Var; 
	vector<int> read_Length;

	if(read_Number<10000)  // update the variance constantly if the read number doesn't reach the threshold
	{
		for(int i=0/*locate_start*/; i<rd_Num/*locate_stop*/; i++)
		{
			// four ends of the paired reads are in the same exon
			if((rd_Vec[i].insert_Len>0)&&(rd_Vec[i].l1>=0)&&(rd_Vec[i].l1==rd_Vec[i].l2)&&(rd_Vec[i].r1==rd_Vec[i].r2)&&(rd_Vec[i].l1==rd_Vec[i].r1))
			{
				len1 += rd_Vec[i].insert_Len;
				valid_Cnt++;
				read_Length.push_back(rd_Vec[i].insert_Len);
			}
		}
	}
	else
	{
		for(int i=0/*locate_start*/; i<rd_Num/*locate_stop*/; i++)
		{
			// four ends of the paired reads are in the same exon
			if((rd_Vec[i].insert_Len>0)&&(rd_Vec[i].l1>=0)&&(rd_Vec[i].l1==rd_Vec[i].l2)&&(rd_Vec[i].r1==rd_Vec[i].r2)&&(rd_Vec[i].l1==rd_Vec[i].r1))
			{
				len1 += rd_Vec[i].insert_Len;
				valid_Cnt++;
			}
		}
	}

	float mean1 = 0;
	if(valid_Cnt>1)
	{
		changed_Flag = true;
		mean1 = len1*1.0/valid_Cnt;		
		float rate = abs(1-mean1/len_Mean);  // relative change of the mean
		float temp_mean = (len1+len_Mean*read_Number)*1.0/(valid_Cnt+read_Number);
		if(rate<0.12||(len_Mean==DEFAULT_MEAN&&rate<0.15))
		{
			mean = temp_mean;
			changed_Flag = true;
		}
#pragma region calculate the variance
		if(read_Number<10000)
		{
			float sum = 0;
			for(int i=0; i<valid_Cnt; i++)
				sum += (read_Length[i]-temp_mean)*(read_Length[i]-temp_mean);
			float new_sq = sum/(valid_Cnt-1);
			float temp_sq = read_Number*(ori_sq+(temp_mean-len_Mean)*(temp_mean-len_Mean))+valid_Cnt*new_sq;
			sigma = sqrt(temp_sq/(read_Number+valid_Cnt));
		}
#pragma endregion
	}

	if(sigma<20||sigma>50)
		sigma = 30;

	read_Number += locate_RdNum;  // count read number
	
	return changed_Flag;
}

// Return the number of reads involved in distribution estimation of fragment length
int Rd_Set::fragLen_Estimation_Initial(vector<int>& read_Length)
{
	double len1 = 0, mean1 = 0;
	int valid_Cnt = 0;
	bool changed_Flag = false;
	rd_Num = rd_Vec.size();  // number of valid reads

	for(int i=0; i<rd_Num; i++)
	{
		// reads are in the same exon
		if((rd_Vec[i].insert_Len>0)&&(rd_Vec[i].l1>=0)&&(rd_Vec[i].l1==rd_Vec[i].l2)&&(rd_Vec[i].r1==rd_Vec[i].r2)&&(rd_Vec[i].l1==rd_Vec[i].r1))
		{
			len1 += rd_Vec[i].insert_Len;
			valid_Cnt++;
			read_Length.push_back(rd_Vec[i].insert_Len);
		}
	}
	return valid_Cnt;
}

// Estimate the probability distribution of fragment length
void Rd_Set::fragLen_Dis(double mu, double sig, int lim1, int lim2)
{
	mean = mu;
	sigma = sig;
	
	int len = lim2 - lim1 + 1;
	fragLen_Pro = new float[len];  // array of the probability distribution of the fragment length

	float den = 2*sigma*sigma;
	float sum = 0;
	for(int i = lim1; i<=lim2; i++)
	{
		float temp = exp(-(i-mean)*(i-mean)/den);
		fragLen_Pro[i-lim1] = temp;
		sum += temp;
	}
	// normalization of the probability
	for(int i = 0; i< len; i++)
	{
		fragLen_Pro[i] = fragLen_Pro[i]/sum;
		// fragLen_F[i] = fragLen_Dis[i]/sum;
	}
}

// Estimate the probability distribution of fragment length based on the sequencing data using discrete intervals
void Rd_Set::fragLen_DisInterval(double mu, double sig, int lim1, int lim2)
{
	fragLen_Dis(mu, sig, lim1, lim2);

	int interval = 10;  // divide the fragment into intervals of length of 10
	int inter_Num = int((lim2-lim1+1)/10+0.99); // number of intervals

	fragLen_Interval = new float[inter_Num];
	for(int i = 0; i<inter_Num; i++)
		fragLen_Interval[i] = 0;

	for(int i = 0; i<inter_Num-1; i++)
	{
		for(int j = 0; j<10; j++)
			fragLen_Interval[i] += fragLen_Pro[i*interval+j];  // accumulate the probablity by the inervals
	}
	int tail = lim2 - (lim1+interval*(inter_Num-1)-1);  // probability in the last interval
	int i = inter_Num-1;
	for(int j = 0; j<tail; j++)
		fragLen_Interval[i] += fragLen_Pro[i*interval+j];

	ofstream outfile;
	outfile.open("D:\\Splicing Graph\\Simulator\\pro_LenInterval.txt",ios::out);
	for(int i = 0; i<inter_Num; i++)
		outfile<<lim1+i*interval<<" "<<lim1+(i+1)*interval<<" "<<fragLen_Interval[i]<<endl;  // output distribution of fragment length to the file

	outfile.close();
}

// Build cost matrix for the junctions
void Rd_Set::junc_MtxBuild_1(Graph_Trans& gph)
{
	int id_Num = idx.size();  // number of valid paths
	junc_Num = gph.junc_Num;  // number of junctions

	junc_Mtx = new float*[junc_Num];
	int rank = 0;
	float temp = 0;

	// initialization
	for(int j = 0; j<junc_Num; j++)
	{
		junc_Mtx[j] = new float[id_Num];
		for(int k = 0; k<id_Num; k++)
		{
			junc_Mtx[j][k] = 0;
		}
	}

	for(int k = 0; k<id_Num; k++)
	{
		int e_Num = path[idx[k]].exon_Number;  // number of exons
		int L_k = path[idx[k]].length - mean;  // valid length of fragment
		for(int j = 0; j<e_Num-1; j++)
		{
			rank = gph.junc_Rank[j][j+1];  // junction serial
			if(rank>0)
			{
				temp = rd_Num*rd_Len/L_k;
				junc_Mtx[rank-1][k] = temp;
			}
		}
	}
}

// Calculate the cost of the first node and the last node
void Rd_Set::cost_Linkage(Graph_Trans& gph)
{
	int id_Num = idx.size();  // number of valid transcripts
	double cost1 = 0, cost2 = 0;
	int e_Num = gph.vtx_Number-2;  // number of exons
	for(int k=0; k<id_Num; k++)
	{
		int i = idx[k];   // real serial of the transcript
		int len = path[i].exon_Number;   // transcript length
		int first = path[i].exon_Node[0];     // the first node
		int last = path[i].exon_Node[len-1];  // the last node
		cost1 = gph.first_NodeCost[first] + gph.last_NodeCost[last];
		if(len<2)
			cost2 = exp(-(double)len/(2.0*e_Num));  // cost of path length
		else
			cost2 = exp(-(double)len/(0.5*e_Num));
		path[i].link_Cost1 = cost1;
		path[i].link_Cost2 = cost2;
	}
}

// Calcuate the entropy of a single read
void Rd_Set::path_Entro1(Graph_Trans& gph, string filename)
{
	int p_Num = path.size();      // number of paths
	int rd_Num = rd_Vec.size();   // number of reads
	int exon_Num = gph.exon_Vec.size();   // number of exons
	Rd cur_Rd;   // current read

	// calcuate the probablity that the read is generated from a given transcript
	rd_Num = locate_stop - locate_start;
	ofstream outfile1;
	outfile1.open("E:\\Splicing Graph\\data\\SRR111873\\explv\\path_Locate.txt",ios::out);
	// int s2 = _min(locate_stop,locate_start+100000);
	for(int i = locate_start; i<locate_stop; i++)
	{
		path_Entro_Locate(gph,rd_Vec[i],rd_Pro[i-locate_start]);
	}
	
	outfile1.close();
	
// output the probability; select the valid transcript; initizalize the posterior probability
#pragma region 
	select_ValidTrans();        //  select valid transcripts

	// calculate the posterior probability
	rd_Pro1 = new r_Type*[rd_Num];
	int id_Num = m_idx.size();  // number of valid transcripts
	
	for(int i = 0; i<rd_Num; i++)
	{
		rd_Pro1[i] = new r_Type[p_Num];
		for(int j = 0; j<p_Num; j++)
		{
			rd_Pro1[i][j] = rd_Pro[i][j];
		}
	}
#pragma endregion

#pragma region calcute the entropies
	double sum = 0, sum1 = 0;
	double eps = 1e-09;
	for(int j = 0; j<p_Num; j++)
		path[j].entro_Abun1 = -10000;  // initialization

	for(int j = 0; j<id_Num; j++)
	{
		sum = 0; sum1 = 0;
		int k = m_idx[j];
		for(int i = 0; i<rd_Num; i++)
		{
			sum += rd_Pro[i][k]*log(rd_Pro[i][k]+eps);
			sum1 += log(rd_Pro[i][k]+eps);
		}
		path[k].entro_proAccu = sum;  // accumulate probability
		path[k].entro_Abun1 = sum1/rd_Num;
	}
#pragma endregion
 	// rd_Vec.clear();

	select_ValidTransScore();  // select valid transcrips based on the generation probability of the read

	// select transcripts when annotations are available
	select_ValidTransAnnoted(gph);
}

// Calcuate the probablity that the read is generated from a given transcript
void Rd_Set::path_Entro_Locate(Graph_Trans& gph, Rd& cur_Rd, r_Type* rd_pro)
{
	int p_Num = path.size();    // number of paths
	double temp = 0;

	int e_Num = gph.exon_Vec.size();   // number of exons
	Rd cur_Rd1;   // µ±Ç°read

	int L1 = cur_Rd.Left1, L2 = cur_Rd.Right1;  // left end
	int R1 = cur_Rd.Left2, R2 = cur_Rd.Right2;  // right end

#pragma region estimate the transcript the read belongs to 
	int k = 0;  // serial of transcript
	bool flag = false;  // flag indicting whether the transcritp contains the left read
	int l1 = cur_Rd.l1, l2 = cur_Rd.l2;
	int r1 = cur_Rd.r1, r2 = cur_Rd.r2;

	vector<int>& left_readIdx = cur_Rd.left_readIdx;
	vector<int>& right_readIdx = cur_Rd.right_readIdx;

	if(left_readIdx.size()==0)  // if the left-index set is empty, add the serial of the boundary exon
	{
		left_readIdx.push_back(l1);
		left_readIdx.push_back(l2);
	}
	if(cur_Rd.right_readIdx.size()==0)  // if the right-index set is empty, add the serial of the boundary exon
	{
		right_readIdx.push_back(r1);
		right_readIdx.push_back(r2);
	}

	int n1 = left_readIdx.size(), n2 = right_readIdx.size();
	// int exon_num1 = _max(n1,n2);
	int exon_num1 = 1+(l2>l1)+(r1>l2)+(r2>r1);

	// search the transcript which probably contains the fragment and calculate the entropy
#pragma region
	vector<int> exon_idx;  // index of the exon

	while(k<p_Num&& path[k].exon_Node[0]<=l1) // search the transcript which contains the left read
	{
		flag = false;
		int subNum = path[k].exon_Number;     // number of exons
		if(subNum>=exon_num1)
		{
			int id1 = -1, id2 = -1, id3 = -1, id4 = -1;  // the indices in the transcript of the exons in which the four points of the read fall
			int ce = left_readIdx[0];   // serial of the first exon of the left read
			int index = 0;
			exon_idx.clear();
			int j1 = 0, j2 = 0;
			for(j1 = 0; j1<subNum; j1++)
			{
				if(path[k].exon_Node[j1]==ce)
				{
					exon_idx.push_back(j1);
					while(index<n1&&left_readIdx[index]==ce)
						index++;
					if(index<n1) ce = left_readIdx[index];
					else break;
				}
			}

			ce = right_readIdx[0];    // serial of the first exon of the right read
			index = 0;
			int n_1 = exon_idx.size();  
			if(n_1>0&&left_readIdx[n1-1]==path[k].exon_Node[exon_idx[n_1-1]])
			{
				id1 = exon_idx[0]; 
				if(exon_idx.size()>1) id2 = exon_idx[n_1-1];
				else id2 = id1;
				if(r1<0)  // single-end read
					rd_pro[k] = rd_Entro(cur_Rd, path[k], gph.exon_Vec,id1,id2,-1,-1); // the transcritp contains the read; calculate the generation probability
				else      // paired-end reads
				{
					for(j2 = 0; j2<subNum; j2++)
					{
						if(path[k].exon_Node[j2]==ce)
						{
							exon_idx.push_back(j2);
							while(index<n2&&right_readIdx[index]==ce)
								index++;
							if(index<n2) ce = right_readIdx[index];
							else break;
						}
					}
					int n_2 = exon_idx.size();
					id4 = exon_idx[n_2-1];
					int r_idx = path[k].exon_Node[id4];

					if(index==n2&&right_readIdx[n2-1]==r_idx){
						id3 = exon_idx[n_1];
						rd_pro[k] = rd_Entro(cur_Rd, path[k], gph.exon_Vec,id1,id2,id3,id4); // the transcritp contains the read; calculate the generation probability	
					}
				}
			}	
		}  // end if compare exon number of the transcript
		k++;
	}//end while
#pragma endregion

#pragma endregion
}

// Select possible transcripts after calculating P(x_i|y_k)
void Rd_Set::select_ValidTrans()
{
	double thresh1 = 1e-05;
	int cnt = (int)rd_Num*0.01;
	int p_Num = path.size();
	m_idx.clear();
	for(int i = 0; i<p_Num; i++)
	{
		int k = 0;
		for(int j = 0; j<rd_Num; j++)
		{
			if(rd_Pro[j][i]>thresh1)
				k++;
		}
		if(k>cnt)
			m_idx.push_back(i);   // the transcript is valid
	}
}

// Select possible transcripts
void Rd_Set::select_ValidTransScore()
{
	vector<int> low_Abun_Id;   // index of low-abundance transcript
	vector<int> high_Abun_Id;  // index of high-abundance transcript
	double maxvalue = -100000, minvalue = 100000, sum = 0;
	int max_Id = 0, min_Id = 0;
	int m_dataDim = m_idx.size();  // number of valid transcripts
	
#pragma region calculate the maximum and minimum generation probability of the read
	for(int k = 0; k<m_dataDim; k++)
	{
		int id = m_idx[k];
		if(path[id].entro_Abun1>maxvalue)
		{
			maxvalue = path[id].entro_Abun1;  // maximum abundance
			max_Id = id;
		}
		if(path[id].entro_Abun1<minvalue)
		{
			minvalue = path[id].entro_Abun1;  // minimum abundance
			min_Id = id;
		}
		sum += path[id].entro_Abun1;
	}
	double meanvalue = sum/m_dataDim;   // avarage abundance
	double midvalue = (maxvalue+minvalue)/2;  // median abundance
#pragma endregion

#pragma region convert generation probability of read to score
	vector<bool> flag(m_dataDim);             
	float* ini_Abun = new float[m_dataDim];   // initial abundance
	float alpha = 0.5;
	float thresh = alpha*meanvalue + (1-alpha)*minvalue;  // threshold
	idx.clear();
	if(m_dataDim==1){
		idx = m_idx;
		path[idx[0]].assemble_Score = 1;
	}
	else{
		for(int k = 0; k<m_dataDim; k++)
		{
			int id = m_idx[k];
			if(path[id].entro_Abun1<thresh)
			{
				low_Abun_Id.push_back(id);    // index of low-abundance transcript
				double temp = (path[id].entro_Abun1-minvalue)/(thresh-minvalue);
				flag[k] = false;
				path[id].assemble_Score = temp*0.5;
			}
			else
			{
				high_Abun_Id.push_back(id);   // index of high-abundance transcript
				double temp = (path[id].entro_Abun1-thresh)/(maxvalue-thresh);
				flag[k] = true;
				path[id].assemble_Score = 0.5 + 0.5*temp;
				high_Abun_Id.push_back(id);   // index of high-abundance transcript
			}
		}

		float thresh1 = 0.4;
		if(path.size()<30)
			thresh1 = 0.4;
		else
			thresh1 = 0.6;
		for(int k = 0; k<m_dataDim; k++)
		{
			int id = m_idx[k];
			if(path[id].assemble_Score>thresh1)
				idx.push_back(id);
		}
	}
#pragma endregion
}

// Select transcripts when annotations are available
void Rd_Set::select_ValidTransAnnoted(Graph_Trans& gph)
{
	int anno_Num = gph.annotated_Idx.size();
	int sel_Num = idx.size();  // number of selected paths
	int cnt = 0;
	vector<bool> Flag(anno_Num,false);
	for(int j=0; j<sel_Num; j++)
	{
		if(cnt<sel_Num)
		{
			for(int i=0; i<anno_Num; i++)
			{
				if(idx[j]==gph.annotated_Idx[i])
				{
					Flag[i] = true;
					cnt++;  // annotated transcript
					if(cnt==anno_Num)
						break;
				}
			}
			if(cnt==anno_Num)
				break;
		}
	}

	if(cnt<anno_Num)  // there are unselected annotated transcripts
	{
		for(int i=0; i<anno_Num; i++)
		{
			if(Flag[i]==false)   // the transcript hasn't been considered
				idx.push_back(gph.annotated_Idx[i]);
		}
	}
}

// Identify the index of the transcript
int Rd_Set::indexTrans(Gene& single_Gene, string trs_Name)
{
	int num = single_Gene.trans.size();
	for(int i = 0; i < num; i++)
	{
		if(single_Gene.trans[i].trsName==trs_Name)
			return i;   // return gene serial
	}

	return -1;  // the gene doesn't exist

}

// Initialize the generation probability of read (p(x_i|t_k)
void Rd_Set::initial_RdPro(int trans_Num)
{
	rd_Num = locate_stop - locate_start;
	rd_Pro = new r_Type*[rd_Num];	

	for(int i = 0; i<rd_Num; i++)
	{
		rd_Pro[i] = new r_Type[trans_Num];
		for(int j = 0; j<trans_Num; j++)
		{
			rd_Pro[i][j] = 0;
		}
	}
}

// Connectivity cost of the path, used to constrain the first exon and the last exon
void Rd_Set::connect_Weight(Graph_Trans& gph)
{
	int id_Num = path.size();  // number of valid paths
	double cost = 0;
	for(int k = 0; k<id_Num; k++)
	{
		int first = path[idx[k]].exon_Node[0];  // the first exon of the path
		path[idx[k]].connect_FirstWeight = gph.first_NodeCost[first];
		int last = path[idx[k]].exon_Node[path.size()-1];  // the last exon of the path
		path[idx[k]].connect_LastWeight = gph.last_NodeCost[last];
	}
}

// Clear the generation probability array of the reads
void Rd_Set::clear_RdPro(int trs_num)
{
	int rd_num = locate_stop - locate_start;
	for(int i = 0; i<rd_num; i++)
	{
		delete[] rd_Pro[i];
		delete[] rd_Pro1[i];
	}
	delete[] rd_Pro;
	delete[] rd_Pro1;
	rd_Pro = NULL;
	rd_Pro1 = NULL; // clear the index array

	// path.clear();
	int path_Num = path.size();  // number of the paths
	for(int i=0; i<path_Num; i++)
	{
		delete[] path[i].exon_Node;  // clear the index array
		if(path[i].bound!=NULL)
		{
			delete[] path[i].bound;
		}
	}
	vector<Path>().swap(path);

	// clear the index array
	// idx.clear();
	vector<int>().swap(idx);
	vector<int>().swap(m_idx);
	vector<int>().swap(sel_Idx);
		
	// clear the probability distribution array of the fragment length
	delete[] fragLen_Pro;
	delete[] fragLen_Interval;
}

// Clear the vector
void Rd_Set::clear_RdVec()
{
	vector<Rd>().swap(rd_Vec);        // clear the array of reads
	vector<Rd>().swap(singlerd_Vec);  // clear the array of single-end reads
}

// Clear the vector
void Rd_Set::clear()
{
	// clear the probability array
	for(int i = 0; i<rd_Num; i++)
	{
		delete[] rd_Pro[i];
	}
	rd_Pro = NULL;
	rd_Pro1 = NULL;

	// clear the read array
	vector<Rd>().swap(rd_Vec);
	// clear the transcript array
	vector<Path>().swap(path);
	// clear the index array
	vector<int>().swap(idx);
}

