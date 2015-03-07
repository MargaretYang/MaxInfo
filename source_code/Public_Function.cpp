
#include "Read_Assignment.h"
using namespace std;

///////////////////////////////////////////////
//// Public functions
// Compare the values of elements with the type PAIR
int cmp(PAIR& x, PAIR& y)
{
	return x.second>y.second;
}

// Compare the values of elements with the type PAIR_INT
int cmp_1(PAIR_INT& x, PAIR_INT& y)
{
	return x.second>y.second;
}

// Compare the orders of the exons
int cmp_Exon(Exon_Node& exon1, Exon_Node& exon2)
{
	return exon1.bound1<exon2.bound1;
}

// Convert to valid filename
string valid_filename(string oriname)
{
	string tPar = ":";
	string tPar1 = ";";
	string::size_type pos, pos0;

	string::size_type pre; 
	pre = 0;

	if(oriname=="aux"||oriname=="com"||oriname=="com2")
		oriname = oriname + "1";
	string outname = oriname;
	while((pos=outname.find(tPar,pre))!=string::npos||
		(pos=outname.find(tPar1,pre))!=string::npos)
	{
		oriname = outname.replace(pos, tPar.length(), "-"); 
		pre = pos + 1; 
		outname = oriname;
	}

	return outname;
}

// Check whether the given gene exists
int isExisting(vector<Gene>& gene, string gene_Name)
{
	int num = gene.size();
	for(int i = 0; i < num; i++)
	{
		if(gene[i].gene_Name==gene_Name)
			return i;   // the gene exists, return gene serial
	}

	return -1;  // the gene doesn't exist
}

// Check whether the given transcript exists
int isExisting_Trans(vector<Trans>& trs, string trs_Name)
{
	int num = trs.size();
	for(int i = 0; i < num; i++)
	{
		if(trs[i].trsName==trs_Name)
			return i;   // the transcript exists, return transcript serial
	}

	return -1;  // the transcript doesn't exist
}

// Load the gene annotations
void load_GTF(char* path, char* filename, char* ext)
{
	ifstream myfile;
	char file1[200];
	sprintf(file1,"%s\\%s.%s",path,filename,ext);
	myfile.open(file1,ios::in);

	ofstream outfile;
	sprintf(file1,"%s\\%s_1.%s",path,filename,ext);
	outfile.open(file1,ios::out);	

	string line, word;
	int serial = 1;
	int cnt = 0;
	int file_Serial = 0;

	while(std::getline(myfile,line)&&cnt<500)
	{
		if(cnt%4000==0){
			file_Serial++;
			outfile.close();
			sprintf(file1,"%s\\%s_%d.%s",path,filename,file_Serial,ext);			
			outfile.open(file1,ios::out);
		}	
		outfile<<line<<endl;
		cout<<line<<endl;
		cnt++;
	}

	cout<<cnt<<endl;
	myfile.close();
	outfile.close();
}