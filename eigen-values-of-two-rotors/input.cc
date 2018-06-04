#include <fstream>
#include <sstream>
#include <iostream>
#include "input.h"
using namespace std;

//----- INPUT PARAMETERS----------------------------

// SYSTEM
const char IO_GRID[] = "GRID";
const char IO_JROT[] = "JROT";
const char IO_SKIP[] = "SKIP";
const char IO_BONDLENGTH[] = "BONDLENGTHS";

int nSize;
int nSizeRot;
int nSkip;
double rCOM;   
double bondLength1;
double bondLength2;

void IOReadParams(const char in_file[])
{
	ifstream inf(in_file,ios::in);
	string params;
	while (inf>>params)
	{
		if (params==IO_GRID)
		{
			inf>>nSize;
		}
		else
		if (params==IO_JROT)
		{
			inf>>nSizeRot;
		}
		else
		if (params==IO_SKIP)
		{
			inf>>nSkip;
		}
		else
		if (params==IO_BONDLENGTH)
		{
			inf>>bondLength1;
			inf>>bondLength2;
			inf>>rCOM;
		}
		else
		{}

		getline(inf,params,'\n');  // skip comments at the end of the line 
	}
	inf.close();
}
