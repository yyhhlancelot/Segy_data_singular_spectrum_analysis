#pragma once
#include "stdafx.h"
#ifndef WRITE_SGY_H
#define WRITE_SGY_H

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <fstream>
#include <math.h>

/*
	description : write segy .h
	function : 1. write reelhead && linehead && trace into file
				2. change binary to ieee/ibc format

	date : 2018-10-23
	editor : yyh
*/

using namespace std;

class write_sgy
{
public:
	write_sgy(void);
	~write_sgy(void);

public:
	//swap bytes
	void swap(unsigned long *tni4);

	void From_Float_32BitFloat(float *input, unsigned long *DataUint32, int SAMPLE);

	void FloatToIeee(unsigned long *pBuf);

	bool write_SGY(char* sourceFile, FILE *fid1, char* destFile, FILE *fid2, fcube cube1, unsigned short inlineNum, unsigned short xlineNum, unsigned short sampleNum, unsigned long traceNum);
	
};

#endif