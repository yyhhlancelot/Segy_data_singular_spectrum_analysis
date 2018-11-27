#include "stdafx.h"
#include "fxymssa.h"
#include <string.h>
#include <io.h>
#include <fstream>
#include "CSegyData.h"
#include <cstring>
#include <math.h>
#include "armadillo"
#include "write_sgy.h"

using namespace std;
using namespace arma;

/*
	description : write segy .cpp
	
	function : 1. write reelhead && linehead && trace into file
				2. change binary to ieee/ibc format

	date : 2018-10-23
	editor : yyh
*/

write_sgy::write_sgy(void)
{
}

write_sgy::~write_sgy(void) 
{
}

void write_sgy::From_Float_32BitFloat(float *input, unsigned long *DataUint32, int SAMPLE)
//  DataUint32:IBM无符号整型的一道地震数据
//  input:十进制型的一道地震数据
//  SAMPLE:每一道的采样点数
{
	long sign;//符号
	long exp, e;//指数
	float input1;//注意，不能使用long input1;
	long fmant;//尾数
	int  j;

	for (j = 0; j < SAMPLE; j++)
	{

		sign = (input[j] < 0 ? 1 : 0);//提取符号

		input[j] = float(input[j] * pow(-1, sign));//abs(input)

		exp = 0;

		input1 = input[j];

		if (input[j] > 0)//非0值才计算
		{
			if ((int)input[j] > 0)
			{
				exp++;
				while ((int)input1 / 16 > 0)
				{
					exp++;
					input1 = input1 / 16;
				}

			}
			else
			{
				while ((int)input1 * 16 == 0)
				{
					exp--;
					input1 = input1 * 16;
				}

				exp++;
			}

			e = (exp + 64);
		}
		else//判断浮点数是否为零
			e = exp;

		fmant = (long)(input[j] * pow(16, -exp)*pow(2, 24));//尾数

		DataUint32[j] = (sign << 31) | (e << 24) | fmant;

		swap(&DataUint32[j]);//高低位字节转换
							 //FloatToIeee(&DataUint32[j]);
	}
}

void write_sgy::swap(unsigned long *tni4)
{
	*tni4 = (((*tni4 >> 24) & 0xff) | ((*tni4 & 0xff) << 24) | ((*tni4 >> 8) & 0xff00) | ((*tni4 & 0xff00) << 8));
}

void write_sgy::FloatToIeee(unsigned long *pBuf)
{
	union tay
	{
		int  a;
		char ch[4];
	} w;
	union tal
	{
		int  b;
		char ch[4];
	} t;
	w.a = *pBuf;
	t.ch[0] = w.ch[3];
	t.ch[1] = w.ch[2];
	t.ch[2] = w.ch[1];
	t.ch[3] = w.ch[0];
	*pBuf = t.b;
}

bool write_sgy::write_SGY(char* sourceFile, FILE *fid1, char* destFile, FILE *fid2, fcube cube1, unsigned short inlineNum, unsigned short xlineNum, unsigned short sampleNum, unsigned long traceNum)
{
	/*
		Earning : cube1 should be at a dimension of (inline x  xline x t)
		
		if not, please use "permute()" to change the order first. "permute()" is at "fxymssa.cpp".

		editor : yyh
	*/

	fid1 = fopen(sourceFile, "rb+");

	if (fid1 == NULL)
	{
		fprintf(stderr, "无法打开此文件!\n");
		return false;
	}

	fid2 = fopen(destFile, "wb");

	if (fid2 == NULL)
	{
		fprintf(stderr, "无法打开此文件!\n");
		return false;
	}

	char *reelHeadBuffer = new char[3200];
	float *lineHeadBuffer = new float[100];


	int readCount_reel = fread(reelHeadBuffer, 1, 3200, fid1);
	int writeCount_reel = fwrite(reelHeadBuffer, 1, 3200, fid2);
	if (writeCount_reel != 3200)
	{
		fprintf(stderr, "文件写入错误!\n");
		return false;
	}

	int readCount_line = fread(lineHeadBuffer, 4, 100, fid1);
	int writeCount_line = fwrite(lineHeadBuffer, 4, 100, fid2);
	if (writeCount_line != 100)
	{
		printf("文件写入错误！\n");
	}

	delete[]reelHeadBuffer;
	delete[]lineHeadBuffer;

	float *headBuffer = new float[60];
	float *dataStore = new float[sampleNum];
	float *dataBuffer = new float[sampleNum];
	unsigned long *Data_ibm = new unsigned long[sampleNum];




	for (int i = 0; i <= traceNum - 1; i++)
	{
		int readCount_head = fread(headBuffer, 4, 60, fid1);

		int writeCount_head = fwrite(headBuffer, 4, 60, fid2);

		if (writeCount_head != 60)
		{
			fprintf(stderr, "文件写入错误!\n");
			return false;
		}

		int readCount_store = fread(dataStore, 4, sampleNum, fid1);

		int in_index = i / xlineNum;
		int x_index = i % xlineNum;

		for (int m = 0; m <= sampleNum - 1; m++)
		{

			dataBuffer[m] = cube1(in_index, x_index, m); //put trace into buffer
		}

		From_Float_32BitFloat(dataBuffer, Data_ibm, sampleNum);

		int writeCountIBM = fwrite(Data_ibm, 4, sampleNum, fid2);
		if (writeCountIBM != sampleNum)
		{
			fprintf(stderr, "文件写入错误!\n");
			return false;
		}

		if (i % 10000 == 0)
		{
			cout << "have written" << i << "traces" << endl;
		}
	}

	delete[] headBuffer;
	delete[] dataBuffer;
	delete[] dataStore;
	delete[] Data_ibm;

	fclose(fid1);
	fclose(fid2);

	return true;
}
