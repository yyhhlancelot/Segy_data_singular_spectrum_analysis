#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include "PostStackData.h"
#include "fxymssa.h"
#include "fxymssa.cpp"
#include "CSegyData.h"

using namespace std;
using namespace ap;
using namespace arma;
/* 
	F-XY domain multichannel singular spectrum analysis (MSSA)matlab仿真的C++改写

	Warning : we solve the border effect by suitablely increasing the value of N which means the SVD's first N singular value.

	Please set N reasonablely.

	Reference:
	[1] Oropeza, V., and M. Sacchi, 2011, Simultaneous seismic data denoising and reconstruction via multichannel singular spectrum analysis, Geophysics, 76, V25-V32.
	[2] Huang, W., R. Wang, M. Zhang, and Y. Chen, 2015, Damped multichannel singular spectrum analysis for 3D random noise attenuation: SEG expanded abstracts: 85th Annual international meeting, 4714鈥�719.
	[3] Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, XX, XX-XX.
	
	input : .sgy
	output : .sgy
	editor : yyh
	first edit date : 2018-10-23
	latest edit : 2018-11-13
*/
int main()
{

	/*test*/
	//const int sLine = 1300;//start inline
	//const int eLine = 1350;//end inline

	//const int sCdp = 750;//start xline
	//const int eCdp = 800;//end xline

	//string sourFile = "J:\\Code\\project\\fxymssa3\\test_data\\YJ_X_750_800_I_1300_1350_T_1200_1340.sgy";
	//string destFile = "J:\\Code\\project\\fxymssa3\\test_data\\YJ_X_750_800_I_1300_1350_T_1200_1340_overlap_5_ave_120_20_20_N_15.sgy";
	////string destFile = "J:\\Code\\project\\fxymssa3\\test_data\\YJ_X_750_800_I_1300_1350_T_1200_1340_directly.sgy";

	////set cut size & overlap size
	//const int t_cut = 120;
	//const int in_cut = 20;
	//const int x_cut = 20;
	//const int overlap = 5;

	////default value
	//const int flow = 1;
	//const int fhigh = 124;
	//const int N = 15;
	//const float sampleInterval = 0.004;
	//const int verb = 0;

	fxymssa<cx_double> instance;
	/*instance.fxymssa_block_main(t_cut, in_cut, x_cut, sCdp, eCdp, sLine, eLine, sourFile, destFile, flow, fhigh, sampleInterval, N, verb);*/
	//instance.fxymssa_slice_main(sCdp, eCdp, sLine, eLine, sourFile, destFile, flow, fhigh, sampleInterval, N, verb);
	//instance.fxymssa_block_overlap_main(t_cut, in_cut, x_cut, overlap, sCdp, eCdp, sLine, eLine, sourFile, destFile, flow, fhigh, sampleInterval, N, verb);
	//instance.fxymssa_block_overlap_ave_main(t_cut, in_cut, x_cut, overlap, sCdp, eCdp, sLine, eLine, sourFile, destFile, flow, fhigh, sampleInterval, N, verb);
	//instance.fxymssa_func_write(sCdp, eCdp, sLine, eLine, sourFile, destFile, flow, fhigh, sampleInterval, N, verb);
	//instance.fxymssa_block_overlap_all_ave_main(t_cut, in_cut, x_cut, overlap, sCdp, eCdp, sLine, eLine, sourFile, destFile, flow, fhigh, sampleInterval, N, verb);

	/************************************** project data ********************************/

	const int sLine = 155;//start inline
	const int eLine = 1450;//end inline

	const int sCdp = 95;//start xline
	const int eCdp = 870;//end xline

    string sourFile = "J:\\Code\\project\\fxymssa3\\data\\yj\\YJ_X_95_870_I_155_1450_T_900_1340.sgy";//原始SGY数据
	string destFile = "J:\\Code\\project\\fxymssa3\\data\\yj\\YJ_X_95_870_I_155_1450_T_900_1340_block_overlap_5_ave_120_20_20.sgy";

	//set cut size & overlap size
	const int t_cut = 120;
	const int in_cut = 20;
	const int x_cut = 20;
	const int overlap = 5;

	//default value
	const int flow = 1;
	const int fhigh = 124;
	const int N = 15;
	const float sampleInterval = 0.004;
	const int verb = 0;

	//fxymssa<cx_double> instance;
	////instance.fxymssa_block_overlap_main(t_cut, in_cut, x_cut, overlap, sCdp, eCdp, sLine, eLine, sourFile, destFile, flow, fhigh, sampleInterval, N, verb);
	//instance.fxymssa_slice_main(sCdp, eCdp, sLine, eLine, sourFile, destFile, flow, fhigh, sampleInterval, N, verb);
	instance.fxymssa_block_overlap_ave_main(t_cut, in_cut, x_cut, overlap, sCdp, eCdp, sLine, eLine, sourFile, destFile, flow, fhigh, sampleInterval, N, verb);

	return 1;
}