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
	该类为F-XY domain multichannel singular spectrum analysis (MSSA)matlab仿真的C++改写
	
	Reference:   
    [1] Oropeza, V., and M. Sacchi, 2011, Simultaneous seismic data denoising and reconstruction via multichannel singular spectrum analysis, Geophysics, 76, V25-V32.
    [2] Huang, W., R. Wang, M. Zhang, and Y. Chen, 2015, Damped multichannel singular spectrum analysis for 3D random noise attenuation: SEG expanded abstracts: 85th Annual international meeting, 4714鈥�719.
    [3] Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, XX, XX-XX.

	editor : yyh
	first-edition : 2018-10-19
*/

//matlab permute
template <typename T>
Cube<T> fxymssa<T>::permute(Cube<T>& cube, const std::tuple<uword, uword, uword>& order)
{
	/*
		discription : "permute()" rearranges the dimensions of the cube so that they 
		are in the order specified by the tuple ORDER.

		the function is the same as Matlab B =  permute(A, ORDER)

		editor : yyh
	*/
	uword idx1 = std::get<0>(order);
	uword idx2 = std::get<1>(order);
	uword idx3 = std::get<2>(order);

	urowvec dimension = shape(cube);

	uword rows = dimension(idx1 - 1);
	uword cols = dimension(idx2 - 1);
	uword slis = dimension(idx3 - 1);

	Cube<T> output;
	output.zeros(rows, cols, slis);

	uword perm = idx1 * 100 + idx2 * 10 + idx3;

	switch (perm)
	{
	case 123:
	{
		output = cube; // identity
	}
	break;
	case 132:
	{
		for (int c = 0; c < cube.n_cols; ++c)
			for (int r = 0; r < cube.n_rows; ++r)
				for (int s = 0; s < cube.n_slices; ++s)
					output(r, s, c) = cube(r, c, s);
	}
	break;
	case 213:
	{
		for (int c = 0; c < cube.n_cols; ++c)
			for (int r = 0; r < cube.n_rows; ++r)
				for (int s = 0; s < cube.n_slices; ++s)
					output(c, r, s) = cube(r, c, s);
	}
	break;
	case 231:
	{
		for (int c = 0; c < cube.n_cols; ++c)
			for (int r = 0; r < cube.n_rows; ++r)
				for (int s = 0; s < cube.n_slices; ++s)
					output(c, s, r) = cube(r, c, s);
	}
	break;
	case 312:
	{
		for (int c = 0; c < cube.n_cols; ++c)
			for (int r = 0; r < cube.n_rows; ++r)
				for (int s = 0; s < cube.n_slices; ++s)
					output(s, r, c) = cube(r, c, s);
	}
	break;
	case 321:
	{
		for (int c = 0; c < cube.n_cols; ++c)
			for (int r = 0; r < cube.n_rows; ++r)
				for (int s = 0; s < cube.n_slices; ++s)
					output(s, c, r) = cube(r, c, s);
	}
	break;
	}

	return output;
}

template <typename T>
fxymssa<T>::fxymssa(void)
{
}

template <typename T>
fxymssa<T>::~fxymssa(void)
{
}

template <typename T>
cx_mat fxymssa<T>::P_H(cx_mat din, int lx, int ly)
{

	/*form block Hankel matrix*/
	
	
	int nx = din.n_rows;
	int ny = din.n_cols;
	int lxx = nx - lx + 1;
	int lyy = ny - ly + 1;
	cx_mat dout(lx * ly, lxx * lyy);

	for (int j = 0; j <= ny - 1; j++) 
	{	
		cx_colvec fi_col = din(span(0, lx - 1), span(j, j));
		cx_rowvec la_vec = din(span(lx - 1, nx - 1), span(j, j)).st();

		cx_mat r = hankel(fi_col, la_vec);
		fi_col.reset();
		la_vec.reset();
		


		if (j < ly - 1)
		{
			for (int id = 0; id <= j; id++)
			{
				int x_start = j * lx - id * lx;
				int x_end = (j + 1) * lx - id * lx - 1;
				int y_start = id * lxx;
				int y_end = lxx + id * lxx - 1;
				dout(span(x_start, x_end), span(y_start, y_end)) = r;
			} 
		}
		else 
		{
			for (int id = 0; id <= ny - j - 1; id++)
			{
				int x_start = (ly - 1) * lx - id * lx;
				int x_end = ly * lx - id * lx - 1;
				int y_start = (j + 1 - ly) * lxx + id * lxx;
				int y_end = (j - ly + 2) * lxx + id * lxx - 1;
				dout(span(x_start, x_end), span(y_start, y_end)) = r;
			}
		}
	}
	din.reset();
	return dout;
}

template <typename T>
cx_mat fxymssa<T>::P_R(cx_mat din, int N)
{	
	
	int nx = din.n_rows;
	int ny = din.n_cols;
	
	cx_mat dout;

	cx_mat U;
	vec s;
	cx_mat V;

	if (arma::svd(U, s, V, din, "std") != 1)
	{
		fprintf(stderr, "C++ implement singular value decomposition fail\n");
		
	}
	cx_vec s_temp = conv_to<cx_vec>::from(s);
	cx_mat S = arma::diagmat(s_temp);

	if (N == 1) //default
	{
		dout = U.col(N - 1) * S(N - 1, N - 1) * V.col(N - 1).t();
	}
	else
	{
		if (N <= S.n_rows)
		{
			dout = U.cols(0, N - 1) * S(span(0, N - 1), span(0, N - 1)) * (V.cols(0, N - 1).t());
		}
		else
		{
			N = S.n_rows;
			dout = U.cols(0, N - 1) * S(span(0, N - 1), span(0, N - 1)) * (V.cols(0, N - 1).t());
		}

	}

	din.reset();
	U.reset();
	s.reset();
	V.reset();

	return dout;
}

template <typename T>
cx_mat fxymssa<T>::P_A(cx_mat d_in, int nx, int ny, int lx, int ly)
{

	// averaging the block Hankel matrix to output the result
	int lxx = nx - lx + 1;
	int lyy = ny - ly + 1;

	cx_mat d_out(nx, ny);

	d_out.zeros(nx, ny);

	for (int j = 0; j <= ny - 1; j++)
	{
		if (j < ly - 1)
		{
			for (int id = 0; id <= j; id++)
			{

				d_out.col(j) = d_out.col(j) + ave_antid(d_in(span(j * lx - id * lx, (j + 1) * lx - id * lx - 1), span(id * lxx, lxx + id * lxx - 1))) / (j + 1);
			}
		}
		else
		{
			for (int id = 0; id <= ny - j - 1; id++)
			{
				d_out.col(j) = d_out.col(j) + ave_antid(d_in(span((ly - 1) * lx - id * lx, ly * lx - id * lx - 1), span((j + 1 - ly) * lxx + id * lxx, (j - ly + 2) * lxx + id * lxx - 1))) / (ny - j);
			}
		}
	}
	d_in.reset();
	
	return d_out;

}

template <typename T>
cx_colvec fxymssa<T>::ave_antid(cx_mat d_in)
{
	int row_din = d_in.n_rows;
	int col_din = d_in.n_cols;

	int n_out = row_din + col_din - 1;
	cx_colvec d_out(n_out);
	d_out.zeros();

	for (int i = 0; i <= n_out - 1; i++)
	{
		if (i < row_din - 1)
		{
			for (int id = 0; id <= i; id++)
			{
				double temp = i + 1;
				d_out(i) = d_out(i) + d_in(i - id, id) / temp;
			}
		}
		else
		{
			for (int id = 0; id <= n_out - 1 - i; id++)
			{
				double temp2 = n_out - i;
				d_out(i) = d_out(i) + d_in(row_din - id - 1, (i + 1 - row_din) + id) / temp2;
			}
		}
	}
	d_in.reset();
	return d_out;
}

template <typename T>
cx_mat fxymssa<T>::hankel(cx_colvec fi_col, cx_rowvec la_row)
{
	/*
	discription : "hankel()" form a hankel matrix whose first column is fi_col 
	and whose last row is la_row.

	editor : yyh
	*/

	int nrow = fi_col.n_elem;
	int ncol = la_row.n_elem;

	cx_mat hankel_out(nrow, ncol);

	hankel_out.col(0) = fi_col;
	hankel_out.row(nrow - 1) = la_row;
	fi_col.reset();
	la_row.reset();
	
	for (int j = 1; j <= ncol - 1; j++)
	{
		for (int i = 0; i <= nrow - 2; i++)
		{
			hankel_out(i, j) = hankel_out(i + 1, j - 1);
		}
	}
	return hankel_out;
}

template <typename T>
mat fxymssa<T>::fxymssa_slice_step(int slice_index, int sampleNum, int inlineNum, cx_mat Di, int flow, int fhigh, float dt, int N, int verb)
{
	/*
		discription : "fxymssa_slice_step()" is used in large-scaled data instance, we deal the data slice-by-slice
		
		FXY_MSSA : F-XY domain multichannel singular spectrum analysis(MSSA)
		
		IN  
		Di: input a post stack segy data whose format is defined in Armadillo
		flow : processing frequency range(lower)
		fhigh : processing frequency range(higher)
		dt : temporal sampling interval
		N : number of singular value to be preserved

		OUT
		D1 : output a segy data

		editor : yyh

	*/
	int xlineNum = 1;

	cx_mat D1(sampleNum, inlineNum);
	D1 = Di;
	

	int nf = next_pow_of_2(sampleNum);
	
	cx_mat DATA_FX(nf, inlineNum);

	DATA_FX = arma::fft(D1, nf);
	

	cx_mat DATA_FX0(nf, inlineNum);
	DATA_FX0.zeros();

	int ilow = floor(flow * dt * nf) + 1;
	int ihigh = floor(fhigh * dt * nf) + 1;

	if (ilow < 1)
	{
		ilow = 1;
	}

	if (ihigh > floor(nf / 2) + 1)
	{
		ihigh = floor(nf / 2) + 1;
	}

	int lx = floor(inlineNum / 2) + 1;
	int lxx = inlineNum - lx + 1;
	int ly = floor(xlineNum / 2) + 1;
	int lyy = xlineNum - ly + 1;

	cx_mat M(lx * ly, lxx * lyy);
	M.zeros();

	//main loop
	for (int k = ilow - 1; k <= ihigh - 1; k++) 
	{
		//if (k % 10 == 0)
		//{
		//	cout << "the " << slice_index << "th slice's" << k << "th loop and total loop is " << ihigh << endl;

		//}

		if (xlineNum == 1)
		{	
			
			cx_mat row_tmp = DATA_FX.row(k);
			M  = P_H(row_tmp.st(), lx, ly);
			
			
		}
		else
		{
			M = P_H(DATA_FX.row(k), lx, ly);

		}

		M = P_R(M, N);

		DATA_FX0.row(k) = P_A(M, inlineNum, xlineNum, lx, ly).st();
		if ((k % 5 == 0) && (verb == 1))
		{
			printf("F %d is done!\n", k);
		}
	}

	for (int k = nf / 2 + 1; k <= nf - 1; k++)
	{
		DATA_FX0.row(k) = conj(DATA_FX0.row(nf - k));
	}

	mat D1_(nf, inlineNum); //real part 

	for (int x_index = 0; x_index <= xlineNum - 1; x_index++) 
	{
		cx_mat Y = arma::ifft(DATA_FX0, nf);
		mat Y_ = real(Y);

		D1_ = Y_;
	}

	mat _D1 = D1_(span(0, sampleNum - 1), span(0, inlineNum - 1));

	return _D1;
}

template <typename T>
bool fxymssa<T>::fxymssa_slice_main(int sCdp, int eCdp, int sLine, int eLine, string _sFileName, string _dFileName, int flow, int fhigh, float dt, int N, int verb)
{
	/*
		discription : main procedure of fxymssa, play a supporting role with "fxymssa_slice_step()"
		
		notion : we have write sgy function on the end.

		FXY_MSSA : F-XY domain multichannel singular spectrum analysis(MSSA)
		
		IN

		sCdp : start Cdp
		eCdp : end Cdp
		sLine : start inline
		eLine : end inline
		_sFileName ; source segy file name
		_dFileName ; target segy file name

		flow : processing frequency range(lower)
		fhigh : processing frequency range(higher)
		dt : temporal sampling interval
		N : number of singular value to be preserved

		OUT
		D1 : output a segy data
		
		editor : yyh
		date : 2018-10-23
	*/


	///////////////////***************   read segy data ****************///////////////////

	CSegyData segyData;
	segyData.OpenSegyFile(_sFileName.c_str(), 0);


	long traceNum = segyData.GetSegyTotalTraceNumber();
	int sampleNum = segyData.GetSegySampleNumber();  //nt

	///////////////////***************   fxymssa   **************/////////////////////

	int inlineNum = eLine - sLine + 1; //nx
	int xlineNum = eCdp - sCdp + 1;; //ny

	cx_cube data_D(sampleNum, inlineNum, xlineNum);
	//cx_cube data_D1 = data_D.zeros();

	TraceHead pTraceHead;
	float *sDataBuffer = new float[sampleNum];

	// put data into data_D (inline x xline x t)
	for (long trace_num = 0; trace_num <= traceNum - 1; trace_num++)
	{
		if (trace_num % 10000 == 0)
		{
			cout << "have read " << trace_num << " traces" << endl;

		}
		segyData.ReadOneTrace(trace_num, &pTraceHead, sDataBuffer);

		int in_index = pTraceHead.RecordNum;
		int x_index = pTraceHead.CdpNum;

		for (int i = 0; i < sampleNum; i++)
		{
			data_D(i, in_index - sLine, x_index - sCdp) = sDataBuffer[i];
		}

	}
	delete[] sDataBuffer;

	//cx_cube VOL1(sampleNum, inlineNum, xlineNum);
	for (int i = 0; i <= data_D.n_slices - 1; i++) //xline direction
	{
		cout << "xline direction is dealing with " << i + 1 << "th slice with " << data_D.n_slices << " in total" << endl;
		mat D1 = fxymssa_slice_step(i, sampleNum, inlineNum, data_D.slice(i), flow, fhigh, dt, N, verb);
		data_D.slice(i) = conv_to<cx_mat>::from(D1); // VOL1.slice
		D1.reset();
	}

	cx_cube D0 = permute(data_D, std::tuple<uword, uword, uword>(1, 3, 2)); //sampleNum, xlineNum, inlineNum
	data_D.reset();
	//cube tmpVOL = permute(VOL1, std::tuple<uword, uword, uword>(2, 3, 1)); // inline x xline x t

	//cx_cube D0(sampleNum, xlineNum, inlineNum);
	

	for (int i = 0; i <= D0.n_slices - 1; i++) //xline direction
	{
		cout << "inline direction is dealing with " << i + 1 << "th slice with " << D0.n_slices << " in total" << endl;
		cx_mat tmp = conv_to<cx_mat>::from(D0.slice(i));
		mat D2 = fxymssa_slice_step(i, sampleNum, xlineNum, tmp, flow, fhigh, dt, N, verb);
		D0.slice(i) = conv_to<cx_mat>::from(D2); //VOL.slice
		D2.reset();
	}
	//D0.reset();
	cx_cube VOL_ =permute(D0, std::tuple<uword, uword, uword>(3, 2, 1));
	D0.reset();
	///////////*********  write segy file ************////////////

	// format change
	const char* const_s = nullptr;
	const_s = _sFileName.c_str();
	char* s = nullptr;
	s = const_cast<char*>(const_s);

	const char* const_d = nullptr;
	const_d = _dFileName.c_str();
	char* d = nullptr;
	d = const_cast<char*>(const_d);


	// write reelhead & linehead & tracehead & tracedata
	write_sgy write_instance;

	FILE *fid1 = NULL;
	FILE *fid2 = NULL;
	fcube VOL_temp = conv_to<fcube>::from(VOL_);

	write_instance.write_SGY(s, fid1, d, fid2, VOL_temp, inlineNum, xlineNum, sampleNum, traceNum);

	return true;
	
}

template <typename T>
bool fxymssa<T>::fxymssa_func_write(const int sCdp, const int eCdp, const int sLine, const int eLine, string _sFileName, string _dFileName, int flow, int fhigh, float dt, int N, int verb)
{
	/*
	discription : this function can only be used directly in small sized test data

	the function is exactly the same as "fxymssa()" in matlab.

	the diffrence is that we have write sgy code attached to it.

	FXY_MSSA : F-XY domain multichannel singular spectrum analysis(MSSA)

	IN

	sCdp : start Cdp
	eCdp : end Cdp
	sLine : start inline
	eLine : end inline
	_sFileName ; source segy file name
	_dFileName ; target segy file name

	flow : processing frequency range(lower)
	fhigh : processing frequency range(higher)
	dt : temporal sampling interval
	N : number of singular value to be preserved

	OUT

	D1 : output a segy data

	editor : yyh
	date : 2018-10-23
	*/

	///////////////////***************   read segy data ****************///////////////////

	CSegyData segyData;
	segyData.OpenSegyFile(_sFileName.c_str(), 0);

	
	long traceNum = segyData.GetSegyTotalTraceNumber();
	int sampleNum = segyData.GetSegySampleNumber();  //nt


	///////////////////***************   fxymssa   **************/////////////////////

	int inlineNum = eLine - sLine + 1; //nx
	int xlineNum = eCdp - sCdp + 1;; //ny

	cx_cube data_D(inlineNum, xlineNum, sampleNum);
	//cx_cube data_D1 = data_D.zeros();
	
	TraceHead pTraceHead;
	float *sDataBuffer = new float[sampleNum];

	// put data into data_D (inline x xline x t)
	for (long trace_num = 0; trace_num <= traceNum - 1; trace_num++)
	{
		if (trace_num % 10000 == 0)
		{
			cout << "have read " << trace_num << " traces" << endl;

		}
		segyData.ReadOneTrace(trace_num, &pTraceHead, sDataBuffer);
		
		int in_index = pTraceHead.RecordNum;
		int x_index = pTraceHead.CdpNum;
		
		for (int i = 0; i < sampleNum; i++)
		{
			data_D(in_index - sLine, x_index - sCdp, i) = sDataBuffer[i];
		}
	
	}
	delete[] sDataBuffer;

	int nf = next_pow_of_2(sampleNum);
	cx_cube data_FX(inlineNum, xlineNum, nf);
	

	for (int x_index = 0; x_index <= xlineNum - 1; x_index++)
	{
		//if (x_index % 10 == 0)
		//{
		//	cout << "have dealt " << x_index << " cdps" << endl;

		//}
		cx_mat temp = data_D.tube(span(0, inlineNum - 1), span(x_index, x_index)); // inline x t
		cx_mat Y = arma::fft(temp.st(), nf); 
		cx_mat Y_temp = Y.st();
		data_FX.tube(span(0, inlineNum - 1), span(x_index, x_index)) = Y_temp;

		temp.reset();
		Y.reset();
		Y_temp.reset();
			
	}

	data_D.reset(); //memory free

	int ilow = floor(flow * dt * nf) + 1;
	int ihigh = floor(fhigh * dt * nf) + 1;

	if (ilow < 1)
	{
		ilow = 1;
	}

	if (ihigh > floor(nf / 2) + 1)
	{
		ihigh = floor(nf / 2) + 1;
	}

	int lx = floor(inlineNum / 2) + 1;
	int lxx = inlineNum - lx + 1;
	int ly = floor(xlineNum / 2) + 1;
	int lyy = xlineNum - ly + 1;

	cx_mat M(lx * ly, lxx * lyy);
	M.zeros();

	cx_cube data_FX0(inlineNum, xlineNum, nf);
	data_FX0.zeros();
	//////////*********** main loop ****************/////////

	for (int k = ilow - 1; k <= ihigh - 1; k++)
	{
		cout << "process is on the " << k << "th waterline with " << ihigh - ilow << " in total" << endl;

		if (xlineNum == 1)
		{
			M = P_H(data_FX.slice(k).st(), lx, ly); // for data_D's dimension is (inlineNum, xlineNum, sampleNum)
		}
		else
		{
			M = P_H(data_FX.slice(k), lx, ly);

		}

		M = P_R(M, N);

		data_FX0.slice(k) = P_A(M, inlineNum, xlineNum, lx, ly);

		if ((k % 5 == 0) && (verb == 1))
		{
			printf("F %d is done!\n", k);
		}
	}
	data_FX.reset(); // memory free
	M.reset();
	
	for (int k = nf / 2 + 1; k <= nf - 1; k++)
	{
		data_FX0.slice(k) = conj(data_FX0.slice(nf - k));
	}

	cube D1_(inlineNum, xlineNum, nf); //real part 

	for (int x_index = 0; x_index <= xlineNum - 1; x_index++)
	{
		cx_mat temp = data_FX0.tube(span(0, inlineNum - 1), span(x_index, x_index));
		cx_mat Y = arma::ifft(temp.st(), nf);
		Y = Y.st();
		mat Y_ = real(Y);
		
		D1_.tube(span(0, inlineNum - 1), span(x_index, x_index)) = Y_;

	}
	data_FX0.reset();
	cube D1_sub(inlineNum, xlineNum, sampleNum); // t: nf -> sampleNum

	D1_sub = D1_.subcube(0, 0, 0, inlineNum - 1, xlineNum - 1, sampleNum - 1);
	D1_.reset(); // memory free

	///////////*********  write segy file ************////////////

	// format change
	const char* const_s = nullptr;
	const_s = _sFileName.c_str();
	char* s = nullptr;
	s = const_cast<char*>(const_s);

	const char* const_d = nullptr;
	const_d = _dFileName.c_str();
	char* d = nullptr;
	d = const_cast<char*>(const_d);


	// write reelhead & linehead & tracehead & tracedata
	write_sgy write_instance;

	FILE *fid1 = NULL;
	FILE *fid2 = NULL;
	fcube D1_sub_temp = conv_to<fcube>::from(D1_sub);

	write_instance.write_SGY(s, fid1, d, fid2, D1_sub_temp, inlineNum, xlineNum, sampleNum, traceNum);

	return true;
}

template <typename T>
bool fxymssa<T>::fxymssa_block_main(const int t_cut, const int in_cut, const int x_cut, const int sCdp, const int eCdp, const int sLine, const int eLine, string _sFileName, string _dFileName, int flow, int fhigh, float dt, int N, int verb)
{
	/*
	discription : deal the segy data block-by-block without overlap, faster than the slice-by-slice method,
	however we have to deal with the problem of boundary effect, which maybe cause a certain amount
	of error.

	FXY_MSSA : F-XY domain multichannel singular spectrum analysis(MSSA)
	
	IN
	t_cut : cut size in time dimension
	in_cut : cut size in inline dimension
	x_cut : cut size in xline dimension
	sCdp : start Cdp
	eCdp : end Cdp
	sLine : start inline
	eLine : end inline
	_sFileName ; source segy file name
	_dFileName ; target segy file name

	flow : processing frequency range(lower)
	fhigh : processing frequency range(higher)
	dt : temporal sampling interval
	N : number of singular value to be preserved

	OUT

	D1 : output a segy data

	editor : yyh

	date : 2018-10-31
	*/

	///////////////////***************   read segy data ****************///////////////////

	CSegyData segyData;
	segyData.OpenSegyFile(_sFileName.c_str(), 0);

	long traceNum = segyData.GetSegyTotalTraceNumber();
	int sampleNum = segyData.GetSegySampleNumber();  //nt

	int inlineNum = eLine - sLine + 1; //nx
	int xlineNum = eCdp - sCdp + 1;; //ny

	cube data_D(sampleNum, inlineNum, xlineNum);
	//cx_cube data_D1 = data_D.zeros();

	TraceHead pTraceHead;
	float *sDataBuffer = new float[sampleNum];

	// put data into data_D (inline x xline x t)
	for (long trace_num = 0; trace_num <= traceNum - 1; trace_num++)
	{
		if (trace_num % 10000 == 0)
		{
			cout << "have read " << trace_num << " traces" << endl;
		}
		segyData.ReadOneTrace(trace_num, &pTraceHead, sDataBuffer);

		int in_index = pTraceHead.RecordNum;
		int x_index = pTraceHead.CdpNum;

		for (int i = 0; i < sampleNum; i++)
		{
			data_D(i, in_index - sLine, x_index - sCdp) = sDataBuffer[i];
		}

	}
	delete[] sDataBuffer;
	
	///////////////////***************   fxymssa   **************/////////////////////


	// CUT METHOD ONLY
	const int t_parts = sampleNum / t_cut; //7
	const int t_remain = sampleNum % t_cut; //21
	const int in_parts = inlineNum / in_cut; //12
	const int in_remain = inlineNum % in_cut; //96
	const int x_parts = xlineNum / x_cut; //9
	const int x_remain = xlineNum % x_cut; //46

	cx_cube block_remain;
	cx_cube block(t_cut, in_cut, x_cut);

	////////////////////////// main_loop
	int counter = 0;

	// CUT METHOD ONLY
	for (int i = 0; i <= in_parts; i++)
	{
		for (int j = 0; j <= x_parts; j++)
		{
			for (int k = 0; k <= t_parts; k++) //special situation, when cut is devided by Num with remaineder
			{
				counter++;

				if ((k == t_parts) || (j == x_parts) || (i == in_parts))
				{
					cout << "the " << counter << "th block is being processed with " << (in_parts + 1) * (x_parts + 1) * (t_parts + 1) << " in total..." << endl;

					if ((k == t_parts) && (i != in_parts) && (j != x_parts)) //the first situation
					{
						block_remain = conv_to<cx_cube>::from(data_D.subcube(k * t_cut, i * in_cut, j * x_cut, k * t_cut + t_remain - 1, (i + 1) * in_cut - 1, (j + 1) * x_cut - 1));

						cube block_dealt = fxymssa_block_func(t_remain, in_cut, x_cut, block_remain, flow, fhigh, dt, N, verb);

						data_D(span(k * t_cut, k * t_cut + t_remain - 1), span(i * in_cut, (i + 1) * in_cut - 1), span(j * x_cut, (j + 1) * x_cut - 1)) = block_dealt;

						block_dealt.reset();
					}

					else if ((j == x_parts) && (k != t_parts) && (i != in_parts)) //the second situation
					{
						block_remain = conv_to<cx_cube>::from(data_D.subcube(k * t_cut, i * in_cut, j * x_cut, (k + 1) * t_cut - 1, (i + 1) * in_cut - 1, j * x_cut + x_remain - 1));

						cube block_dealt = fxymssa_block_func(t_cut, in_cut, x_remain, block_remain, flow, fhigh, dt, N, verb);

						data_D(span(k * t_cut, (k + 1) * t_cut - 1), span(i * in_cut, (i + 1) * in_cut - 1), span(j * x_cut, j * x_cut + x_remain - 1)) = block_dealt;

						block_dealt.reset();
					}

					else if ((k == t_parts) && (j == x_parts) && (i != in_parts)) //the third situation
					{
						block_remain = conv_to<cx_cube>::from(data_D.subcube(k * t_cut, i * in_cut, j * x_cut, k * t_cut + t_remain - 1, (i + 1) * in_cut - 1, j * x_cut + x_remain - 1));

						cube block_dealt = fxymssa_block_func(t_remain, in_cut, x_remain, block_remain, flow, fhigh, dt, N, verb);

						data_D(span(k * t_cut, k * t_cut + t_remain - 1), span(i * in_cut, (i + 1) * in_cut - 1), span(j * x_cut, j * x_cut + x_remain - 1)) = block_dealt;

						block_dealt.reset();
					}

					else if ((i == in_parts) && (k != t_parts) && (j != x_parts)) //the forth situation
					{
						block_remain = conv_to<cx_cube>::from(data_D.subcube(k * t_cut, i * in_cut, j * x_cut, (k + 1) * t_cut - 1, i * in_cut + in_remain - 1, (j + 1) * x_cut - 1));

						cube block_dealt = fxymssa_block_func(t_cut, in_remain, x_cut, block_remain, flow, fhigh, dt, N, verb);

						data_D(span(k * t_cut, (k + 1) * t_cut - 1), span(i * in_cut, i * in_cut + in_remain - 1), span(j * x_cut, (j + 1) * x_cut - 1)) = block_dealt;

						block_dealt.reset();
					}

					else if ((i == in_parts) && (k == t_parts) && (j != x_parts)) //the fifth situation
					{
						block_remain = conv_to<cx_cube>::from(data_D.subcube(k * t_cut, i * in_cut, j * x_cut, k * t_cut + t_remain - 1, i * in_cut + in_remain - 1, (j + 1) * x_cut - 1));

						cube block_dealt = fxymssa_block_func(t_remain, in_remain, x_cut, block_remain, flow, fhigh, dt, N, verb);

						data_D(span(k * t_cut, k * t_cut + t_remain - 1), span(i * in_cut, i * in_cut + in_remain - 1), span(j * x_cut, (j + 1) * x_cut - 1)) = block_dealt;

						block_dealt.reset();
					}

					else if ((i == in_parts) && (j == x_parts) && (k != t_parts)) //the sixth situation
					{
						block_remain = conv_to<cx_cube>::from(data_D.subcube(k * t_cut, i * in_cut, j * x_cut, (k + 1) * t_cut - 1, i * in_cut + in_remain - 1, j * x_cut + x_remain - 1));

						cube block_dealt = fxymssa_block_func(t_cut, in_remain, x_remain, block_remain, flow, fhigh, dt, N, verb);

						data_D(span(k * t_cut, (k + 1) * t_cut - 1), span(i * in_cut, i * in_cut + in_remain - 1), span(j * x_cut, j * x_cut + x_remain - 1)) = block_dealt;

						block_dealt.reset();
					}

					else if ((i == in_parts) && (j == x_parts) && (k == t_parts)) //the seventh situation
					{
						block_remain = conv_to<cx_cube>::from(data_D.subcube(k * t_cut, i * in_cut, j * x_cut, k * t_cut + t_remain - 1, i * in_cut + in_remain - 1, j * x_cut + x_remain - 1));

						cube block_dealt = fxymssa_block_func(t_remain, in_remain, x_remain, block_remain, flow, fhigh, dt, N, verb);

						data_D(span(k * t_cut, k * t_cut + t_remain - 1), span(i * in_cut, i * in_cut + in_remain - 1), span(j * x_cut, j * x_cut + x_remain - 1)) = block_dealt;

						block_dealt.reset();
					}

				}

				else // common situation
				{
					cout << "the " << counter << "th block is being processed with " << (in_parts + 1) * (x_parts + 1) * (t_parts + 1) << " in total..." << endl;

					block = conv_to<cx_cube>::from(data_D.subcube(k * t_cut, i * in_cut, j * x_cut, (k + 1) * t_cut - 1, (i + 1) * in_cut - 1, (j + 1) * x_cut - 1));

					cube block_dealt = fxymssa_block_func(t_cut, in_cut, x_cut, block, flow, fhigh, dt, N, verb);

					data_D(span(k * t_cut, (k + 1) * t_cut - 1), span(i * in_cut, (i + 1) * in_cut - 1), span(j * x_cut, (j + 1) * x_cut - 1)) = block_dealt;

					block_dealt.reset();
				}

			}
		}
	}

	block.reset();

	//remember to permute the cube, for you have to put a (inline x xline x t) size cube into the write segy function. 
	fxymssa<double> change;
	fcube data_D_permute = conv_to<fcube>::from(change.permute(data_D, std::tuple<uword, uword, uword>(2, 3, 1)));

	data_D.reset();

	///////////*********  write segy file ************////////////

	// format change
	const char* const_s = nullptr;
	const_s = _sFileName.c_str();
	char* s = nullptr;
	s = const_cast<char*>(const_s);

	const char* const_d = nullptr;
	const_d = _dFileName.c_str();
	char* d = nullptr;
	d = const_cast<char*>(const_d);


	// write reelhead & linehead & tracehead & tracedata
	write_sgy write_instance;

	FILE *fid1 = NULL;
	FILE *fid2 = NULL;
	//fcube data_D_permute_temp = conv_to<fcube>::from(data_D_permute);

	write_instance.write_SGY(s, fid1, d, fid2, data_D_permute, inlineNum, xlineNum, sampleNum, traceNum);

	return true;
}

template <typename T>
bool fxymssa<T>::fxymssa_block_overlap_ave_main(const int t_cut, const int in_cut, const int x_cut, const int overlap, const int sCdp, const int eCdp, const int sLine, const int eLine, string _sFileName, string _dFileName, int flow, int fhigh, float dt, int N, int verb)
{
	/*
		discription : deal the segy data block-by-block with overlap, faster than the slice-by-slice method,

		by computing the average of the overlapped part to make the spectrum more continuous.

		FXY_MSSA : F-XY domain multichannel singular spectrum analysis(MSSA)
	
		IN

		t_cut : cut size in time dimension
		in_cut : cut size in inline dimension
		x_cut : cut size in xline dimension
		sCdp : start Cdp
		eCdp : end Cdp
		sLine : start inline
		eLine : end inline
		_sFileName ; source segy file name
		_dFileName ; target segy file name

		flow : processing frequency range(lower)
		fhigh : processing frequency range(higher)
		dt : temporal sampling interval
		N : number of singular value to be preserved

		OUT

		D1 : output a segy data
		
		editor : yyh

		date : 2018-11-06
	*/

	///////////////////***************   read segy data ****************///////////////////

	CSegyData segyData;
	segyData.OpenSegyFile(_sFileName.c_str(), 0);

	long traceNum = segyData.GetSegyTotalTraceNumber();
	int sampleNum = segyData.GetSegySampleNumber();  //nt

	int inlineNum = eLine - sLine + 1; //nx
	int xlineNum = eCdp - sCdp + 1;; //ny

	cube data_D(sampleNum, inlineNum, xlineNum);
	//cx_cube data_D1 = data_D.zeros();

	TraceHead pTraceHead;
	float *sDataBuffer = new float[sampleNum];

	// put data into data_D (inline x xline x t)
	for (long trace_num = 0; trace_num <= traceNum - 1; trace_num++)
	{
		if (trace_num % 10000 == 0)
		{
			cout << "have read " << trace_num << " traces" << endl;
		}
		segyData.ReadOneTrace(trace_num, &pTraceHead, sDataBuffer);

		int in_index = pTraceHead.RecordNum;
		int x_index = pTraceHead.CdpNum;

		for (int i = 0; i < sampleNum; i++)
		{
			data_D(i, in_index - sLine, x_index - sCdp) = sDataBuffer[i];
		}

	}
	delete[] sDataBuffer;

	//////////////////////////********************   fxymssa   **********************/////////////////////////////

	//////////////// block-by-block

	// CUT & OVRELAP & AVERAGE METHOD
	const int in_parts = inlineNum / (in_cut - overlap);
	const int in_remain = inlineNum % (in_cut - overlap);
	const int x_parts = xlineNum / (x_cut - overlap);
	const int x_remain = xlineNum % (x_cut - overlap);
	const int t_parts = sampleNum / (t_cut - overlap);
	const int t_remain = sampleNum % (t_cut - overlap);

	cx_cube block_old;
	cx_cube block_new;

	///////// main_loop
	int counter = 0;

	// CUT & OVRELAP & AVERAGE METHOD

	// coordinate
	int t_s_old, t_e_old, t_s_new, t_e_new; //t_s_old means the old block 's start point in dimension of time, and t_e_new means the new block's end point in dimension of time

	int in_s_old, in_e_old, in_s_new, in_e_new;

	int x_s_old, x_e_old, x_s_new, x_e_new;

	// initialize in_line coordinate
	in_s_old = 0; in_e_old = in_cut - 1;
	in_s_new = 0; in_e_new = in_cut - 1;

	unsigned short flag_in = 0;

	while (in_e_new != inlineNum - 1)
	{
		if (in_e_new > inlineNum - 1)
		{
			in_e_new = inlineNum - 1;
		}
		
		// initialize x_line coordinates
		x_s_old = 0; x_e_old = x_cut - 1;
		x_s_new = 0; x_e_new = x_cut - 1;

		unsigned short flag_x = 0;

		while (x_e_new != xlineNum - 1)
		{
			if (x_e_new > xlineNum - 1)
			{
				x_e_new = xlineNum - 1;
			}

			// initialize t coordinate
			t_s_old = 0; t_e_old = t_cut - 1;
			t_s_new = 0; t_e_new = t_cut - 1;    


			unsigned short flag_t = 0; 

			while (t_e_new != sampleNum - 1)
			{
				if (t_e_new > sampleNum - 1)
				{
					t_e_new = sampleNum - 1;
				}

				counter++;

				cout << "the " << counter << "th block is being processed with " << (in_parts + 1) * (x_parts + 1) * (t_parts + 1) << " in total..." << endl;

				if ((flag_t == 0) && (flag_x == 0) && (flag_in == 0)) // at the begining of three dimension, not overlapped yet
				{
					cx_cube block = conv_to<cx_cube>::from(data_D.subcube(t_s_new, in_s_new, x_s_new, t_e_new, in_e_new, x_e_new));

					cube block_dealt = fxymssa_block_func(t_e_new - t_s_new + 1, in_e_new - in_s_new + 1, x_e_new - x_s_new + 1,  block, flow, fhigh, dt, N, verb);

					data_D(span(t_s_new, t_e_new), span(in_s_new, in_e_new), span(x_s_new, x_e_new)) = block_dealt;

					block.reset();

					block_dealt.reset();
					
				}

				else // overlapped situation
				{
					// extract and store overlapped old data 
					cx_cube block_extract_overlap_old = conv_to<cx_cube>::from(data_D.subcube(t_s_new, in_s_new, x_s_new, t_e_old, in_e_old, x_e_old));

					// update new block
					cx_cube block_new = conv_to<cx_cube>::from(data_D.subcube(t_s_new, in_s_new, x_s_new, t_e_new, in_e_new, x_e_new));

					cube block_dealt_new = fxymssa_block_func(t_e_new - t_s_new + 1, in_e_new - in_s_new + 1, x_e_new - x_s_new + 1, block_new, flow, fhigh, dt, N, verb);
					
					data_D(span(t_s_new, t_e_new), span(in_s_new, in_e_new), span(x_s_new, x_e_new)) = block_dealt_new; 

					// extract and store overlapped new data 
					cx_cube block_extract_overlap_new = conv_to<cx_cube>::from(data_D.subcube(t_s_new, in_s_new, x_s_new, t_e_old, in_e_old, x_e_old))
						;
					// compute the average of the old part and the new part
					cube block_overlap_final = (conv_to<cube>::from(block_extract_overlap_old) + conv_to<cube>::from(block_extract_overlap_new)) / 2;

					// update overlapped part
					data_D(span(t_s_new, t_e_old), span(in_s_new, in_e_old), span(x_s_new, x_e_old)) = block_overlap_final;

					block_dealt_new.reset();
					block_extract_overlap_old.reset();
					block_extract_overlap_new.reset();
					block_overlap_final.reset();

				}

				if (t_e_new == sampleNum - 1)
				{
					break;
				}

				// update time new coordinate

					t_s_old = t_s_new; t_e_old = t_e_new;

					t_s_new = t_e_old - overlap + 1;

					t_e_new = t_s_new + t_cut - 1;
				
					flag_t++;
			}

			if (x_e_new == xlineNum - 1)
			{
				break;
			}

			//update xline new coordinate
			x_s_old = x_s_new; x_e_old = x_e_new;

			x_s_new = x_e_old - overlap + 1;

			x_e_new = x_s_new + x_cut - 1;

			flag_x++;
		}
		
		if (in_e_new == inlineNum - 1)
		{
			break;
		}

		// update inline new coordinate
		in_s_old = in_s_new; in_e_old = in_e_new;

		in_s_new = in_e_old - overlap + 1;

		in_e_new = in_s_new + in_cut - 1;

		flag_in++;

	}

	//remember to permute the cube, for you have to put a (inline x xline x t) size cube into the write segy function. 
	fxymssa<double> change;
	fcube data_D_permute = conv_to<fcube>::from(change.permute(data_D, std::tuple<uword, uword, uword>(2, 3, 1)));

	data_D.reset();

	///////////*********  write segy file ************////////////

	// format change
	const char* const_s = nullptr;
	const_s = _sFileName.c_str();
	char* s = nullptr;
	s = const_cast<char*>(const_s);

	const char* const_d = nullptr;
	const_d = _dFileName.c_str();
	char* d = nullptr;
	d = const_cast<char*>(const_d);


	// write reelhead & linehead & tracehead & tracedata
	write_sgy write_instance;

	FILE *fid1 = NULL;
	FILE *fid2 = NULL;
	//fcube data_D_permute_temp = conv_to<fcube>::from(data_D_permute);

	write_instance.write_SGY(s, fid1, d, fid2, data_D_permute, inlineNum, xlineNum, sampleNum, traceNum);

	return true;

}



template <typename T>
bool fxymssa<T>::fxymssa_block_overlap_main(const int t_cut, const int in_cut, const int x_cut, const int overlap, const int sCdp, const int eCdp, const int sLine, const int eLine, string _sFileName, string _dFileName, int flow, int fhigh, float dt, int N, int verb)
{
	/*
		discription : deal the segy data block-by-block with overlap, faster than the slice-by-slice method,

		weakness ; still exists the problem of boundary effct.

		FXY_MSSA : F-XY domain multichannel singular spectrum analysis(MSSA)
		
		IN

		t_cut : cut size in time dimension
		in_cut : cut size in inline dimension
		x_cut : cut size in xline dimension
		sCdp : start Cdp
		eCdp : end Cdp
		sLine : start inline
		eLine : end inline
		_sFileName ; source segy file name
		_dFileName ; target segy file name

		flow : processing frequency range(lower)
		fhigh : processing frequency range(higher)
		dt : temporal sampling interval
		N : number of singular value to be preserved

		OUT

		D1 : output a segy data

		editor : yyh

		date : 2018-10-31
	*/

	///////////////////***************   read segy data ****************///////////////////

	CSegyData segyData;
	segyData.OpenSegyFile(_sFileName.c_str(), 0);

	long traceNum = segyData.GetSegyTotalTraceNumber();
	int sampleNum = segyData.GetSegySampleNumber();  //nt

	int inlineNum = eLine - sLine + 1; //nx
	int xlineNum = eCdp - sCdp + 1;; //ny

	cube data_D(sampleNum, inlineNum, xlineNum);
	//cx_cube data_D1 = data_D.zeros();

	TraceHead pTraceHead;
	float *sDataBuffer = new float[sampleNum];

	// put data into data_D (inline x xline x t)
	for (long trace_num = 0; trace_num <= traceNum - 1; trace_num++)
	{
		if (trace_num % 10000 == 0)
		{
			cout << "have read " << trace_num << " traces" << endl;
		}
		segyData.ReadOneTrace(trace_num, &pTraceHead, sDataBuffer);

		int in_index = pTraceHead.RecordNum;
		int x_index = pTraceHead.CdpNum;

		for (int i = 0; i < sampleNum; i++)
		{
			data_D(i, in_index - sLine, x_index - sCdp) = sDataBuffer[i];
		}

	}
	delete[] sDataBuffer;

	///////////////////***************   fxymssa   **************/////////////////////
	
	///////// block-by-block 
	
	// CUT & OVERLAP METHOD
	const int in_parts = inlineNum / (in_cut - overlap);
	const int in_remain = inlineNum % (in_cut - overlap);
	const int x_parts = xlineNum / (x_cut - overlap);
	const int x_remain = xlineNum % (x_cut - overlap);
	const int t_parts = sampleNum / (t_cut - overlap);
	const int t_remain = sampleNum % (t_cut - overlap);

	cx_cube block;

	///////// main_loop
	int counter = 0;

	// CUT & OVERLAP METHOD

	// coordinate
	int t_fa1, t_fa2, t_re1, t_re2; //t_fa means t_fake whose data will be overlapped, t_re means t_real whose data will be maintained
	
	int in_fa1, in_fa2, in_re1, in_re2;
	
	int x_fa1, x_fa2, x_re1, x_re2;
	

	// initialize in_line coordinate
	in_fa1 = 0; in_fa2 = in_cut - 1;

	in_re1 = in_fa1; in_re2 = in_re1 + (in_cut - overlap);


	while (in_fa2 != inlineNum - 1)
	{
		if (in_fa2 > inlineNum - 1)
		{
			in_fa2 = inlineNum - 1;
		}

		// initialize x_line coordinate
		x_fa1 = 0; x_fa2 = x_cut - 1;

		x_re1 = x_fa1; x_re2 = x_re1 + (x_cut - overlap);

		while (x_fa2 != xlineNum -1)
		{
			if (x_fa2 > xlineNum - 1)
			{
				x_fa2 = xlineNum - 1;
			}

			// initialize t coordinate
			t_fa1 = 0; t_fa2 = t_cut - 1;

			t_re1 = t_fa1; t_re2 = t_re1 + (t_cut - overlap);

			while (t_fa2 != sampleNum - 1)
			{
				if (t_fa2 > sampleNum - 1)
				{
					t_fa2 = sampleNum - 1;
				}

				counter++;

				cout << "the " << counter << "th block is being processed with " << (in_parts + 1) * (x_parts + 1) * (t_parts + 1) << " in total..." << endl;
				
				// process of fxymssa
				block = conv_to<cx_cube>::from(data_D.subcube(t_fa1, in_fa1, x_fa1, t_fa2, in_fa2, x_fa2));

				cube block_dealt = fxymssa_block_func(t_fa2 - t_fa1 + 1, in_fa2 - in_fa1 + 1, x_fa2 - x_fa1 + 1, block, flow, fhigh, dt, N, verb);

				data_D(span(t_fa1, t_fa2), span(in_fa1, in_fa2), span(x_fa1, x_fa2)) = block_dealt;

				block_dealt.reset();

				block.reset();

				if(t_fa2 == sampleNum - 1)
				{
					break;
				}

				// update t coordinate
				t_re1 = t_re2;

				t_fa1 = t_re2;

				t_re2 = t_re2 + (t_cut - overlap);

				t_fa2 = t_fa1 + t_cut - 1;

			}

			if (x_fa2 == xlineNum - 1)
			{
				break;
			}

			//update x_line coordinate
			x_re1 = x_re2;

			x_fa1 = x_re2;

			x_re2 = x_re2 + (x_cut - overlap);

			x_fa2 = x_fa1 + x_cut - 1;


		}

		if (in_fa2 == inlineNum - 1)
		{
			break;
		}

		// update in_line coordinate

		in_re1 = in_re2;

		in_fa1 = in_re2;

		in_re2 = in_re2 + (in_cut - overlap);

		in_fa2 = in_fa1 + in_cut - 1;
	}




	//remember to permute the cube, for you have to put a (inline x xline x t) size cube into the write segy function. 
	fxymssa<double> change;
	fcube data_D_permute = conv_to<fcube>::from(change.permute(data_D, std::tuple<uword, uword, uword>(2, 3, 1)));

	data_D.reset();

	///////////*********  write segy file ************////////////

	// format change
	const char* const_s = nullptr;
	const_s = _sFileName.c_str();
	char* s = nullptr;
	s = const_cast<char*>(const_s);

	const char* const_d = nullptr;
	const_d = _dFileName.c_str();
	char* d = nullptr;
	d = const_cast<char*>(const_d);


	// write reelhead & linehead & tracehead & tracedata
	write_sgy write_instance;

	FILE *fid1 = NULL;
	FILE *fid2 = NULL;
	//fcube data_D_permute_temp = conv_to<fcube>::from(data_D_permute);

	write_instance.write_SGY(s, fid1, d, fid2, data_D_permute, inlineNum, xlineNum, sampleNum, traceNum);

	return true;

}

template <typename T>
bool fxymssa<T>::fxymssa_block_overlap_all_ave_main(const int t_cut, const int in_cut, const int x_cut, const int overlap, const int sCdp, const int eCdp, const int sLine, const int eLine, string _sFileName, string _dFileName, int flow, int fhigh, float dt, int N, int verb)
{
	/*
	discription : deal the segy data block-by-block with overlap, faster than the slice-by-slice method,

	after the main loop of fxymssa, compute the average of neighboring slice slice-by-slice and finally 
	
	wash the data and get the total average.

	FXY_MSSA : F-XY domain multichannel singular spectrum analysis(MSSA)

	IN

	t_cut : cut size in time dimension
	in_cut : cut size in inline dimension
	x_cut : cut size in xline dimension
	sCdp : start Cdp
	eCdp : end Cdp
	sLine : start inline
	eLine : end inline
	_sFileName ; source segy file name
	_dFileName ; target segy file name

	flow : processing frequency range(lower)
	fhigh : processing frequency range(higher)
	dt : temporal sampling interval
	N : number of singular value to be preserved

	OUT

	D1 : output a segy data

	editor : yyh

	date : 2018-11-09
	*/

	///////////////////***************   read segy data ****************///////////////////

	CSegyData segyData;
	segyData.OpenSegyFile(_sFileName.c_str(), 0);

	long traceNum = segyData.GetSegyTotalTraceNumber();
	int sampleNum = segyData.GetSegySampleNumber();  //nt

	int inlineNum = eLine - sLine + 1; //nx
	int xlineNum = eCdp - sCdp + 1;; //ny

	cube data_D(sampleNum, inlineNum, xlineNum);
	//cx_cube data_D1 = data_D.zeros();

	TraceHead pTraceHead;
	float *sDataBuffer = new float[sampleNum];

	// put data into data_D (inline x xline x t)
	for (long trace_num = 0; trace_num <= traceNum - 1; trace_num++)
	{
		if (trace_num % 10000 == 0)
		{
			cout << "have read " << trace_num << " traces" << endl;
		}
		segyData.ReadOneTrace(trace_num, &pTraceHead, sDataBuffer);

		int in_index = pTraceHead.RecordNum;
		int x_index = pTraceHead.CdpNum;

		for (int i = 0; i < sampleNum; i++)
		{
			data_D(i, in_index - sLine, x_index - sCdp) = sDataBuffer[i];
		}

	}
	delete[] sDataBuffer;

	///////////////////***************   fxymssa   **************/////////////////////

	///////// block-by-block 

	// CUT & OVERLAP METHOD
	const int in_parts = inlineNum / (in_cut - overlap);
	const int in_remain = inlineNum % (in_cut - overlap);
	const int x_parts = xlineNum / (x_cut - overlap);
	const int x_remain = xlineNum % (x_cut - overlap);
	const int t_parts = sampleNum / (t_cut - overlap);
	const int t_remain = sampleNum % (t_cut - overlap);

	cx_cube block;

	///////// main_loop
	int counter = 0;

	// CUT & OVERLAP METHOD

	// coordinate
	int t_fa1, t_fa2, t_re1, t_re2; //t_fa means t_fake whose data will be overlapped, t_re means t_real whose data will be maintained

	int in_fa1, in_fa2, in_re1, in_re2;

	int x_fa1, x_fa2, x_re1, x_re2;


	// initialize in_line coordinate
	in_fa1 = 0; in_fa2 = in_cut - 1;

	in_re1 = in_fa1; in_re2 = in_re1 + (in_cut - overlap);


	while (in_fa2 != inlineNum - 1)
	{
		if (in_fa2 > inlineNum - 1)
		{
			in_fa2 = inlineNum - 1;
		}

		// initialize x_line coordinate
		x_fa1 = 0; x_fa2 = x_cut - 1;

		x_re1 = x_fa1; x_re2 = x_re1 + (x_cut - overlap);

		while (x_fa2 != xlineNum - 1)
		{
			if (x_fa2 > xlineNum - 1)
			{
				x_fa2 = xlineNum - 1;
			}

			// initialize t coordinate
			t_fa1 = 0; t_fa2 = t_cut - 1;

			t_re1 = t_fa1; t_re2 = t_re1 + (t_cut - overlap);

			while (t_fa2 != sampleNum - 1)
			{
				if (t_fa2 > sampleNum - 1)
				{
					t_fa2 = sampleNum - 1;
				}

				counter++;

				cout << "the " << counter << "th block is being processed with " << (in_parts + 1) * (x_parts + 1) * (t_parts + 1) << " in total..." << endl;

				// process of fxymssa
				block = conv_to<cx_cube>::from(data_D.subcube(t_fa1, in_fa1, x_fa1, t_fa2, in_fa2, x_fa2));

				cube block_dealt = fxymssa_block_func(t_fa2 - t_fa1 + 1, in_fa2 - in_fa1 + 1, x_fa2 - x_fa1 + 1, block, flow, fhigh, dt, N, verb);

				data_D(span(t_fa1, t_fa2), span(in_fa1, in_fa2), span(x_fa1, x_fa2)) = block_dealt;

				block_dealt.reset();

				block.reset();

				if (t_fa2 == sampleNum - 1)
				{
					break;
				}

				// update t coordinate
				t_re1 = t_re2;

				t_fa1 = t_re2;

				t_re2 = t_re2 + (t_cut - overlap);

				t_fa2 = t_fa1 + t_cut - 1;

			}

			if (x_fa2 == xlineNum - 1)
			{
				break;
			}

			//update x_line coordinate
			x_re1 = x_re2;

			x_fa1 = x_re2;

			x_re2 = x_re2 + (x_cut - overlap);

			x_fa2 = x_fa1 + x_cut - 1;


		}

		if (in_fa2 == inlineNum - 1)
		{
			break;
		}

		// update in_line coordinate

		in_re1 = in_re2;

		in_fa1 = in_re2;

		in_re2 = in_re2 + (in_cut - overlap);

		in_fa2 = in_fa1 + in_cut - 1;
	}

	///////////////////////  average all
	int in_cut_re = in_cut - overlap;
	int x_cut_re = x_cut - overlap;

	int k = 0;
	// xline direction
	while (k <= xlineNum - 2)
	{
		//int index_x = x_cut_re * k;
		//mat total_temp = data_D.slice(index_x) + data_D.slice(index_x + 1) + data_D.slice(index_x + 2) + data_D.slice(index_x - 1) + data_D.slice(index_x - 2);
		//mat late_ave = total_temp / 5;
		////update
		//data_D.slice(index_x) = late_ave;
		//data_D.slice(index_x - 1) = (late_ave + data_D.slice(index_x - 1)) / 2;
		//data_D.slice(index_x + 1) = (late_ave + data_D.slice(index_x + 1)) / 2;
		//k++;

		//update
		data_D.slice(k) = (data_D.slice(k) + data_D.slice(k + 1)) / 2;
		k++;

		//total_temp.reset();
		//late_ave.reset();
	}

	// inline direction
	k = 0;
	while (k <= inlineNum - 2)
	{
		//int index_in = in_cut_re * k;
		//mat total_temp = data_D.slice(index_in) + data_D.slice(index_in + 1) + data_D.slice(index_in + 2) + data_D.slice(index_in - 1) + data_D.slice(index_in - 2);
		//mat late_ave = total_temp / 5;
		////update
		//data_D.slice(index_in) = late_ave;
		//data_D.slice(index_in - 1) = (late_ave + data_D.slice(index_in - 1)) / 2;
		//data_D.slice(index_in + 1) = (late_ave + data_D.slice(index_in + 1)) / 2;
		//

		data_D.tube(span(0, sampleNum - 1), span(k, k)) = (data_D.tube(span(0, sampleNum - 1), span(k, k)) + data_D.tube(span(0, sampleNum - 1), span(k + 1, k + 1))) / 2;
		k++;
		//total_temp.reset();
		//late_ave.reset();
	}

	//remember to permute the cube, for you have to put a (inline x xline x t) size cube into the write segy function. 
	fxymssa<double> change;
	fcube data_D_permute = conv_to<fcube>::from(change.permute(data_D, std::tuple<uword, uword, uword>(2, 3, 1)));

	data_D.reset();

	///////////*********  write segy file ************////////////

	// format change
	const char* const_s = nullptr;
	const_s = _sFileName.c_str();
	char* s = nullptr;
	s = const_cast<char*>(const_s);

	const char* const_d = nullptr;
	const_d = _dFileName.c_str();
	char* d = nullptr;
	d = const_cast<char*>(const_d);


	// write reelhead & linehead & tracehead & tracedata
	write_sgy write_instance;

	FILE *fid1 = NULL;
	FILE *fid2 = NULL;
	//fcube data_D_permute_temp = conv_to<fcube>::from(data_D_permute);

	write_instance.write_SGY(s, fid1, d, fid2, data_D_permute, inlineNum, xlineNum, sampleNum, traceNum);

	return true;

}


template <typename T>
cube fxymssa<T>::fxymssa_block_func(int nt, int nx, int ny, Cube<T>& data_block, int flow, int fhigh, float dt, int N, int verb)
{
	/*
	discription : supports the "fxymssa_block_main" and "fxymssa_block_overlap_main", deal with the block data 

	and return with the same size.

	FXY_MSSA : F-XY domain multichannel singular spectrum analysis(MSSA)
	IN
	data_block: a block data which was cutted from the original data, whose data format is Armadillo Cube;
	nt : size of time dimension of the data_block
	nx ; size of inline dimension of the data_block
	ny : size of xline dimension of the data_block

	flow : processing frequency range(lower)
	fhigh : processing frequency range(higher)
	dt : temporal sampling interval
	N : number of singular value to be preserved

	OUT
	D1 : output a cube which has same size with the input data_block

	editor : yyh

	date : 2018-10-31
	*/

	//cube data_block(nt, nx, ny);

	int nf = next_pow_of_2(nt);

	cx_cube data_FX_b(nf, nx, ny);

	for (int x_index = 0; x_index <= ny - 1; x_index++)
	{
		if (x_index % 10 == 0)
		{
			cout << "have dealt " << x_index << " cdps" << endl;

		}
		
		cx_mat Y = arma::fft(data_block.slice(x_index), nf);
		
		data_FX_b.slice(x_index) = Y;

		
		Y.reset();
		

	}

	data_block.reset();

	int ilow = floor(flow * dt * nf) + 1;
	int ihigh = floor(fhigh * dt * nf) + 1;

	if (ilow < 1)
	{
		ilow = 1;
	}

	if (ihigh > floor(nf / 2) + 1)
	{
		ihigh = floor(nf / 2) + 1;
	}

	int lx = floor(nx / 2) + 1;
	int lxx = nx - lx + 1;
	int ly = floor(ny / 2) + 1;
	int lyy = ny - ly + 1;

	cx_mat M(lx * ly, lxx * lyy);
	M.zeros();

	cx_cube data_FX0_b(nf, nx, ny);

	data_FX0_b.zeros();

	//////////*********** main loop ****************/////////

	for (int k = ilow - 1; k <= ihigh - 1; k++)
	{
		if (ny == 1)
		{
			cx_mat temp = data_FX_b.tube(span(k, k), span(0, nx - 1));
			temp = temp.st();
			M = P_H(temp, lx, ly);
			temp.reset();
		}
		else
		{
			cx_mat temp = data_FX_b.tube(span(k, k), span(0, nx - 1));
			M = P_H(temp, lx, ly);
			temp.reset();

		}

		M = P_R(M, N);

		data_FX0_b.tube(span(k, k), span(0, nx - 1)) = P_A(M, nx, ny, lx, ly);

		if ((k % 5 == 0) && (verb == 1))
		{
			printf("F %d is done!\n", k);
		}
	}
	data_FX_b.reset(); // memory free
	M.reset();

	for (int k = nf / 2 + 1; k <= nf - 1; k++)
	{
		data_FX0_b.tube(span(k, k), span(0, nx - 1)) = conj(data_FX0_b.tube(span(nf - k, nf - k), span(0, nx - 1)));
	}

	cube D1_b(nf, nx, ny);
	for (int x_index = 0; x_index <= ny - 1; x_index++)
	{
		cx_mat Y = arma::ifft(data_FX0_b.slice(x_index), nf);
		mat Y_ = real(Y);

		D1_b.slice(x_index) = Y_;

	}

	data_FX0_b.reset();

	cube D1_b_sub = D1_b.subcube(0, 0, 0, nt - 1, nx - 1, ny - 1);
	D1_b.reset();
	return D1_b_sub;

}