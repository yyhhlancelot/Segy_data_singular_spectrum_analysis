#pragma once
#include "stdafx.h"
#ifndef FXYMSSA_H
#define FXYMSSA_H

#include <math.h>
#include "CSegyData.h"
#include <iostream>
#include <string.h>
#include <cstring>
#include <vector>
#include <complex>
#include "fcntl.h"
#include <time.h>
#include <armadillo>

using namespace std;
using namespace arma;


template <typename T> 
class fxymssa
{

public:
	fxymssa(void);
	~fxymssa(void);

public:

	cx_mat P_H(cx_mat din, int lx, int ly);
	cx_mat P_R(cx_mat din, int N);
	cx_mat P_A(cx_mat din, int nx, int ny, int lx, int ly);
	cx_colvec ave_antid(cx_mat din);
	inline int next_pow_of_2(int v);
	cx_mat hankel(cx_colvec fir_col, cx_rowvec la_row);
	bool fxymssa_func_write(const int sCdp, const int eCdp, const int sLine, const int eLine, string _sFileName, string _dFileName, int flow, int fhigh, float dt, int N, int verb);
	mat fxymssa_slice_step(int slice_index, int sampleNum, int inlineNum, cx_mat Di, int flow, int fhigh, float dt, int N, int verb);
	bool fxymssa_slice_main(const int sCdp, const int eCdp, const int sLine, const int eLine, string _sFileName, string _dFileName, int flow, int fhigh, float dt, int N, int verb);
	inline urowvec fxymssa::shape(const Cube<T>& x);
	Cube<T> permute(Cube<T>& cube, const std::tuple<uword, uword, uword>& order);
	bool fxymssa::fxymssa_block_overlap_main(const int t_cut, const int in_cut, const int x_cut, const int overlap, const int sCdp, const int eCdp, const int sLine, const int eLine, string _sFileName, string _dFileName, const int flow, const int fhigh, const float dt, const int N, const int verb);
	bool fxymssa::fxymssa_block_overlap_all_ave_main(const int t_cut, const int in_cut, const int x_cut, const int overlap, const int sCdp, const int eCdp, const int sLine, const int eLine, string _sFileName, string _dFileName, const int flow, const int fhigh, const float dt, const int N, const int verb);
	bool fxymssa::fxymssa_block_overlap_ave_main(const int t_cut, const int in_cut, const int x_cut, const int overlap, const int sCdp, const int eCdp, const int sLine, const int eLine, string _sFileName, string _dFileName, const int flow, const int fhigh, const float dt, const int N, const int verb);
	
	bool fxymssa::fxymssa_block_main(const int t_cut, const int in_cut, const int x_cut, const int sCdp, const int eCdp, const int sLine, const int eLine, string _sFileName, string _dFileName, int flow, int fhigh, float dt, int N, int verb);
	cube fxymssa::fxymssa_block_func(int nt, int nx, int ny, Cube<T>& data_block, int flow, int fhigh, float dt, int N, int verb);
};

template <typename T>
inline int fxymssa<T>::next_pow_of_2(int v)
{
	int next;
	int out;
	next = ceil(log(v) / log(2.0));
	out = pow(2, next);
	return out;
}

template <typename T>
inline urowvec fxymssa<T>::shape(const Cube<T>& x)
{
	return { x.n_rows, x.n_cols, x.n_slices };
}

#endif