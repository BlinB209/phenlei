//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      vector.hpp
//! @brief     The function of array operations.
//! @author    Peng Zerui, Xu Pingyu.

#include <vector>
#ifndef VECTOR_TWO_IMENSION_H_FILE
#define VECTOR_TWO_IMENSION_H_FILE

using namespace std;

vector<int> creatmatrix_1D_int(int, int);
vector<double> creatmatrix_1D(int, double);
vector<vector<int>> creatmatrix_int(int, int, int value = 0);
vector<vector<double>> creatmatrix_double(int, int, double value = 0);
vector<vector<double>> creatDiagMatrix_double(const vector<double>&);
vector<double> add1D_num(const vector<double>&, double);
vector<double> multiply1D_num(const vector<double>&, double);
vector<vector<double>> multiply2D_num(const vector<vector<double>>&, double);
vector<double> add1D_matrix(const vector<double>&, const vector<double>&);
vector<double> minus1D_matrix(const vector<double>&, const vector<double>&);
vector<vector<int>> trans_int(const vector<vector<int>>&);
vector<vector<double>> trans_double(const vector<vector<double>>&);
vector<vector<double>> multiply_double(const vector<vector<double>>&, const vector<vector<double>>&);
vector<double> multiply_two_one(const vector<vector<double>>& A, const vector<double>& B);
void show_matrix_int(const vector<vector<int>>&);
void show_matrix_double(const vector<vector<double>>&);
vector<vector<double>> tile_copy_double(vector<double>, int);
vector<double> sum_Rows(vector<vector<double>>);
#endif