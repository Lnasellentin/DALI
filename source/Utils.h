/*
 * Ln(a) Sellentin
 * Universität Heidelberg
 * Header for the utilities
 */

#ifndef __UTILS_H_INCLUDED__
#define __UTILS_H_INCLUDED__

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm> //can sort a c++ vector
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
using namespace std;



void printxyz(double x, double y, double z);

double vec1Matvec2(gsl_vector* vec1, gsl_matrix* mat, gsl_vector* vec2, int dim);

void printmatrix(gsl_matrix* mat, int zeilen, int spalten);

void invmatrix(gsl_matrix* tobeinverted, unsigned int dim, gsl_matrix* inverted);

double vMv(vector<double> vec1,gsl_matrix* M, vector<double> vec2,int dim );
void retSubFish(int p1, int p2, int dimbigfish, gsl_matrix* bigfish, gsl_matrix* subfish);

double trace(gsl_matrix* m, int dim );
void mult4mat(gsl_matrix* m1, gsl_matrix* m2, gsl_matrix* m3, gsl_matrix* m4, int dim, gsl_matrix* result);
void mult3mat(gsl_matrix* m1, gsl_matrix* m2, gsl_matrix* m3, int dim, gsl_matrix* result);
void mult2mat(gsl_matrix* m1, gsl_matrix* m2, int dim, gsl_matrix* result);
void mult2mat(gsl_matrix* m1, gsl_matrix* m2, int dim, gsl_matrix* result, bool alldiag);


double determinant(gsl_matrix* f, unsigned int dim) ;
void Gnuplotconfidencecontours(string specifier, vector<double> testprob);

void pythoncontours(vector<double> testprob, int xcount, int ycount, string executablename,string datafile,string xaxis, string yaxis);

void pythoncontours_w(vector<double> testprob, int xcount, int ycount, string executablename,string datafile,string xaxis, string yaxis);

#endif
