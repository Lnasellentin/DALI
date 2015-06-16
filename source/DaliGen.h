/*
 DaliGen Header
 Ln(a) Sellentin
 Universität Heidelberg
 2015
 */


#include "DaliBase.h"
using namespace std;

class DaliGen : public DaliBase
{

private:
    // flags to check what the code has done. Touch at your own risk!
    bool twiddlecheck;
    bool checkmatderivs;
    bool checkmatdoublederivs;
    bool checkmattriplederivs;
    bool lowmem; //enable a low-memory execution
    
    //how many derivatives shall maximally be taken
    int  derivorder;
    
   //finite difference steps for covariance matrix    
   double fstep_C, dstep_C, tstep_C;

    
    //steps for matrix differentiation, depending on parameter, not yet fully implemented
    // TODO: implement 
    vector<double> M_step;
    vector<double> M_dstep;
    vector<double> M_tstep;
    
public:

    /*A member function pointer used to point to the function DataCovMat*/
    typedef double (DaliGen::*matrixfiller)(vector <double> eval, int i, int j);

    vector<gsl_matrix*> DCov; /*Stores first derivatives of the data covariance matrix*/
    vector<vector<gsl_matrix*> > DDCov; /*Stores second derivatives of the data covariance matrix*/
    vector<vector<vector<gsl_matrix*> > > DDDCov; /*Stores third derivatives of the data covariance matrix*/

    vector<gsl_matrix* > Ctwiddle; /*Shorthand for C_0^{-1}C,_alpha*/
    vector<vector<gsl_matrix* > > CDoubletwiddle; /*Shorthand for C_0^{-1}C,_alphabeta*/
    vector<vector<vector<gsl_matrix* > > > CTripletwiddle; /*Shorthand for C_0^{-1}C,_alphabetagamma*/


    DaliGen(int theo, int indep, int noof, vector<double> fidin, int derivorder_in, string outdir_in);

    
    //old version that subtracts matrix entries
    double secOneMat(matrixfiller p, vector<double> eval, int par, int i, int j, double step=2e-2);
     
    //new version that subtracts matrices
    void DMatrix( matrixfiller p, vector<double> eval, int par, int matrixdim, gsl_matrix* result, double step = 0.01);
    double secOneMat(matrixfiller p, vector<double> eval, int par,int matrixdim,gsl_matrix* result, double step=0.01);
    void DDMatrix(matrixfiller p, vector<double> eval, int par1, int par2, int matrixdim, gsl_matrix* result, double step = 0.01);
    void DDDMatrix(matrixfiller p, vector<double> eval, int par1, int par2, int par3, int matrixdim, gsl_matrix* result, double step = 0.01);
    
    void SetSteps_Cov(double f_step_C, double d_step_C, double t_step_C);
    
 
    void Calc_DCs(); //calculates all first derivs of the covariance matrix
    void Calc_DDCs(); //calculates all second derivs of the covariance matrix
    void Calc_DDDCs(); //calculates all third derivs of the covariance matrix
    void FillCtwiddle();
    void DoubleTwiddle();

    void CalcFisher();
    void Calc_BeyondFish();

    void Free();
    void LowMemory();

};