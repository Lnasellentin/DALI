/*
 DALI - Base class
 Ln(a) Sellentin
 Universität Heidelberg
 2015
 */

#ifndef __DALIBASE_H_INCLUDED__
#define __DALIBASE_H_INCLUDED__


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <gsl/gsl_matrix.h>
#include "Utils.h"
#include <sys/stat.h>  //needed to create folders, if they don't exist

using namespace std;

class DaliBase
{

private:
    const static double tiny = 1e-10;

public:
    string outdir;

    bool modeldep; //does the model depend on parameters?
    bool covmatdep; //does the covariance matrix depend on parameters?
    bool is_diag;   //whether the covariance matrix is diagonal
    int theodim;  // number of theoretical parameters (Om_m, Om_L, w_0, w_a...)
    int indepdim; // number of independent variables (z, k, theta, phi...)
    
    /**The number of data points.**/
    int noofdat;

    double fstep; //step for first derivs of model
    double dstep; //step for second derivs model
    double tstep;  //step for third derivs model

    int verbosity;
    int RecommendedChecks;


    /**
     *Storage room for the data. The structure <vector<vector<double > > works like a C-array ar[i][j], and is also called like this: datavaribs[i][j] is the jth independent variable if the ith datapoint. I.e. each datavaribs[i][j] is of type double. Each datavaribs[i] is a vector of doubles.
     * Each of the vectors datavaribs[i] will later be inserted into the
     * PhysMod, as a point where derivatives are taken. Make sure to correctly distribute your data into this array. If you have 50 data points, the index i of datavaribs[i][j] will run from 0 to 49. If you have only one independent variable, only datavaribs[i][0] will exist. So a cosmologist could for example set datavaribs[i][0] = redshift[i]. The code resizes datavaribs for you: it will have dimensions noofdat (the number of datapoints) times indep (the number of independent variables).  Do not resize this vector to your convenience! Also don't use push_back().
     **/
    vector<vector< double > > datavaribs;
    
    /**
     A gsl_matrix* that allows you to multiply data points with a weight. The weights must be constant; i.e. they cannot depend on parameters. Abbreviating the datapoint_weights matrix by W, the loglikelihood is given by 
     L = -1/2 ( (\vec{d} - \vec{m}) W C^{-1} (\vec{d} - \vec{m} ) ). 
     **/
    gsl_matrix* datapoint_weights;

    gsl_matrix* fisher; //this will be the fisher-matrix, which is symmetric

    //The higher order Dali-tensors
    vector<vector<vector<double> > > Sabg; 
    vector<vector<vector<vector<double> > > >  Qabgd; 
    vector<vector<vector<vector<vector<double> > > > > Pabgde;
    vector<vector<vector<vector<vector<vector<double> > > > > > Habgdef ;
    
    //some flags in order to check which entries have already been calculated
    vector<vector<vector <int> > > Sabg_compute_flags;
    vector<vector<vector<vector<int> > > > Qabgd_compute_flags;
    vector<vector<vector< vector<vector<int> > > > > Pabgde_compute_flags;
    vector<vector<vector< vector< vector<vector<int> > > > > > Habgdef_compute_flags;
    
    
    gsl_matrix* DataCov; //the data covariance matrix
    gsl_matrix* InvDataCov; //the inverse of it

    /**A flag that must be set by the function DeclareData(),
    after the user has stored the data in the vector datavaribs**/
    bool datavaribsdeclared;

    bool checkd, checkdd, checkddd; //whether derivatives have been calculated
    bool covcheck;

    /**Stores the values of the parameters at the maximum likelihood point (for forecasting, this is also called 'the fiducial'. At this point, DALI will evaluate the derivatives. The resulting DALI-tensors will depend on the choice of the fiducial point! If you change the fiducial, you will have to recalculate your tensors.**/
    vector<double> fid;

    vector<gsl_vector* > dmod; // first derivatives of the model
    vector<vector<gsl_vector* > > ddmod; //second derivatives of the model
    vector<vector<vector<gsl_vector* > > > dddmod; // third order derivatives

    /*member function pointer: looks crazy but it's very handy :) */
    typedef double (DaliBase::*MemFunPoint)(vector <double> x);


    /*Member functions are documented in the .cpp files, apart from the pure virtuals.*/
    DaliBase(int theo, int indep, int noof, vector<double> fidin, string outdir_in);  
    void BaseFree(); 
    
    double first(MemFunPoint p, vector<double> eval, int par, double step = 1e-3);
    double secOne(MemFunPoint p, vector<double> eval, int par, double step = 1e-2);
    double second(MemFunPoint p, vector<double> eval, int par1, int par2, double step = 1e-2);
    double third(MemFunPoint p, vector<double> eval, int par1, int par2, int par3, double step = 1e-3);


    /**A pure virtual function that must be overwritten by the user. The vector pars_and_indep has length theodim+indep. Its first theodim entries are reserved for the parameters, its last indep entries are reserved for independent variables. The function shall return the modelprediction as a function of the theoretical parameters and the independent variables.
     The ordering is important: The code will only take derivatives of this function with respect to its first theodim parameters.**/
    virtual double PhysMod(vector<double> pars_and_indep) = 0;

    /**A pure virtual function that must be overwritten by the user. See the Tutorial on how to achieve this.
     - pars: is a vector of length theodim that holds the theoretical parameters.
     - i and j: the entries/elements of the data covariance matrix. If your data covariance matrix depends on the data point datavaribs[i][j], use i and j to access datavaribs[i][j]. The function shall return the (i,j)-element of the Data covariance matrix. **/
    virtual double DataCovMatrix(vector <double> pars, int i, int j) = 0;

    /**A pure virtual function that must be overwritten by the user. The function must set the datavaribs[i][j], where i runs in [0,noofdat-1] and j runs in [0,indep-1]. At constant i, datavaribs contains the complete ith measurement. The jth independent variable of the ith measurement is then datavaribs[i][j]. The function must end with "datavaribsdeclared=1;" This tells the code that it now knows the independent variables.**/
    virtual void DeclareData() = 0;


    void CalcFirstDmod(); // calculates all first derivs of the model
    void CalcScndDmod();  // calculates all second derivs of the model
    void CalcThrdDmod();  // calculates all third derivs of the model

    void PrintFirstDmod();
    void PrintScndDmod();
    void PrintThrdDmod();
    void PlotMatrix(gsl_matrix* ToBePlotted, int dim, string filesuffix);

    void FillDataCovMat();
    void DaliCheckDataDecl();
    void DaliSave(string filesuffix);

    void Parameterdependence(bool model, bool covmat);
    void SetSteps(double f_step, double d_step, double t_step);
    void CheckDataSize();
    void SetVerbosity(int v);
    void SetDiagonal();
    void ConductRecommendedChecks(int a);

    
};



#endif
