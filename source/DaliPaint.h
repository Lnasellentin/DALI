/*
 DaliPaint Header
 Ln(a) Sellentin
 Universität Heidelberg
 2015
 */


#ifndef __DALIPAINT_H_INCLUDED__
#define __DALIPAINT_H_INCLUDED__

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include <cmath>
#include <algorithm>
#include "Utils.h"
#include <sys/stat.h>
#include <gsl/gsl_histogram2d.h>
using namespace std;


class DaliPaint
{

private:
 //MCMC-things that Metropolis and Hamiltonian need
    int seed_p;
    gsl_rng * r;
    int startflag;
    double pratio;  
    double plast;
    
    stringstream seedstr;
    string seedstring;
    string ChainName;
    FILE* Chainout;
    int flush;
    
      
    //things that Metropolis MCMC needs to run
    double pcurrent;
    gsl_matrix* Cholesky;
    gsl_vector* jump;
    gsl_vector* draw;   
    bool Cholesky_jump;

    vector <double> testx;
    vector <double> curx;
    
    //storing last samples in case the new one gets rejected
    vector <double> lastx;  
    double lastp;


  /*----DALI-Specific things---------*/

    /*vectors to store the minimum and maximum of all sampling-coordinates
     (needed for automatic generation of plots)*/
    vector<double>minx;
    vector<double>maxx;

    double tiny;
    int fisher_alloced;
    bool loaded_chain;
    

public:
   
    vector<vector <double> > chain;
    /*this stores the probablities at the samples of the chain*/
    vector<double> chain_y;
      
    //These vector elements can be filled with Latex, such that python prints nice axis lables
    vector<string> paramnames;
    string outdir;

    vector<double> fid;
    int theodim;
    vector<double> upperbounds;
    vector<double> lowerbounds;
    bool boundsgiven;
    int derivorder;

    /*Since these are public members, they can be asessed directly
     * from the main. Do so, if you want to add priors and stuff*/
    gsl_matrix* fisher;

    //introduce typedefs here, since we'll need these
    //things for adding/combining probes/block-diagonal matrices
    typedef vector<vector<vector<double> > > flex;
    typedef vector<vector<vector<vector<double> > > > quarx;
    typedef vector<vector<vector<vector<vector<double> > > > > pent;
    typedef vector<vector<vector<vector<vector<vector<double> > > > > > hex;

    flex Sabg;
    quarx Qabgd;
    pent Pabgde;
    hex Habgdef ;


    DaliPaint(int theo, vector<double> fidin, string outdir_name, int derivorder_in);
    void Free();
    void SetBounds(vector<double> lowerin, vector<double> upperin);


    void Accept();
    void Reject();
    void sampling(int repetitions,vector<double> sigmas,vector<double> startx, int seed);
    void print_chain();
    void best_fit();

    void gridevaluation(string outname,double binwidth);
    void marginDown(int p1, int p2, string outname, double binwidth = 0.03);

    void DaliFlashBack_Fish(string filename, int minus);
    void DaliFlashBack_Sabg(string filename, int minus);
    void DaliFlashBack_Qabgd(string filename, int minus);
    void DaliFlashBack_Pabgde(string filename, int minus);
    void DaliFlashBack_Habgdef(string filename, int minus);

    virtual double ToBeSampled(vector<double> x);

    double LogLike(vector<double> peval);
    void FisherMargins(gsl_matrix* fisherin, int dim, double binwidths, int ex);

    void WriteChainToFile();
    void LoadChainFromFile(string filename);

    void PrepareCholesky();
    
    void AddFishers(gsl_matrix* fisher_in);
    void AddSabgs(flex Sabg_in);
    void AddQabgds(quarx Qabgd_in);
    void AddPabgdes(pent Pabgdef_in);
    void AddHabgdefs(hex Habgdef_in);

};


#endif