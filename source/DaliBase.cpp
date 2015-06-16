/*
 DaliBase class source code.
 Ln(a) Sellentin
 Universität Heidelberg
 */

#include "DaliBase.h"

using namespace std;




/**
 Constructor of the DaliBase class. This constructor is automatically called by the code, and does not need to be acessed by the user. 
 It initializes all the Dali tensors to zero, allocates memory and sets default values for various parameters.
 - theo: the number of theoretical parameters that the user's application depends on
 - indep: the number of independent variables that each data point contains. Maybe you have a case where a datapoint is a 3-dimensional direction, then you would like to store three coordinates per data point, and consequently you would put indep =3.
 -noof: the number of data points, where each data point is considered as one measurement that simultaneously had indep variables under control
 -fidin: a vector of length theo. It stores the fiducial values of your theoretical parameters/the maximum likelihood values of your theoretical parameters
 -outdir_in: specifies the folder in which DALI will store your results. If the folder does not exist, DALI creates it
 **/
DaliBase::DaliBase(int theo, int indep, int noof, vector<double> fidin, string outdir_in):theodim(theo),
    indepdim(indep), noofdat(noof),fid(fidin),outdir(outdir_in)
{

    // This must be changed to 1 by the function "DeclareData()", if DaliCalc shall execute
    datavaribsdeclared = 0;
    datapoint_weights = gsl_matrix_calloc(noofdat,noofdat);
    gsl_matrix_set_identity(datapoint_weights); // if the user does not provide weights, the weighting matrix is the identity matrix

    checkd = 0;
    checkdd = 0;
    checkddd = 0;
    covcheck = 0;
    is_diag = 0; //assume covariance matrix not to be diagonal by default

    //default values for the finite differences-stepwidths
    fstep = 1e-3;
    dstep = 1e-3;
    tstep = 1e-3;
    verbosity = 0;

    fisher     = gsl_matrix_calloc(theo, theo);

    InvDataCov = gsl_matrix_calloc(noof, noof); // This is the inv. covariance of the data, and therefore noof times noof
    DataCov = gsl_matrix_calloc(noof, noof);


    /*Resizing the tensors*/
    Sabg.resize(theo);
    Qabgd.resize(theo);
    Pabgde.resize(theo);
    Habgdef.resize(theo);
    
    
    Sabg_compute_flags.resize(theo);
    Qabgd_compute_flags.resize(theo);
    Pabgde_compute_flags.resize(theo);
    Habgdef_compute_flags.resize(theo);

    for(int i = 0; i < theo; i++)
    {
        Sabg[i].resize(theo);
        Qabgd[i].resize(theo);
        Pabgde[i].resize(theo);
        Habgdef[i].resize(theo);
	
	Sabg_compute_flags[i].resize(theo);
        Qabgd_compute_flags[i].resize(theo);
        Pabgde_compute_flags[i].resize(theo);
        Habgdef_compute_flags[i].resize(theo);
	

        for(int j = 0; j < theo; j++)
        {
            Sabg[i][j].resize(theo);
            Qabgd[i][j].resize(theo);
            Pabgde[i][j].resize(theo);
            Habgdef[i][j].resize(theo);
	    
	    Sabg_compute_flags[i][j].resize(theo);
            Qabgd_compute_flags[i][j].resize(theo);
            Pabgde_compute_flags[i][j].resize(theo);
            Habgdef_compute_flags[i][j].resize(theo);

            for(int k = 0; k < theo; k++)
            {
                Sabg[i][j][k] = 0.0;
		Sabg_compute_flags[i][j][k] = 0;
		
                Qabgd[i][j][k].resize(theo);
                Pabgde[i][j][k].resize(theo);
                Habgdef[i][j][k].resize(theo);
		
		Qabgd_compute_flags[i][j][k].resize(theo);
                Pabgde_compute_flags[i][j][k].resize(theo);
                Habgdef_compute_flags[i][j][k].resize(theo);

                for(int l = 0; l < theo; l++)
                {
                    Qabgd[i][j][k][l] = 0.0;
		    Qabgd_compute_flags[i][j][k][l] = 0;
		    
                    Pabgde[i][j][k][l].resize(theo);
                    Habgdef[i][j][k][l].resize(theo);
		    
		    Pabgde_compute_flags[i][j][k][l].resize(theo);
                    Habgdef_compute_flags[i][j][k][l].resize(theo);

                    for(int m = 0; m < theo; m++)
                    {
                        Pabgde[i][j][k][l][m] = 0.0;
			Pabgde_compute_flags[i][j][k][l][m] = 0;
			
                        Habgdef[i][j][k][l][m].resize(theo);
			Habgdef_compute_flags[i][j][k][l][m].resize(theo);

                        for(int n = 0; n < theo; n++)
                        {
                            Habgdef[i][j][k][l][m][n] = 0.0;
			    Habgdef_compute_flags[i][j][k][l][m][n] = 0;
                        }

                    }
                }

            }
        }
    }



    /*resizing the array that shall store the data*/
    datavaribs.resize(noofdat); //noofdat data points must be stored
    for(int i = 0; i < noofdat; i++)
    {
        datavaribs.at(i).resize(indepdim); //and each data point has indep coordinates/entries
    }


    cout << endl;
    cout << "--------DALI---------" << endl;
    cout << "Memory for data allocated: " << noofdat << " data points can be stored, each point is a vector of " << indep << " entries." << endl;


    if(fidin.size() != theodim)
    {
        cout << "You provided a fiducial parameter set of more or less values than specified theoretical parameters" << endl;
	abort();
    }

    /*Determine whether Out-directory exists, and if not, create it*/
    struct stat dircheck;

    /*Function returns zero, if directory exists*/
    if(stat(outdir.c_str(),&dircheck) == 0 && S_ISDIR(dircheck.st_mode))
    {
        cout << "Output directory: " << outdir << endl;
    }

    else /*Make directory, if not existent*/
    {

        string createdirectory = (string)"mkdir " + outdir;

        system(createdirectory.c_str());

        cout << "Newly created output directory: " << outdir << endl;
    }


    dmod.resize(theo);
    for(int i = 0; i < theo; i++)
        dmod.at(i) = gsl_vector_calloc(noof); // because at each data pint a derivative will be taken



    /*The following creates a rank-2 tensor (we want double derivatives, hence rank-2).
     *Each tensor element is a gsl_vector with the length of the data set
     * calloc puts all vector elements to zero, which we will later use
     * for checking whether higher order derivatives are needed.
     */
    ddmod.resize(theo);
    dddmod.resize(theo);
    for(int i = 0; i < theo; i++)
    {
        ddmod.at(i).resize(theo);
        dddmod.at(i).resize(theo);
        for(int j = 0; j < theo; j++)
        {
            ddmod.at(i).at(j) = gsl_vector_calloc(noof);
            dddmod.at(i).at(j).resize(theo);

            for(int k = 0; k < theo; k++)
            {
                dddmod.at(i).at(j).at(k) = gsl_vector_calloc(noof);
            }

        }
    }



    //by default, the model and the covariance matrix are assumed to depend on parameters
    modeldep = 1;
    covmatdep = 1;
    RecommendedChecks = 0;

}


/**Telling DALI whether the model and/or the covariance matrix are parameter dependend. If model == 0, then the model is treated as parameter indepenend, if covmat == 0, then the covariance matrix is treated as independend, and no derivatives of the corresponding quantities will be calculated. This speeds up DALI.**/
void DaliBase::Parameterdependence(bool model, bool covmat)
{
    modeldep = model;
    covmatdep = covmat;
}

/**If the data covariance matrix is diagonal, the user can call this member function to speed up the execution of the code.**/
void DaliBase::SetDiagonal()
{
  cout << "Assuming diagonal covariance matrix." << endl;
  is_diag =1;
}


/**If a is greater than zero, the code will produce output files of the derivatives that the user has to check himself.**/
void DaliBase::ConductRecommendedChecks(int a)
{
  RecommendedChecks = a;
  
}





/**Allows the user to change the finite difference steps, that are used for the numerical evaluation of the derivatives for the MODEL. A corresponding function exists for the finite difference steps for the covariance matrix.
 -f_step: the finite difference step for the first derivatives of the model
 -d_step: for the second derivatives of the model
 -t_step: for the thrid derivatives of the model**/
void DaliBase::SetSteps(double f_step, double d_step, double t_step)
{
    fstep = f_step;
    dstep = d_step;
    tstep = t_step;
}

/**If v > 1, the code will print more messages to the screen.**/
void DaliBase::SetVerbosity(int v)
{
    verbosity = v;
}



/**Deallocates the memory of all pointer-like members which DALI uses**/
void DaliBase::BaseFree()
{
    gsl_matrix_free(fisher);
    gsl_matrix_free(InvDataCov);
    gsl_matrix_free(DataCov);
    gsl_matrix_free(datapoint_weights);

    for(int i = 0; i < theodim; i++)
    {
        gsl_vector_free(dmod.at(i));
        for(int j = 0; j < theodim; j++)
        {
            gsl_vector_free(ddmod.at(i).at(j));
            for(int k = 0; k < theodim; k++)
            {
                gsl_vector_free(dddmod.at(i).at(j).at(k));
            }

        }
    }

    cout << "Freed memory associated to the class DaliBase" << endl;
}




/*Takes the first derivative of the scalar function (i.e. one that returns a scalar quantity), to which the pointer MemFunPoint p points.
 -eval: the evaluation point of the derivatives
 -par: the number of the parameter with respect to which the function derives
 -step: the width of the finite difference step*/
double DaliBase::first(MemFunPoint p, vector<double> eval, int par, double step)
{

    if(par > theodim)
    {
        cout << "Taking derivative with respect to non-Existing parameter." << endl;
        abort();
    }

    if(eval.size() != theodim+indepdim)
    {
        cout << "Inserted vector does not contain the right amount of theoretical parameters." << endl;
        abort();
    }

    vector<double> plusdel(eval);
    vector<double> minusdel(eval);

    plusdel.at(par) += step;
    minusdel.at(par) -= step;

    /*2-point rule for finite difference step*/
    return 0.5*( (this->*p)(plusdel) - (this->*p)(minusdel)   )/step;

}



/*second deriv wrt to one parameter at the point eval*/
double DaliBase::secOne(MemFunPoint p, vector<double> eval, int par, double step)
{

    if(par > theodim)
    {
        cout << "Taking derivative with respect to non-Existing parameter" << endl;
        abort();
    }

    if(eval.size() != theodim + indepdim)
    {
        cout << "Inserted vector does not contain the right amount of theoretical parameters." << endl;
        abort();
    }

    vector<double> plusdel(eval);
    vector<double> minusdel(eval);

    plusdel.at(par) += step;
    minusdel.at(par)-= step;

    //2nd Order accuracy
    double hilf = ((this->*p)(plusdel) + (this->*p)(minusdel) -2.0*(this->*p)(eval));

    /*if hilf is smaller than tiny, floting point division by step will go wrong and yield a finite derivative,
     although it should be pretty much zero. So better return zero directly.*/
    if(fabs(hilf) < tiny)
        return 0.0;

    else
        return hilf/step/step; // returns the derivative

}


/*second derivative wrt two parameters at the point eval*/
double DaliBase::second(MemFunPoint p, vector<double> eval, int par1, int par2, double step)
{


    if( (par1 > theodim) || (par2 > theodim))
    {
        cout << "Taking derivative with respect to non-Existing parameter" << endl;
        abort();
    }


    if(eval.size() != theodim + indepdim)
    {
        cout << "Inserted vector does not contain the right amount of theoretical parameters." << endl;
        abort();
    }



    if(par1==par2)
        return secOne(p,eval,par1,step);


    vector<double> plusdel(eval);
    vector<double> minusdel(eval);

    plusdel.at(par1) += step;
    plusdel.at(par2) += step;
    minusdel.at(par1)-= step;
    minusdel.at(par2)-= step;


    return 0.5*(  ( (this->*p)(plusdel) + (this->*p)(minusdel) - 2.0*(this->*p)(eval) ) /step/step
                  - secOne(p,eval,par1) - secOne(p,eval,par2));

}


/*third derivative, using again the finite difference of second derivatives. This is numerically not the best way to do it - but the easiest way in multiple dimensions*/
double DaliBase::third(MemFunPoint p, vector<double> eval, int par1, int par2, int par3, double step)
{

    if( (par1 > theodim) || (par2 > theodim) || (par3 > theodim))
    {
        cout << "Taking derivative with respect to non-Existing parameter" << endl;
        abort();
    }

    if(eval.size() != theodim + indepdim)
    {
        cout << "Inserted vector does not contain the right amount of theoretical parameters." << endl;
        abort();
    }

    vector<double> plusdel(eval);
    vector<double> minusdel(eval);

    plusdel.at(par3) += step; //par3 here
    minusdel.at(par3)-= step;

    double hilf = (second(p,plusdel,par1,par2,step)-second(p,minusdel,par1,par2,step));

    if(fabs(hilf) < tiny)
        return 0.0;

    else
        return 0.5*hilf/step;

}


/*Checks whether the user has resized the datavaribs - which would lead to erroneous calculations*/
void DaliBase::CheckDataSize()
{
    if(datavaribs.size() != noofdat)
    {
        cout << "It seems at some point, you have resized your datavaribs in the first dimension: it does not have as many data points anymore, as you specified. Abort();" << endl;
        abort();
    }

    for(int i = 0; i < noofdat; i++)
    {

        if(datavaribs[i].size() != indepdim)
        {
            cout << "It seems at some point, you have resized your datavaribs in the second dimension: at least one entry does not store indep variables anymore. Abort();" << endl;
            abort();
        }

    }

}


/*calculates all first order derivatives of the model*/
void DaliBase::CalcFirstDmod()
{
    DaliCheckDataDecl();
    CheckDataSize();

    if(modeldep != 0)
    {
        cout << "Calculating first derivatives of the model." << endl;
        for(int alpha = 0; alpha < theodim; alpha ++) //Parameter-loop
        {
            if(verbosity > 1)
            {
                cout << "Deriving wrt parameter " << alpha << endl;
            }

            for(int i = 0; i < noofdat; i ++) //data-loop
            {
                vector<double> eval (fid);
                eval.insert(eval.end(), datavaribs.at(i).begin(), datavaribs.at(i).end()); // like this, a vector of length theodim + indepdim should be created
                gsl_vector_set( dmod[alpha], i, first(&DaliBase::PhysMod, eval, alpha, fstep) );
            }
        }
    }

    //set flag also, if model is param independent, to avoid infinite loop
    checkd =1;
}


/*calculates all second order derivatives of the model*/
void DaliBase::CalcScndDmod()
{
    DaliCheckDataDecl();
    CheckDataSize();

    if(modeldep != 0)
    {
        cout << "Calculating second derivatives of the model." << endl;
        if(!checkd)
            CalcFirstDmod();


        for(int alpha = 0; alpha < theodim; alpha++)
        {
            /*calculating only the upper triangle, and setting the lower the same*/
            for(int beta = alpha; beta < theodim; beta++)
            {

                if(verbosity > 1)
                {
                    cout << "Deriving wrt parameter " << alpha << " and " << beta << endl;
                }

                for(int i = 0; i < noofdat; i++)
                {
                    vector<double> eval (fid);
                    eval.insert(eval.end(), datavaribs.at(i).begin(), datavaribs.at(i).end());

                    double hilf = second(&DaliBase::PhysMod, eval, alpha, beta,dstep); // store the result, before setting all permutations

                    gsl_vector_set(ddmod.at(alpha).at(beta),i,hilf);
                    gsl_vector_set(ddmod.at(beta).at(alpha),i,hilf);


                }
            }
        }

    }
    checkdd = 1;
}


/*calculates all third order derivatives of the model*/
void DaliBase::CalcThrdDmod()
{
    DaliCheckDataDecl();
    CheckDataSize();

    if(modeldep != 0)
    {
        cout << "Calculating third derivatives of the model." << endl;
        for(int alpha = 0; alpha < theodim; alpha++)
        {
            for(int beta = alpha; beta < theodim; beta++)
            {
                for(int gamma = beta; gamma < theodim; gamma++)
                {

                    if(verbosity > 1)
                    {
                        cout << "Deriving wrt parameter " << alpha << ", " << beta << ", " << gamma << endl;
                    }


                    for(int i = 0; i < noofdat; i++)
                    {

                        vector<double> eval (fid);
                        eval.insert(eval.end(), datavaribs.at(i).begin(), datavaribs.at(i).end());

                        // store the result, before setting all permutations
                        double hilf = third(&DaliBase::PhysMod, eval, alpha, beta, gamma,tstep);

                        /*all six permutations of the three indices*/
                        gsl_vector_set(dddmod.at(alpha).at(beta).at(gamma),i,hilf);
                        gsl_vector_set(dddmod.at(beta).at(gamma).at(alpha),i,hilf);
                        gsl_vector_set(dddmod.at(gamma).at(alpha).at(beta),i,hilf);

                        gsl_vector_set(dddmod.at(alpha).at(gamma).at(beta),i,hilf);
                        gsl_vector_set(dddmod.at(gamma).at(beta).at(alpha),i,hilf);
                        gsl_vector_set(dddmod.at(beta).at(alpha).at(gamma),i,hilf);


                    }

                }
            }
        }

    }
    checkddd = 1;
}




/**Prints the first derivatives of the model to a file. Plotting them allows to check whether the finite difference has been adjusted appropriately.**/
void DaliBase::PrintFirstDmod()
{
    if(!checkd)
    {
        cout << "First Derivatives of the model were not yet calculated." << endl;
        abort();
    }

    string filename = outdir + (string)"/FirstDerivs";

    FILE * firsts;
    firsts = fopen(filename.c_str(), "w");

    for(int i = 0; i < noofdat; i ++)
    {
        for(int dat = 0; dat < indepdim; dat++)
        {
            fprintf(firsts, " %.8e " ,datavaribs[i][dat]);
        }

        for(int alpha = 0; alpha < theodim; alpha ++)
        {
            fprintf(firsts, " %.8e " ,gsl_vector_get( dmod[alpha], i));
        }
        fprintf(firsts,"\n");
    }

    fclose(firsts);
}



/**Writes the second derivatives of the model to the file outdir/SecondDerivs**/
void DaliBase::PrintScndDmod()
{
    if(!checkdd)
    {
        cout << "Second Derivatives were not yet calculated." << endl;
        abort();
    }

    FILE * scnds;

    string filename = outdir + (string)"/SecondDerivs";
    scnds = fopen(filename.c_str(), "w");

    for(int alpha = 0; alpha < theodim; alpha ++)
    {
        for(int beta = alpha; beta < theodim; beta++)
        {
            fprintf(scnds, "#       %d%d     " ,alpha,beta);
        }
    }
    fprintf(scnds,"\n");


    for(int i = 0; i < noofdat; i ++)
    {
        for(int dat = 0; dat < indepdim; dat++)
        {
            fprintf(scnds, " %.8e " ,datavaribs[i][dat]);
        }
        for(int alpha = 0; alpha < theodim; alpha ++)
        {
            for(int beta = alpha; beta < theodim; beta++)
            {
                fprintf(scnds, " %.8e " ,gsl_vector_get( ddmod[alpha][beta], i));
            }
        }
        fprintf(scnds,"\n"); // new line for new element from the data set
    }
    fclose(scnds);
}


/**Writes the third derivatives of the model to the file outdir/ThirdDerivs**/
void DaliBase::PrintThrdDmod()
{
    if(!checkddd)
    {
        cout << "Third Derivatives were not yet calculated." << endl;
        abort();
    }

    FILE * thrds;

    string filename = outdir + (string)"/ThirdDerivs";

    thrds = fopen(filename.c_str(), "w");

    for(int alpha = 0; alpha < theodim; alpha ++)
    {
        for(int beta = alpha; beta < theodim; beta++)
        {
            for(int gamma = beta; gamma < theodim; gamma++)
            {
                /*should be indented such, that the title stands over correct column, even if number of indep variables is increased/decreased*/
                fprintf(thrds, " #    %d%d%d       " ,alpha,beta, gamma);
            }
        }
    }
    fprintf(thrds,"\n");



    for(int i = 0; i < noofdat; i ++)
    {
        for(int dat = 0; dat < indepdim; dat++)
        {
            fprintf(thrds, " %.8e " ,datavaribs[i][dat]);
        }
        for(int alpha = 0; alpha < theodim; alpha ++)
        {
            for(int beta = alpha; beta < theodim; beta++)
            {
                for(int gamma = beta; gamma < theodim; gamma++)
                {
                    fprintf(thrds, " %.8e " ,gsl_vector_get( dddmod[alpha][beta][gamma], i));
                }
            }

        }
        fprintf(thrds,"\n"); // new line for new element from the data set
    }

    fclose(thrds);
}


/*Calculates the Data Covariance matrix*/
void DaliBase::FillDataCovMat()
{
    if(!datavaribsdeclared)
        DeclareData();

    cout << "Calculating data covariance matrix............";
    for(int i = 0; i < noofdat; i++)
    {
        /*Calculating only the upper triangle and setting the
        *lower triangle equal.*/
        for(int j = i; j < noofdat; j++)
        {
            double hilf = DataCovMatrix(fid, i, j);
            gsl_matrix_set(DataCov, i, j, hilf);
            gsl_matrix_set(DataCov, j, i, hilf);

        }
    }


    cout << "done." << endl;
    //printmatrix(DataCov,noofdat,noofdat);

    cout << "Calculating inverse data covariance matrix......";
    if(!is_diag)
    invmatrix(DataCov, noofdat, InvDataCov);
    
    if(is_diag)
    {
      cout << "covariance matrix assumed to be diagonal....calculating....";
      
      gsl_matrix_set_identity(InvDataCov); //only to get the offdiagonal zeros
      
      //overwriting the diagonal:
      for(int i = 0; i < noofdat; i++)
      {
       gsl_matrix_set(InvDataCov,i,i, 1./gsl_matrix_get(DataCov,i,i));
      }
      
      
    }
  

    cout << "done." << endl;
    //printmatrix(InvDataCov,noofdat,noofdat);

    covcheck = 1;
}




void DaliBase::DaliCheckDataDecl()
{
    /*First, this function checks whether the Data have been declared. If not, it tries to execute the DeclareData() function. If the data have been already declared in the main(), this step is skipped, and the function proceeds happily to calculate the Dali-tensors.*/
    if(datavaribsdeclared == 0)
    {
        cout << "Trying to declare data....";
        DeclareData();
        cout << "data declared." << endl;
    }

    if(datavaribsdeclared == 0)
    {
        cout << "You probably forgot to put 'datadeclared =1;' at the end of your DeclareData() function." << endl;
        abort();
    }

}




/**Saves the Fisher matrix and the Dali tensors to separate files in the outdir.
 * The formatting of the files is
 * "index index index Sabg_Tensor entry"
 * in case of the Sabg-tensor. For the Fisher matrix and the higher order tensors
 * the same formatting is chosen, just with more or less indices, separated by spaces.
 * Keep in mind which index number refers to which parameter.
 * A common question is: "Why are these tensors not invariant under all permutations
 * of the indices? Why is my Sabg[0][0][1] not the same as my Sabg[0][1][0]?"
 * The reason is that the tensors here used are those of Eq.(14) in arxiv 1401.6892_v2.
 -filesuffix: you may provide some suffix in order to identify the files.**/
void DaliBase::DaliSave(string filesuffix)
{
    cout << "Save Tensors to file." << endl;

    string FishName    = outdir + (string)"/Fish"    + filesuffix.c_str();
    string SabgName    = outdir + (string)"/Sabg"    + filesuffix.c_str();
    string QabgdName   = outdir + (string)"/Qabgd"   + filesuffix.c_str();
    string PabgdeName  = outdir + (string)"/Pabgde"  + filesuffix.c_str();
    string HabgdefName = outdir + (string)"/Habgdef" + filesuffix.c_str();

    FILE* Fishout = fopen(FishName.c_str(), "w");

    FILE* Sabgout = fopen(SabgName.c_str(),"w");
    FILE* Qabgdout = fopen(QabgdName.c_str(), "w");
    FILE* Pabgdeout = fopen(PabgdeName.c_str(),"w");
    FILE* Habgdefout = fopen(HabgdefName.c_str(),"w");

    for(int i = 0; i < theodim; i++)
    {
        for(int j = 0; j < theodim; j++)
        {
            fprintf(Fishout, "%d %d %.18e\n",i, j, gsl_matrix_get(fisher,i,j));
            for(int k = 0; k < theodim; k++)
            {
                fprintf(Sabgout, "%d %d %d %.18e\n",i,j,k, Sabg[i][j][k]);
                for(int l = 0; l < theodim; l++)
                {
                    fprintf(Qabgdout, "%d %d %d %d %.18e\n",i,j,k,l,Qabgd[i][j][k][l]);
                    for(int m = 0; m < theodim; m++)
                    {
                        fprintf(Pabgdeout, "%d %d %d %d %d %.18e\n",i,j,k,l,m,Pabgde[i][j][k][l][m]);
                        for(int n = 0; n < theodim; n++)
                        {
                            fprintf(Habgdefout, "%d %d %d %d %d %d %.18e\n",
                                    i,j,k,l,m,n,Habgdef[i][j][k][l][m][n]);
                        }
                    }

                }
            }
        }
    }


    fclose(Fishout);
    fclose(Sabgout);
    fclose(Qabgdout);
    fclose(Pabgdeout);
    fclose(Habgdefout);
}


/**Creates a file that allows to plot a square matrix. This can be used in order to check whether the derivatives of the covariance matrix seem okay.**/
void DaliBase::PlotMatrix(gsl_matrix* ToBePlotted, int dim, string filesuffix)
{
    string MatrixName    = outdir + (string)"/Matrix"    + filesuffix.c_str();
    FILE* Matrixout = fopen(MatrixName.c_str(), "w");

    //if(is_diag)
    {
      for(int i = 0; i < dim; i++)
      {
        fprintf(Matrixout, "%d %d %.12e\n",i,i,gsl_matrix_get(ToBePlotted,i,i));       
      }
      
    }
    
  /*  else
    {
    
    for(int i = 0; i < dim; i++)
    {
        for(int j = 0; j < dim; j++)
        {
            fprintf(Matrixout, "%d %d %.12e\n",i,j,gsl_matrix_get(ToBePlotted,i,j));
        }
        fprintf(Matrixout,"\n");
    }
    }*/

    fclose(Matrixout);


}

