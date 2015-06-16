/*
 * DaliGen source code
 * Ln(a) Sellentin
 * Universität Heidelberg
 * 
 *Source code for Dalilizing a likelihood which contains a parameter dependent model, as well as a parameter dependent covariance matrix. Functions for deriving matrices with respect to their parameters can be found here. The Dali-tensors are also calculated here.
 */

#include "DaliGen.h"
using namespace std;


/**
  Constructor of a general dali object. The vectors which hold your parameters, the covariance matrix and its derivatives are allocated. Boolean flags are set to zero, indicating that nothing has been calculated yet.
  This constructor needs to be given the following parameters:
     - int theo: the number of theoretical parameters, i.e. the dimensionality of your likelihood in parameter space
     - int indep: the number of independent variables, that your data depend on. So basically everything that is not a theoretical parameter but still goes into the data. In cosmology this could be e.g. the redshift or the mass of some object. indep can be larger than 1, if your data depend on more than one independent variable, for example on time AND temperature, as in case of the fishes that are described in the manual.
     - int noof: the number of data points in your data set
     - vector <double> fidin: a vector of length theo, that contains the values of your fiducial parameters. The order of the parameters must be kept the same throughout the code! So please remember which parameter you store at position fidin[i].
     - int derivorder: if set to one DALI will calculate third derivatives, if zero, DALI will stop at second order derivatives
     - string outdir_in: tells the code where to store its results. If the folder outdir_in does not exist, DALI will create it. If it exists, it will not be overwritten.
**/
DaliGen::DaliGen(int theo, int indep, int noof, vector<double> fidin, int derivorder_in, string outdir_in):DaliBase(theo,indep,noof,fidin, outdir_in)
{

    cout << endl;
    cout << "Constructor: Creating a general Dali-Object." << endl;


    if(fidin.size() != theo)
    {
        cout << "The vector containing your fiducial values does not have the same dimension as the number of parameters that you specified. Abort()." << endl;
        abort();
    }

    int set = 0;
    for(int i = 0; i < theodim; i++)
    {
        if(fabs(fidin[i]) > 1e-10)
            set = 1;
    }

    if(!set)
    {
        cout << "WARNING: All your fiducial parameters are zero. Did you forget to assign them values?" << endl;
    }


    if((derivorder_in != 2 ) && (derivorder_in != 3))
    {
        cout << "Requested derivorder is: " << derivorder_in << endl;

        cout << "You can only choose between second order derivatives or third order derivatives. If you try to use DALI as a Fisher matrix code by having wished to only use first order derivatives: Why :) ? If you want to draw Fisher ellipses for comparison, then you can simply avoid passing on the DALI tensors to the Dalipainter. " << endl;
        abort();
    }

    derivorder = derivorder_in;

    DCov.resize(theo);
    DDCov.resize(theo);
    DDDCov.resize(theo);

    // Ctwiddle = the inverse covariance matrix times a once derived covariance matrix. This combination appears very often throughout the code, such that an abbreviation is handy
    Ctwiddle.resize(theo);

    //Same as Ctwiddle, but with second derivatives
    CDoubletwiddle.resize(theo);
    //Same as Ctwiddle, but with third derivatives
    CTripletwiddle.resize(theo);

    for(int p = 0; p < theo; p++)
    {
        DCov.at(p) =  gsl_matrix_alloc(noofdat, noofdat);
        Ctwiddle.at(p) = gsl_matrix_alloc(noofdat, noofdat);

        DDCov.at(p).resize(theo);
        DDDCov.at(p).resize(theo);
        CDoubletwiddle.at(p).resize(theo);
        CTripletwiddle.at(p).resize(theo);
        for(int q = 0; q < theo; q++)
        {

            DDDCov.at(p).at(q).resize(theo); //results in a theo /times theo /times theo cube to store the 3rd derivs
            CTripletwiddle.at(p).at(q).resize(theo);
        }
    }




    for(int p = 0; p < theo; p++)
    {
        for(int q = 0; q <= p; q++)
        {
            DDCov.at(p).at(q) = gsl_matrix_alloc(noofdat, noofdat);
            CDoubletwiddle.at(p).at(q) = gsl_matrix_alloc(noofdat, noofdat);

            for(int rr = 0; rr <= q; rr++)
            {
                /*only allocate memory at one position. Later, the pointers
                * of the other permutations are pointing to this *single* memory */
                DDDCov.at(p).at(q).at(rr) = gsl_matrix_alloc(noofdat, noofdat);
                CTripletwiddle.at(p).at(q).at(rr) = gsl_matrix_alloc(noofdat, noofdat);
            }
        }
    }


    /*defaultsteps for matrix differentiation: not yet ready
    for(int i = 0; i < theodim; i++)
    {
      M_step.push_back(1e-3);  //first derivative
      M_dstep.push_back(1e-2); //second derivatives
      M_tstep.push_back(1e-2); //third derivatives
    }

    */
    
    

    //setting flags to zero
    checkmatderivs=0;
    checkmatdoublederivs = 0;
    checkmattriplederivs = 0;
    twiddlecheck = 0;
    lowmem = 0;
    
    //default values for the finite difference steps of the covariancematrix
    fstep_C = 0.01;
    dstep_C = 0.01;
    tstep_C = 0.01;


}


/**If this function is called before the DaliTensors are calculated, then the code will not keep the derivatives of the covariance matrix in memory. The derivative matrices then cannot be plotted. The DaliTensors can still be calculated (because the code uses the Ctwiddles for that.)**/
void DaliGen::LowMemory()
{
  lowmem = 1; 
}


/**Allows the user to set the step widths of the finite difference steps when taking derivatives of the covariance matrix. A similar function exists for the derivatives of the physical model.**/
void DaliGen::SetSteps_Cov(double f_step_C, double d_step_C, double t_step_C)
{
    fstep_C = f_step_C;
    dstep_C = d_step_C;
    tstep_C = t_step_C;
}



/**Frees the memory associated with the pointers in the DaliGen class.**/
void DaliGen::Free()
{

    BaseFree();


    for(int p = 0; p < theodim; p++)
    {
        gsl_matrix_free(Ctwiddle.at(p));
	
        if(!lowmem) gsl_matrix_free(DCov.at(p));
	
        for(int q = 0; q <= p; q++)
        {
            if(p==q)
            {
                //need this to make double-triple-quadruple-sure that the same memory is not freed twice
                gsl_matrix_free(CDoubletwiddle.at(p).at(p));
                
                if(!lowmem) gsl_matrix_free( DDCov.at(p).at(p));

            }

            if(p!= q)
            {
                //don't need to free the other permutation, since the pointers point to the same memory
               if(!lowmem)  gsl_matrix_free( DDCov.at(p).at(q));
	       
               gsl_matrix_free(CDoubletwiddle.at(p).at(q));

            }
            for(int rr = 0; rr <= q; rr++)
            {
                if((p==q) && (q==rr))
                {
                    gsl_matrix_free(CTripletwiddle.at(p).at(p).at(p));
		    
                    if(!lowmem) gsl_matrix_free(DDDCov.at(p).at(p).at(p));		                   
                }

                else
                {
                    /*same story as above: the other permutations
                    * point to the same memory space, and therefore
                    * don't need to be freed separately*/
		    
		    gsl_matrix_free(CTripletwiddle.at(p).at(q).at(rr));
		    
                    if(!lowmem) gsl_matrix_free(DDDCov.at(p).at(q).at(rr));

                }

            }


        }
    }

    cout << "Freed memory associated to the class DaliGen" << endl;
}



/*A faster function to evaluate the derivatives of the matrix: subtracts matrices from each other*/
void DaliGen::DMatrix( matrixfiller p, vector<double> eval, int par, int matrixdim, gsl_matrix* result, double step)
{

    DaliCheckDataDecl();


    if(par >= eval.size())
    {
        cout << "Deriving with respect to non existing parameter." << endl;
        abort();
    }

    //the matrix result must already be allocated

    vector<double> plusdel(eval);
    vector<double> minusdel(eval);

    plusdel.at(par) += step;
    minusdel.at(par) -= step;

    gsl_matrix* plusmat = gsl_matrix_calloc(matrixdim, matrixdim);
    gsl_matrix* minusmat = gsl_matrix_calloc(matrixdim, matrixdim);


    //first calculating the complete plus mat, since typically, it will be faster to compute one setup first, and then all entries, instead of computing the setup for each entry
    for(int i = 0; i < matrixdim; i ++)
    {
        for(int j = i; j < matrixdim; j++)
        {
            //since matrix is symmetric
            double hilf = (this->*p)(plusdel, i, j);
            gsl_matrix_set(plusmat, i, j, hilf);
            gsl_matrix_set(plusmat, j, i, hilf);
        }
    }


    for(int i = 0; i < matrixdim; i ++)
    {
        for(int j = i; j < matrixdim; j++)
        {
            //since matrix is assumed to be symmetric
            double hilf = (this->*p)(minusdel, i, j);
            gsl_matrix_set(minusmat, i, j, hilf);
            gsl_matrix_set(minusmat, j, i, hilf);
        }
    }

    gsl_matrix_sub (plusmat, minusmat); //subtracts minusmat from plusmat and stores the result in plusmat
    gsl_matrix_scale (plusmat, 0.5/step); // multiplies the matrix difference, in order to get the derivative

    gsl_matrix_memcpy(result,plusmat); //copies the derived matrix to the result matrix

    gsl_matrix_free(plusmat);
    gsl_matrix_free(minusmat);


}



/*Takes double derivative of one matrix, by subtracting matrices, speed-up version: Only works for the covariance matrix now.*/
double DaliGen::secOneMat(matrixfiller p, vector<double> eval, int par, int matrixdim, gsl_matrix* result, double step)
{
    if(par > eval.size())
    {
        cout << "Taking derivative with respect to non-Existing parameter" << endl;
        abort();
    }

    vector<double> plusdel(eval);
    vector<double> minusdel(eval);

    plusdel.at(par) += step;
    minusdel.at(par)-= step;

    gsl_matrix* plusmat = gsl_matrix_calloc(noofdat,noofdat);
    gsl_matrix* minusmat = gsl_matrix_calloc(noofdat, noofdat);
    gsl_matrix* mat = gsl_matrix_calloc(noofdat,noofdat);
    
    gsl_matrix_memcpy(mat,DataCov);
    
    if(!is_diag)
    {
    for(int i = 0; i < noofdat; i++)
    {
        for(int j = i; j < noofdat; j++)
        {
            double hilf = (this->*p)(plusdel,i,j);
            gsl_matrix_set(plusmat, i, j, hilf);
            gsl_matrix_set(plusmat, j, i, hilf);
        }
    }

    for(int i = 0; i < noofdat; i++)
    {
        for(int j = i; j < noofdat; j++)
        {
            double hilf = (this->*p)(minusdel,i,j);
            gsl_matrix_set(minusmat, i, j, hilf );
            gsl_matrix_set(minusmat, j, i, hilf );
        }
    }

    }
    else
    {
      
      for(int i = 0; i < noofdat; i++)
      {
            double hilf = (this->*p)(plusdel,i,i);
            gsl_matrix_set(plusmat, i, i, hilf);        
      }

      for(int i = 0; i < noofdat; i++)
      {
            double hilf = (this->*p)(minusdel,i,i);
            gsl_matrix_set(minusmat, i, i, hilf );       
      }      
           
    }

    /*the following lines calculate derivatives, equivalently to the following element-wise case:*/
    //return ((this->*p)(plusdel,i,j) + (this->*p)(minusdel,i,j) -2.0*(this->*p)(eval,i,j))/step/step;

    gsl_matrix_add (plusmat, minusmat); //adds minusmat to plusmat and stores the result in plusmat
    gsl_matrix_scale (mat, 2.0);
    gsl_matrix_sub (plusmat, mat);  //subtracts the scaled mat from plus mat and stores in plusmat

    gsl_matrix_scale (plusmat, pow(step, -2.0)); // multiplies the matrix difference, in order to get the derivative

    gsl_matrix_memcpy(result,plusmat); //copies the derived matrix to the result matrix

    gsl_matrix_free(plusmat);
    gsl_matrix_free(minusmat);
    gsl_matrix_free(mat);
}


/*Second derivatives of a matrix, if p points to the function that calculates the entries of the matrix. Result is the second derivative. Subtracts whole matrices from each other. Speed-up version: Only works for the covariance matrix now. Will make errors if p does not point to DataCovMat*/
void DaliGen::DDMatrix(matrixfiller p, vector<double> eval, int par1, int par2, int matrixdim, gsl_matrix* result, double step)
{

    DaliCheckDataDecl();

    /*The matrix result needs to be already allocated*/

    if( (par1 > eval.size()) || (par2 > eval.size()))
    {
        cout << "Taking derivative with respect to non-Existing parameter" << endl;
        abort();
    }

    /*calculating double deriv wrt par 1, always need this*/
    gsl_matrix* doublederivs1 = gsl_matrix_calloc(noofdat, noofdat);
    secOneMat(p,eval,par1,noofdat,doublederivs1, step);


    if(par1==par2)
    {
        gsl_matrix_memcpy(result,doublederivs1);
        gsl_matrix_free(doublederivs1);
    }


    else
    {

        gsl_matrix* doublederivs2 = gsl_matrix_calloc(noofdat, noofdat);
        secOneMat(p,eval,par2,noofdat,doublederivs2, step);

	//matrices are initialized to be zero everywhere
        gsl_matrix* plusmat = gsl_matrix_calloc(noofdat, noofdat);
        gsl_matrix* minusmat = gsl_matrix_calloc(noofdat, noofdat); 	
	
        gsl_matrix* mat = gsl_matrix_calloc(noofdat, noofdat);	
        gsl_matrix_memcpy(mat,DataCov);


        vector<double> plusdel(eval);
        vector<double> minusdel(eval);

        plusdel.at(par1) += step;
        plusdel.at(par2) += step;
        minusdel.at(par1)-= step;
        minusdel.at(par2)-= step;


	if(!is_diag)
	{
          for(int i = 0; i < noofdat; i ++)
          {
            for(int j = i; j < noofdat; j++)
            {
                double hilf =  (this->*p)(plusdel,i,j) ;
                gsl_matrix_set(plusmat, i, j, hilf);
                gsl_matrix_set(plusmat, j, i, hilf);
            }
          }
          
          for(int i = 0; i < noofdat; i ++)
          {
            for(int j = i; j < noofdat; j++)
            {

                double hilf =  (this->*p)(minusdel,i,j) ;
                gsl_matrix_set(minusmat, i, j, hilf);
                gsl_matrix_set(minusmat, j, i, hilf);
            }
        }
        
	}
	
	else
	{
          for(int i = 0; i < noofdat; i ++)
          {
                 double hilf =  (this->*p)(plusdel,i,i) ;
                 gsl_matrix_set(plusmat, i, i, hilf);
           //the non-diagonal entries were initialized to be zero and stay zero if diagonal
          }
          
          for(int i = 0; i < noofdat; i ++)
          {
                double hilf =  (this->*p)(minusdel,i,i) ;
                gsl_matrix_set(minusmat, i, i, hilf);
             //the non-diagonal entries were initialized to be zero and stay zero if diagonal
            
           }
          
          
          
	  
	}
        
        



        gsl_matrix_add(plusmat,minusmat);
        gsl_matrix_scale(mat, 2.0);
        gsl_matrix_sub(plusmat,mat);
        gsl_matrix_scale(plusmat, pow(step,-2.0));

        gsl_matrix_sub(plusmat, doublederivs1);
        gsl_matrix_sub(plusmat, doublederivs2);
        gsl_matrix_scale(plusmat,0.5);

        gsl_matrix_memcpy(result,plusmat);

        gsl_matrix_free(plusmat);
        gsl_matrix_free(minusmat);
        gsl_matrix_free(mat);
        gsl_matrix_free(doublederivs1);
        gsl_matrix_free(doublederivs2);

    }

}



/*Third derivative of a matrix. You might want to adjust the parameter "step" which adjusts the size of the finite difference step.*/
/*void DaliGen::DDDMatrix(matrixfiller p, vector<double> eval, int par1, int par2, int par3, int matrixdim, gsl_matrix* result, double step)
{

  
  
    DaliCheckDataDecl();

    //the matrix result needs to be already allocated

    if( (par1 > eval.size()) || (par2 > eval.size()) || (par3 > eval.size()))
    {
        cout << "Taking derivative with respect to non-Existing parameter" << endl;
        abort();
    }

    
     if((par1 == par2) && (par2 == par3))
     {
       
       
       
      vector<double> plusdel(eval);
      vector<double> minusdel(eval);

      vector<double> plus2del(eval);
      vector<double> minus2del(eval);
      
      plusdel.at(par3) += step;
      minusdel.at(par3)-= step; 
      
      plus2del.at(par3) += 2.0*step;
      minus2del.at(par3)-= 2.0*step;
  
      gsl_matrix* plusmat = gsl_matrix_calloc(matrixdim, matrixdim);
      gsl_matrix* minusmat = gsl_matrix_calloc(matrixdim, matrixdim);
      gsl_matrix* plus2mat = gsl_matrix_calloc(matrixdim, matrixdim);
      gsl_matrix* minus2mat = gsl_matrix_calloc(matrixdim, matrixdim);  
      
      
      for(int i = 0; i < matrixdim; i ++)
      {
        for(int j = i; j < matrixdim; j++)
        {
            //since matrix is symmetric
            double hilf = (this->*p)(plusdel, i, j);
            gsl_matrix_set(plusmat, i, j, hilf);
            gsl_matrix_set(plusmat, j, i, hilf);
        }
      }

      for(int i = 0; i < matrixdim; i ++)
      {
        for(int j = i; j < matrixdim; j++)
        {
            //since matrix is symmetric
            double hilf = (this->*p)(minusdel, i, j);
            gsl_matrix_set(minusmat, i, j, hilf);
            gsl_matrix_set(minusmat, j, i, hilf);
        }
      }

      
      for(int i = 0; i < matrixdim; i ++)
      {
        for(int j = i; j < matrixdim; j++)
        {
            //since matrix is symmetric
            double hilf = (this->*p)(plus2del, i, j);
            gsl_matrix_set(plus2mat, i, j, hilf);
            gsl_matrix_set(plus2mat, j, i, hilf);
        }
      }

      for(int i = 0; i < matrixdim; i ++)
      {
        for(int j = i; j < matrixdim; j++)
        {
            //since matrix is symmetric
            double hilf = (this->*p)(minus2del, i, j);
            gsl_matrix_set(minus2mat, i, j, hilf);
            gsl_matrix_set(minus2mat, j, i, hilf);
        }
      }
      
            
      
      gsl_matrix_scale(plus2mat, 0.5);
      gsl_matrix_scale(minus2mat, 0.5);
      
      gsl_matrix_add (minusmat, plus2mat); //adds arg2 to arg1 and stores the result in arg1
      gsl_matrix_sub (minusmat, plusmat);
      gsl_matrix_sub (minusmat, minus2mat);
      
      
      cout << "Check stability" << endl;
      //Check for stability befor dividing by the tinynumber step**3 
      for(int i = 0; i < matrixdim; i ++)
      {
        for(int j = i; j < matrixdim; j++)
        {
	  double triple_tiny = 1e-10*gsl_matrix_get(DataCov, i,j);
	  
           if(fabs(gsl_matrix_get(minusmat, i, j)) < fabs(triple_tiny))
	   {
	     gsl_matrix_set(minusmat,i,j,0.0);
	  }

        }
      }
      
      gsl_matrix_scale(minusmat, pow(step, -3.0));
      
      gsl_matrix_memcpy(result,minusmat);
      
      gsl_matrix_free(plusmat);
      gsl_matrix_free(minusmat);
      gsl_matrix_free(plus2mat);
      gsl_matrix_free(minus2mat);
      
     }
     
     
     else
     {
    
    
    
    vector<double> plusdel(eval);
    vector<double> minusdel(eval);

    plusdel.at(par3) += step;
    minusdel.at(par3)-= step;

    gsl_matrix* plusmat = gsl_matrix_calloc(matrixdim, matrixdim);
    gsl_matrix* minusmat = gsl_matrix_calloc(matrixdim, matrixdim);

    //maybe better use the default step 2e-2 here. After all, the accuracy with which the sec. deriv can be estimated should be exploited. Instead of using the crude third deriv step
    DDMatrix(p, plusdel, par1, par2, matrixdim, plusmat,5e-3);
    DDMatrix(p, minusdel, par1, par2, matrixdim, minusmat, 5e-3);


    gsl_matrix_sub (plusmat, minusmat); //subtracts minusmat from plusmat and stores the result in plusmat
    
       cout << "Check stability" << endl;
      //Check for stability befor dividing by the tinynumber step 
      for(int i = 0; i < matrixdim; i ++)
      {
        for(int j = i; j < matrixdim; j++)
        {
	  double triple_tiny = 1e-10*gsl_matrix_get(DataCov, i,j);
	  
           if(fabs(gsl_matrix_get(minusmat, i, j)) < fabs(triple_tiny))
	   {
	     gsl_matrix_set(minusmat,i,j,0.0);
	  }

        }
      }
    
    
    gsl_matrix_scale (plusmat, 0.5/step); // multiplies the matrix difference, in order to get the derivative

    gsl_matrix_memcpy(result,plusmat);

    gsl_matrix_free(plusmat);
    gsl_matrix_free(minusmat);

    
     }


}*/


/*more accurate, although numerically more costly evaluation of the third derivatives.*/
void DaliGen::DDDMatrix(matrixfiller p, vector<double> eval, int par1, int par2, int par3, int matrixdim, gsl_matrix* result, double step)
{
  
  double h = step; //different finite difference steps
  double H = step;
  double cH = step; // ch = 'curly H', also a finite difference step
  
  
  /*the naming of these vectors goes:
   ph: all steps h are added
   pH: all steps H are added
   pcH: all steps cH are added
   mh: all steps h are subtracted
   mH: all steps H are subtracted
   mcH: all steps cH are subtracted*/
  vector<double> ph_pH_pcH(eval);
  ph_pH_pcH[par1] += h;
  ph_pH_pcH[par2] += H;
  ph_pH_pcH[par3] += cH;
  
  gsl_matrix* ph_pH_pcH_mat = gsl_matrix_calloc(matrixdim, matrixdim);
  
  for(int i = 0; i < matrixdim; i++)
  {
    for(int j = i; j < matrixdim; j++)
    {
      double hilf = (this->*p)(ph_pH_pcH, i, j);
      gsl_matrix_set(ph_pH_pcH_mat,i,j,hilf);
      gsl_matrix_set(ph_pH_pcH_mat,j,i,hilf);
    }
  }
  
  
  
  vector<double> ph_pH_mcH(eval);
  ph_pH_mcH[par1] += h;
  ph_pH_mcH[par2] += H;
  ph_pH_mcH[par3] -= cH; 
  
  gsl_matrix* ph_pH_mcH_mat = gsl_matrix_calloc(matrixdim, matrixdim);
  
  for(int i = 0; i < matrixdim; i++)
  {
    for(int j = i; j < matrixdim; j++)
    {
      double hilf = (this->*p)(ph_pH_mcH, i, j);
      gsl_matrix_set(ph_pH_mcH_mat,i,j,hilf);
      gsl_matrix_set(ph_pH_mcH_mat,j,i,hilf);
    }
  }
  
  
  gsl_matrix_sub (ph_pH_pcH_mat, ph_pH_mcH_mat); //subtracts ph_pH_mcH_mat from ph_pH_pcH_mat and stores the result in ph_pH_pcH_mat
  
  
  
  vector<double> ph_mH_pcH(eval);
  ph_mH_pcH[par1] += h;
  ph_mH_pcH[par2] -= H;
  ph_mH_pcH[par3] += cH;
  
  gsl_matrix* ph_mH_pcH_mat = gsl_matrix_calloc(matrixdim, matrixdim);
  
  for(int i = 0; i < matrixdim; i++)
  {
    for(int j = i; j < matrixdim; j++)
    {
      double hilf = (this->*p)(ph_mH_pcH, i, j);
      gsl_matrix_set(ph_mH_pcH_mat,i,j,hilf);
      gsl_matrix_set(ph_mH_pcH_mat,j,i,hilf);
    }
  }
  
  
  gsl_matrix_sub(ph_pH_pcH_mat, ph_mH_pcH_mat);
  
  
  vector<double> ph_mH_mcH(eval);
  ph_mH_mcH[par1] += h;
  ph_mH_mcH[par2] -= H;
  ph_mH_mcH[par3] -= cH;
  
  gsl_matrix* ph_mH_mcH_mat = gsl_matrix_calloc(matrixdim, matrixdim);
  
  for(int i = 0; i < matrixdim; i++)
  {
    for(int j = i; j < matrixdim; j++)
    {
      double hilf = (this->*p)(ph_mH_mcH, i, j);
      gsl_matrix_set(ph_mH_mcH_mat,i,j,hilf);
      gsl_matrix_set(ph_mH_mcH_mat,j,i,hilf);
    }
  }
  
  
  gsl_matrix_add(ph_pH_pcH_mat, ph_mH_mcH_mat);
  
  
  
  vector<double> mh_pH_pcH(eval);
  mh_pH_pcH[par1] -= h;
  mh_pH_pcH[par2] += H;
  mh_pH_pcH[par3] += cH;
  
  gsl_matrix* mh_pH_pcH_mat = gsl_matrix_calloc(matrixdim, matrixdim);
  
  for(int i = 0; i < matrixdim; i++)
  {
    for(int j = i; j < matrixdim; j++)
    {
      double hilf = (this->*p)(mh_pH_pcH, i, j);
      gsl_matrix_set(mh_pH_pcH_mat,i,j,hilf);
      gsl_matrix_set(mh_pH_pcH_mat,j,i,hilf);
    }
  }
  
  
  gsl_matrix_sub(ph_pH_pcH_mat, mh_pH_pcH_mat);
  
  
  vector<double> mh_pH_mcH(eval);
  mh_pH_mcH[par1] -= h;
  mh_pH_mcH[par2] += H;
  mh_pH_mcH[par3] -= cH;
  
  gsl_matrix* mh_pH_mcH_mat = gsl_matrix_calloc(matrixdim, matrixdim);
  
  for(int i = 0; i < matrixdim; i++)
  {
    for(int j = i; j < matrixdim; j++)
    {
      double hilf = (this->*p)(mh_pH_mcH, i, j);
      gsl_matrix_set(mh_pH_mcH_mat,i,j,hilf);
      gsl_matrix_set(mh_pH_mcH_mat,j,i,hilf);
    }
  }
  
  gsl_matrix_add(ph_pH_pcH_mat, mh_pH_mcH_mat);
  
  
  vector<double> mh_mH_mcH(eval);
  mh_mH_mcH[par1] -= h;
  mh_mH_mcH[par2] -= H;
  mh_mH_mcH[par3] -= cH;
  
  gsl_matrix* mh_mH_mcH_mat = gsl_matrix_calloc(matrixdim, matrixdim);
  
  for(int i = 0; i < matrixdim; i++)
  {
    for(int j = i; j < matrixdim; j++)
    {
      double hilf = (this->*p)(mh_mH_mcH, i, j);
      gsl_matrix_set(mh_mH_mcH_mat,i,j,hilf);
      gsl_matrix_set(mh_mH_mcH_mat,j,i,hilf);
    }
  }
  
  
  gsl_matrix_sub(ph_pH_pcH_mat, mh_mH_mcH_mat);
  
  vector<double> mh_mH_pcH(eval);
  mh_mH_pcH[par1] -= h;
  mh_mH_pcH[par2] -= H;
  mh_mH_pcH[par3] += cH;
  
  gsl_matrix* mh_mH_pcH_mat = gsl_matrix_calloc(matrixdim, matrixdim);
  
  for(int i = 0; i < matrixdim; i++)
  {
    for(int j = i; j < matrixdim; j++)
    {
      double hilf = (this->*p)(mh_mH_pcH, i, j);
      gsl_matrix_set(mh_mH_pcH_mat,i,j,hilf);
      gsl_matrix_set(mh_mH_pcH_mat,j,i,hilf);
    }
  } 
  
  
 gsl_matrix_add(ph_pH_pcH_mat, mh_mH_pcH_mat); 
  

      //Check for stability befor dividing by the tinynumber step 
      for(int i = 0; i < matrixdim; i ++)
      {
        for(int j = i; j < matrixdim; j++)
        {
	  double triple_tiny = 1e-10*gsl_matrix_get(DataCov, i,j);
	  
           if(fabs(gsl_matrix_get(ph_pH_pcH_mat, i, j)) < fabs(triple_tiny))
	   {
	     gsl_matrix_set(ph_pH_pcH_mat,i,j,0.0);
	  }

        }
      }
 
     gsl_matrix_scale (ph_pH_pcH_mat, 1./(8.0*h*H*cH   )); // multiplies the matrix difference, in order to get the derivative

    gsl_matrix_memcpy(result,ph_pH_pcH_mat);
    
    
    gsl_matrix_free(ph_pH_pcH_mat);
    gsl_matrix_free(ph_pH_mcH_mat);
    gsl_matrix_free(ph_mH_pcH_mat);
    gsl_matrix_free(ph_mH_mcH_mat);
    
    gsl_matrix_free(mh_pH_pcH_mat);
    gsl_matrix_free(mh_pH_mcH_mat);
    gsl_matrix_free(mh_mH_pcH_mat);
    gsl_matrix_free(mh_mH_mcH_mat);
    
  
}




/*Calculating some shorthand matrices: Ctwiddle*/
void DaliGen::FillCtwiddle()
{

    CheckDataSize();

    if(!covcheck)
        FillDataCovMat();
    if(!checkmatderivs)
        Calc_DCs();

    for(int par = 0; par < theodim; par++)
    {
        //multiplies the inverse data cov matrix and DCov, and stores it in Ctwiddle
        mult2mat(InvDataCov, DCov.at(par), noofdat,Ctwiddle.at(par),is_diag );
	
	if(lowmem) gsl_matrix_free(DCov.at(par));
	
    }


    //set flag that these matrices were calculated
    twiddlecheck = 1;
}






/*Calculate all the derivatives of the data covariance matrix*/
void DaliGen::Calc_DCs()
{

    CheckDataSize();

    /*checking whether the data covariance matrix has been calculated.
      * If not, then DALI calculates it now.*/
    if(!covcheck)
        FillDataCovMat();

    if(covmatdep != 0)
    {
        cout << "Calculating 1st derivatives of the Covariance matrix." << endl;

        for(int par = 0; par < theodim; par++)
        {
            if(verbosity > 1)
            {
                cout << "Deriving wrt parameter " << par << endl;
            }

            // the matrices where the derivatives are stored were in the constructor
            //calculating the derivatives
            DMatrix(&DaliGen::DataCovMatrix, fid, par, noofdat, DCov.at(par), fstep_C);
	    
	    if(RecommendedChecks)
	    {
	       stringstream suff;  
               suff << par;
               string suffstring = (string)"DCov" + suff.str();     
               PlotMatrix(DCov[par], noofdat, suffstring); 
	    }
	    

        }
    }
    //set a flag that tells other functions that the derivatives have been calculated
    checkmatderivs = 1;
}


/*Calculate all the second derivatives of the data covariance matrix*/
void DaliGen::Calc_DDCs()
{
    CheckDataSize();

    //checking whether the data covariance matrix has been calculated. If not, then DALI calculates it now.
    if(!covcheck)
        FillDataCovMat();

    if(covmatdep != 0)
    {
        cout << "Calculating 2nd derivatives of the Covariance matrix." << endl;



        //loop only over a triangle of parameter combinations, since the other triangle will be the same
        for(int par = 0; par < theodim; par++)
        {
            for(int par2 = 0; par2 <= par; par2++)
            {
                if(verbosity > 1)
                {
                    cout << "Deriving wrt parameter " << par << ", " << par2 << endl;
                }

                //calculating the matrices
                DDMatrix(&DaliGen::DataCovMatrix, fid, par, par2, noofdat, DDCov.at(par).at(par2),dstep_C);

                //multiplies two matrices, stores result in CDoubletwiddle
                mult2mat(InvDataCov, DDCov.at(par).at(par2), noofdat, CDoubletwiddle.at(par).at(par2),is_diag);

                //setting the matrices, that are the same, and don't need to be calculated a second time
                DDCov.at(par2).at(par) = DDCov.at(par).at(par2);
                CDoubletwiddle.at(par2).at(par) = CDoubletwiddle.at(par).at(par2);
		
		
	       if(RecommendedChecks)
	       {
	        stringstream suff;  
                suff << par << par2;
                string suffstring = (string)"DDCov" + suff.str();     
                PlotMatrix(DDCov[par][par2], noofdat, suffstring); 
	       }

		if(lowmem) gsl_matrix_free(DDCov.at(par).at(par2));
		
            }
        }

    }
    //set flag to let other functions know these things were calculated
    checkmatdoublederivs = 1;
}


/*calculate the third derivatives of the covariance matrix*/
void DaliGen::Calc_DDDCs()
{
    CheckDataSize();
    //checking whether the data covariance matrix has been calculated. If not, then DALI calculates it now.
    if(!covcheck)
        FillDataCovMat();

    if(covmatdep)
    {
        if(derivorder==3)
        {
            cout << "Calculating 3rd derivatives of the Covariance matrix." << endl;


            for(int par = 0; par < theodim; par++)
            {
                for(int par2 = 0; par2 <= par; par2++)
                {
                    for(int par3 = 0; par3 <= par2; par3++)
                    {
                        if(verbosity > 1)
                        {
                            cout << "Deriving wrt parameter " << par << ", " << par2 << ", " << par3 << endl;
                        }
                        DDDMatrix(&DaliGen::DataCovMatrix, fid, par, par2, par3, noofdat, DDDCov.at(par).at(par2).at(par3),tstep_C);

                        //multiplies two matrices, stores result in CTripletwiddle
                        mult2mat(InvDataCov, DDDCov.at(par).at(par2).at(par3), noofdat, CTripletwiddle.at(par).at(par2).at(par3),is_diag);
                        /*exploiting pointer arithmetic in order to set the other permutations of the 3 derivatives*/


                        DDDCov.at(par2).at(par3).at(par)          = DDDCov.at(par).at(par2).at(par3);
                        CTripletwiddle.at(par2).at(par3).at(par)  = CTripletwiddle.at(par).at(par2).at(par3);

                        DDDCov.at(par3).at(par).at(par2)          = DDDCov.at(par).at(par2).at(par3);
                        CTripletwiddle.at(par3).at(par).at(par2)  = CTripletwiddle.at(par).at(par2).at(par3);

                        DDDCov.at(par3).at(par2).at(par)          = DDDCov.at(par).at(par2).at(par3);
                        CTripletwiddle.at(par3).at(par2).at(par)  = CTripletwiddle.at(par).at(par2).at(par3);

                        DDDCov.at(par2).at(par).at(par3)          = DDDCov.at(par).at(par2).at(par3);
                        CTripletwiddle.at(par2).at(par).at(par3)  = CTripletwiddle.at(par).at(par2).at(par3);

                        DDDCov.at(par).at(par3).at(par2)          = DDDCov.at(par).at(par2).at(par3);
                        CTripletwiddle.at(par).at(par3).at(par2)  = CTripletwiddle.at(par).at(par2).at(par3);

			if(RecommendedChecks)
	                {
	                 stringstream suff;  
                         suff << par << par2 << par3;
                         string suffstring = (string)"DDDCov" + suff.str();     
                         PlotMatrix(DDDCov[par][par2][par3], noofdat, suffstring); 
	                }
			
			
                        if(lowmem) gsl_matrix_free(DDDCov.at(par2).at(par3).at(par));
			
                        /*  printmatrix(DDDCov.at(par).at(par2).at(par3),noofdat, noofdat);
                          cout << endl;
                         printmatrix(DDDCov.at(par3).at(par2).at(par),noofdat, noofdat);  */
                    }
                }
            }

        }
    }
//setting flag to let the other functions know this stuff was calculated
    checkmattriplederivs = 1;

}



/** Calculates the Fisher matrix and stores it in the gsl_matrix* DaliGen_object.fisher .
 *  The code will automatically check whether everything has been computed that it needs for the fisher matrix. If not, it will do so. If the Fisher matrix is singular, it doesn't matter much for the DALI contours of the likelihood. The code will however warn you that you shall not use the fisher matrix as a proposal distribution for the MCMC sampler.
   **/
void DaliGen::CalcFisher()
{
    CheckDataSize();

    //check that everything is there, that is needed
    DaliCheckDataDecl();
    if(!checkd)
        CalcFirstDmod();
    if(!covcheck)
        FillDataCovMat();
    if(!checkmatderivs)
        Calc_DCs();
    if(!twiddlecheck)
        FillCtwiddle();

    for(int alpha = 0; alpha < theodim; alpha ++)
    {
        for(int beta = 0; beta < theodim; beta ++)
        {
            double a = vec1Matvec2( dmod[alpha], InvDataCov,dmod[beta], noofdat);

            //create an intermediate matrix of which the trace will be taken
            gsl_matrix* fortrace = gsl_matrix_alloc(noofdat,noofdat);
            gsl_matrix* fortrace_w = gsl_matrix_alloc(noofdat,noofdat);
            mult2mat(Ctwiddle[alpha], Ctwiddle[beta], noofdat, fortrace,is_diag);

            //including data point weights:
            mult2mat(fortrace, datapoint_weights,noofdat,fortrace_w,is_diag);

            double b =  0.5* trace(fortrace_w,noofdat);

            double hilf = a + b;

            //setting the fisher matrix
            gsl_matrix_set(fisher, alpha, beta, hilf);

            gsl_matrix_free(fortrace);
            gsl_matrix_free(fortrace_w);
        }
    }

    if(theodim < 15)
    {
        cout << "The Fisher matrix..........................." << endl;
        printmatrix(fisher,theodim,theodim);
        cout << endl;
    }

    if(theodim >= 15)
    {
        cout << "Not displaying fisher matrix: it would be too large." << endl;
    }

    double det = determinant(fisher,theodim);
    if(fabs(det) < 1e-10)
    {
        cout << "WARNING: Fisher matrix is singular. Don't use it as an input for the MCMC sampler." << endl;
        cout << "It's determinant is det F = " << det << endl;
    }


}


/**Calculates the higher order Dali-Tensors that approximate the likelihood beyond what the fisher matrix can do. It checks whether everything has been calculated that is needed. If not, it will do so. The resulting Dali tensors are stored in DaliGen_object.Sabg, DaliGen_object.Qabgd, DaliGen_object.Pabgde, DaliGen_object.Habgdef, where DaliGen_object is the name of the object that you created in your main()
 ATTENTION: depending on whether you requested second or third order derivatives for your model and/or your likelihood, the tensor Qabgd will change, depending on the number of derivatives.  **/
void DaliGen::Calc_BeyondFish()
{

    //check that everything is there that is needed
    DaliCheckDataDecl();
    if(!checkd)
        CalcFirstDmod();
    if(!checkdd)
        CalcScndDmod();
    if(!checkddd)
        CalcThrdDmod();
    if(!covcheck)
        FillDataCovMat();
    if(!checkmatderivs)
        Calc_DCs();
    if(!twiddlecheck)
        FillCtwiddle();
    if(!checkmatdoublederivs)
        Calc_DDCs();
    if(!checkmattriplederivs)
        Calc_DDDCs();
    
    


    //looping over parameter space, calculating all derivatives wrt parameters
    for(int alpha = 0; alpha < theodim; alpha ++)
    {
        for(int beta = 0; beta < theodim; beta ++)
        {
            for(int gamma = 0; gamma < theodim; gamma++)
            {   /*The order of the indices doesn't matter. Any permutation can be chosen.
                  It only needs to be kept constant during the execution of the code*/

                  if(verbosity > 1)
		  {
		    cout << "Calculating S_abg tensor entries " << alpha << " " << beta
		         << " " << gamma << endl;
		    
		  }
                  
                 if(!Sabg_compute_flags[alpha][beta][gamma])
		 {
                 
                  if(covmatdep)
                  {
                      gsl_matrix* Shelper = gsl_matrix_alloc(noofdat,noofdat);
                      mult2mat(Ctwiddle[gamma],CDoubletwiddle[alpha][beta],noofdat,Shelper,is_diag);

                      gsl_matrix* Shelper_w = gsl_matrix_alloc(noofdat,noofdat);
                      mult2mat(Shelper, datapoint_weights,noofdat,Shelper_w,is_diag);

                      Sabg[alpha][beta][gamma] += 1.5*trace( Shelper_w,noofdat);
		      

                      gsl_matrix_free(Shelper_w);
                      gsl_matrix_free(Shelper);
                  }

                  if(modeldep)
                  {
                      Sabg[alpha][beta][gamma] += 3.0 * vec1Matvec2( ddmod[alpha][beta], InvDataCov,dmod[gamma], noofdat);
                  } 
     
     
                   Sabg_compute_flags[alpha][beta][gamma] = 1;
     
		   //exploiting the alpha-beta symmetry
                   Sabg[beta][alpha][gamma]= Sabg[alpha][beta][gamma];
		   Sabg_compute_flags[beta][alpha][gamma] = 1;
	        } 
                
	    }
	    
	}
    }
                

                
   for(int alpha = 0; alpha < theodim; alpha ++)
    {
        for(int beta = 0; beta < theodim; beta ++)
        {
            for(int gamma = 0; gamma < theodim; gamma++)
            {   /*The order of the indices doesn't matter. Any permutation can be chosen.
                  It only needs to be kept constant during the execution of the code*/             

                for(int delta = 0; delta < theodim; delta++)
                {
		  
		  if(verbosity > 1)
		  {
		    cout << "Calculating Q_abgd tensor entries and higher (if requested)" << alpha << " " << beta << " " << gamma << " " << delta << endl;	    
		  }
		  
		  /*compute_flags for Qabgd must still be implemented*/
                    if(covmatdep)
                    {
                        gsl_matrix* Qhelper = gsl_matrix_alloc(noofdat,noofdat);
                        gsl_matrix* Qhelper_w = gsl_matrix_alloc(noofdat,noofdat);

                        mult2mat(CDoubletwiddle[gamma][delta],
                                 CDoubletwiddle[alpha][beta], noofdat, Qhelper,is_diag);

                        mult2mat(Qhelper, datapoint_weights,noofdat,Qhelper_w,is_diag);

                        Qabgd[alpha][beta][gamma][delta] +=1.5*trace(Qhelper_w,noofdat);

                        gsl_matrix_free(Qhelper);
                        gsl_matrix_free(Qhelper_w);
                    }


                    if(modeldep)
                    {
                        Qabgd[alpha][beta][gamma][delta] +=
                            3.0* vec1Matvec2(
                                ddmod[alpha][beta], InvDataCov,ddmod[delta][gamma]
                                , noofdat);
                    }


                    if(derivorder == 3)
                    {

                        if(covmatdep)
                        {
                            gsl_matrix* Q2helper = gsl_matrix_alloc(noofdat, noofdat);
                            gsl_matrix* Q2helper_w = gsl_matrix_alloc(noofdat, noofdat);
                            mult2mat(CTripletwiddle[beta][gamma][delta],
                                     Ctwiddle[alpha] ,noofdat, Q2helper,is_diag);

                            mult2mat(Q2helper, datapoint_weights,noofdat,Q2helper_w,is_diag);

                            Qabgd[alpha][beta][gamma][delta] += 2.0*trace(Q2helper_w,noofdat);

                            gsl_matrix_free(Q2helper_w);
                            gsl_matrix_free(Q2helper);
                        }

                        if(modeldep)
                        {
                            Qabgd[alpha][beta][gamma][delta] +=
                                4.0*vec1Matvec2( dddmod[alpha][beta][gamma], InvDataCov,dmod[delta], noofdat);
                        }



                        for(int epsil = 0; epsil < theodim; epsil++)
                        {
			  
			  if(!Pabgde_compute_flags[alpha][beta][gamma][delta][epsil])
			  {
                            if(covmatdep)
                            {
                                gsl_matrix* Phelper = gsl_matrix_alloc(noofdat,noofdat);
                                gsl_matrix* Phelper_w = gsl_matrix_alloc(noofdat,noofdat);

                                mult2mat(CDoubletwiddle[delta][epsil],
                                         CTripletwiddle[alpha][beta][gamma],noofdat,Phelper,is_diag);
                                mult2mat(Phelper, datapoint_weights, noofdat,Phelper_w,is_diag);

                                Pabgde[alpha][beta][gamma][delta][epsil] += 5.0*trace(Phelper_w,noofdat);

                                gsl_matrix_free(Phelper_w);
                                gsl_matrix_free(Phelper);
                            }

                            if(modeldep)
                            {
                                Pabgde[alpha][beta][gamma][delta][epsil] +=
                                    10.0*vec1Matvec2( dddmod[alpha][beta][gamma], InvDataCov,ddmod[delta][epsil], noofdat);
                            }
                            
                          Pabgde_compute_flags[alpha][beta][gamma][delta][epsil] = 1;
  
			  Pabgde[alpha][beta][gamma][epsil][delta] = Pabgde[alpha][beta][gamma][delta][epsil];			  
			  Pabgde_compute_flags[alpha][beta][gamma][epsil][delta] = 1;
			  
			 Pabgde[alpha][gamma][beta][epsil][delta] 
			        = Pabgde[alpha][gamma][beta][delta][epsil] 
			        = Pabgde[alpha][beta][gamma][delta][epsil];
			 
			  Pabgde_compute_flags[alpha][gamma][beta][delta][epsil] = 1;
			  Pabgde_compute_flags[alpha][gamma][beta][epsil][delta] = 1;
			  
			  
			  Pabgde[beta][alpha][gamma][delta][epsil] 
			  = Pabgde[beta][alpha][gamma][epsil][delta]
			  = Pabgde[alpha][beta][gamma][delta][epsil];
			  
			  Pabgde_compute_flags[beta][alpha][gamma][delta][epsil] = 1;
			  Pabgde_compute_flags[beta][alpha][gamma][epsil][delta] = 1;
			  
			  Pabgde[gamma][beta][alpha][delta][epsil]
			  = Pabgde[gamma][beta][alpha][epsil][delta]
			  = Pabgde[alpha][beta][gamma][delta][epsil];
			  
			  Pabgde_compute_flags[gamma][beta][alpha][delta][epsil] = 1;
			  Pabgde_compute_flags[gamma][beta][alpha][epsil][delta] = 1;
			  
			  
			   Pabgde[beta][gamma][alpha][delta][epsil]
			  = Pabgde[beta][gamma][alpha][epsil][delta]
			  = Pabgde[alpha][beta][gamma][delta][epsil];
			  
			  Pabgde_compute_flags[beta][gamma][alpha][delta][epsil] = 1;
			  Pabgde_compute_flags[beta][gamma][alpha][epsil][delta] = 1;
                            

			  }

			  
			  			  

                            for(int phi = 0; phi < theodim; phi++)
                            {
                               if(!Habgdef_compute_flags[alpha][beta][gamma][delta][epsil][phi])
			       { 
				 
				 
                                if(covmatdep)
                                {
                                    gsl_matrix* Hhelper = gsl_matrix_alloc(noofdat,noofdat);
                                    gsl_matrix* Hhelper_w = gsl_matrix_alloc(noofdat,noofdat);

                                    mult2mat(CTripletwiddle[alpha][beta][gamma],
                                             CTripletwiddle[delta][epsil][phi],noofdat,Hhelper,is_diag);

                                    mult2mat(Hhelper,datapoint_weights,noofdat,Hhelper_w,is_diag);

                                    Habgdef[alpha][beta][gamma][delta][epsil][phi] += 5.0*trace(Hhelper_w,noofdat);

                                    gsl_matrix_free(Hhelper_w);
                                    gsl_matrix_free(Hhelper);
                                }

                                if(modeldep)
                                {
                                    Habgdef[alpha][beta][gamma][delta][epsil][phi] +=
                                        10.0*vec1Matvec2( dddmod[alpha][beta][gamma],
                                                          InvDataCov,dddmod[delta][epsil][phi],
                                                          noofdat);

                                }
                                
                                Habgdef_compute_flags[alpha][beta][gamma][delta][epsil][phi] = 1;
				
				
				/*Muss nur die alpha. beta gamma permutieren*/
				/*1*/
				
				Habgdef[alpha][beta][gamma][phi][epsil][delta]
				= Habgdef[alpha][beta][gamma][phi][delta][epsil]
				= Habgdef[alpha][beta][gamma][epsil][phi][delta]
				= Habgdef[alpha][beta][gamma][epsil][delta][phi]
				= Habgdef[alpha][beta][gamma][delta][phi][epsil]
				
				= Habgdef[alpha][beta][gamma][delta][epsil][phi];
				
				
				Habgdef_compute_flags[alpha][beta][gamma][phi][epsil][delta] = 1;
				Habgdef_compute_flags[alpha][beta][gamma][phi][delta][epsil] = 1;
				Habgdef_compute_flags[alpha][beta][gamma][epsil][phi][delta] = 1;
				Habgdef_compute_flags[alpha][beta][gamma][epsil][delta][phi] = 1;
				Habgdef_compute_flags[alpha][beta][gamma][delta][phi][epsil] = 1;
                              
                               /*2*/
			       
			       	Habgdef[alpha][gamma][beta][phi][epsil][delta]
				= Habgdef[alpha][gamma][beta][phi][delta][epsil]
				= Habgdef[alpha][gamma][beta][epsil][phi][delta]
				= Habgdef[alpha][gamma][beta][epsil][delta][phi]
				= Habgdef[alpha][gamma][beta][delta][phi][epsil]
				= Habgdef[alpha][gamma][beta][delta][epsil][phi]
				
				= Habgdef[alpha][beta][gamma][delta][epsil][phi];
				
				
				Habgdef_compute_flags[alpha][gamma][beta][phi][epsil][delta] = 1;
				Habgdef_compute_flags[alpha][gamma][beta][phi][delta][epsil] = 1;
				Habgdef_compute_flags[alpha][gamma][beta][epsil][phi][delta] = 1;
				Habgdef_compute_flags[alpha][gamma][beta][epsil][delta][phi] = 1;
				Habgdef_compute_flags[alpha][gamma][beta][delta][phi][epsil] = 1;
				Habgdef_compute_flags[alpha][gamma][beta][delta][epsil][phi] = 1;

				
				/*3*/
			       	Habgdef[gamma][alpha][beta][phi][epsil][delta]
				= Habgdef[gamma][alpha][beta][phi][delta][epsil]
				= Habgdef[gamma][alpha][beta][epsil][phi][delta]
				= Habgdef[gamma][alpha][beta][epsil][delta][phi]
				= Habgdef[gamma][alpha][beta][delta][phi][epsil]
				= Habgdef[gamma][alpha][beta][delta][epsil][phi]
				
				= Habgdef[alpha][beta][gamma][delta][epsil][phi];
				
				
				Habgdef_compute_flags[gamma][alpha][beta][phi][epsil][delta] = 1;
				Habgdef_compute_flags[gamma][alpha][beta][phi][delta][epsil] = 1;
				Habgdef_compute_flags[gamma][alpha][beta][epsil][phi][delta] = 1;
				Habgdef_compute_flags[gamma][alpha][beta][epsil][delta][phi] = 1;
				Habgdef_compute_flags[gamma][alpha][beta][delta][phi][epsil] = 1;
				Habgdef_compute_flags[gamma][alpha][beta][delta][epsil][phi] = 1;
				
                              /*4*/
				Habgdef[beta][alpha][gamma][phi][epsil][delta]
				= Habgdef[beta][alpha][gamma][phi][delta][epsil]
				= Habgdef[beta][alpha][gamma][epsil][phi][delta]
				= Habgdef[beta][alpha][gamma][epsil][delta][phi]
				= Habgdef[beta][alpha][gamma][delta][phi][epsil]
				= Habgdef[beta][alpha][gamma][delta][epsil][phi]
				
				= Habgdef[alpha][beta][gamma][delta][epsil][phi];
				
				
				Habgdef_compute_flags[beta][alpha][gamma][phi][epsil][delta] = 1;
				Habgdef_compute_flags[beta][alpha][gamma][phi][delta][epsil] = 1;
				Habgdef_compute_flags[beta][alpha][gamma][epsil][phi][delta] = 1;
				Habgdef_compute_flags[beta][alpha][gamma][epsil][delta][phi] = 1;
				Habgdef_compute_flags[beta][alpha][gamma][delta][phi][epsil] = 1;
				Habgdef_compute_flags[beta][alpha][gamma][delta][epsil][phi] = 1;
				
				
				/*5*/
				Habgdef[beta][gamma][alpha][phi][epsil][delta]
				= Habgdef[beta][gamma][alpha][phi][delta][epsil]
				= Habgdef[beta][gamma][alpha][epsil][phi][delta]
				= Habgdef[beta][gamma][alpha][epsil][delta][phi]
				= Habgdef[beta][gamma][alpha][delta][phi][epsil]
				= Habgdef[beta][gamma][alpha][delta][epsil][phi]
				
				= Habgdef[alpha][beta][gamma][delta][epsil][phi];
				
				
				Habgdef_compute_flags[beta][gamma][alpha][phi][epsil][delta] = 1;
				Habgdef_compute_flags[beta][gamma][alpha][phi][delta][epsil] = 1;
				Habgdef_compute_flags[beta][gamma][alpha][epsil][phi][delta] = 1;
				Habgdef_compute_flags[beta][gamma][alpha][epsil][delta][phi] = 1;
				Habgdef_compute_flags[beta][gamma][alpha][delta][phi][epsil] = 1;
				Habgdef_compute_flags[beta][gamma][alpha][delta][epsil][phi] = 1;
				
			      /*6*/
			      	Habgdef[gamma][beta][alpha][phi][epsil][delta]
				= Habgdef[gamma][beta][alpha][phi][delta][epsil]
				= Habgdef[gamma][beta][alpha][epsil][phi][delta]
				= Habgdef[gamma][beta][alpha][epsil][delta][phi]
				= Habgdef[gamma][beta][alpha][delta][phi][epsil]
				= Habgdef[gamma][beta][alpha][delta][epsil][phi]
				
				= Habgdef[alpha][beta][gamma][delta][epsil][phi];
				
				
				Habgdef_compute_flags[gamma][beta][alpha][phi][epsil][delta] = 1;
				Habgdef_compute_flags[gamma][beta][alpha][phi][delta][epsil] = 1;
				Habgdef_compute_flags[gamma][beta][alpha][epsil][phi][delta] = 1;
				Habgdef_compute_flags[gamma][beta][alpha][epsil][delta][phi] = 1;
				Habgdef_compute_flags[gamma][beta][alpha][delta][phi][epsil] = 1;
				Habgdef_compute_flags[gamma][beta][alpha][delta][epsil][phi] = 1;
                              
			    }
                                
                                
                            }//phi
                        }//epsil
                    }//derivorder
                }//delta
            }//gamma
        }//beta
    }//alpha

}






