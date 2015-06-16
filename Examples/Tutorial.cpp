#include "../source/DaliGen.h"
#include "../source/DaliPaint.h"
using namespace std;


/*derives a class "TutorialApplication" from the class DaliGen. The User can do with TutorialApplication whatever is needed in order to customize DALI to the physical problem at hand.*/
class TutorialApplication : public DaliGen
{
public:
  /*public member variables of the class*/
  //none needed
  
  /*public member functions of the class*/  
  TutorialApplication(int theo, int indep, int noof, vector<double> fidin, int order_in, string outdir_name):DaliGen(theo,indep,noof,fidin, order_in, outdir_name)
  {
    
   vector<double> pars;
   pars.resize(theo);
    
   DeclareData();
   PhysMod(pars);
   DataCovMatrix(pars,0,0);
   
   cout << "Once you have conducted these steps, decomment the commented section in the main()." << endl;
   
  };
  
 void DeclareData()
  {
  cout << "Data shall be declared: Either read in from a file, or calculated by some code to be specified here. Aim of this function is to fill the 2-dimensional array 'datavaribs[number-of-data-points][number-of-indep-variables]'." << endl;
  cout << "The entries datavaribs[some data point][i] hold the i-th independent variable of the 'some data point'" << endl;
  cout << "Example: Imagine a data-set of only 3 points, where the number N of fishes was measured, as a function of time t and temperature T. This example is further detailed in the manual. The datavarib[][]-array would then be:" << endl;
  
  cout << "data[0][0] = T_1 " << endl;
  cout << "data[0][1] = t_1 " << endl;
  cout << endl;
  
  cout << "data[1][0] = T_2" << endl;
  cout << "data[1][1] = t_2 " << endl;
  cout << endl;

  cout << "data[2][0] = T_3 " << endl;
  cout << "data[2][1] = t_3 " << endl;
  
  cout << endl;
  cout << "What Dali does NOT need, is the actual measurement at the above specified points. In our case that would be the number N of fishes, and it does not need to be specified in this array. (Just as with the Fisher matrix: also there you don't need to provide the actual result of the measurements.) So don't wonder why there is no second grid that stores" << endl;
  cout << "dataresult[0] = N(T_1, t_1) " << endl;
  cout << "dataresult[1] = N(T_2, t_2) " << endl;
  cout << "dataresult[2] = N(T_3, t_3) " << endl;
   cout << "The code will later calculate these quantitites N(T_3, t_3) itself. That is, what the function PhysMod does." << endl;
    
  cout << "end this function with 'datavaribsdeclared=1; If not, then the code will always complain that you never specified how you measure.'" << endl;
  datavaribsdeclared = 1;
  }
  
  double PhysMod(vector <double> pars_and_indep)
  {
    cout << "If you already have a Fisher-Matrix Code in C/C++, then replace the following return 0 with your physical model function of which your Fisher-Code takes derivatives." << endl;
    cout << "In general, this function shall compute the physical model's prediction as a function of the parameters-to-be-measured, and the independent variables. This function could e.g. be the N(t,T,n_mg, n_K) from the fish-example in the manual." << endl;
    
    return 0;
    
  }
  
  double DataCovMatrix(vector <double> pars, int i, int j)
  {
    cout << "If you already have a Fisher-Matrix Code in C/C++, then replace the following return 0 with the function that returns your inverse covariance matrix therefrom. For a constant inverse covarariance matrix, this function could e.g. simply read in from a data-file." << endl
    << endl;
    
    cout << "In general, this function shall return the i-j-element of the data covariance matrix. Note: NOT the inverse data covariance matrix. The code does the inversion by itself." << endl;
    
 
    return 0;
    
  }
   
};

int main()
{
  cout << endl << endl << "DALI-TUTORIAL" << endl << endl;
  
  cout << "If you get any segmentation faults during the tutorial, don't worry. It will be that you use more or less parameters than the example code here. DALI really depends crucially on always having the right number of parameters or datapoint everywhere, so the segmentation faults really teach you: 'Aha, I forgot to adapt something.' So each time you decomment a new section, and you might forget to resize something, it will first segfault. But the problems are always evident and repaired quickly. Still, if the code's behaviour is cryptic, put it into the debugger gdb. It will tell you at which line exactly the code stopped. (Ignore that gdb usually also displays where inside DALI the segfault happened. You'll only need the first gdb messge that tells you in which line of your main() the segmentation fault occured.)" << endl;
  cout << endl << "Decomment the next section of the main() for the next tutorial step, or run the Testgen.cpp if you want to see an immediately working example: make gen && ./gen" << endl;
  
  
  /** If you get a segmentation fault while trying this tutorial,
   * a likely reason is that you forgot to adjust the vectors fidin, sigmas or paramnames.
   * Please don't write to vector elements, which the vector does not have. **/
  
  
  
  /**Tutorial Step 1: Specify your application**/
  
  /*specify some numbers of your physical problem, that are 
   *needed for constructing an object of your problem
   */
  
/*  
  int howmanyparameters = 2; // these numbers were invented...
  int howmanyindependentvariables = 1;
  int numberofdatapoints = 40;
  
  vector<double> fidin;
  fidin.resize(howmanyparameters);
  
  //caution: you cannot acces vector elements beyond [resize-1]
  fidin[0] = 1.11;
  fidin[1] = 2.718;
  //fidin[2] = .... if you have more than 2 parameters
  
  
  //Here, the DALI-object for your problem is constructed. TA is its name, and it is an object of class //"Tutorialapplication". Above the main() you can see what this object was told to do.
   
  TutorialApplication
  TA(howmanyparameters,howmanyindependentvariables,numberofdatapoints,fidin,2, "Tutorial");

*/

  /**Tutorial Step 2: Executing DALI**/
  /*decomment the code below ( -> and recompile <- ) once you
   *overwrote the functions of the Tutorialapplication with
   *either the examples from Examples.cpp or your own
   *application.*/
/*   
  
   //Tell DALI whether the model and/or the covariance matrix depend on parameters
   // first "1": yes, model depends on parameters (put to zero else)
   // second "1": yes, covariance matrix depends on parameters (put to zero else)
   TA.Parameterdependence(1,1);
   //Tell DALI how the measurements were distributed
   TA.DeclareData(); 
   //Calculates the Fisher matrix
   TA.CalcFisher();  
   //Calculates DALI-Tensors
   TA.Calc_BeyondFish(); 
*/
  
  /**Tutorial Step 3: Pass the DALI Tensors on **/
  /*The likelihood approximation is assembled in an object
   * of the class DaliPaint. Here, we pass on the Dali
   * tensors to that class. */

/*  
  DaliPaint Painter(howmanyparameters, fidin, "Tutorial",2);
  
  //the fisher matrix is a pointer. It can cause segmentation faults.
  //Often, this also produces a GSL-error message
  Painter.fisher  = TA.fisher;
  
  //If you want to plot ellipses, don't give the painter the following tensors.
  Painter.Sabg    = TA.Sabg;
  
  //Attention: Qabgd will change depending on whether you calculate third order derivatives, or only second order.
  // So make sure, to really pass to the painter the Qabgd that was calculated with second derivatives only, if
  // you want to have a second-derivatives-only plot. It is not enough to simply leave out the Pabgde and Habgdef 
  // tensors
  Painter.Qabgd   = TA.Qabgd;
  
  //These tensors must only be passed on, if third derivatives were calculated.
  Painter.Pabgde  = TA.Pabgde;
  Painter.Habgdef = TA.Habgdef;

*/
  
  /**Tutorial Step 4: Customize the Plot creation**/
/*  
  
  //set some parameter names, if you want to. Can be omitted.
  Painter.paramnames[0] = "$\\Omega_m$";
  Painter.paramnames[1] = "$w_0$";
  
  
  //give the MCMC sampler a proposal distribution
  vector<double> sigmas;
  double fact = 0.3;
  
  //The Fisher matrix will approximately have the extent
  //of the likelihood, so if it isn't singular, one can use
  //it as a guess for the proposal distribution
  // Else, one has to come up with guesses
  sigmas.push_back(fact/sqrt(gsl_matrix_get(Painter.fisher,0,0)));
  sigmas.push_back(fact/sqrt(gsl_matrix_get(Painter.fisher,1,1)));

*/
 
  /**Tutorial Step 5: Sample and make plots**/

/*
  //MCMC sampling
  Painter.sampling(2e6,sigmas,fidin, 103275); 
  //marginalizing and giving names for the python-executable
  Painter.marginDown(0,1,"Omwo");
  
  //this might or might not work, depending on your operating system
  system("python Tutorial/Omwo01.py");
  
*/  
 
 /*Cleaning up, to prevent memory leaks*/

/*
 
  TA.Free();
  Painter.Free();
 
*/

}