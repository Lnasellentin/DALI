/*An Example that creates the DALI-Approximation to the boxy ring.*/

#include "../source/DaliBase.h"
#include "../source/DaliGen.h"
#include "../source/DaliPaint.h"
#include "../source/Utils.h"
using namespace std;

class TestGen : public DaliGen
{
public:
  TestGen(int theo, int indep, int noof, vector<double> fidin, int derivorder_in, string folder):DaliGen(theo, indep, noof, fidin, derivorder_in, folder)
  {
   /*don't need any funny things inside the constructor apart from the data declaration*/ 
   
   DeclareData();
  }
  
  void DeclareData()
  {
    for(int i = 0; i < noofdat; i++)
    {
      /*it doesn't matter what we declare here because PhysMod
      returns always zero. Just the number of samples matters*/
      datavaribs.at(i).at(0) = 0.0; 
    }
    
    /*tells Dali that the Data have been declared.
     * So the code knows which independent variables
     * must be held fixed when taking derivatives
     * with respect to parameters
     */
  datavaribsdeclared = 1;  
  }
  
  
  double PhysMod(vector <double> pars_and_indep)
  {
   return 0.0; //looking at noise only, mean is always zero    
  }
  
  
  /*This tells DALI how to create the Data covariance matrix. i and j are the indices of the matrix. */
  double DataCovMatrix(vector <double> pars, int i, int j)
  {
    /*makes a diagonal matrix*/
    if(i==j)
    return (pow(pars[0]*pars[0],2)+pow(pars[1]*pars[1],2) ); //this sums up to be the noise^2
    
    else
      return 0.0;
    
  }
  

  
}; //End of the class.



int main()
{
  /*We have two parameters, so the fiducial point needs two entries*/
  vector<double> fiduc;
  fiduc.resize(2);
  
  /*Sets the values of the fiducial: these must be maximum likelihood points. For the boxy ring, the below points to maximize the likelihood.*/
  fiduc[0] = pow(0.5,0.25);
  fiduc[1] = pow(0.5,0.25);

  
  /*Creates an object of the TestGen class. From left to right the arguments specify the object to have two theoretical parameters, one independent variable, five samples shall be drawn, derivatives shall be evaluated at the point specified by fiduc, and 3rd derivatives shall not be taken. Only second. If you change the last 2 to a 3, then third derivatives will be taken.*/
  TestGen Gen(2,1,50, fiduc,3, "Gen_Folder");
  Gen.Parameterdependence(0,1);
  /*This initializes the data. In the class TestGen, this is already done in the constructor, so here this is a repeated execution of DeclareData(). If you don't know what a constructor is, then you will want to execute this function in your main().*/
  Gen.DeclareData();
  
  /*Makes the Code write out the derivatives of the covariance matrix.
   Plotting them is not particularly exciting in the case at hand.
   */
  Gen.ConductRecommendedChecks(1);
  
  /*Calculates the Fisher matrix*/
  Gen.CalcFisher();
 
  
  
  /*Calculates DALI-Tensors*/
  Gen.Calc_BeyondFish(); 
  
  /*Prints the derivatives of the model to a file. In this case, the result will be a file full of zeros, since the model does not return anything but zero.*/
  Gen.PrintFirstDmod();
  Gen.PrintScndDmod();
  Gen.PrintThrdDmod();
  
  
  Gen.DaliSave("Gen_run");
  
  
  /*Now, all Dali Tensors were calculated, and we start to make plots. This is done by the DALI paint class. If you wonder why another class is doing this: Because as the code will grow, DaliPaint will be able to add different DaliTensors together and so on.*/
  
 
  /*Create a DaliPaint object, telling it that it has 2 parameters and the fiducial is fiduc*/
  DaliPaint BM(2, fiduc, "Plots",3);
  

  
  /*Can be omitted: Sets some parameter names in python-latex. 
   *These will be displayed on the axes of the plots.
   If you leave this out, no names will appear, but the code runs fine.
   If you forget the $$, python will complain.*/
  BM.paramnames[0] = "$xxxxxx$";
  BM.paramnames[1] = "$yyyyyy$";
  
  
  /*in multiple dimensions, an MCMC sampler runs through the DALI-approximated likelihood. You may wish to give it upper and lower bounds, in case you have not yet nicely found out which proposal distribution to give to the sampler.*/
  vector<double> lowerbounds;
  vector<double> upperbounds;
  lowerbounds.push_back(-3);
  lowerbounds.push_back(-3);
  upperbounds.push_back(3);
  upperbounds.push_back(3);
  
  BM.SetBounds(lowerbounds,upperbounds);
  
  
  /*pass on the Fisher matrix and the DaliTensors to the DaliPaint object "BM". */
  BM.fisher = Gen.fisher;
  BM.Sabg   = Gen.Sabg;
  BM.Qabgd  = Gen.Qabgd;
  BM.Pabgde = Gen.Pabgde;
  BM.Habgdef = Gen.Habgdef;
  
  /*you might want to see your fishermatrix on screen. 
   *2,2 tells it that it is a square twodimensional matrix*/
  //printmatrix(BM.fisher,2,2);
  
  /*Tell the MCMC sampler where to start, and with what proposal distribution: Startx is where it draws the first sample, which you will want to set to the fiducial evaluation point, typically. This avoids burnin.*/
  vector<double> startx;
  vector<double> sigmas;
  
  startx.push_back(fiduc[0]);
  startx.push_back(fiduc[1]);

  sigmas.push_back(0.08);
  sigmas.push_back(0.08);

  /*Doesn't use the MCMC sampler, but plots the DALI-contours on a grid instead. The MCMC sampler is demonstrated below.*/
 BM.gridevaluation("Square",0.005);
  

 
 /*Tells linux to also execute the python executable that was created by bridevaluation. This makes the plot*/
 
 system("python Plots/SquareGrid.py");
  
 int seeed = 111;
 BM.sampling(1e6,sigmas,startx, seeed);
  
  //This prints the Chain to your screen.
  //You'll only want to do this for cheking tiny little short chains for debugging
 // BM.print_chain();
  
  /*Plotting with the MCMC-chain: 0.01 is the relative binwidth*/
BM.marginDown(0,1,"Square", 0.01);
  
system("python Plots/Square01.py");
  
 
 /*Clean up*/
  Gen.Free();
  BM.Free();
}