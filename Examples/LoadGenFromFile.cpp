 

#include "../source/DaliPaint.h"
using namespace std;

int main()
{
  
  cout << "You can only read in if you already ran 'make gen && ./gen' " << endl;
  
  vector<double> fidin(2);
  fidin[0] = pow(0.5,0.25);
  fidin[1] = pow(0.5,0.25);

  /*Before trying this, make sure the TestGen.cpp has written to the folder Gen_folder. Else, the input files will not be found.*/
  DaliPaint Painter(2, fidin, "Gen_Folder",3);
  
  Painter.DaliFlashBack_Fish("./Gen_Folder/FishGen_run",0);
  Painter.DaliFlashBack_Sabg("./Gen_Folder/SabgGen_run",0);
  Painter.DaliFlashBack_Qabgd("./Gen_Folder/QabgdGen_run",0);
  Painter.DaliFlashBack_Pabgde("./Gen_Folder/PabgdeGen_run",0);
  Painter.DaliFlashBack_Habgdef("./Gen_Folder/HabgdefGen_run",0);
  
  
  Painter.paramnames[0] = "$1$";
  Painter.paramnames[1] = "$2$";

   
  vector<double> sigmas;

  sigmas.push_back(0.08);
  sigmas.push_back(0.08);
  
  vector<double> lowerbounds;
  vector<double> upperbounds;
  lowerbounds.push_back(-3);
  lowerbounds.push_back(-3);
  upperbounds.push_back(3);
  upperbounds.push_back(3);
  
  Painter.SetBounds(lowerbounds,upperbounds);
  
  Painter.gridevaluation("Square_loaded",0.005);
  
  system("python ./Gen_Folder/Square_loadedGrid.py");
 
 
}
