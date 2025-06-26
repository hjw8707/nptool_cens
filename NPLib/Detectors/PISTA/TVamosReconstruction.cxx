#include "TVamosReconstruction.h"

// STL
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


using namespace std;

ClassImp(TVamosReconstruction);

//////////////////////////////////////////////////////////
TVamosReconstruction::TVamosReconstruction(){

}

//////////////////////////////////////////////////////////
TVamosReconstruction::~TVamosReconstruction(){
}

/*
//////////////////////////////////////////////////////////
void TVamosReconstruction::RandomInit()
{  
  Array = Ptr = new float[255];
  for(unsigned int i=0;i<255;i++)
    *(Ptr++) = ((float) i)/254.0 - 0.5;
  Ptr = Array;
}

//////////////////////////////////////////////////////////
float TVamosReconstruction::GetRandom()
{
  float value = *Ptr;
 
  if(Ptr < Array+254)
    Ptr++;
  else
    Ptr = Array;
 
  return value;
}
*/


//////////////////////////////////////////////////////////
void TVamosReconstruction::ReadMatrix(string FName)
{
  int Len=255;
  int i,j;
  char Line[100];
  stringstream *InOut;
  ifstream File;

   ifstream IF;
   IF.open(FName.c_str());
   
   while(IF.getline(Line,Len))
     {
       InOut = new stringstream();
       *InOut << Line;
       *InOut >> i;
       *InOut >> j;
       *InOut >> MRec[i][j][0];
       *InOut >> MRec[i][j][1];
       *InOut >> MRec[i][j][2];
       delete InOut;
     }

   cout << "*** Reading the VAMOS magnetic field Matrix ***" << endl; 
}

 
 
//////////////////////////////////////////////////////////
void TVamosReconstruction::CalculateReconstruction(double Xf, double Tf, double BrhoRef, double&Brho, double&Theta, double&Path)
{
  double Brhot,Thetat,Patht;
  int VecI[2];
    
  Brho = Theta = Path = -500;                  

  Brhot  = 0.;
  Thetat = 0.;
  Patht  = 0.;

  VecI[0] = (int) (Xf+600.);
  VecI[1] = (int) (Tf+200.);
  if(
     (VecI[0] >=0 && VecI[0] < 1100)
     &&
     (VecI[1] >=0 && VecI[1] < 550)
     )
    {
      
      Brhot = (((double) MRec[VecI[0]][VecI[1]][0]) + rand.Uniform(-0.5,0.5))/1000.;
      Thetat = (((double) MRec[VecI[0]][VecI[1]][1]) + rand.Uniform(-0.5,0.5)) - 200.;
      Thetat *=-1;
      Patht = (((double) MRec[VecI[0]][VecI[1]][2]) + rand.Uniform(-0.5,0.5)) /10.;
    }

  if(Brhot >0.001 && Thetat > -300. && Thetat < 300.
                                                && Patht >0 && Patht < 2000.)
    {
      
      Brho = BrhoRef*((float) Brhot);
      Theta = (float) Thetat*-1;
      Path = (float) Patht;
    }
}

