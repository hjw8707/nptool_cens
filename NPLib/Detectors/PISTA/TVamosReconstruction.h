#ifndef TVamosReconstruction_h
#define TVamosReconstruction_h

// STL
using namespace std;

// ROOT
#include "TObject.h"
#include "TRandom3.h"

class TVamosReconstruction{
  public:
    TVamosReconstruction();
    ~TVamosReconstruction();

  public:
    void ReadMatrix(string FName);
    void RandomInit();
    void CalculateReconstruction(double,double,double,double&,double&,double&);
  
  private:
    float GetRandom();

  private:
    int MRec[1100][550][3];
    TRandom3 rand;

    ClassDef(TVamosReconstruction,1)
};

#endif

