#ifndef __NeuLANDDATA__
#define __NeuLANDDATA__
/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : December 2019                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold NeuLAND Raw data                                          *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// STL
#include <vector>
using namespace std;

// ROOT
#include "TObject.h"

class TNeuLANDData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    // 0 // 
    vector<int>   fNeuLAND_0_ID;
    vector<int>   fNeuLAND_0_Sam;
    vector<int>   fNeuLAND_0_Gtb;
    vector<int>   fNeuLAND_0_Module;
    vector<int>   fNeuLAND_0_Channel;
    vector<int>   fNeuLAND_0_Cycle;
    vector<double>   fNeuLAND_0_ADC;
    vector<double>   fNeuLAND_0_TAC;
    
    // 1 // 
    vector<int>   fNeuLAND_1_ID;
    vector<int>   fNeuLAND_1_Sam;
    vector<int>   fNeuLAND_1_Gtb;
    vector<int>   fNeuLAND_1_Module;
    vector<int>   fNeuLAND_1_Channel;
    vector<int>   fNeuLAND_1_Cycle;
    vector<double>   fNeuLAND_1_ADC;
    vector<double>   fNeuLAND_1_TAC;

    // Veto // 
    // Charge down
    vector<int>   fNeuLAND_QdVETO_ID;
    vector<double>   fNeuLAND_QdVETO_Charge;
    
    // Time down
    vector<int>   fNeuLAND_TdVETO_ID;
    vector<double>   fNeuLAND_TdVETO_Time;
    
    // Charge up
    vector<int>   fNeuLAND_QuVETO_ID;
    vector<double>   fNeuLAND_QuVETO_Charge;
    
    // Time up
    vector<int>   fNeuLAND_TuVETO_ID;
    vector<double>   fNeuLAND_TuVETO_Time;

   //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TNeuLANDData();
    ~TNeuLANDData();
    

  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public:
    void Clear();
    void Clear(const Option_t*) {};
    void Dump() const;


  //////////////////////////////////////////////////////////////
  // Getters and Setters
  // Prefer inline declaration to avoid unnecessary called of 
  // frequently used methods
  // add //! to avoid ROOT creating dictionnary for the methods
  public:
    //////////////////////    SETTERS    ////////////////////////
    // 0 // 
    inline void Set0(const int& ID, const int& Sam, const int& Gtb, const int& Module, const int& Channel, const int& Cycle, const double& ADC, const double& TAC){
    fNeuLAND_0_ID.push_back(ID);
    fNeuLAND_0_Sam.push_back(Sam);
    fNeuLAND_0_Gtb.push_back(Gtb);
    fNeuLAND_0_Module.push_back(Module);
    fNeuLAND_0_Channel.push_back(Channel);
    fNeuLAND_0_Cycle.push_back(Cycle);
    fNeuLAND_0_ADC.push_back(ADC);
    fNeuLAND_0_TAC.push_back(TAC);
    };//!

    // 1 // 
    inline void Set1(const int& ID, const int& Sam, const int& Gtb, const int& Module, const int& Channel, const int& Cycle, const double& ADC, const double& TAC){
    fNeuLAND_1_ID.push_back(ID);
    fNeuLAND_1_Sam.push_back(Sam);
    fNeuLAND_1_Gtb.push_back(Gtb);
    fNeuLAND_1_Module.push_back(Module);
    fNeuLAND_1_Channel.push_back(Channel);
    fNeuLAND_1_Cycle.push_back(Cycle);
    fNeuLAND_1_ADC.push_back(ADC);
    fNeuLAND_1_TAC.push_back(TAC);
    };//!
    // Veto //
    // UP // 
    // Charge
    inline void SetChargeVetoUp(const Double_t& ID, const Double_t& Charge){
    fNeuLAND_QuVETO_ID.push_back(ID);
    fNeuLAND_QuVETO_Charge.push_back(Charge);
    };//!

    // Time
    inline void SetTimeVetoUp(const Double_t& ID, const Double_t& Time){
    fNeuLAND_TuVETO_ID.push_back(ID);
    fNeuLAND_TuVETO_Time.push_back(Time);
    };//!

    // DOWN // 
    // Charge
    inline void SetChargeVetoDown(const Double_t& ID, const Double_t& Charge){
    fNeuLAND_QdVETO_ID.push_back(ID);
    fNeuLAND_QdVETO_Charge.push_back(Charge);
    };//!

    // Time
    inline void SetTimeVetoDown(const Double_t& ID, const Double_t& Time){
    fNeuLAND_TdVETO_ID.push_back(ID);
    fNeuLAND_TdVETO_Time.push_back(Time);
    };//!


    //////////////////////    GETTERS    ////////////////////////
    // MULT //
    // 0 
    inline unsigned int Get0Mult() const
      {return fNeuLAND_0_ID.size();};
    // 1  
    inline unsigned int Get1Mult() const
      {return fNeuLAND_1_ID.size();};
    // Veto
    // ChargeVETO 
    inline unsigned int GetChargeVETOUpMult() const
      {return fNeuLAND_QuVETO_ID.size();};
    // TimeVETO
    inline unsigned int GetTimeVETOUpMult() const
      {return fNeuLAND_TuVETO_ID.size();};
    // ChargeVETO
    inline unsigned int GetChargeVETODownMult() const
      {return fNeuLAND_QdVETO_ID.size();};
    // TimeVETO
    inline unsigned int GetTimeVETODownMult() const
      {return fNeuLAND_TdVETO_ID.size();};


    
    // Value // 
    // 0 //
    inline int Get0ID(unsigned int& i){ return   fNeuLAND_0_ID[i];};
    inline int Get0Sam(unsigned int& i){ return   fNeuLAND_0_Sam[i];};
    inline int Get0Gtb(unsigned int& i){ return   fNeuLAND_0_Gtb[i];};
    inline int Get0Module(unsigned int& i){ return   fNeuLAND_0_Module[i];};
    inline int Get0Channel(unsigned int& i){ return   fNeuLAND_0_Channel[i];};
    inline int Get0Cycle(unsigned int& i){ return   fNeuLAND_0_Cycle[i];};
    inline double Get0ADC(unsigned int& i){ return   fNeuLAND_0_ADC[i];};
    inline double Get0TAC(unsigned int& i){ return   fNeuLAND_0_TAC[i];};

    // 1 //
    inline int Get1ID(unsigned int& i){ return   fNeuLAND_1_ID[i];};
    inline int Get1Sam(unsigned int& i){ return   fNeuLAND_1_Sam[i];};
    inline int Get1Gtb(unsigned int& i){ return   fNeuLAND_1_Gtb[i];};
    inline int Get1Module(unsigned int& i){ return   fNeuLAND_1_Module[i];};
    inline int Get1Channel(unsigned int& i){ return   fNeuLAND_1_Channel[i];};
    inline int Get1Cycle(unsigned int& i){ return   fNeuLAND_1_Cycle[i];};
    inline double Get1ADC(unsigned int& i){ return   fNeuLAND_1_ADC[i];};
    inline double Get1TAC(unsigned int& i){ return   fNeuLAND_1_TAC[i];};
 // MULT //
     // Value // 
    // ChargeVETO 
    inline UShort_t GetChargeVETOUpID(unsigned int& i) const
      {return fNeuLAND_QuVETO_ID[i];};
    inline double GetChargeVETOUp(unsigned int& i) const
      {return fNeuLAND_QuVETO_Charge[i];};
    // TimeVETO 
    inline UShort_t GetTimeVETOUpID(unsigned int& i) const
      {return fNeuLAND_TuVETO_ID[i];};
    inline double GetTimeVETOUp(unsigned int& i) const
      {return fNeuLAND_TuVETO_Time[i];};
    // ChargeVETO 
    inline UShort_t GetChargeVETODownID(unsigned int& i) const
      {return fNeuLAND_QdVETO_ID[i];};
    inline double GetChargeVETODown(unsigned int& i) const
      {return fNeuLAND_QdVETO_Charge[i];};
    // TimeVETO 
    inline UShort_t GetTimeVETODownID(unsigned int& i) const
      {return fNeuLAND_TdVETO_ID[i];};
    inline double GetTimeVETODown(unsigned int& i) const
      {return fNeuLAND_TdVETO_Time[i];};


   
  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TNeuLANDData,1)  // NeuLANDData structure
};

#endif
