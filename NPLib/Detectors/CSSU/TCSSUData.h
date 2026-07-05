#ifndef __CSSUDATA__
#define __CSSUDATA__

#include <vector>

#include "TObject.h"
using namespace std;

class TCSSUData : public TObject {
 private:
  vector<short> fCSSU_E_DetectorNbr;
  vector<double> fCSSU_Energy;
  vector<double> fCSSU_Time;

  vector<short> fCSSU_PMT_DetectorNbr;
  vector<short> fCSSU_PMT_Number;
  vector<int> fCSSU_PMT_PhotonCount;
  vector<double> fCSSU_PMT_FirstTime;

 public:
  TCSSUData();
  virtual ~TCSSUData();

  void Clear();
  void Clear(const Option_t*) {};
  void Dump() const;

  inline unsigned int GetEnergyMult() const { return fCSSU_Energy.size(); }
  inline int GetEnergyDetectorNbr(const unsigned int& i) const { return fCSSU_E_DetectorNbr[i]; }
  inline double GetEnergy(const unsigned int& i) const { return fCSSU_Energy[i]; }
  inline double GetTime(const unsigned int& i) const { return fCSSU_Time[i]; }

  inline unsigned int GetPMTMult() const { return fCSSU_PMT_PhotonCount.size(); }
  inline int GetPMTDetectorNbr(const unsigned int& i) const { return fCSSU_PMT_DetectorNbr[i]; }
  inline int GetPMTNumber(const unsigned int& i) const { return fCSSU_PMT_Number[i]; }
  inline int GetPMTPhotonCount(const unsigned int& i) const { return fCSSU_PMT_PhotonCount[i]; }
  inline double GetPMTFirstTime(const unsigned int& i) const { return fCSSU_PMT_FirstTime[i]; }

  inline void SetEnergy(const int& det, const double& energy, const double& time) {
    fCSSU_E_DetectorNbr.push_back(det);
    fCSSU_Energy.push_back(energy);
    fCSSU_Time.push_back(time);
  }

  inline void SetPMT(const int& det, const int& pmt, const int& photons, const double& firstTime) {
    fCSSU_PMT_DetectorNbr.push_back(det);
    fCSSU_PMT_Number.push_back(pmt);
    fCSSU_PMT_PhotonCount.push_back(photons);
    fCSSU_PMT_FirstTime.push_back(firstTime);
  }

  ClassDef(TCSSUData, 1)
};

#endif
