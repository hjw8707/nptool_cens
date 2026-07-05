#include <iostream>

#include "NPDetectorFactory.h"
#include "TCSSUData.h"

ClassImp(TCSSUData)

TCSSUData::TCSSUData() {}

TCSSUData::~TCSSUData() {}

void TCSSUData::Clear() {
  fCSSU_E_DetectorNbr.clear();
  fCSSU_Energy.clear();
  fCSSU_Time.clear();
  fCSSU_PMT_DetectorNbr.clear();
  fCSSU_PMT_Number.clear();
  fCSSU_PMT_PhotonCount.clear();
  fCSSU_PMT_FirstTime.clear();
}

void TCSSUData::Dump() const {
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New CSSU Event XXXXXXXXXXXXXXXXX" << endl;

  for (unsigned int i = 0; i < fCSSU_Energy.size(); ++i) {
    cout << "CSSU detector " << fCSSU_E_DetectorNbr[i] << " Energy: " << fCSSU_Energy[i]
         << " Time: " << fCSSU_Time[i] << endl;
  }

  for (unsigned int i = 0; i < fCSSU_PMT_PhotonCount.size(); ++i) {
    cout << "CSSU detector " << fCSSU_PMT_DetectorNbr[i] << " PMT " << fCSSU_PMT_Number[i]
         << " Photons: " << fCSSU_PMT_PhotonCount[i] << " FirstTime: " << fCSSU_PMT_FirstTime[i] << endl;
  }
}

extern "C" {
class proxy_npl_CSSU {
 public:
  proxy_npl_CSSU() { NPL::DetectorFactory::getInstance()->AddToken("CSSU", "CSSU"); }
};

proxy_npl_CSSU p_npl_CSSU;
}
