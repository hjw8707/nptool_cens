#include "TExogamGeo.h"
#include <string>


// Position of flanges assuming target always in (0,0,0). TargetZ to check consistence with HitPosition
// Target X and Y are not important (if we are sure that the beam hits the target), only (X,Y) position of the beam matters.
TExogamGeo::TExogamGeo()
:HitPosition(Vector3D()), Beam(Vector3D()), Recoil(Vector3D()),
CloverNumber(-1), CrystalNumber(-1), OuterNumber(-1), EGamma(-1.)
{
  ReadPhotonCS();
}

double TExogamGeo::p_function(const double alpha, Vector3D& SegmentPos,const Vector3D& TargetPos, Vector3D& UnitVector) {
    Vector3D GammaDir = SegmentPos - TargetPos;
    // GammaDir.print();
    return cos(alpha)*(GammaDir.norm());
}
// Définir la fonction f(theta)

double TExogamGeo::equation(const double alpha, Vector3D& SegmentPos,Vector3D& TargetPos, Vector3D& UnitVector) {
    double p_f = p_function(alpha,SegmentPos,TargetPos,UnitVector);
    // std::cout <<   p_f - (SegmentPos - TargetPos).scalarp(UnitVector) << std::endl;
    return p_f - (SegmentPos - TargetPos).scalarp(UnitVector);
}

// Définir la dérivée approximative de la fonction f(theta) par rapport à theta (différence finie)
double TExogamGeo::equation_derivative(const double alpha, Vector3D& SegmentPos, Vector3D& TargetPos, Vector3D& UnitVector) {
    double h = 1e-6;  // Petite valeur pour l'approximation numérique de la dérivée
    return (equation(alpha+h,SegmentPos,TargetPos, UnitVector) - equation(alpha,SegmentPos,TargetPos, UnitVector)) / h;
}




void TExogamGeo::SetGammaInteractionPoint(const short& CloverNumber_,const short& CrystalNumber_,
    const short& OuterNumber_,const double& EGamma_){
      SetCloverNbr(CloverNumber_);
      SetCrystalNbr(CrystalNumber_);
      SetOuterNbr(OuterNumber_);
      SetEGamma(EGamma_);
      SetInteractionDepth(EGamma_);
}

void TExogamGeo::SetCloverNbr(const short& CloverNumber_){
  CloverNumber = CloverNumber_;
}

void TExogamGeo::SetCrystalNbr(const short& CrystalNumber_){
  CrystalNumber = CrystalNumber_;
}

void TExogamGeo::SetOuterNbr(const short& OuterNumber_){
  OuterNumber = OuterNumber_;
}

void TExogamGeo::SetEGamma(const double& EGamma_){
  EGamma = EGamma_;
}

void TExogamGeo::SetInteractionDepth(const double& EGamma_){
  InteractionDepth = interaction_depth(EGamma_);
}

double TExogamGeo::interaction_depth(double Energy){
  auto b = Map_PhotonCS.lower_bound(Energy);
  auto a = prev(b);
  if(b == Map_PhotonCS.begin()){
    a = b;
    b++;
  }
  else if(b == Map_PhotonCS.end() && b!= Map_PhotonCS.begin()){
    b--;
    a = prev(b);
  }
  double coeff = (Energy - a->first)/(b->first - a->first);

  double PhotonCrossSection = a->second + coeff*(b->second - a->second); // mm2/g

  return 1./(GeDensity*PhotonCrossSection);
}
void TExogamGeo::ReadPhotonCS(){
  // path to photon cross section
  std::string NPPath = getenv("NPTOOL");
  std::string CSFilename = NPPath + "/Inputs/PhotonCrossSection/CoherentGe.xcom";
  // std::string CSFilename = "./CoherentGe.xcom";
  std::string LineBuffer;

  std::ifstream CSFile;
  CSFile.open(CSFilename.c_str());

  if (!CSFile.is_open()) {
    std::cout << " No CS file found "
      << CSFilename << std::endl;
    return;
  }
  while(CSFile.good()){
    double gammaE, CrossSection;
    getline(CSFile, LineBuffer);
    std::istringstream ss(LineBuffer);
    ss >> gammaE >> CrossSection; // E in MeV, converted to keV, CrossSection in cm2/g
    CrossSection *= 100.;
    Map_PhotonCS[gammaE] = CrossSection;
  }
}

Vector3D TExogamGeo::getouterpos(const double alpha, const double InteractionDepth){
  double a,b,c,f;
  double l = std::cos(alpha)*InteractionDepth; // length normal to the crystal
  switch(OuterNumber){
    case 0:
      if(l >=30.){
        a = 0.;
        b = -39.;
        c = 0.;
        f = -39.;
      }
      else{
        a = -0.273;
        b = -30.81;
        c = -0.273;
        f = -30.81;
      }
      break;
    case 1:
      if(l >=30.){
        a = 0.;
        b = -12.25;
        c = 0.;
        f = -39.;
      }
      else{
        a = -0.066;
        b = -10.27;
        c = -0.273;
        f = -30.81;
      }
      break;
    case 2:
      if(l >=30.){
        a = 0.;
        b = -12.25;
        c = 0.;
        f = -12.25;
      }
      else{
        a = -0.066;
        b = -10.27;
        c = -0.066;
        f = -10.27;
      }
      break;
    case 3:
      if(l >=30.){
        a = 0.;
        b = -39.;
        c = 0.;
        f = -12.25;
      }
      else{
        a = -0.273;
        b = -30.81;
        c = -0.066;
        f = -10.27;
      }
      break;
  }
  return Vector3D(a*std::cos(alpha)*InteractionDepth + b, 0., c*std::cos(alpha)*InteractionDepth + f);
}



void TExogamGeo::ResetGamma(){
  GammaDir.clear();
  Beam.clear();
  HitPosition.clear();
  Recoil.clear();
  CloverNumber=-1;
  CrystalNumber=-1;
  OuterNumber=-1;
  EGamma=0.;
  InteractionDepth=0.;
}

Vector3D TExogamGeo::newton_raphson(double initial_guess,
    Vector3D& TargetPos, double tolerance = 1e-6, int max_iterations = 100) {
    double alpha = initial_guess;
    Vector3D unit(0.,1.,0.);
    Vector3D fail ;
    setclover(&unit);
    // std::cout << "new NewtonRaphson" << std::endl;
    for (int i = 0; i < max_iterations; ++i) {
        Vector3D outerpos = getouterpos(alpha,InteractionDepth);
        Vector3D SegmentPos(outerpos.X(), getdistance() + std::cos(alpha) * InteractionDepth, outerpos.Z());
        // Vector3D SegmentPos(outerpos.X(), getdistance() Distances[CloverNumber] + std::cos(alpha) * InteractionDepth, outerpos.Z());
        setcrystal(&SegmentPos);
        setclover(&SegmentPos);


        double f_alpha = equation(alpha,SegmentPos,TargetPos,unit);
        double f_prime_alpha = equation_derivative(alpha,SegmentPos,TargetPos,unit);

        // Mise à jour de theta
        double alpha_new = alpha - f_alpha / f_prime_alpha;

        // Si la différence est inférieure à la tolérance, on considère la solution trouvée
        if (std::fabs(alpha_new - alpha) < tolerance) {
            return SegmentPos - TargetPos;
        }
        // std::cout << alpha << " " << alpha_new << std::endl;
        alpha = alpha_new;
        fail = SegmentPos-TargetPos ;
    }

    // If failed convergence is wrong by more than 100% of initial guess just return initial guess
    if ((std::fabs(alpha - initial_guess) > (initial_guess *2)) ){
      Vector3D outerpos = getouterpos(initial_guess,InteractionDepth);
      Vector3D SegmentPos(outerpos.X(), getdistance() + std::cos(initial_guess) * InteractionDepth, outerpos.Z());
      fail = SegmentPos - TargetPos;
    }
    return fail;  // Retourne la dernière valeur trouvée si la convergence n'est pas atteinte
}

void TExogamStructure::setclover(Vector3D* SegmentPos){
  switch (CloverNumber){
    case 1:
      (*SegmentPos).rotateVector(PI/4,0.);
      break;
    case 2:
      (*SegmentPos).rotateVector(-PI/2,-PI/4);
      break;
    case 3:
      (*SegmentPos).rotateVector(-PI/2,PI/4);
      break;
    case 4:
      (*SegmentPos).rotateVector(-PI/4,0.);
      break;
    case 5:
      (*SegmentPos).rotateVector(-3*PI/4,0.);
      break;
    case 6:
      // (*SegmentPos).rotateVector(-PI/2,-PI/2);
      (*SegmentPos).rotateVector(-PI/2,PI/2);
      break;
    case 7:
      // (*SegmentPos).rotateVector(-3*PI/4,-PI/2);
      (*SegmentPos).rotateVector(-3*PI/4,PI/2);
      break;
    case 8:
      (*SegmentPos).rotateVector(-PI,0.);
      break;
    case 9:
      // (*SegmentPos).rotateVector(-3*PI/4,PI/2);
      (*SegmentPos).rotateVector(-3*PI/4,-PI/2);
      break;
    case 10:
      // (*SegmentPos).rotateVector(-PI/2,PI/2);
      (*SegmentPos).rotateVector(-PI/2,-PI/2);
      break;
    case 11:
      // (*SegmentPos).rotateVector(-PI/4,PI/2);
      (*SegmentPos).rotateVector(-PI/4,-PI/2);
      break;
    case 12: // Reference flange
      break;
    case 13:
      // (*SegmentPos).rotateVector(-PI/4,-PI/2);
      (*SegmentPos).rotateVector(-PI/4,PI/2);
      break;
    case 14:
      (*SegmentPos).rotateVector(3*PI/4,0.);
      break;
    case 15:
      (*SegmentPos).rotateVector(-PI/2,-3*PI/4);
      break;
    case 16:
      (*SegmentPos).rotateVector(-PI/2,3*PI/4);
      break;
    case 17:
      (*SegmentPos).rotateVector(-PI/2,0.); // Beam in
      break;
    case 18:
      (*SegmentPos).rotateVector(PI/2,0.); // Beam out
      break;


  }
}

void TExogamStructure::setcrystal(Vector3D* SegmentPos){
  return (*SegmentPos).rotateAroundY(-CrystalNumber*PI/2); // Crystal A= 0, Crystal B= 1 etc.
}

TExogamStructure::TExogamStructure(std::map<unsigned int, double> Distances_)
:TExogamGeo(),Distances(Distances_){}

TExogamCartesian::TExogamCartesian(std::map<unsigned int, std::map<unsigned int,  std::map<unsigned int, Vector3D>>> Coordinates_)
:TExogamGeo(),Coordinates(Coordinates_){}

Vector3D TExogamCartesian::GetCloverCenter(){
  Vector3D MidAC = 0.5*(GetCrystalCenter(0) + GetCrystalCenter(2));
  Vector3D MidBD = 0.5*(GetCrystalCenter(1) + GetCrystalCenter(3));
  if((MidAC - MidBD).norm() > tolerance){ // tolerance on a measurement error to warn the user
    std::cout << "WARNING: the middle of Clover " << CloverNumber << "is not consistent using the outer measures" << std::endl;
    std::cout << "The norm diff between the 2 vectors is " << (MidAC - MidBD).norm() << " mm, mean of the 2 vecs is assumed to be the good one" << std::endl;
    return 0.5*(MidAC+MidBD);
  }
  else
    return MidAC;
}

Vector3D TExogamCartesian::GetCrystalCenter(const unsigned int Crystal){
  Vector3D Mid02 = 0.5*(Coordinates[CloverNumber][Crystal][0] + Coordinates[CloverNumber][Crystal][2]);
  Vector3D Mid13 = 0.5*(Coordinates[CloverNumber][Crystal][1] + Coordinates[CloverNumber][Crystal][3]);
  if((Mid02 - Mid13).norm() > tolerance){ // tolerance on a measurement error to warn the user
    std::cout << "WARNING: the middle of Clover " << CloverNumber << " Crystal " << Crystal << "is not consistent using the outer measures" << std::endl;
    std::cout << "The norm diff between the 2 vectors is " << (Mid02 - Mid13).norm() << " mm, mean of the 2 vecs is assumed to be the good one" << std::endl;
    return 0.5*(Mid02+Mid13);
  }
  else // In this case, the 2 vectors are assumed to be the same
    return Mid02;

}

void TExogamCartesian::setcrystal(Vector3D* SegmentPos){
  Vector3D RefCrystalA(-1.,0.,-1.);
  RefCrystalA.unit();
  setclover(&RefCrystalA);
  double RotaClover =  RefCrystalA.angle(GetCrystalCenter(0) - GetCloverCenter());
  return (*SegmentPos).rotateAroundY(-CrystalNumber*PI/2 + RotaClover); // Crystal A= 0, Crystal B= 1 etc.
}

void TExogamCartesian::setclover(Vector3D* SegmentPos){
  Vector3D RefRotAroundX(0.,1.,0.);
  Vector3D RefRotAroundY(0.,0.,1.);
  Vector3D CloverCenter = GetCloverCenter();
  
  Vector3D ProjCoverCenterYZ(0.,CloverCenter.Y(),CloverCenter.Z());
  Vector3D ProjCoverCenterXZ(CloverCenter.X(),0.,CloverCenter.Z());
  ProjCoverCenterXZ.unit();
  ProjCoverCenterYZ.unit();
  CloverCenter.unit();


  return (*SegmentPos).rotateVector(ProjCoverCenterYZ.angle(RefRotAroundX),-ProjCoverCenterXZ.angle(RefRotAroundY));
}
double TExogamCartesian::getdistance(){
  return GetCloverCenter().norm();
}
