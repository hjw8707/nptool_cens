// NPL
#include "NPCalibrationSource.h"

// Root
#include"TDirectory.h"

////////////////////////////////////////////////////////////////////////////////
NPL::CalibrationSource::CalibrationSource(){
m_SourceSignal = 0 ;
m_FitFunctionType = 1 ;
}

////////////////////////////////////////////////////////////////////////////////
NPL::CalibrationSource::~CalibrationSource(){
}

////////////////////////////////////////////////////////////////////////////////
void NPL::CalibrationSource::ComputeSourceSignal(){


 if (m_FitFunctionType == 0){ //3-alpha source golobal (8 gaussians)
     m_SourceSignal = new TF1("AlphaSourceGauss",AlphaSourceGaus,0,1e9,23);
     m_SourceSignal->SetNpx(20000);

     unsigned int arg_number = 0;
     m_SourceSignal->SetParameter(arg_number++,10);

     for(unsigned int i = 0 ; i < m_Energies.size() ; i++){
         m_SourceSignal->SetParameter(arg_number++, 1000); //main peak amp
         m_SourceSignal->SetParameter(arg_number++, 5000+i*1000); //main peak mean
         for(unsigned int j = 0 ; j < m_Energies[i].size() ; j++){
             m_SourceSignal->SetParameter(arg_number++,m_BranchingRatio[i][j]/m_BranchingRatio[i][0]);
             m_SourceSignal->SetParameter(arg_number++,m_Energies[i][j]/m_Energies[i][0]);
         }
     }

 } else if(m_FitFunctionType==1){// triplet of alpha peaks (3 gaussian with tail)

     //m_SourceSignal = new TF1("AlphaTriplet",AlphaTriplet,0,1e9,9);
     m_SourceSignal = new TF1("AlphaTripletWithTail",AlphaTripletWithTwoTail,0,1e9,12);
     m_SourceSignal->SetNpx(20000);

     unsigned int arg_number = 0;

     m_SourceSignal->SetParameter(arg_number++,1000);
     m_SourceSignal->SetParameter(arg_number++,5000+1*1000);
     m_SourceSignal->SetParameter(arg_number++,20);
     m_SourceSignal->SetParameter(arg_number++,20);
     m_SourceSignal->FixParameter(arg_number++,m_BranchingRatio[1][0]/m_BranchingRatio[1][0]);
     m_SourceSignal->FixParameter(arg_number++,m_BranchingRatio[1][1]/m_BranchingRatio[1][0]);
     m_SourceSignal->SetParameter(arg_number++,m_BranchingRatio[1][2]/m_BranchingRatio[1][0]);
     m_SourceSignal->SetParameter(arg_number++,m_Energies[1][0]/m_Energies[1][0]);
     m_SourceSignal->SetParameter(arg_number++,m_Energies[1][1]/m_Energies[1][0]);
     m_SourceSignal->SetParameter(arg_number++,m_Energies[1][2]/m_Energies[1][0]);
     m_SourceSignal->SetParLimits(9,0.95*m_Energies[1][2]/m_Energies[1][0],1.05*m_Energies[1][2]/m_Energies[1][0]);
     m_SourceSignal->SetParameter(arg_number++,100);
     //m_SourceSignal->SetParLimits(11,10,1000);
     m_SourceSignal->FixParameter(arg_number++,0.4);
     //m_SourceSignal->SetParLimits(11,0.35,0.45);
 }
 else if (m_FitFunctionType == 2){//old way, inline i definition, 8 gaussians.
     TString contribution;
     unsigned int arg_number = 0;

     for(unsigned int i = 0 ; i < m_Energies.size() ; i++){
         for(unsigned int j = 0 ; j < m_Energies[i].size() ; j++){
             contribution +=
                 Form("(%f)*[%i]*exp(-(x-[%i]*%f)*(x-[%i]*%f)/(2*[%i]*[%i]))",
                         m_BranchingRatio[i][j]/m_BranchingRatio[i][0], arg_number,
                         arg_number+1, m_Energies[i][j]/m_Energies[i][0],
                         arg_number+1, m_Energies[i][j]/m_Energies[i][0],
                         arg_number+2,arg_number+2);
             if(j!=m_Energies[i].size()-1) contribution+="+";
         }

         arg_number+=3;
         if(i!=m_Energies.size()-1) contribution+="+";
     }

     m_SourceSignal = new TF1("np_source_signal",contribution,0,1e9);
     m_SourceSignal->SetNpx(20000);
     arg_number = 0;
     for(unsigned int i = 0 ; i < m_Energies.size() ; i++){
         m_SourceSignal->SetParameter(arg_number++,1000);
         m_SourceSignal->SetParameter(arg_number++,5000+i*1000);
         m_SourceSignal->SetParameter(arg_number++,20);
     }
 }

}

////////////////////////////////////////////////////////////////////////////////
double AlphaSourceGaus(double *x,double *p)
{  
        // eight gaussians: (3 for 239Pu, 3 for 241Am, 2 for Cm244)
        // same sigma for all gaussians
        
        //p0 =  sigma of all peaks
        
        //239 Pu
        //p1 =  amplitude of most intense peak
        //p2 =  mean of most intense peak
        //p3 =  relative scaling in amplitude 
        //p4 =  relative scaling in mean for main peak
        //p5 =  relative scaling in amplitude for second peak
        //p6 =  relative scaling in mean for second peak
        //p7 =  relative scaling in amplitude for third peak
        //p8 =  relative scaling in mean for third peak
        //241 Am
        //p9 =  amplitude of most intense peak
        //p10 =  mean of most intense peak
        //p2 =  sigma of most intense peak
        //p11 =  relative scaling in amplitude 
        //p12 =  relative scaling in mean for main peak
        //p13 =  relative scaling in amplitude for second peak
        //p14 =  relative scaling in mean for second peak
        //p15 =  relative scaling in amplitude for third peak
        //p16 =  relative scaling in mean for third peak
        //244 Cm
        //p17 =  amplitude of most intense peak
        //p18 =  mean of most intense peak
        //p2 =  sigma of most intense peak
        //p19 =  relative scaling in amplitude 
        //p20 =  relative scaling in mean for main peak
        //p21 =  relative scaling in amplitude for second peak
        //p22 =  relative scaling in mean for second peak
        double arg = 0;
        double fitval = 0;
        int ngaus = 8;
        int k = 0;
        int rest = 0;
        
        if (p[2]==0) {cout<<"sigma=0 in AlphaTripletGaus"<<endl; return fitval;}

        for(int i=0; i<ngaus; i++){
            rest = i%3;
            if (i!=0 && rest==0) k++;
            arg = (x[0] - p[2+k*ngaus]*p[4+k*ngaus+2*rest])/p[0];
            fitval += p[1+k*ngaus]*p[3+k*ngaus+2*rest]*TMath::Exp(-0.5*arg*arg);
        }

        return fitval;
}
////////////////////////////////////////////////////////////////////////////////
double AlphaTripletGaus(double *x,double *p)
{
        // eight gaussians, same sigma 
        //p0 =  amplitude of most intense peak
        //p1 =  mean of most intense peak
        //p2 =  sigma most intense peak
        //p3 =  relative scaling in amplitude 
        //p4 =  relative scaling in mean for main peak
        //p5 =  relative scaling in amplitude for second peak
        //p6 =  relative scaling in mean for second peak
        //p7 =  relative scaling in amplitude for third peak
        //p8 =  relative scaling in mean for third peak
        double arg = 0;
        if (p[2]!=0) arg = (x[0] - p[1]*p[4])/p[2];
        double fitval = p[0]*p[3]*TMath::Exp(-0.5*arg*arg);
        if (p[2]!=0) arg = (x[0] - p[1]*p[6])/p[2];
              fitval += p[0]*p[5]*TMath::Exp(-0.5*arg*arg);
        if (p[2]!=0) arg = (x[0] - p[1]*p[8])/p[2];
              fitval += p[0]*p[7]*TMath::Exp(-0.5*arg*arg);
        return fitval;
}
////////////////////////////////////////////////////////////////////////////////
double AlphaTripletWithTail(double *x,double *p)
{
        //p0 =  area 
        //p1 =  mean
        //p2 =  sigma
        //p3 =  tail
        //p4 =  amplitude scaling first peak
        //p5 =  amplitude scaling second peak
        //p6 =  amplitude scaling third peak
        //p7 =  energy scaling first peak
        //p8 =  energy scaling second peak
        //p9 =  energy scaling third peak
        double arg1 = 0, arg2 = 0, arg3=0, arg4=0, arg5=0, arg6=0;
        if (p[2]!=0 && p[3]!=0) {
            arg1 = (x[0] - p[1]*p[7])/p[3] + p[2]*p[2]/(2*p[3]*p[3]); 
            arg2 = ( (x[0] - p[1]*p[7])/p[2] + p[2]/p[3] )/TMath::Sqrt(2); 

            arg3 = (x[0] - p[1]*p[8])/p[3] + p[2]*p[2]/(2*p[3]*p[3]); 
            arg4 = ( (x[0] - p[1]*p[8])/p[2] + p[2]/p[3] )/TMath::Sqrt(2); 

            arg5 = (x[0] - p[1]*p[9])/p[3] + p[2]*p[2]/(2*p[3]*p[3]); 
            arg6 = ( (x[0] - p[1]*p[9])/p[2] + p[2]/p[3] )/TMath::Sqrt(2); 
        }
        double fitval = p[0]/(2*p[3])*p[4]*TMath::Exp(arg1)*TMath::Erfc(arg2)
                      + p[0]/(2*p[3])*p[5]*TMath::Exp(arg3)*TMath::Erfc(arg4)
                      + p[0]/(2*p[3])*p[6]*TMath::Exp(arg5)*TMath::Erfc(arg6);
        return fitval;
}
////////////////////////////////////////////////////////////////////////////////
double AlphaTripletWithTwoTail(double *x,double *p)
{
        //p0 =  area 
        //p1 =  mean
        //p2 =  sigma
        //p3 =  tail
        //p4 =  amplitude scaling first peak
        //p5 =  amplitude scaling second peak
        //p6 =  amplitude scaling third peak
        //p7 =  energy scaling first peak
        //p8 =  energy scaling second peak
        //p9 =  energy scaling third peak
        //p10 =  tail2
        //p11 =  eta
        double arg1 = 0, arg2 = 0, arg3=0, arg4=0, arg5=0, arg6=0;
        double arg1b = 0, arg2b = 0, arg3b=0, arg4b=0, arg5b=0, arg6b=0;
        if (p[2]!=0 && p[3]!=0) {
            arg1 = (x[0] - p[1]*p[7])/p[3] + p[2]*p[2]/(2*p[3]*p[3]); 
            arg2 = ( (x[0] - p[1]*p[7])/p[2] + p[2]/p[3] )/TMath::Sqrt(2); 
            arg1b = (x[0] - p[1]*p[7])/p[10] + p[2]*p[2]/(2*p[10]*p[10]); 
            arg2b = ( (x[0] - p[1]*p[7])/p[2] + p[2]/p[10] )/TMath::Sqrt(2); 

            arg3 = (x[0] - p[1]*p[8])/p[3] + p[2]*p[2]/(2*p[3]*p[3]); 
            arg4 = ( (x[0] - p[1]*p[8])/p[2] + p[2]/p[3] )/TMath::Sqrt(2); 
            arg3b = (x[0] - p[1]*p[8])/p[10] + p[2]*p[2]/(2*p[10]*p[10]); 
            arg4b = ( (x[0] - p[1]*p[8])/p[2] + p[2]/p[10] )/TMath::Sqrt(2); 

            arg5 = (x[0] - p[1]*p[9])/p[3] + p[2]*p[2]/(2*p[3]*p[3]); 
            arg6 = ( (x[0] - p[1]*p[9])/p[2] + p[2]/p[3] )/TMath::Sqrt(2); 
            arg5b = (x[0] - p[1]*p[9])/p[10] + p[2]*p[2]/(2*p[10]*p[10]); 
            arg6b = ( (x[0] - p[1]*p[9])/p[2] + p[2]/p[10] )/TMath::Sqrt(2); 
        }
        double fitval = 
            p[0]/2*p[4]*((1-p[11])/p[3]*TMath::Exp(arg1)*TMath::Erfc(arg2)
                         +  p[11]/p[10]*TMath::Exp(arg1b)*TMath::Erfc(arg2b))
          + p[0]/2*p[5]*((1-p[11])/p[3]*TMath::Exp(arg3)*TMath::Erfc(arg4)
                         +  p[11]/p[10]*TMath::Exp(arg3b)*TMath::Erfc(arg4b))
          + p[0]/2*p[6]*((1-p[11])/p[3]*TMath::Exp(arg5)*TMath::Erfc(arg6)
                         +  p[11]/p[10]*TMath::Exp(arg5b)*TMath::Erfc(arg6b));
        return fitval;
}
////////////////////////////////////////////////////////////////////////////////
TF1* NPL::CalibrationSource::GetSourceSignal(){
  if(!m_SourceSignal)
    ComputeSourceSignal();

  return m_SourceSignal;
}

////////////////////////////////////////////////////////////////////////////////
void NPL::CalibrationSource::AddParticleContribution(vector<double> Energies ,  vector<double> ErrEnergies, vector<double> BranchingRatio){
  m_Energies.push_back(Energies);
  m_ErrEnergies.push_back(ErrEnergies);
  m_BranchingRatio.push_back(BranchingRatio);
}

////////////////////////////////////////////////////////////////////////////////
void NPL::CalibrationSource::Set_ThreeAlphaSource(){
 Set_239Pu();
 Set_241Am();
 Set_244Cm();
}

////////////////////////////////////////////////////////////////////////////////
void NPL::CalibrationSource::Set_239Pu(){
  vector<double> Energies;
  vector<double> ErrEnergies;
  vector<double> BranchingRatio;
  
  Energies.push_back(5.15659);
  Energies.push_back(5.1443);
  Energies.push_back(5.1055);

  ErrEnergies.push_back(1.4e-4);
  ErrEnergies.push_back(8e-5);
  ErrEnergies.push_back(8e-4);
  
  BranchingRatio.push_back(70.77);
  BranchingRatio.push_back(17.11);
  BranchingRatio.push_back(11.94);

  AddParticleContribution(Energies,ErrEnergies,BranchingRatio);

}

////////////////////////////////////////////////////////////////////////////////
void NPL::CalibrationSource::Set_241Am(){
  vector<double> Energies;
  vector<double> ErrEnergies;
  vector<double> BranchingRatio;
  
  Energies.push_back(5.48556);
  Energies.push_back(5.44280);
  Energies.push_back(5.388);

  ErrEnergies.push_back(1.2e-4);
  ErrEnergies.push_back(1.3e-4);
  ErrEnergies.push_back(1e-3); // to be checked
  
  BranchingRatio.push_back(84.8);
  BranchingRatio.push_back(13.1);
  BranchingRatio.push_back(1.66);

  AddParticleContribution(Energies,ErrEnergies,BranchingRatio);
}

////////////////////////////////////////////////////////////////////////////////
void NPL::CalibrationSource::Set_244Cm(){
  vector<double> Energies;
  vector<double> ErrEnergies;
  vector<double> BranchingRatio;
  
  Energies.push_back(5.80477);
  Energies.push_back(5.76264);

  ErrEnergies.push_back(5e-5);
  ErrEnergies.push_back(3e-5);
  
  BranchingRatio.push_back(76.40);
  BranchingRatio.push_back(23.60);

  AddParticleContribution(Energies,ErrEnergies,BranchingRatio);
}


