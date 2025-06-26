// NPL
#include"NPDetectorManager.h"
#include"NPOptionManager.h"
#include"NPVAnalysis.h"
#include"NPAnalysisFactory.h"
#include"RootInput.h"
#include"RootHistogramsCalib.h"

// STL
#include<iostream>
#include<stdio.h>
#include<dlfcn.h>
#include<stdlib.h>
#include<unistd.h>
#include<sys/time.h>

// Root
#include"TTreeReader.h"
#include"TKey.h"
#include"TEnv.h" 
#include"TROOT.h"
void ProgressDisplay(clock_t&,clock_t&,unsigned long&, unsigned long&, unsigned long&, double&, unsigned long&, int&, int&);
void ProgressDisplay(struct timeval& begin, struct timeval& end, unsigned long& treated,unsigned long& inter,unsigned long& total,double& mean_rate,unsigned long& displayed, int& current_tree, int& total_tree);
////////////////////////////////////////////////////////////////////////////////
int main(int argc , char** argv){
  // command line parsing
  NPOptionManager* myOptionManager = NPOptionManager::getInstance(argc,argv);
  std::string docalibrationfileName    = myOptionManager->GetDoCalibrationFile();
  std::string OutputfileName      = myOptionManager->GetOutputFile();
  if(!myOptionManager->IsTestCalibration()){
    myOptionManager->SetIsCalibration();
  }
  else{
    std::string TreeName="TestCalibration";
    RootOutput::getInstance(OutputfileName,TreeName);
  }
  myOptionManager->SetReader(true);
  std::string inputfilename = myOptionManager->GetRunToReadFile();

  if (myOptionManager->IsDefault("DoCalibration") && myOptionManager->IsCalibration()) {
    std::cout << "Please use a valid DoCalibration File" << std::endl;
    return 0;
  }

  // get input files from NPOptionManager


  // Instantiate the detector using a file
  NPL::DetectorManager* myDetector = new NPL::DetectorManager();
  if (myOptionManager->IsCalibration()){
    RootHistogramsCalib::getInstance(OutputfileName);
    myDetector->ReadDoCalibrationFile(docalibrationfileName);
  }
  else{
    std::string runToReadfileName = NPOptionManager::getInstance()->GetRunToReadFile();
    RootInput::getInstance(runToReadfileName);
    myDetector->ReadConfigurationFile(myOptionManager->GetDetectorFile());
  }
  myDetector->InitializeRootInput();
  TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
  myDetector->SetTreeReader(inputTreeReader);
  
  
  std::cout << std::endl << "///////// Starting Calibration ///////// "<< std::endl;
  TChain* Chain = RootInput:: getInstance()->GetChain();
  myOptionManager->GetNumberOfEntryToAnalyse();

	unsigned long first_entry = myOptionManager->GetFirstEntryToAnalyse(); // defaults to zero
  unsigned long nentries = Chain->GetEntries();
  if(nentries> myOptionManager->GetNumberOfEntryToAnalyse() && myOptionManager->GetNumberOfEntryToAnalyse()>0)
    nentries = myOptionManager->GetNumberOfEntryToAnalyse() ;
	if(nentries + first_entry > Chain->GetEntries()) {nentries = first_entry+Chain->GetEntries();}
	

  TString ChainName = Chain->GetName();
  std::cout << " Number of Event to be treated : " << nentries << " on chain " << ChainName << std::endl;

  unsigned long inter = 0;
  unsigned long treated = 0;
  double mean_rate =0;
  unsigned long displayed=0;
  clock_t end;
  clock_t begin = clock();
  struct timeval tv_end;
  struct timeval tv_begin; gettimeofday(&tv_begin,NULL);
  unsigned long new_nentries = 0 ;
  int current_tree = 0 ;
  int total_tree = Chain->GetNtrees();
  int entry_max = NPOptionManager::getInstance()->GetNumberOfEntryToAnalyse();

  bool IsCalibration = myOptionManager->IsCalibration();
  if(IsCalibration){
    myDetector->InitializeRootHistogramsCalib(); 
    while(inputTreeReader->Next()){
      myDetector->FillHistogramsCalib();
      current_tree = Chain->GetTreeNumber()+1;
      ProgressDisplay(tv_begin,tv_end,treated,inter,nentries,mean_rate,displayed,current_tree,total_tree);
      if(entry_max >= 0 && treated> entry_max)
        break;
    }
    myDetector->DoCalibration();
    myDetector->WriteHistogramsCalib();
  }
  else if(myOptionManager->IsTestCalibration()){
    myDetector->InitializeRootOutput();
    while(inputTreeReader->Next()){
      myDetector->TestCalibration();
      myDetector->FillOutputTree();
      current_tree = Chain->GetTreeNumber()+1;
      ProgressDisplay(tv_begin,tv_end,treated,inter,nentries,mean_rate,displayed,current_tree,total_tree);
      if(entry_max >= 0 && treated> entry_max)
        break;
    }

  }
  else{
    std::cout << "Error, IsCalibration Option not well configured !\n"; 
  }

#if __cplusplus > 199711L && NPMULTITHREADING
  myDetector->StopThread();
#endif

  current_tree = Chain->GetTreeNumber()+1;
  //ProgressDisplay(begin,end,treated,inter,nentries,mean_rate,displayed,current_tree,total_tree);

  ProgressDisplay(tv_begin,tv_end,treated,inter,nentries,mean_rate,displayed,current_tree,total_tree);

  RootInput::Destroy();
  RootOutput::Destroy();
  RootHistogramsCalib::Destroy();
  return 0;
}
////////////////////////////////////////////////////////////////////////////////
void ProgressDisplay(struct timeval& begin, struct timeval& end, unsigned long& treated,unsigned long& inter,unsigned long& total,double& mean_rate,unsigned long& displayed, int& current_tree, int& total_tree){
  gettimeofday(&end,NULL);
  long double elapsed= end.tv_sec-begin.tv_sec+(end.tv_usec-begin.tv_usec)*1e-6;
  if(elapsed>1||treated>=total ){
    displayed++;
    double event_rate = inter/elapsed;
    mean_rate += (event_rate-mean_rate)/(displayed);
    double percent = 100*treated/total;
    double remain = (total-treated)/mean_rate;

    char* timer;
    double check;
    int minutes = remain/60.;
    int seconds = remain -60*minutes;
    if(remain>60)
      check=asprintf(&timer,"%dmin%ds",minutes,seconds);
    else
      check=asprintf(&timer,"%ds",(int)(remain));
    
    static char star[10];
    if(displayed%2==0 || treated==total)
      sprintf(star,"*******");
    else
      sprintf(star,"-------");

    if(treated!=total)
      printf("\r \033[1;31m %s Progress: \033[1;36m%.1f%% \033[1;31m| Rate: %.1fk evt/s | Remain: %s | Tree: %d/%d %s   \033[0m         ", star,percent,mean_rate/1000.,timer, current_tree,total_tree,star);

    else{
      printf("\r \033[1;32m %s Progress: %.1f%% | Rate: %.1fk evt/s | Remain: %s | Tree: %d/%d %s   \033[0m         ", star,percent,mean_rate/1000.,timer, current_tree, total_tree,star);
    }
    fflush(stdout);
    inter=0;
    gettimeofday(&begin,NULL);
  }

  treated++;
  inter++;
}
////////////////////////////////////////////////////////////////////////////////
void ProgressDisplay(clock_t& begin, clock_t& end, unsigned long& treated,unsigned long& inter,unsigned long& total,double& mean_rate,unsigned long& displayed, int& current_tree, int& total_tree){
  end = clock();
  if((end-begin)>CLOCKS_PER_SEC||treated>=total ){
    displayed++;
    long double elapsed =(long double) (end-begin)/CLOCKS_PER_SEC;
    double event_rate = inter/elapsed;
    mean_rate += (event_rate-mean_rate)/(displayed);
    double percent = 100*treated/total;
    double remain = (total-treated)/mean_rate;

    char* timer;
    double check;
    int minutes = remain/60.;
    int seconds = remain -60*minutes;
    if(remain>60)
      check=asprintf(&timer,"%dmin%ds",minutes,seconds);
    else
      check=asprintf(&timer,"%ds",(int)(remain));
    
    static char star[10];
    if(displayed%2==0 || treated==total)
      sprintf(star,"*******");
    else
      sprintf(star,"-------");

    if(treated!=total)
      printf("\r \033[1;31m %s Progress: \033[1;36m%.1f%% \033[1;31m| Rate: %.1fk evt/s | Remain: %s | Tree: %d/%d %s   \033[0m         ", star,percent,mean_rate/1000.,timer, current_tree,total_tree,star);

    else{
      printf("\r \033[1;32m %s Progress: %.1f%% | Rate: %.1fk evt/s | Remain: %s | Tree: %d/%d %s   \033[0m         ", star,percent,mean_rate/1000.,timer, current_tree, total_tree,star);
    }
    fflush(stdout);
    inter=0;
    begin = clock() ;
  }

  treated++;
  inter++;
}
