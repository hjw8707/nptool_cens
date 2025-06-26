#ifndef EXOGAMGEO_H
#define EXOGAMGEO_H
/************************************************************************************************


The first function return the structur "Clover_struc" for the flange number given in argument (from 1 to 16 in EXOGAM)

	
	
@ The convention :										
												
 theta = scattering angle by respect to the beam axis										
 phi = azimutal angle										
 phi convention : 0 degree for flange 12 (top-theta=90 deg) and we go in the clockwise looking in the beam direction				
 
 BEWARE that the output phi is the angle of the PROJECTION on the (0,x,y) plan of 		
 the gamma vector !!!!!!!!!!!!!!! (see ROOT manual Chapter 18)										
 
@ Output are in rad
												   											   *
@ The distance between the clover and the target has to be defined by the user	(see #define)		
												
@ The routine works only with the ROOT Lib
								
@ To print the result : switch COMMENTFLAG to 1	


************************************************************************************************

The second function make the Doppler correction of an input energy according to the scattering
(see function at the end of the file for input)		

*************************************************************************************************


	E. Clment CEA-Saclay/SPhN
	June 2006
	
	E553 update E.Clement GANIL
	Sept 2008						   											   
	E530 update N. de Sereville IPNO (target can now be at any position)
	April 2009
	
	E553 update E.Clement GANIL
	May 2009: *more correct crystal size and can thickness
		  *segment mean angle depends on the interaction depth due to the crystal shape
						   											   
************************************************************************************************																
************************************************************************************************/																
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TChain.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TMath.h"
#include "TRint.h"
#include "TString.h"
#include "TCutG.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TF1.h"
#include "TLatex.h"
#include "TDatime.h"
#include "TSystem.h"



struct Clover_struc {
	
	float Theta_Clover;  	// clover center scattering angle by respect to beam axe
	float Phi_Clover;	// clover center phi angle of the PROJECTION on the (0,x,y) plan of the gamma vector 
	float X_Clover,Y_Clover,Z_Clover;
	float Theta_Crystal[4]; //crystalA=0 crystaB=1 etc ....
	float Phi_Crystal[4];
	float X_Crystal[4],Y_Crystal[4],Z_Crystal[4];
	float Theta_Crystal_Seg[4][4]; //CrystalA Segment1 =0,0 CrystalA Segment2=0,1 etc ...
	float Phi_Crystal_Seg[4][4];
	float X_Crystal_Seg[4][4],Y_Crystal_Seg[4][4],Z_Crystal_Seg[4][4];
	std::map<unsigned int, double> D_CloveFlange_Targ_Map; // key is clove nb, value is distance to target (in mm)
};

#define COMMENTFLAG	0


//#define D_CloveFlange_Targ1 	(140.0 +7.0)  //mm + distance capot-crystal (np18-22-03)
//#define D_CloveFlange_Targ2 	(140.0 +7.0)  //mm
//#define D_CloveFlange_Targ3 	(140.0 +7.0)  //mm
//#define D_CloveFlange_Targ4 	(140.0 +7.0)  //mm
//#define D_CloveFlange_Targ5 	(140.0 +7.0)  //mm
//#define D_CloveFlange_Targ6 	(140.0 +7.0)  //mm
//#define D_CloveFlange_Targ7 	(140.0 +7.0)  //mm
//#define D_CloveFlange_Targ8 	(140.0 +7.0)  //mm
//#define D_CloveFlange_Targ9 	(140.0 +7.0)  //mm
//#define D_CloveFlange_Targ10	(140.0 +7.0)  //mm
//#define D_CloveFlange_Targ11	(140.0 +7.0)  //mm
//#define D_CloveFlange_Targ12	(140.0 +7.0)  //mm
//#define D_CloveFlange_Targ13	(140.0 +7.0)  //mm
//#define D_CloveFlange_Targ14	(140.0 +7.0)  //mm
//#define D_CloveFlange_Targ15	(140.0 +7.0)  //mm
//#define D_CloveFlange_Targ16	(140.0 +7.0)  //mm 7mm distance capot cryostat
//#define D_CloveFlange_Targ17	(140.0 +7.0)  //mm

#define	TARGET_POSITION_X	0.	// mm
#define	TARGET_POSITION_Y	0.	// mm
// #define	TARGET_POSITION_Z	0.	// mm
#define InteractionDepth	20.	// mm

inline struct Clover_struc Ask_For_Angles(int flange, double InterDepth, double Position, double TargetZ){ // Position in mm	

struct Clover_struc Result;
int i,j;
Result.D_CloveFlange_Targ_Map[flange] = Position;
// float Real_position[17] ={D_CloveFlange_Targ1, D_CloveFlange_Targ2, D_CloveFlange_Targ3, D_CloveFlange_Targ4, D_CloveFlange_Targ5, D_CloveFlange_Targ6,
// 				D_CloveFlange_Targ7, D_CloveFlange_Targ8, D_CloveFlange_Targ9, D_CloveFlange_Targ10, D_CloveFlange_Targ11,
// 				D_CloveFlange_Targ12, D_CloveFlange_Targ13, D_CloveFlange_Targ14, D_CloveFlange_Targ15, D_CloveFlange_Targ16,D_CloveFlange_Targ17};


//Elements definition
// TVector3 flange12(D_CloveFlange_Targ12,0,0);  	//clover vector
TVector3 flange12(Result.D_CloveFlange_Targ_Map[flange],0,0);  	//clover vector
TVector3 flange12Crist[4];			//crystals vectors
TVector3 flange12CristSeg[4][4];		//segments vectors

//Initialisation by default
TVector3 v1(1,0,0);
TVector3 v2(0,1,0);
TVector3 v3(0,0,1);

	for(i=0;i<4;i++){
		flange12Crist[i].SetX(1.);	
		flange12Crist[i].SetY(1.);	
		flange12Crist[i].SetZ(1.);
		flange12Crist[i].SetTheta(1.);	
		flange12Crist[i].SetPhi(1.);	
		flange12Crist[i].SetMag(1.);	
		for(j=0;j<4;j++){
			flange12CristSeg[i][j].SetX(1.);
			flange12CristSeg[i][j].SetY(1.);
			flange12CristSeg[i][j].SetZ(1.);
			flange12CristSeg[i][j].SetTheta(1.);
			flange12CristSeg[i][j].SetPhi(1.);
			flange12CristSeg[i][j].SetMag(1.);
		}
	}
/*Initialisation of all the vectors to the flange 12 BUT with the correct distance for the asked clover
 between target and germanium crystal*/
	
	//Clover position
	flange12.SetTheta(90.0*TMath::Pi()/180.0);
	flange12.SetPhi(0.0*TMath::Pi()/180.0);
	flange12.SetMag(Result.D_CloveFlange_Targ_Map[flange]); 


//segment mean angle depends on the interaction depth due to the crystal shape; this creates a linear function with InteractionDepth as input
// see plan np18-22-05 && np18-22-04
TF1 *ShapeC = new TF1("ShapeC","0.132*x+20.54");  // Provient du bizautage des clobvers EXOGAM dans la directoon de la cible
TF1 *ShapeS1 = new TF1("ShapeS1","0.273*x+30.81");  //gd
TF1 *ShapeS2 = new TF1("ShapeS2","0.066*x+10.27"); //pt Ancien clover EXOGAM, n'existe pas dans les nouveaux EXOGAM


float EXOGAM_Crystal_Center; 
float EXOGAM_Segment_Pos1,EXOGAM_Segment_Pos2 ;

	// Original code by E. Clément using fixed interaction depth
	//if(InteractionDepth>=30){EXOGAM_Crystal_Center= 24.5;EXOGAM_Segment_Pos1=39.; EXOGAM_Segment_Pos2=12.25; }
	//else{EXOGAM_Crystal_Center=ShapeC->Eval(InteractionDepth);
	//     EXOGAM_Segment_Pos1=ShapeS1->Eval(InteractionDepth); //gd
	//     EXOGAM_Segment_Pos2=ShapeS2->Eval(InteractionDepth); //pt
	//}
	// Modified code by H. Jacob using input interaction depth 30 mm sortie du bizautage après 30 mm 
	if(InterDepth>=30){EXOGAM_Crystal_Center= 24.5;EXOGAM_Segment_Pos1=39.; EXOGAM_Segment_Pos2=12.25; }
	else{EXOGAM_Crystal_Center=ShapeC->Eval(InterDepth);
	     EXOGAM_Segment_Pos1=ShapeS1->Eval(InterDepth); //gd
	     EXOGAM_Segment_Pos2=ShapeS2->Eval(InterDepth); //pt
	}
	


	//Crystal A
	flange12Crist[0].SetY(flange12.Y()+(EXOGAM_Crystal_Center));
	flange12Crist[0].SetZ(flange12.Z()-(EXOGAM_Crystal_Center));
	flange12Crist[0].SetX(Result.D_CloveFlange_Targ_Map[flange]);
		//segment1
	flange12CristSeg[0][0].SetY(flange12.Y()+(EXOGAM_Segment_Pos1)); //gd
	flange12CristSeg[0][0].SetZ(flange12.Z()-(EXOGAM_Segment_Pos1)); //gd
	flange12CristSeg[0][0].SetX(Result.D_CloveFlange_Targ_Map[flange]);
		//segment2
	flange12CristSeg[0][1].SetY(flange12.Y()+(EXOGAM_Segment_Pos2));
	flange12CristSeg[0][1].SetZ(flange12.Z()-(EXOGAM_Segment_Pos1));
	flange12CristSeg[0][1].SetX(Result.D_CloveFlange_Targ_Map[flange]);
		//segment3
	flange12CristSeg[0][2].SetY(flange12.Y()+(EXOGAM_Segment_Pos2));
	flange12CristSeg[0][2].SetZ(flange12.Z()-(EXOGAM_Segment_Pos2));
	flange12CristSeg[0][2].SetX(Result.D_CloveFlange_Targ_Map[flange]);
		//segment4
	flange12CristSeg[0][3].SetY(flange12.Y()+(EXOGAM_Segment_Pos1));
	flange12CristSeg[0][3].SetZ(flange12.Z()-(EXOGAM_Segment_Pos2));
	flange12CristSeg[0][3].SetX(Result.D_CloveFlange_Targ_Map[flange]);

	//Crystal B
	flange12Crist[1].SetY(flange12.Y()-(EXOGAM_Crystal_Center));
	flange12Crist[1].SetZ(flange12.Z()-(EXOGAM_Crystal_Center));
	flange12Crist[1].SetX(Result.D_CloveFlange_Targ_Map[flange]);
		//segment1
	flange12CristSeg[1][0].SetY(flange12.Y()-(EXOGAM_Segment_Pos1));
	flange12CristSeg[1][0].SetZ(flange12.Z()-(EXOGAM_Segment_Pos1));
	flange12CristSeg[1][0].SetX(Result.D_CloveFlange_Targ_Map[flange]);
		//segment2
	flange12CristSeg[1][1].SetY(flange12.Y()-(EXOGAM_Segment_Pos1));
	flange12CristSeg[1][1].SetZ(flange12.Z()-(EXOGAM_Segment_Pos2));
	flange12CristSeg[1][1].SetX(Result.D_CloveFlange_Targ_Map[flange]);
		//segment3
	flange12CristSeg[1][2].SetY(flange12.Y()-(EXOGAM_Segment_Pos2));
	flange12CristSeg[1][2].SetZ(flange12.Z()-(EXOGAM_Segment_Pos2));
	flange12CristSeg[1][2].SetX(Result.D_CloveFlange_Targ_Map[flange]);
		//segment4
	flange12CristSeg[1][3].SetY(flange12.Y()-(EXOGAM_Segment_Pos2));
	flange12CristSeg[1][3].SetZ(flange12.Z()-(EXOGAM_Segment_Pos1));
	flange12CristSeg[1][3].SetX(Result.D_CloveFlange_Targ_Map[flange]);
	
	//Crystal C
	flange12Crist[2].SetY(flange12.Y()-(EXOGAM_Crystal_Center));
	flange12Crist[2].SetZ(flange12.Z()+(EXOGAM_Crystal_Center));
	flange12Crist[2].SetX(Result.D_CloveFlange_Targ_Map[flange]);
		//segment1
	flange12CristSeg[2][0].SetY(flange12.Y()-(EXOGAM_Segment_Pos1));
	flange12CristSeg[2][0].SetZ(flange12.Z()+(EXOGAM_Segment_Pos1));
	flange12CristSeg[2][0].SetX(Result.D_CloveFlange_Targ_Map[flange]);
		//segment2
	flange12CristSeg[2][1].SetY(flange12.Y()-(EXOGAM_Segment_Pos2));
	flange12CristSeg[2][1].SetZ(flange12.Z()+(EXOGAM_Segment_Pos1));
	flange12CristSeg[2][1].SetX(Result.D_CloveFlange_Targ_Map[flange]);
		//segment3
	flange12CristSeg[2][2].SetY(flange12.Y()-(EXOGAM_Segment_Pos2));
	flange12CristSeg[2][2].SetZ(flange12.Z()+(EXOGAM_Segment_Pos2));
	flange12CristSeg[2][2].SetX(Result.D_CloveFlange_Targ_Map[flange]);
		//segment4
	flange12CristSeg[2][3].SetY(flange12.Y()-(EXOGAM_Segment_Pos1));
	flange12CristSeg[2][3].SetZ(flange12.Z()+(EXOGAM_Segment_Pos2));
	flange12CristSeg[2][3].SetX(Result.D_CloveFlange_Targ_Map[flange]);
	
	//Crystal D
	flange12Crist[3].SetY(flange12.Y()+(EXOGAM_Crystal_Center));
	flange12Crist[3].SetZ(flange12.Z()+(EXOGAM_Crystal_Center));
	flange12Crist[3].SetX(Result.D_CloveFlange_Targ_Map[flange]);
		//segment1
	flange12CristSeg[3][0].SetY(flange12.Y()+(EXOGAM_Segment_Pos1));
	flange12CristSeg[3][0].SetZ(flange12.Z()+(EXOGAM_Segment_Pos1));
	flange12CristSeg[3][0].SetX(Result.D_CloveFlange_Targ_Map[flange]);
		//segment2
	flange12CristSeg[3][1].SetY(flange12.Y()+(EXOGAM_Segment_Pos1));
	flange12CristSeg[3][1].SetZ(flange12.Z()+(EXOGAM_Segment_Pos2));
	flange12CristSeg[3][1].SetX(Result.D_CloveFlange_Targ_Map[flange]);
		//segment3
	flange12CristSeg[3][2].SetY(flange12.Y()+(EXOGAM_Segment_Pos2));
	flange12CristSeg[3][2].SetZ(flange12.Z()+(EXOGAM_Segment_Pos2));
	flange12CristSeg[3][2].SetX(Result.D_CloveFlange_Targ_Map[flange]);
		//segment4
	flange12CristSeg[3][3].SetY(flange12.Y()+(EXOGAM_Segment_Pos2));
	flange12CristSeg[3][3].SetZ(flange12.Z()+(EXOGAM_Segment_Pos1));
	flange12CristSeg[3][3].SetX(Result.D_CloveFlange_Targ_Map[flange]);

	if (flange >=1 && flange <=17){
	
			switch(flange){ //which flange ??
			
			case 1: 
				flange12.RotateY(-45.0*TMath::Pi()/180.0);
				for(i=0;i<4;i++){
					flange12Crist[i].RotateY(-45.0*TMath::Pi()/180.0);
					for(j=0;j<4;j++){
						flange12CristSeg[i][j].RotateY(-45.0*TMath::Pi()/180.0);
					}
				}
				break;
				
				
				
			case 2: 
				flange12.RotateY(90.0*TMath::Pi()/180.0);
				flange12.RotateX(-45.0*TMath::Pi()/180.0);
				for(i=0;i<4;i++){
					flange12Crist[i].RotateY(90.0*TMath::Pi()/180.0);
					flange12Crist[i].RotateX(-45.0*TMath::Pi()/180.0);
					for(j=0;j<4;j++){
						flange12CristSeg[i][j].RotateY(90.0*TMath::Pi()/180.0);
						flange12CristSeg[i][j].RotateX(-45.0*TMath::Pi()/180.0);
					}
				}
				break;
				
				
				
			case 3: 
				flange12.RotateY(90.0*TMath::Pi()/180.0);
				flange12.RotateX(45.0*TMath::Pi()/180.0);
				for(i=0;i<4;i++){
					flange12Crist[i].RotateY(90.0*TMath::Pi()/180.0);
					flange12Crist[i].RotateX(45.0*TMath::Pi()/180.0);
					for(j=0;j<4;j++){
						flange12CristSeg[i][j].RotateY(90.0*TMath::Pi()/180.0);
						flange12CristSeg[i][j].RotateX(45.0*TMath::Pi()/180.0);
					}
				}
				break;
				
				
				
			case 4:
				flange12.RotateY(45.0*TMath::Pi()/180.0);
				for(i=0;i<4;i++){
					flange12Crist[i].RotateY(45.0*TMath::Pi()/180.0);
					for(j=0;j<4;j++){
						flange12CristSeg[i][j].RotateY(45.0*TMath::Pi()/180.0);
					}
				}
				break;
				
				
				
			case 5: 
				flange12.RotateY(135.0*TMath::Pi()/180.0);
				for(i=0;i<4;i++){
					flange12Crist[i].RotateY(135.0*TMath::Pi()/180.0);
					for(j=0;j<4;j++){
						flange12CristSeg[i][j].RotateY(135.0*TMath::Pi()/180.0);
					}
				}
				break;
				
				
			case 6: 
				flange12.RotateY(90.0*TMath::Pi()/180.0);
				flange12.RotateX(90.0*TMath::Pi()/180.0);
				for(i=0;i<4;i++){
					flange12Crist[i].RotateY(90.0*TMath::Pi()/180.0);
					flange12Crist[i].RotateX(90.0*TMath::Pi()/180.0);
					for(j=0;j<4;j++){
						flange12CristSeg[i][j].RotateY(90.0*TMath::Pi()/180.0);
						flange12CristSeg[i][j].RotateX(90.0*TMath::Pi()/180.0);
					}
				}
				break;
				
			case 7:  
				flange12.RotateY(135.0*TMath::Pi()/180.0);
				flange12.RotateX(90.0*TMath::Pi()/180.0);
				for(i=0;i<4;i++){
					flange12Crist[i].RotateY(135.0*TMath::Pi()/180.0);
					flange12Crist[i].RotateX(90.0*TMath::Pi()/180.0);
					for(j=0;j<4;j++){
						flange12CristSeg[i][j].RotateY(135.0*TMath::Pi()/180.0);
						flange12CristSeg[i][j].RotateX(90.0*TMath::Pi()/180.0);
					}
				}
				break;
				
				
			case 8:   
				flange12.RotateY(180.0*TMath::Pi()/180.0);
				for(i=0;i<4;i++){
					flange12Crist[i].RotateY(180.0*TMath::Pi()/180.0);
					for(j=0;j<4;j++){
						flange12CristSeg[i][j].RotateY(180.0*TMath::Pi()/180.0);
					}
				}
				break;
				
				
			case 9:   
				flange12.RotateY(135.0*TMath::Pi()/180.0);
				flange12.RotateX(-90.0*TMath::Pi()/180.0);
				for(i=0;i<4;i++){
					flange12Crist[i].RotateY(135.0*TMath::Pi()/180.0);
					flange12Crist[i].RotateX(-90.0*TMath::Pi()/180.0);
					for(j=0;j<4;j++){
						flange12CristSeg[i][j].RotateY(135.0*TMath::Pi()/180.0);
						flange12CristSeg[i][j].RotateX(-90.0*TMath::Pi()/180.0);
					}
				}
				break;
				
				
			case 10:   
				flange12.RotateY(90.0*TMath::Pi()/180.0);
				flange12.RotateX(-90.0*TMath::Pi()/180.0);
				for(i=0;i<4;i++){
					flange12Crist[i].RotateY(90.0*TMath::Pi()/180.0);
					flange12Crist[i].RotateX(-90.0*TMath::Pi()/180.0);
					for(j=0;j<4;j++){
						flange12CristSeg[i][j].RotateY(90.0*TMath::Pi()/180.0);
						flange12CristSeg[i][j].RotateX(-90.0*TMath::Pi()/180.0);
					}
				}
				break;
				
				
			case 11:    
				flange12.RotateY(45.0*TMath::Pi()/180.0);
				flange12.RotateX(-90.0*TMath::Pi()/180.0);
				for(i=0;i<4;i++){
					flange12Crist[i].RotateY(45.0*TMath::Pi()/180.0);
					flange12Crist[i].RotateX(-90.0*TMath::Pi()/180.0);
					for(j=0;j<4;j++){
						flange12CristSeg[i][j].RotateY(45.0*TMath::Pi()/180.0);
						flange12CristSeg[i][j].RotateX(-90.0*TMath::Pi()/180.0);
					}
				}
				break;
				
				
			case 12: //Clover of initialisation --> not move
				break;
				
				
			case 13:    
				flange12.RotateY(45.0*TMath::Pi()/180.0);
				flange12.RotateX(90.0*TMath::Pi()/180.0);
				for(i=0;i<4;i++){
					flange12Crist[i].RotateY(45.0*TMath::Pi()/180.0);
					flange12Crist[i].RotateX(90.0*TMath::Pi()/180.0);
					for(j=0;j<4;j++){
						flange12CristSeg[i][j].RotateY(45.0*TMath::Pi()/180.0);
						flange12CristSeg[i][j].RotateX(90.0*TMath::Pi()/180.0);
					}
				}
				break;
				
				
			case 14:    
				flange12.RotateY(-135.0*TMath::Pi()/180.0);
				for(i=0;i<4;i++){
					flange12Crist[i].RotateY(-135.0*TMath::Pi()/180.0);
					for(j=0;j<4;j++){
						flange12CristSeg[i][j].RotateY(-135.0*TMath::Pi()/180.0);
					}
				}
				break;
				
				
			case 15:     
				flange12.RotateY(90.0*TMath::Pi()/180.0);
				flange12.RotateX(-135.0*TMath::Pi()/180.0);
				for(i=0;i<4;i++){
					flange12Crist[i].RotateY(90.0*TMath::Pi()/180.0);
					flange12Crist[i].RotateX(-135.0*TMath::Pi()/180.0);
					for(j=0;j<4;j++){
						flange12CristSeg[i][j].RotateY(90.0*TMath::Pi()/180.0);
						flange12CristSeg[i][j].RotateX(-135.0*TMath::Pi()/180.0);
					}
				}
				break;	
				
			case 16:    
				flange12.RotateY(90.0*TMath::Pi()/180.0);
				flange12.RotateX(135.0*TMath::Pi()/180.0);
				for(i=0;i<4;i++){
					flange12Crist[i].RotateY(90.0*TMath::Pi()/180.0);
					flange12Crist[i].RotateX(135.0*TMath::Pi()/180.0);
					for(j=0;j<4;j++){
						flange12CristSeg[i][j].RotateY(90.0*TMath::Pi()/180.0);
						flange12CristSeg[i][j].RotateX(135.0*TMath::Pi()/180.0);
					}
				}
			
				break;
				
			case 17:
				flange12.RotateY(90.0*TMath::Pi()/180.0);
				for(i=0;i<4;i++){
					flange12Crist[i].RotateY(90.0*TMath::Pi()/180.0);
					for(j=0;j<4;j++){
						flange12CristSeg[i][j].RotateY(90.0*TMath::Pi()/180.0);
					}
				}
				break;	
				
			default : 
				break;
			}
		
		// Take into account the case of a target which is not at the center of the EXOGAM array
                // Target position
                TVector3 targetPos(TARGET_POSITION_X, TARGET_POSITION_Y, TargetZ);

		// case of the selected EXOGAM detector
/*		for (Int_t k = 0; k < 3; k++) {
		   cout << "flange12 avant : " << flange12(k) << endl;
		   cout << "targetPos      : " << targetPos(k) << endl;
		}*/
		flange12 -= targetPos;
/*		for (Int_t k = 0; k < 3; k++) {
		   cout << "flange12 apres : " << flange12(k) << endl;
		}*/
		// loop on cristals
		for (Int_t ii = 0; ii < 4; ii++) {
		   flange12Crist[ii] -= targetPos;
		   // loop on segments
		   for (Int_t jj = 0; jj < 4; jj++) {
		      flange12CristSeg[ii][jj] -= targetPos;
		   }
		}
			
		#if COMMENTFLAG	
		printf(" flange %d theta %f   phi Proj %f \n",flange,flange12.Theta()*180.0/(TMath::Pi()),flange12.Phi()*180.0/(TMath::Pi()));
		for(i=0;i<4;i++){
			printf(" flange %d  cristal %d theta %f   phi Proj %f \n",
			flange,i+1,flange12Crist[i].Theta()*180.0/(TMath::Pi()),flange12Crist[i].Phi()*180.0/(TMath::Pi()));
			for(j=0;j<4;j++){
				printf(" flange %d  cristal %d  seg %d theta %f   phi Proj %f \n",
				flange,i+1,j+1,flange12CristSeg[i][j].Theta()*180.0/(TMath::Pi()),flange12CristSeg[i][j].Phi()*180.0/(TMath::Pi()));		
			}
		}
		 printf("#########################\n"); 
		#endif
		
		
		//Output result in rad in the "Clover_struc" structur
		Result.Theta_Clover=flange12.Theta();
		Result.Phi_Clover=flange12.Phi();
		Result.X_Clover=flange12.X();
		Result.Y_Clover=flange12.Y();
		Result.Z_Clover=flange12.Z();
			
		for(i=0;i<4;i++){
			Result.Theta_Crystal[i]=flange12Crist[i].Theta();
			Result.Phi_Crystal[i]=flange12Crist[i].Phi();
			Result.X_Crystal[i]=flange12Crist[i].X();
			Result.Y_Crystal[i]=flange12Crist[i].Y();
			Result.Z_Crystal[i]=flange12Crist[i].Z();
			for(j=0;j<4;j++){
				Result.Theta_Crystal_Seg[i][j]=flange12CristSeg[i][j].Theta();
				Result.Phi_Crystal_Seg[i][j]=flange12CristSeg[i][j].Phi();
				Result.X_Crystal_Seg[i][j]=flange12CristSeg[i][j].X();
				Result.Y_Crystal_Seg[i][j]=flange12CristSeg[i][j].Y();
				Result.Z_Crystal_Seg[i][j]=flange12CristSeg[i][j].Z();	
						
			}
		}
		  
		 

	}
	
	else {printf("Bad flange number, flange %d  doesn't exist in EXOGAM !!!!!! \n",flange);}

	delete ShapeC;	
	delete ShapeS1;	
	delete ShapeS2;	
  
  
  return Result;
}


// Routine of doppler correction
inline float Doppler_Correction(float Theta_Gamma, float Phi_Gamma, float Theta_Part, float Phi_Part, float Beta_Part, float energie_Mes){  //rad, v/c
float energievraie,cosinusPSI;
			  
		  cosinusPSI =TMath::Sin(Theta_Part)*TMath::Cos(Phi_Part)*TMath::Sin(Theta_Gamma)*TMath::Cos(Phi_Gamma)+
		  	      TMath::Sin(Theta_Part)*TMath::Sin(Phi_Part)*TMath::Sin(Theta_Gamma)*TMath::Sin(Phi_Gamma)+
			      TMath::Cos(Theta_Part)*TMath::Cos(Theta_Gamma);
			      
	energievraie = energie_Mes*(1.-Beta_Part*cosinusPSI)/sqrt(1.-Beta_Part*Beta_Part);
	
	return energievraie;
}

#endif