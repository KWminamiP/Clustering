#include "TROOT.h"
#include "TStyle.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TMath.h"
#include "TH2F.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>

//#include "../tree/tree.C"
#include "../tree.C"
#include "../GetSectorID.C"
//#include "GetPhiID.C"
void ClusteringResolution(){

  // Input                                                                                 

  TChain *chain = new TChain("tree");
  /////chain -> Add("./findbin_check.root");//findbin_newselection.root
  //chain -> Add("../checkCluster.root");//findbin_newselection.root
  //chain->Add("../ClusterSppedy.root");
  chain->Add("../NewClusterWeight_ver2.root");
  tree tree(chain);

  // Output

  TCanvas *c1 = new TCanvas("deltatheta_distoribution");
  TH1F *histN = new TH1F("count","count",60,-0.03,0.03);

  TH1F *histN5 = new TH1F("5entry","5entry",60,-0.03,0.03);
  TH1F *histN6 = new TH1F("5entry","5entry",60,-0.03,0.03);
  TH1F *histN7 = new TH1F("5entry","5entry",60,-0.03,0.03);
  TH1F *histN8 = new TH1F("5entry","5entry",60,-0.03,0.03);

  TH1F *histSize = new TH1F("answer","answer",501,0,500);
  TH1F *histF = new TH1F("all","all",501,0,500);
  TH1F *histCluster = new TH1F("Cluster","Cluster",501,0,500);
  
  Double_t MDTr=0.;
  Double_t MDTx=0.;
  Double_t MDTy=0.;
  Double_t MDTz=0.;
  Double_t MDTthe=0.;
  Double_t MDTeta=0.;
  Double_t MDTphi=0.;
  Double_t MDT_pr=0.;
  Double_t MDT_px=0.;
  Double_t MDT_py=0.;
  Double_t MDT_pz=0.;
  Double_t MDT_pthe=0.;
  Double_t MDT_rho=0.;
  Double_t Houghthe=0.;
  Double_t Houghrho=0.;
  Double_t TGCthe=0.;
  Double_t TGCeta=0.;
  Double_t TGCphi=0.;
  Double_t deleta=0.;
  Double_t delphi=0.;
  Double_t delthe=0.;
  Double_t delrho=0.;
  Double_t deltathe=0.;
  Double_t h=0.;
  Double_t R=0.;
  Double_t R2=0.;
  Int_t Match=0;
  Int_t N_Match=0;
  Int_t N_Offline=0;
  Int_t Nall=0;
  Int_t N_Match_all=0;
  Int_t N1=0;
  Int_t N2=0;
  Int_t A=0;
  Int_t Answersize=0;
  Int_t AnswerLabel=0;
  Int_t LabelSize=0;
  Int_t Label=0;
  std::vector<int> test;
  std::vector<int> Size;
  std::vector<int> Max;
  std::vector<double> Hough;
  std::vector<double> deltaR;
  // Loop                                                                                                   

  std::cout<<"Start"<<endl;

  for(int i = 0; i < chain->GetEntries(); i++){//Event_Loop
    tree.GetEntry(i);
    //if(i>1000) break;                                
    //if(tree.EventNumber!=56) continue;
    //std::cout<<"event"<<tree.EventNumber<<std::endl;
    test.clear();
    Hough.clear();
    Size.clear();
    deltaR.clear();
    Max.clear();
    //cout<<"mu_muonTyope "<<tree.mu_muonType->at(0)<<endl;
    
    N_Offline = 0;
    N_Match =0;
    if(tree.museg_x->size()==0) continue;
    for(int j = 0; j < tree.museg_x->size(); j++){//MDT_Loop
      //cout<<"MDT Number = "<<j<<endl;
      //N1+=1;
      
      MDTr = sqrt(tree.museg_x->at(j)*tree.museg_x->at(j)+tree.museg_y->at(j)*tree.museg_y->at(j));
      MDTx = tree.museg_x->at(j);
      MDTy = tree.museg_y->at(j);
      MDTz = tree.museg_z->at(j);
      MDTthe = fabs(TMath::ATan(MDTr/MDTz));
      MDTeta = -1*TMath::Log(TMath::Tan(MDTthe/2));
      MDTphi = TMath::ATan2(MDTy,MDTx);

      
      MDT_pr = sqrt(tree.museg_px->at(j)*tree.museg_px->at(j)+tree.museg_py->at(j)*tree.museg_py->at(j));
      MDT_pz = tree.museg_pz->at(j);
      MDT_pthe = fabs(TMath::ATan(MDT_pr/MDT_pz));
      MDT_rho = MDTr*TMath::Cos(MDT_pthe)-fabs(MDTz)*TMath::Sin(MDT_pthe);

      if(fabs(MDTz)<13500 || fabs(MDTz)>15000) continue;
      if(fabs(MDTeta)<1.05 || fabs(MDTeta)>2.4) continue;
      if(tree.museg_ndof<3) continue;
      //cout<<"MDTr ="<<MDTr<<" : MDTz ="<<MDTz<<endl;
      //N2+=1;
      N_Offline = N_Offline + 1;
      //cout<<"ETA "<<MDTeta<<" :z "<<MDTz<<endl;

      //      cout<<"MDT_the"<<MDT_pthe<<" : MDT_rho = "<<MDT_rho<<endl;//CHECK
	    
      if(tree.Hough_rad->size()==0) continue;
      Match=0;
      //N3+=1;

      //if(tree.Hough_label->size()==32)  std::cout<<"event"<<tree.EventNumber<<std::endl;
      //cout<<"Hough_rad_size"<<endl;
      delthe = 100;
      R2 = 200;

      for(int k = 0; k < tree.Hough_rad->size(); k++){//TGC_Loop                                            
	//histF->Fill(tree.label_size->at(k));
        Houghthe = tree.Hough_rad->at(k);
        Houghrho = tree.Hough_rho->at(k);
	//cout<<"Hough_rad:"<<Houghthe<<" : Houghrho;"<<Houghrho<<endl;
	h = Houghrho/TMath::Cos(Houghthe)+15234.9*TMath::Sin(Houghthe)/TMath::Cos(Houghthe);
	//cout<<"h = "<<h<<endl;
        TGCthe = TMath::ATan(h/15234.9);
        TGCeta =  -1*TMath::Log(TMath::Tan(TGCthe/2));
        TGCphi = tree.Hough_phi->at(k);
	//	cout<<"Hough_phi  "<<TGCphi<<"  MDTphi "<<MDTphi<<endl;	//CHECk
        deleta = MDTeta-TGCeta;
	//cout<<"deltaeta = "<<deleta<<endl;
        delphi = MDTphi-TGCphi;
	deltathe = MDT_pthe-Houghthe;
	//cout<<"MDTthe"<<MDT_pthe<<"   delthe"<<deltathe<<endl;
	delrho = MDT_rho-Houghrho;
	AnswerLabel = tree.Hough_label->at(k);
        R = sqrt((deleta*deleta)+(delphi*delphi));
	Label =  tree.Hough_label->at(k);	
	LabelSize = tree.label_size->at(k);

	Hough.push_back(deltathe);	
	Size.push_back(LabelSize);	
	deltaR.push_back(R);
	Max.push_back(tree.max->at(k));
	//if(R2<0.2 && LabelSize> 100 ) histN->Fill(deltathe);//resolution
	//if(R2<0.2 && LabelSize>100 ) Match=Match+1;//efficiency 
	//cout<<"max"<<tree.max->at(k)<<endl;

	if(R>0.2) continue;//add 2/26
	if( LabelSize<20 && tree.max->at(k) == 5) continue;
	else if( LabelSize<70 && tree.max->at(k) == 6) continue;
	else if( LabelSize<110 && tree.max->at(k) == 7) continue;
	else if( LabelSize<150 && tree.max->at(k) > 7) continue;
	//histN->Fill(deltathe);

	//	cout <<"eventnumber"<<tree.EventNumber<<" : k "<< k <<" : "<< tree.label_size->at(k) <<" bin"<<" : delthe(k=0) = "<< delthe <<" : deltathe ="<< deltathe <<" : R="<<R<<" : the ="<<Houghthe<<" : rho = "<<Houghrho<<"delphi"<<delphi<<endl;//CHECK

        if( k == 0 ){
          delthe = deltathe;
	  R2 = R;
	  AnswerLabel = Label;
	  Answersize = LabelSize;
	  Number = k;
	  // cout<<"k=0 : R2= "<<R2<<" : delthe= "<<delthe<<" : Label= "<<AnswerLabel<<" : size= "<<Answersize<<endl;
        }else if( k > 0 && fabs(delthe) > fabs(deltathe) ){
          delthe = deltathe;
	  R2 = R;
	  AnswerLabel = Label;
	  Answersize = LabelSize;
	  Number = k;
	  //cout<<"k="<<k<<" : R2= "<<R2<<" : delthe= "<<delthe<<" : Label= "<<AnswerLabel<<" : size= "<<Answersize<<endl;
	}

      }//TGC_Loop
      if(R2>0.2) continue;
      test.push_back(Number);
      //cout<<"deltathe = "<<delthe<<endl;

      //if(N_Offline == 1) histN->Fill(delthe);N_Offline==1

      histN->Fill(delthe);//N_offline !=1

      if(fabs(delthe)>0.02){
	//	cout<<"R_"<<R<<endl;
		std::cout<<"eventNumber_RMS ="<<tree.EventNumber<<" :delthe ="<<delthe<<std::endl;
		std::cout<<"    MDTthe ="<<MDT_pthe<<" : MDTrho = "<<MDT_rho<<std::endl;
		std::cout<<"    Houghthe ="<<Houghthe<<" : Houghrho = "<<Houghrho<<std::endl;
	}
      /*
	if( tree.max->at(0) == 5) histN5->Fill(delthe);
	else if( tree.max->at(0) == 6) histN6->Fill(delthe);
	else if( tree.max->at(0) == 7) histN7->Fill(delthe);
	else if( tree.max->at(0) > 7) histN8->Fill(delthe);
      */
      if(Match>0) N_Match=N_Match+1;

      //cout<<"LAST!! R2= "<<R2<<" : delthe= "<<delthe<<" : Label= "<<AnswerLabel<<" : size= "<<Answersize<<" : Number "<<Number<<endl;
    }//MDT_Loop

    //cout<<"Number of MDT segment "<<N_Offline<<endl;
   
    Nall= Nall + N_Offline;
    N_Match_all= N_Match_all + N_Match;
    
    //cout<<"labelSize"<<tree.Hough_label->size()<<" : HoughSize"<<Hough.size()<<endl;
    if(Hough.size()== 0 ) continue;
    for(int l = 0; l < tree.Hough_label->size(); l++){//TGC_Loop
      A=0;
      //cout<<"A: "<<A<<endl;
      for(int m =0; m < test.size(); m++){
	if(test.at(m) == l) A=A+1;
	//cout<<"m: "<<test.at(m)<<" l: "<<l<<endl;
      }
      //cout<<"deltathe = "<<Hough.at(l)<<" : A = "<<A<<" : l= "<<l<<" : Labelsize= "<<tree.label_size->at(l)<<endl;
      //if(N_Offline ==1)histN->Fill(delthe);
      if(A==1 && N_Offline == 1) histSize->Fill(tree.label_size->at(l));
      //if(A==1 && N_Offline == 1) cout<<"Answer Size ="<<tree.label_size->at(l)<<endl;

      //if(N_Offline == 1) histF->Fill(tree.label_size->at(l));      
      // if(tree.label_size->at(l) > 10 && tree.label_size->at(l) < 15) cout<<"EventNumber = "<<tree.EventNumber<<endl;
    }

    
    for(int l =0; l < Hough.size(); l++){
      //if(N_Offline ==1 )cout<<"SIZE:"<<Size.at(l)<<endl;      
      //cout<<"deltaR = "<<deltaR.at(l)<<endl;
      if(N_Offline ==1 ) histF->Fill(Size.at(l));
     
    }

    //if(N_Offline > 1) cout<<"EventNumber"<<tree.EventNumber<<endl;

  }//Event_Loop
  
  //      c1->SetLogy();
      histN->Draw();
  
  /*
    histN5->SetLineColor(3);
    histN5->Draw("same");
    histN6->SetLineColor(4);
    histN6->Draw("same");
    histN7->SetLineColor(2);
    histN7->Draw("same");
    histN8->SetLineColor(5);
    histN8->Draw("same");  
  */
  //gStyle->SetOptStat(0);
  //histCluster->Draw();
  //cout<<"RMSerror"<<histN->GetRMSError()<<endl;
  /*
    TCanvas *c2 = new TCanvas("answer");
    //c2->SetLogy(1);
    histSize->Draw();
    //cout<<"A : "<<histSize->GetBinContent(1)+histSize->GetBinContent(2)+histSize->GetBinContent(3)+histSize->GetBinContent(4)<<endl;
    TCanvas *c3 = new TCanvas("fail");
    c3->SetLogy(1);
    histF->Draw();
  */


  //cout<<"F : "<<histF->Intergral(1,3)<<endl;
  //cout<<"F : "<<histF->GetBinContent(1)+histF->GetBinContent(2)+histF->GetBinContent(3)+histF->GetBinContent(4)<<endl;

  cout<<"N_Match = "<<N_Match_all<<endl;
  cout<<"N_all = "<<Nall<<endl;
  return;
}
