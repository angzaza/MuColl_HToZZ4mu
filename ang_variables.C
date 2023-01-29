//macro to select good muons and events 
//Pt resolution of selected muons is also evaluated

//   good muons:
//   eta<2.5
//   Pt>5 GeV
//   D0<2 mm
//   Z0 <10 mm


#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "Math/Vector4D.h"
#include "TRotation.h"

void ang_variables(){

// declaration of file paths
	TChain* c_sig = new TChain("c_sig");
  c_sig->Add("ntuple_tae_4muvevt_200k_3tev.root/MyLCTuple");
  
  
  TVector3 boostX,boost_z1,boost_z2,a_1;
  TLorentzRotation l,l_z1,l_z2;
  TLorentzVector mup_z1, mum_z1, mup_z2, mum_z2;
  TLorentzVector mu_final[4];
  TVector3 mup_z1_tr, mum_z1_tr, mup_z2_tr, mum_z2_tr,z2_tr;
  double cosTheta_star, phi, phi_1,cosTheta_1,cosTheta_2;
  TVector3 n_1, n_2, n_sc, n_z;
  n_z.SetXYZ(0,0,1);
  
  
  
// cross section and luminosity
  double xs_sig=9.29e-3;
  //double xs_bkg=0.1375;
  double cm_energy=1.5;
  double lumi=500 ;
  
// VARIABLES DECLARATION
  int n_events_sig=0; // number of generated events 
  int n_events_bkg=0;
  
  double bibw_muon[4]={0};
  
  double w1=0.898715;
  double w2=0.92004;
  double w3=0.919569;
  double w4=0.922237;
  double w5=0.930556;
  double w6=0.956689; 
  
  int n=0;
  int n_rc_muon=0; //number of reco muons for each event
  int counter_eta_mu=0, counter_pt_mu=0, counter_d0_mu=0, counter_z0_mu=0; //counter for good muons
  int count_eta_mu_evt=0, count_pt_mu_evt=0, count_d0_mu_evt=0, count_z0_mu_evt=0, count_mup_evt=0, count_mum_evt=0; //counter for the number of muons in each events passing the selections
  int counter_eta_evt=0, counter_pt_evt=0, counter_d0_evt=0, counter_z0_evt=0, counter_OS_evt=0; //counter for events with at least 4 muons passing the corresponing selections
  int counter_0mu_evt=0, counter_1mu_evt=0, counter_2mu_evt=0, counter_3mu_evt=0, counter_4mu_evt=0, counter_5mu_evt=0, counter_6mu_evt=0, counter_6mup_evt=0; //counter for events with 0,1,2,3,4,5,6,>6 reco muons
  
  
  
  int counter_dR=0, counter_pt10=0, counter_pt20=0;
  int condition_4mu_sel=0;
  int condition_dR=0, condition_pt=0;
  int condition_Z=0, condition_Z1=0, condition_ZZ=0, condition_mass_4mu=0;
  int counter_h=0;
  
  int counter_dR_evt=0, counter_2pt_evt=0, counter_ZZ_evt=0, counter_Z1_evt=0, counter_mass_evt=0;
  
  int za_index=0,zb_index=0,z1_index=0,z2_index=0;
  
  TLorentzVector rcmu,mcmu,z1,z2,higgs;
  TLorentzVector mup[4];
  TLorentzVector mum[4];
  TLorentzVector mup_sel[2];
  TLorentzVector mum_sel[2];
  
  TLorentzVector z1_vec;
  TLorentzVector z2_vec;
  TLorentzVector H_vec;
  
  
  double mass_z=91.1876; 
  double deltaR;
  double delta_m1=0,delta_m2=0,delta_m=0;
  double delta_min=500, delta_min_2=500, pt_max=-1;
  double inv_mass_z1, inv_mass_z2, inv_mass_H, pt_H;
  double inv_mass[4]={0};
  double inv_mass_4mu=0;
  double pt_sum_z2=0;
// VARIABLES FROM BRANCHES
   Int_t n_mcp=0;
   Float_t *mc_px;
  Float_t *mc_py;
  Float_t *mc_pz;
  Float_t *mc_e;
  Float_t *mc_vtx;
  Float_t *mc_vty;
  Float_t *mc_vtz;
  Int_t *mc_pdg;
  Int_t *mc_pa0;
  Int_t *mc_gst;
  Int_t *match_trk;
  
  int num=800000;
  mc_px = (float *) malloc(sizeof(float) * num);
  mc_py = (float *) malloc(sizeof(float) * num);
  mc_pz = (float *) malloc(sizeof(float) * num);
  mc_e = (float *) malloc(sizeof(float) * num);
  mc_vtx = (float *) malloc(sizeof(float) * num);
  mc_vty = (float *) malloc(sizeof(float) * num);
  mc_vtz = (float *) malloc(sizeof(float) * num);
  mc_pdg = (int *) malloc(sizeof(int) * num);
  mc_pa0 = (int *) malloc(sizeof(int) * num);
  mc_gst = (int *) malloc(sizeof(int) * num);
  match_trk = (int *) malloc(sizeof(int) * num);
  
  Float_t trk_d0[30000]={0};
  Float_t trk_curv[30000]={0.};
  Float_t trk_phi[30000]={0};
  Float_t trk_z0[30000]={0};
  Float_t trk_tanl[30000];
  Int_t trk_atIP[30000]={0};
  Int_t tr_fts[30000]={0};
  Int_t trk_n=0;
  
  Int_t n_rec;
  Int_t rc_typ[30000];
  Int_t rc_ntr[30000];
  Int_t rc_ftr[30000];
  Float_t rc_mox[30000];
  Float_t rc_moy[30000];
  Float_t rc_moz[30000];
  Float_t rc_ene[30000];
  Float_t rc_cha[1000];
  
    Int_t r2m_nrel;
  Int_t r2m_t[1000];
  Int_t r2m_f[1000];
  Float_t r2m_w[1000];

  
  /// reco tracks branches/////////////////////////////
  c_sig->SetBranchAddress("trsip", trk_atIP);
  c_sig->SetBranchAddress("tsdze", trk_d0);
  c_sig->SetBranchAddress("tsphi", trk_phi);
  c_sig->SetBranchAddress("tsome", trk_curv);
  c_sig->SetBranchAddress("tszze", trk_z0);
  c_sig->SetBranchAddress("tstnl", trk_tanl);
  c_sig->SetBranchAddress("trfts", tr_fts);
  c_sig->SetBranchAddress("ntrk", &trk_n);
  
  //reco particles branches///////////////
  c_sig->SetBranchAddress("nrec", &n_rec);
  c_sig->SetBranchAddress("rctyp", rc_typ);
  c_sig->SetBranchAddress("rcntr", rc_ntr);
  c_sig->SetBranchAddress("rcmox", rc_mox);
  c_sig->SetBranchAddress("rcmoy", rc_moy);
  c_sig->SetBranchAddress("rcmoz", rc_moz);
  c_sig->SetBranchAddress("rcene", rc_ene);
  c_sig->SetBranchAddress("rcftr", rc_ftr);
  c_sig->SetBranchAddress("rccha", rc_cha);
  
  
   TLorentzVector mc_mu,rc_mu;
  
   double deltaPhi=0,deltaPt=0;
  // int i=0,k=0,n=0;
  
  double mz1,mz2,mH, weight;
  TString filename = "MELA_variables_bkg_200k_3Tev_DLF.root";
  TFile *hfile =0;
  hfile= TFile::Open(filename,"RECREATE");
  TTree *tree = new TTree("T","MELA_variables");
  tree->Branch("mZ1",&mz1,"mz1/D");
  tree->Branch("mZ2",&mz2,"mz2/D");
  tree->Branch("cosTheta_star",&cosTheta_star,"cosTheta_star/D");
  tree->Branch("cosTheta_1",&cosTheta_1,"cosTheta_1/D");
  tree->Branch("cosTheta_2",&cosTheta_2,"cosTheta_2/D");
  tree->Branch("phi",&phi,"phi/D");
  tree->Branch("phi_1",&phi_1,"phi_1/D");
  tree->Branch("mH",&mH,"mH/D");
  tree->Branch("weight",&weight,"weight/D");

   ///////////////////////////////////////////////////////////////////////////////////////////////////////
   
   //////////// HISTOGRAMS //////////////////////////////////////////
  TH1F* invM_Z1=new TH1F("Z1mass_sig","Z1mass_sig",100,10,130);
  TH1F* invM_Z2=new TH1F("Z2mass_sig","Z2mass_sig",100,10,130);
  TH1F* angle_z1=new TH1F("angle_z1","angle_z1",50,0,3.15);
  TH1F* angle_z2=new TH1F("angle_z2","angle_z2",50,0,3.15);
  
  TH1F* theta_star_hist= new TH1F("CosTheta_star_sig","CosTheta_star_sig",20,-1,1);
  TH1F* theta1_hist= new TH1F("CosTheta1_sig","CosTheta1_sig",20,-1,1);
  TH1F* theta2_hist= new TH1F("CosTheta2_sig","CosTheta2_sig",20,-1,1);
  TH1F* phi_hist= new TH1F("phi_sig","phi_sig",20,-3.15,3.15);
  TH1F* phi1_hist= new TH1F("phi1_sig","phi1_sig",20,-3.15,3.15);

  
  /////// CUTS //////////////////////////////////////////////
  // for muons//////
  double eta_cut=2.5;
  double pt_cut=5.;
  double d0_cut=2.;
  double z0_cut=10.;
  
  //for events ////
  double dR_cut=0.02;
  double Pt_1=10.;
  double Pt_2=20.;
  double inv_mass_Z_low=12.;
  double inv_mass_Z_high=120.;
  double inv_mass_Z1_cut=40.;
  double inv_mass_4mu_cut=70.;
  ////////////////////////////////////////////////////////////////////////
  

  
  //cycle on events
// for(unsigned int ientry=0; ientry<5000; ++ientry){
 for(unsigned int ientry=0; ientry<c_sig->GetEntries(); ++ientry){
  	c_sig->GetEntry(ientry);
   	n_rc_muon=0;
   	condition_4mu_sel=0; //this variable is switched to one if the event contains at least 4 muons passing the selection
   	
   	// selection of muons /////////////////////////////////
   	
   	count_eta_mu_evt=0; count_pt_mu_evt=0; count_d0_mu_evt=0; count_z0_mu_evt=0; //counters for the number of muons that pass the selections for each event -> they are set to 0 for each event
   	count_mup_evt=0; count_mum_evt=0;
   	
   	//cycle on reco particles
   	for(int i_rc=0;i_rc<n_rec;i_rc++){
   		//condition on muons with an associated track
   		if(abs(rc_typ[i_rc])==13 && rc_ntr[i_rc]==1){
   			n_rc_muon++;
   			rcmu.SetPxPyPzE(rc_mox[i_rc],rc_moy[i_rc],rc_moz[i_rc],rc_ene[i_rc]);
   			
   			if(abs(rcmu.Eta())<eta_cut){
   				counter_eta_mu++;
   				count_eta_mu_evt++;
   				
   				if(rcmu.Pt()>pt_cut){
   					counter_pt_mu++;
   					count_pt_mu_evt++;
   					
   					if(abs(trk_d0[tr_fts[rc_ftr[i_rc]]])<d0_cut){
   						counter_d0_mu++;
   						count_d0_mu_evt++;
   						
   						if(abs(trk_z0[tr_fts[rc_ftr[i_rc]]])<z0_cut){
   							counter_z0_mu++;
   							count_z0_mu_evt++;
   							
     						
     						// muons TLorentz vectors will be stored in two different arrays according to the charge sign
     						if(rc_cha[i_rc]>0){
     							mup[count_mup_evt].SetPxPyPzE(rc_mox[i_rc],rc_moy[i_rc],rc_moz[i_rc],rc_ene[i_rc]);
     							count_mup_evt++;
     						}
   							else if(rc_cha[i_rc]<0){
   								mum[count_mum_evt].SetPxPyPzE(rc_mox[i_rc],rc_moy[i_rc],rc_moz[i_rc],rc_ene[i_rc]);
   								count_mum_evt++;
   							}
   							
   						} // end condition on z0
   					} // end condition on d0
   			  } //end condition on pt
   	  	} //end condition on eta
   	} //end selection of muons
   } // end cycle on reco particles
   
  
   
   // selection of events //////////////////////////////////////
if(count_eta_mu_evt>=4){
  	counter_eta_evt++;
  	if(count_pt_mu_evt>=4){
  		counter_pt_evt++;
  		if(count_d0_mu_evt>=4){
  			counter_d0_evt++;
  			if(count_z0_mu_evt>=4){
  				counter_z0_evt++;
  				if(count_mup_evt>=2 && count_mum_evt>=2){
   					counter_OS_evt++;	
 			 			condition_4mu_sel=1;
 					}
  			}
  		}
  	}
  }
 
  if(condition_4mu_sel==1){
  
  counter_dR=0; counter_pt10=0; counter_pt20=0;
  condition_dR=0; condition_pt=0; condition_Z=0; condition_Z1=0; condition_mass_4mu=0; condition_ZZ=0;
  
  delta_min_2=500.;
  pt_max=-1.;
  
  //in general we have count_mup_evt positive muons and count_mum_evt negative muons,  with count_mup_evt>=4 and count_mum_evt>=4
  //in order to reconstruct the ZZ candidates 4 muons combinations are constructed with the following cycles in i,j,k,l
  for(int i=0;i<count_mup_evt;i++){
  	for (int j=0;j<count_mum_evt;j++){
    	for(int k=i;k<count_mup_evt;k++){
     		for (int l=0;l<count_mum_evt;l++){	
     			if(k!=i && l!=j){	
     				counter_dR=0; counter_pt10=0; counter_pt20=0;
     				
     				mup_sel[0]=mup[i];
     				mup_sel[1]=mup[k];
     				mum_sel[0]=mum[j];
     				mum_sel[1]=mum[l];	
     				
     				//ZZ candidates are required to have deltaR between each of the four leptons > 0.02
     				for(int a=0;a<2;a++){
     					//deltaR between same sign muons
     					for(int b=a+1;b<2;b++){
     						//deltaR=sqrt((mum_sel[a].Eta()-mum_sel[b].Eta())*(mum_sel[a].Eta()-mum_sel[b].Eta())+(mum_sel[a].Phi()-mum_sel[b].Phi())*(mum_sel[a].Phi()-mum_sel[b].Phi()));
     						deltaR=mum_sel[a].DeltaR(mum_sel[b]);
     						if(deltaR>dR_cut){ //condition on deltaR between negative sign muons
   								counter_dR++;
   							}
   							//deltaR=sqrt((mup_sel[a].Eta()-mup_sel[b].Eta())*(mup_sel[a].Eta()-mup_sel[b].Eta())+(mup_sel[a].Phi()-mup_sel[b].Phi())*(mup_sel[a].Phi()-mup_sel[b].Phi()));
   							deltaR=mup_sel[a].DeltaR(mup_sel[b]);
   							if(deltaR>dR_cut){ //condition on deltaR between positive sign muons
   								counter_dR++;
   							}
     					}
     					//deltaR between opposite sign muons
     					for(int b=0;b<2;b++){
     						//deltaR=sqrt((mum_sel[a].Eta()-mup_sel[b].Eta())*(mum_sel[a].Eta()-mup_sel[b].Eta())+(mum_sel[a].Phi()-mup_sel[b].Phi())*(mum_sel[a].Phi()-mup_sel[b].Phi()));
     						deltaR=mum_sel[a].DeltaR(mup_sel[b]);
   							if(deltaR>dR_cut){ //condition on deltaR between opposite sign muons
   								counter_dR++;
   							}
     					}
     				}
     				
     				if(counter_dR==6){
     					condition_dR++;
     					
     					for(int a=0;a<2;a++){
   							if((mup_sel[a].Pt())>Pt_2){counter_pt20++;} //20
   							if((mup_sel[a].Pt())>Pt_1){counter_pt10++;}
   							if((mum_sel[a].Pt())>Pt_1){counter_pt10++;}
   							if((mum_sel[a].Pt())>Pt_2){counter_pt20++;} //20
   						}
   							
   							//selected events should have at least 2 muons with Pt,i>10 GeV and Pt,j>20 GeV
   							if(counter_pt10>=2 && counter_pt20>=1){
   								condition_pt++;
   								
   								for(int a=0;a<2;a++){
   									for(int b=0;b<2;b++){	
   										inv_mass[b+2*a]=(mum_sel[a]+mup_sel[b]).M(); //invariant mass of Z candidates
   									}
   								}
   								//non-overlapping ZZ invariant mass candidates are inv_mass[0]&inv_mass[3]  and inv_mass[1]&inv_mass[2
   								
   								delta_min=500.;
   								condition_mass_4mu=0;
   								//cycle on non-overlapping ZZ candidates
   								for(n=0;n<2;n++){
   									//Z candidates are required to have 12<invMass<120 GeV
   									if((inv_mass[n]>inv_mass_Z_low && inv_mass[n]<inv_mass_Z_high) && (inv_mass[3-n]>inv_mass_Z_low && inv_mass[3-n]<inv_mass_Z_high)){  //selection ZZ candidates 
   										condition_Z++; //condition_Z is increased by 1 if both Z candidates have mass between 12 and 120 GeV
   										delta_m1=abs(inv_mass[n]-mass_z);
   										delta_m2=abs(inv_mass[3-n]-mass_z);
   										za_index=n;
   										zb_index=3-n;
   										delta_m=delta_m1;
   										if(delta_m2<delta_m1){
   											za_index=3-n;
   											zb_index=n;
   											delta_m=delta_m2;
   										}
   										//Z1 candidates are required to have an invariant mass larger than 40 GeV
   										if(inv_mass[za_index]>inv_mass_Z1_cut){
   											condition_Z1++; //condition_Z1 is increased by 1 if the Z1 candidate has a mass greater than 40 GeV
   											if(delta_m<delta_min){
   												delta_min=delta_m;
   												z1_index=za_index;
   												z2_index=zb_index;
   												pt_sum_z2=mum[(int)(z2_index/2)].Pt()+mup[(int)(z2_index%2)].Pt();
   												//ZZcandidates are required to have 4mu invariant mass > 70 GeV
   												if(((mum_sel[0]+mum_sel[1]+mup_sel[0]+mup_sel[1]).M())>70){
   													condition_mass_4mu=1; //condition_mass_4mu is switched to 1 if the invariant mass of the selected 4 muons is larger than 70 GeV
   												}
   											}
   										}
   									}
   								} // end cycle on the ZZ combinations 
   								
   						
   								if(condition_mass_4mu==1){
   									condition_ZZ++; //condition_ZZ is increased by 1 if there is at least one ZZ candidate passing the selection
   									if(fabs(inv_mass[z1_index]-mass_z)<delta_min_2){
   										delta_min_2=fabs(inv_mass[z1_index]-mass_z);
   										inv_mass_z1=inv_mass[z1_index];
   										inv_mass_z2=inv_mass[z2_index];
   										inv_mass_H=(mum_sel[0]+mum_sel[1]+mup_sel[0]+mup_sel[1]).M();
   										pt_H=(mum_sel[0]+mum_sel[1]+mup_sel[0]+mup_sel[1]).Pt();
   										H_vec=mum_sel[0]+mum_sel[1]+mup_sel[0]+mup_sel[1];
   										z1_vec=mum_sel[(int)(z1_index/2)]+mup_sel[(int)(z1_index%2)];
   										z2_vec=mum_sel[(int)(z2_index/2)]+mup_sel[(int)(z2_index%2)];
   										mup_z1=mup_sel[(int)(z1_index%2)];
    									mum_z1=mum_sel[(int)(z1_index/2)];
    									mup_z2=mup_sel[(int)(z2_index%2)];
    									mum_z2=mum_sel[(int)(z2_index/2)];
									pt_max=pt_sum_z2;
   												
   									}
   									else if(fabs(inv_mass[z1_index]-mass_z)==delta_min_2 && pt_sum_z2>pt_max){
   										delta_min_2=fabs(inv_mass[z1_index]-mass_z);
   										pt_max=pt_sum_z2;
   										inv_mass_z1=inv_mass[z1_index];
   										inv_mass_z2=inv_mass[z2_index];
   										inv_mass_H=(mum_sel[0]+mum_sel[1]+mup_sel[0]+mup_sel[1]).M();
   										pt_H=(mum_sel[0]+mum_sel[1]+mup_sel[0]+mup_sel[1]).Pt();
   										H_vec=mum_sel[0]+mum_sel[1]+mup_sel[0]+mup_sel[1];
   										z1_vec=mum_sel[(int)(z1_index/2)]+mup_sel[(int)(z1_index%2)];
   										z2_vec=mum_sel[(int)(z2_index/2)]+mup_sel[(int)(z2_index%2)];
   										mup_z1=mup_sel[(int)(z1_index%2)];
    									mum_z1=mum_sel[(int)(z1_index/2)];
    									mup_z2=mup_sel[(int)(z2_index%2)];
    									mum_z2=mum_sel[(int)(z2_index/2)];
    							
   									}
   								}
   						
   								
   								
   							}// end condition on pt
     					
     				}//end condition on dR
     				
     				
     				
     				
     			}
     		}
     	}
    }
  }
  
  if(condition_ZZ>=1){
  	//invM_Z1->Fill(inv_mass_z1);
	//if(inv_mass_H>=105 && inv_mass_H<=140){
		mH=inv_mass_H;
  		mz1=(mup_z1+mum_z1).M();
  		mz2=(mup_z2+mum_z2).M();
  		invM_Z1->Fill((mup_z1+mum_z1).M());
   		invM_Z2->Fill((mup_z2+mum_z2).M());
    		angle_z1->Fill(fabs(H_vec.Angle(z1_vec.Vect())));
    		angle_z2->Fill(fabs(H_vec.Angle(z2_vec.Vect())));
    		mu_final[0]=mup_z1;
    		mu_final[1]=mum_z1;
    		mu_final[2]=mup_z2;
    		mu_final[3]=mum_z2;
    //weights
    for(int i_mu=0;i_mu<4;i_mu++){
			//mup_z1
			if(mu_final[i_mu].Pt()>=5 && mu_final[i_mu].Pt()<20){
				bibw_muon[i_mu]=w1;
      }
			if(mu_final[i_mu].Pt()>=20 && mu_final[i_mu].Pt()<30){
        bibw_muon[i_mu]=w2;
      }
      if(mu_final[i_mu].Pt()>=30 && mu_final[i_mu].Pt()<40){
        bibw_muon[i_mu]=w3;
      }
		 if(mu_final[i_mu].Pt()>=40 && mu_final[i_mu].Pt()<60){
        bibw_muon[i_mu]=w4;
     }
	   if(mu_final[i_mu].Pt()>=60 && mu_final[i_mu].Pt()<100){
        bibw_muon[i_mu]=w5;
     }
     if(mu_final[i_mu].Pt()>=100 && mu_final[i_mu].Pt()<200){
        bibw_muon[i_mu]=w6;
     }
    }
    
    weight=bibw_muon[0]*bibw_muon[1]*bibw_muon[2]*bibw_muon[3];
    		
    		
    
   // boost=H_vec.BoostVector();
    //l.Boost(boost);
    
    //cout<<"z1: "<<z1_vec.Px()<<" "<<z1_vec.Py()<<" "<<z1_vec.Pz()<<" "<<z1_vec.E()<<" mass: "<<z1_vec.M()<<endl;
    //cout<<"boost "<<boost.X()<<" "<<boost.Y()<<" "<<boost.Z()<<endl;
    
    // HIGGS REST FRAME ///////////////////
    boostX=-(H_vec.BoostVector());
    
    TLorentzVector z1_X(z1_vec);
    z1_X.Boost(boostX);
    TVector3 z1_X_p3= TVector3(z1_X.X(),z1_X.Y(),z1_X.Z());
    //cout<<"z1_tr_p3: "<<z1_tr_p3.X()<<" "<<z1_tr_p3.Y()<<" "<<z1_tr_p3.Z()<<endl;
    
    TLorentzVector z2_X(z2_vec);
    z2_X.Boost(boostX);
    TVector3 z2_X_p3= TVector3(z2_X.X(),z2_X.Y(),z2_X.Z());
    //cout<<"z2_tr: "<<z2_tr.X()<<" "<<z2_tr.Y()<<" "<<z2_tr.Z()<<endl;
    
    TLorentzVector mupZ1_X(mup_z1);
    mupZ1_X.Boost(boostX);
    TVector3 mupZ1_X_p3= TVector3(mupZ1_X.X(),mupZ1_X.Y(),mupZ1_X.Z());
    
    TLorentzVector mumZ1_X(mum_z1);
    mumZ1_X.Boost(boostX);
    TVector3 mumZ1_X_p3= TVector3(mumZ1_X.X(),mumZ1_X.Y(),mumZ1_X.Z());
    
    TLorentzVector mupZ2_X(mup_z2);
    mupZ2_X.Boost(boostX);
    TVector3 mupZ2_X_p3= TVector3(mupZ2_X.X(),mupZ2_X.Y(),mupZ2_X.Z());
    
    TLorentzVector mumZ2_X(mum_z2);
    mumZ2_X.Boost(boostX);
    TVector3 mumZ2_X_p3= TVector3(mumZ2_X.X(),mumZ2_X.Y(),mumZ2_X.Z());
    
    
    n_1=(mumZ1_X_p3.Cross(mupZ1_X_p3)).Unit();
    //cout<<"n_1 : "<<n_1.X()<<" "<<n_1.Y()<<" "<<n_1.Z()<<endl;
    n_2=(mumZ2_X_p3.Cross(mupZ2_X_p3)).Unit();
    //cout<<"n_2 : "<<n_2.X()<<" "<<n_2.Y()<<" "<<n_2.Z()<<endl;
    n_sc=(n_z.Cross(z1_X_p3)).Unit();
    //cout<<"n_sc : "<<n_sc.X()<<" "<<n_sc.Y()<<" "<<n_sc.Z()<<endl;
    
    cosTheta_star=(z1_X_p3.Unit()).Z();
    phi=((z1_X_p3.Dot(n_1.Cross(n_2)))/fabs((z1_X_p3.Dot(n_1.Cross(n_2)))))*TMath::ACos(-n_1.Dot(n_2));
    phi_1=((z1_X_p3.Dot(n_1.Cross(n_sc)))/fabs(z1_X_p3.Dot(n_1.Cross(n_sc))))*TMath::ACos(n_1.Dot(n_sc));
    //cout<<"phi: "<<phi<<endl;
    
    // Z_i REST FRAME //////////////////////////////////////
    
    boost_z1=-(z1_vec.BoostVector());
    boost_z2=-(z2_vec.BoostVector());
    //leptons in Z parent reference frame
    
    TLorentzVector mupZ1_Z1(mup_z1);
    mupZ1_Z1.Boost(boost_z1);
    TVector3 mupZ1_Z1_p3= TVector3(mupZ1_Z1.X(),mupZ1_Z1.Y(),mupZ1_Z1.Z());
    
    TLorentzVector mumZ1_Z1(mum_z1);
    mumZ1_Z1.Boost(boost_z1);
    TVector3 mumZ1_Z1_p3= TVector3(mumZ1_Z1.X(),mumZ1_Z1.Y(),mumZ1_Z1.Z());
    
    TLorentzVector mupZ2_Z2(mup_z2);
    mupZ2_Z2.Boost(boost_z2);
    TVector3 mupZ2_Z2_p3= TVector3(mupZ2_Z2.X(),mupZ2_Z2.Y(),mupZ2_Z2.Z());
    
    TLorentzVector mumZ2_Z2(mum_z2);
    mumZ2_Z2.Boost(boost_z2);
    TVector3 mumZ2_Z2_p3= TVector3(mumZ2_Z2.X(),mumZ2_Z2.Y(),mumZ2_Z2.Z());
    
    // Z bosons in reciprocal reference frame
    TLorentzVector z1_Z2(z1_vec);
    z1_Z2.Boost(boost_z2);
    TVector3 z1_Z2_p3= TVector3(z1_Z2.X(),z1_Z2.Y(),z1_Z2.Z());
    
    TLorentzVector z2_Z1(z2_vec);
    z2_Z1.Boost(boost_z1);
    TVector3 z2_Z1_p3= TVector3(z2_Z1.X(),z2_Z1.Y(),z2_Z1.Z());
    
    
    n_1=(mumZ1_Z1_p3.Cross(mupZ1_Z1_p3)).Unit();
    n_2=(mumZ2_Z2_p3.Cross(mupZ2_Z2_p3)).Unit();
    n_sc=(n_z.Cross(z1_Z2_p3)).Unit();
    
    
    //cout<<"z1_tr.Z: "<<z1_tr.Z()<<endl;
    
    cosTheta_1=-(z2_Z1_p3.Unit()).Dot(mumZ1_Z1_p3.Unit());
    cosTheta_2=-(z1_Z2_p3.Unit()).Dot(mumZ2_Z2_p3.Unit());
    //cout<<"cosTheta_1: "<<cosTheta_1<<endl;
    
    
    tree->Fill();
    
    //cout<<"cosTheta: "<<cosTheta_star<<endl;
    theta_star_hist->Fill(cosTheta_star);
    phi_hist->Fill(phi);
    phi1_hist->Fill(phi_1);
    theta1_hist->Fill(cosTheta_1);
    theta2_hist->Fill(cosTheta_2);
   			
		if(inv_mass_H>105 && inv_mass_H<140){
			counter_h++;
		}
	  //}
	}
  
  
  if(condition_dR>=1){
  	counter_dR_evt++;
  	if(condition_pt>=1){
  		counter_2pt_evt++;
  		if(condition_Z>=1){
  			counter_ZZ_evt++;
  			if(condition_Z1>=1){
  				counter_Z1_evt++;
  				if(condition_ZZ>=1){
  					counter_mass_evt++;
  				}
  			}
  	}
  }
  
  }
  
  
  	
  }
  
 
  

}//end cycle on events


tree->Print();
tree->Write();

TFile *rootFile = new TFile("angular_variables_bkg_3tev_200k_plots.root","RECREATE");
invM_Z1->Write();
invM_Z2->Write();
theta_star_hist->Write();
phi_hist->Write();
phi1_hist->Write();
theta1_hist->Write();
theta2_hist->Write();




}
