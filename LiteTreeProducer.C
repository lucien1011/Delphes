/*
This macro is to pre-process the delphes tree to make lite ntuples from delphes
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "TFile.h"
#include "TTree.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TChain.h"
#include "TSystem.h"
#include <TLorentzVector.h>
#include <iostream>
#else
class ExRootTreeReader;
#endif

//------------------------------------------------------------------------------

using namespace std;

bool debug = false;

float m4lHighCut = 9999999.;
float m4lLowCut = 70.0;
float mZ2HighCut=120.0;
float mZ2LowCut=12.0;
float mZ1HighCut=120.0;
float mZ1LowCut=40.0;
float leadingPtCut = 20.;
float subleadingPtCut = 10.;
float deltaRCut = 0.02;
float isoCutMu = 0.35;
float isoCutEl = 0.35;
float minMllCut=4.0;

const float Zmass = 91.1876;
const float muon_mass = 0.10566;
const float electron_mass = 0.5110e-3;

std::vector<float> *lep_id = 0;
std::vector<float> *lep_pt = 0; 
std::vector<float> *lep_eta = 0; 
std::vector<float> *lep_phi = 0;
std::vector<float> *lep_mass = 0;
std::vector<float> *lep_RelIso = 0;
int lep_Hindex[4];

bool passedFullSelection;
bool passedZXCRSelection;
float mass4l,mass4mu,mass4e,mass2e2mu;
float massZ1;
float massZ2;
float pTL1, etaL1;
float pTL2, etaL2;
float pTL3, etaL3;
float pTL4, etaL4;
float phiL1;
float phiL2;
float phiL3;
float phiL4;
int idL1, idL2, idL3, idL4;
float cosTheta1, cosTheta2, cosThetaStar, Phi, Phi1;
int nZXCRFailedLeptons;

void computeAngles(TLorentzVector p4M11, int Z1_lept1Id,
             TLorentzVector p4M12, int Z1_lept2Id,
             TLorentzVector p4M21, int Z2_lept1Id,
             TLorentzVector p4M22, int Z2_lept2Id,
             float & costhetastar, 
             float & costheta1, 
             float & costheta2, 
             float & Phi, 
             float & Phi1){

    //build Z 4-vectors
    TLorentzVector p4Z1 = p4M11 + p4M12;
    TLorentzVector p4Z2 = p4M21 + p4M22;
    
    // Sort Z1 leptons so that:
    if ( (Z1_lept1Id*Z1_lept2Id<0 && Z1_lept1Id<0) || // for OS pairs: lep1 must be the negative one
         (Z1_lept1Id*Z1_lept2Id>0 && p4M11.Phi()<=p4M12.Phi()) //for SS pairs: use random deterministic convention
         ) {
      swap(p4M11, p4M12);
    }
    
    // Same for Z2 leptons
    if ( (Z2_lept1Id*Z2_lept2Id<0 && Z2_lept1Id<0) ||
         (Z2_lept1Id*Z2_lept2Id>0 && p4M21.Phi()<=p4M22.Phi()) 
         ) {
      swap(p4M21, p4M22);
    }
    
    
    //build H 4-vectors
    TLorentzVector p4H = p4Z1 + p4Z2; 
    
    // -----------------------------------
    
    //// costhetastar
    TVector3 boostX = -(p4H.BoostVector());
    TLorentzVector thep4Z1inXFrame( p4Z1 );
    TLorentzVector thep4Z2inXFrame( p4Z2 );
    thep4Z1inXFrame.Boost( boostX );
    thep4Z2inXFrame.Boost( boostX );
    TVector3 theZ1X_p3 = TVector3( thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z() );
    TVector3 theZ2X_p3 = TVector3( thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z() );    
    costhetastar = theZ1X_p3.CosTheta();
    
    //// --------------------------- costheta1
    TVector3 boostV1 = -(p4Z1.BoostVector());
    TLorentzVector p4M11_BV1( p4M11 );
    TLorentzVector p4M12_BV1( p4M12 );
    TLorentzVector p4M21_BV1( p4M21 );
    TLorentzVector p4M22_BV1( p4M22 );
    p4M11_BV1.Boost( boostV1 );
    p4M12_BV1.Boost( boostV1 );
    p4M21_BV1.Boost( boostV1 );
    p4M22_BV1.Boost( boostV1 );
      
    TLorentzVector p4V2_BV1 = p4M21_BV1 + p4M22_BV1;
    //// costheta1
    costheta1 = -p4V2_BV1.Vect().Dot( p4M11_BV1.Vect() )/p4V2_BV1.Vect().Mag()/p4M11_BV1.Vect().Mag();
    
    //// --------------------------- costheta2
    TVector3 boostV2 = -(p4Z2.BoostVector());
    if (boostV2.Mag()>=1.) {
      cout << "Warning: Mela::computeAngles: Z2 boost with beta=1, scaling down" << endl;
      boostV2*=0.9999;
    }
    TLorentzVector p4M11_BV2( p4M11 );
    TLorentzVector p4M12_BV2( p4M12 );
    TLorentzVector p4M21_BV2( p4M21 );
    TLorentzVector p4M22_BV2( p4M22 );
    p4M11_BV2.Boost( boostV2 );
    p4M12_BV2.Boost( boostV2 );
    p4M21_BV2.Boost( boostV2 );
    p4M22_BV2.Boost( boostV2 );
      
    TLorentzVector p4V1_BV2 = p4M11_BV2 + p4M12_BV2;
    //// costheta2
    costheta2 = -p4V1_BV2.Vect().Dot( p4M21_BV2.Vect() )/p4V1_BV2.Vect().Mag()/p4M21_BV2.Vect().Mag();
      
    //// --------------------------- Phi and Phi1 (old phistar1 - azimuthal production angle)
    //    TVector3 boostX = -(thep4H.BoostVector());
    TLorentzVector p4M11_BX( p4M11 );
    TLorentzVector p4M12_BX( p4M12 );
    TLorentzVector p4M21_BX( p4M21 );
    TLorentzVector p4M22_BX( p4M22 );
      
    p4M11_BX.Boost( boostX );
    p4M12_BX.Boost( boostX );
    p4M21_BX.Boost( boostX );
    p4M22_BX.Boost( boostX );
      
    TVector3 tmp1 = p4M11_BX.Vect().Cross( p4M12_BX.Vect() );
    TVector3 tmp2 = p4M21_BX.Vect().Cross( p4M22_BX.Vect() );    
      
    TVector3 normal1_BX( tmp1.X()/tmp1.Mag(), tmp1.Y()/tmp1.Mag(), tmp1.Z()/tmp1.Mag() ); 
    TVector3 normal2_BX( tmp2.X()/tmp2.Mag(), tmp2.Y()/tmp2.Mag(), tmp2.Z()/tmp2.Mag() ); 
    
    //// Phi
    TLorentzVector p4Z1_BX = p4M11_BX + p4M12_BX;    
    float tmpSgnPhi = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normal2_BX) );
    float sgnPhi = tmpSgnPhi/fabs(tmpSgnPhi);
    Phi = sgnPhi * acos( -1.*normal1_BX.Dot( normal2_BX) );
      
      
    //////////////
      
    TVector3 beamAxis(0,0,1);
    TVector3 tmp3 = (p4M11_BX + p4M12_BX).Vect();
      
    TVector3 p3V1_BX( tmp3.X()/tmp3.Mag(), tmp3.Y()/tmp3.Mag(), tmp3.Z()/tmp3.Mag() );
    TVector3 tmp4 = beamAxis.Cross( p3V1_BX );
    TVector3 normalSC_BX( tmp4.X()/tmp4.Mag(), tmp4.Y()/tmp4.Mag(), tmp4.Z()/tmp4.Mag() );
          
    //// Phi1
    float tmpSgnPhi1 = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normalSC_BX) );
    float sgnPhi1 = tmpSgnPhi1/fabs(tmpSgnPhi1);    
    Phi1 = sgnPhi1 * acos( normal1_BX.Dot( normalSC_BX) );    

}


//------------------------------------------------------------------------------

void AnalyseEvents(ExRootTreeReader* treeReader, TTree* outTree) {

    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchTrack = treeReader->UseBranch("Track");
    TClonesArray *branchMuon = treeReader->UseBranch("Muon");
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
    
    Long64_t allEntries = treeReader->GetEntries();
    
    std::cout << "** Chain contains " << allEntries << " events" << std::endl;
    
    Long64_t entry;
    
    Int_t i, j, k, pdgCode;
    
    Track *track;
    Muon *muon;
    Electron *electron;
    GenParticle *particle, *gentrack;

    for(entry = 0; entry < allEntries; ++entry) {
        treeReader->ReadEntry(entry);
        
        if(entry %1000 == 0) std::cout<<entry<< std::endl;

        unsigned int Nele = branchElectron->GetEntriesFast();
        unsigned int Nmu = branchMuon->GetEntriesFast();
        unsigned int Nlep = Nele + Nmu;

        bool foundHiggsCandidate=false;

        if(Nlep < 4) continue;

        bool properLep_ID = false;
        int Nmm = 0; int Nmp = 0; int Nem = 0; int Nep = 0;
        for(unsigned int i =0; i<Nmu; i++) {
            muon = (Muon*) branchMuon->At(i);

            lep_pt->push_back(muon->PT);
            lep_eta->push_back(muon->Eta);
            lep_phi->push_back(muon->Phi);
            lep_id->push_back(muon->Charge*13);
            lep_mass->push_back(muon_mass);
            lep_RelIso->push_back(muon->SumPt/muon->PT);
            
            if(muon->Charge<0.) Nmm = Nmm+1;
            if(muon->Charge<0.) Nmp = Nmp+1;
        };

        for(unsigned int i =0; i<Nele; i++) {
            electron = (Electron*) branchElectron->At(i); 
            
            lep_pt->push_back(electron->PT);
            lep_eta->push_back(electron->Eta);
            lep_phi->push_back(electron->Phi);
            lep_id->push_back(electron->Charge*11);
            lep_mass->push_back(electron_mass);
            lep_RelIso->push_back(electron->SumPt/electron->PT);

            if(electron->Charge<0.) Nem = Nem+1;
            if(electron->Charge>0.) Nep = Nep+1;
        };

        if(Nmm>=2 && Nmp>=2) properLep_ID = true; //4mu
        if(Nem>=2 && Nep>=2) properLep_ID = true; //4e
        if(Nmm>0 && Nmp>0 && Nem>0 && Nep>0) properLep_ID = true; //2e2mu 

        if (!properLep_ID) continue;

        int n_Zs=0;
        std::vector<int> Z_lepindex1;
        std::vector<int> Z_lepindex2;
        std::vector<float> Z_pt, Z_eta, Z_phi, Z_mass;
 
        for(unsigned int i=0; i<Nlep; i++){
            for(unsigned int j=i+1; j<Nlep; j++){

                if (((*lep_id)[i]+(*lep_id)[j])!=0) continue;

                TLorentzVector li, lj;
                li.SetPtEtaPhiM((*lep_pt)[i],(*lep_eta)[i],(*lep_phi)[i],(*lep_mass)[i]);
                lj.SetPtEtaPhiM((*lep_pt)[j],(*lep_eta)[j],(*lep_phi)[j],(*lep_mass)[j]);
                
                TLorentzVector Z = li+lj;
                
                if (debug) cout<<"this Z mass: "<<Z.M()<<" mZ2Low: "<<mZ2LowCut<<endl;
                
                if (Z.M()>0.0) {
                    n_Zs++;
                    Z_pt.push_back(Z.Pt());
                    Z_eta.push_back(Z.Eta());
                    Z_phi.push_back(Z.Phi());
                    Z_mass.push_back(Z.M());
                    Z_lepindex1.push_back(i);
                    Z_lepindex2.push_back(j);
                    if (debug) cout<<" add Z_lepindex1: "<<i<<" Z_lepindex2: "<<j<<endl;
                }
                
            } // lep i
        } // lep j

        // Consider all ZZ candidates
        TLorentzVector Z1Vec, Z2Vec, HVec;
        double minZ1DeltaM_SR=9999.9; double minZ1DeltaM_CR=99999.9;
        double maxZ2SumPt_SR=0.0; double maxZ2SumPt_CR=0.0;
        bool foundSRCandidate=false;

        passedFullSelection=false;

        std::vector<int> Z_Hindex;
        for (int i=0; i<4; i++) {
            if (i<2) Z_Hindex.push_back(-1);
            lep_Hindex[i]=-1;
        }

        for (int i=0; i<n_Zs; i++) {
            for (int j=i+1; j<n_Zs; j++) {
                
                int i1 = Z_lepindex1[i]; int i2 = Z_lepindex2[i];                            
                int j1 = Z_lepindex1[j]; int j2 = Z_lepindex2[j];                            
                
                if (i1 == j1 || i1 == j2 || i2 == j1 || i2 == j2) continue;
                
                TLorentzVector lep_i1, lep_i2, lep_j1, lep_j2;
                lep_i1.SetPtEtaPhiM((*lep_pt)[i1],(*lep_eta)[i1],(*lep_phi)[i1],(*lep_mass)[i1]);
                lep_i2.SetPtEtaPhiM((*lep_pt)[i2],(*lep_eta)[i2],(*lep_phi)[i2],(*lep_mass)[i2]);
                lep_j1.SetPtEtaPhiM((*lep_pt)[j1],(*lep_eta)[j1],(*lep_phi)[j1],(*lep_mass)[j1]);
                lep_j2.SetPtEtaPhiM((*lep_pt)[j2],(*lep_eta)[j2],(*lep_phi)[j2],(*lep_mass)[j2]);
                
                TLorentzVector Zi, Zj;
                Zi.SetPtEtaPhiM(Z_pt[i],Z_eta[i],Z_phi[i],Z_mass[i]);
                Zj.SetPtEtaPhiM(Z_pt[j],Z_eta[j],Z_phi[j],Z_mass[j]);
                
                if (debug) {cout<<"ZZ candidate Zi->M() "<<Zi.M()<<" Zj->M() "<<Zj.M()<<endl;}
                
                TLorentzVector Z1, Z2;
                int Z1index, Z2index;
                int Z1_lepindex[2] = {0,0};
                int Z2_lepindex[2] = {0,0};
                double Z1DeltaM, Z2SumPt;
                
                if (abs(Zi.M()-Zmass)<abs(Zj.M()-Zmass)) { 
                    Z1index = i; Z2index = j;
                    Z1 = Zi; Z2 = Zj;                 
                    if (lep_i1.Pt()>lep_i2.Pt()) { Z1_lepindex[0] = i1;  Z1_lepindex[1] = i2; }
                    else { Z1_lepindex[0] = i2;  Z1_lepindex[1] = i1; }                
                    if (lep_j1.Pt()>lep_j2.Pt()) { Z2_lepindex[0] = j1;  Z2_lepindex[1] = j2; } 
                    else { Z2_lepindex[0] = j2;  Z2_lepindex[1] = j1; }                
                    Z1DeltaM = abs(Zi.M()-Zmass); 
                    Z2SumPt = lep_j1.Pt()+lep_j2.Pt();
                }
                else { 
                    Z1index = j; Z2index = i;
                    Z1 = Zj; Z2 = Zi; 
                    if (lep_j1.Pt()>lep_j2.Pt()) { Z1_lepindex[0] = j1;  Z1_lepindex[1] = j2; }
                    else { Z1_lepindex[0] = j2;  Z1_lepindex[1] = j1; }
                    if (lep_i1.Pt()>lep_i2.Pt()) { Z2_lepindex[0] = i1;  Z2_lepindex[1] = i2; }
                    else { Z2_lepindex[0] = i2;  Z2_lepindex[1] = i1; }
                    Z1DeltaM = abs(Zj.M()-Zmass); 
                    Z2SumPt = lep_i1.Pt()+lep_i2.Pt();
                }
                
                // Check isolation cut (without FSR ) for Z1 leptons
                if (debug) {cout << "RelIso NoFSR Z1: "<< (*lep_RelIso)[Z1_lepindex[0]] << " " << (*lep_RelIso)[Z1_lepindex[1]] << endl;}
                if ((*lep_RelIso)[Z1_lepindex[0]]>((abs((*lep_id)[Z1_lepindex[0]])==11) ? isoCutEl : isoCutMu)) continue;
                if ((*lep_RelIso)[Z1_lepindex[1]]>((abs((*lep_id)[Z1_lepindex[1]])==11) ? isoCutEl : isoCutMu)) continue;
                
                // Check Leading and Subleading pt Cut
                std::vector<double> allPt;
                allPt.push_back(lep_i1.Pt()); allPt.push_back(lep_i2.Pt());
                allPt.push_back(lep_j1.Pt()); allPt.push_back(lep_j2.Pt());
                std::sort(allPt.begin(), allPt.end());
                if (debug) cout<<" leading pt: "<<allPt[3]<<" cut: "<<leadingPtCut
                               <<" subleadingPt: "<<allPt[2]<<" cut: "<<subleadingPtCut<<endl;
                if (allPt[3]<leadingPtCut || allPt[2]<subleadingPtCut ) continue;
                
                // Check dR(li,lj)>0.02 for any i,j
                vector<double> alldR;
                alldR.push_back(lep_i1.DeltaR(lep_i2));
                alldR.push_back(lep_i1.DeltaR(lep_j1));
                alldR.push_back(lep_i1.DeltaR(lep_j2));
                alldR.push_back(lep_i2.DeltaR(lep_j1));
                alldR.push_back(lep_i2.DeltaR(lep_j2));
                alldR.push_back(lep_j1.DeltaR(lep_j2));
                if (debug) cout<<" minDr: "<<*min_element(alldR.begin(),alldR.end())<<endl;
                if (*min_element(alldR.begin(),alldR.end())<deltaRCut) continue;
                
                // Check M(l+,l-)>4.0 GeV for any OS pair
                vector<double> allM;
                TLorentzVector i1i2;
                i1i2 = (lep_i1)+(lep_i2); allM.push_back(i1i2.M());
                TLorentzVector j1j2;
                j1j2 = (lep_j1)+(lep_j2); allM.push_back(j1j2.M());

                if ((*lep_id)[i1]*(*lep_id)[j1]<0) {
                    TLorentzVector i1j1;
                    i1j1 = (lep_i1)+(lep_j1); allM.push_back(i1j1.M());
                    TLorentzVector i2j2;
                    i2j2 = (lep_i2)+(lep_j2); allM.push_back(i2j2.M());
                } else {
                    TLorentzVector i1j2;
                    i1j2 = (lep_i1)+(lep_j2); allM.push_back(i1j2.M());
                    TLorentzVector i2j1;
                    i2j1 = (lep_i2)+(lep_j1); allM.push_back(i2j1.M());
                }
                if (debug) cout<<" min m(l+l-): "<<*min_element(allM.begin(),allM.end())<<endl;
                if (*min_element(allM.begin(),allM.end())<minMllCut) { continue;}

                // Do not include FSR photons
                //  if (*min_element(allM.begin(),allM.end())<0.1) { continue;}
                // Check the "smart cut": !( |mZa-mZ| < |mZ1-mZ| && mZb<12)
                // only for 4mu or 4e ZZ candidates
               // 

                bool passSmartCut=true;
                if ( abs((*lep_id)[i1])==abs((*lep_id)[j1])) {
                    TLorentzVector Za, Zb;
                    if ((*lep_id)[i1]==(*lep_id)[j1]) {                  
                        Za = (lep_i1)+(lep_j2);
                        Zb = (lep_i2)+(lep_j1);                    
                    } else {
                        Za = (lep_i1)+(lep_j1);
                        Zb = (lep_i2)+(lep_j2);
                    }                
                    if ( abs(Za.M()-Zmass)<abs(Zb.M()-Zmass) ) {
                        if (debug) cout<<"abs(Za.M()-Zmass)-abs(Z1.M()-Zmass): "
                                       <<abs(Za.M()-Zmass)-abs(Z1.M()-Zmass)<<" Zb.M(): "<<Zb.M()<<endl;
                        if ( abs(Za.M()-Zmass)<abs(Z1.M()-Zmass) && Zb.M()<mZ2LowCut ) passSmartCut=false;
                    }
                    else {
                        if (debug) cout<<"abs(Zb.M()-Zmass)-abs(Z1.M()-Zmass): "
                                       <<abs(Zb.M()-Zmass)-abs(Z1.M()-Zmass)<<" Za.M(): "<<Za.M()<<endl;
                        if ( abs(Zb.M()-Zmass)<abs(Z1.M()-Zmass) && Za.M()<mZ2LowCut ) passSmartCut=false;
                    }
                    
                }
                
                if (!passSmartCut) continue; 
                if (debug) cout<<" massZ1: "<<Z1.M()<<" massZ2: "<<Z2.M()<<endl;
                if (Z1.M() < mZ1LowCut) continue;
                if (Z1.M() > mZ1HighCut) continue;
                if (Z2.M() < mZ2LowCut) continue;
                if (Z2.M() > mZ2HighCut) continue; 
                if (debug) cout<<" pass Z mass cuts"<<endl;
                
                
                // Signal region if Z2 leptons are both tight ID Iso
                bool signalRegion=true;
                if ((*lep_RelIso)[Z2_lepindex[0]]>((abs((*lep_id)[Z2_lepindex[0]])==11) ? isoCutEl : isoCutMu)) signalRegion=false;
                if ((*lep_RelIso)[Z2_lepindex[1]]>((abs((*lep_id)[Z2_lepindex[1]])==11) ? isoCutEl : isoCutMu)) signalRegion=false;
                //if (!((*lep_tightId)[Z2_lepindex[0]])) signalRegion=false; // checking tight lepton ID
                //if (!((*lep_tightId)[Z2_lepindex[1]])) signalRegion=false; // checking tight lepton ID
                
                // Check if this candidate has the highest D_bkg_kin
                std::vector<TLorentzVector> P4s;
                P4s.clear();
                std::vector<int> tmpIDs;
                tmpIDs.clear();
                
                if (Z1_lepindex[0] == i1) {
                    P4s.push_back(lep_i1); P4s.push_back(lep_i2);
                    if (Z2_lepindex[0] == j1) {
                        P4s.push_back(lep_j1); P4s.push_back(lep_j2);
                    } else {
                        P4s.push_back(lep_j2); P4s.push_back(lep_j1);
                    }
                } else if (Z1_lepindex[0] == i2) {
                    P4s.push_back(lep_i2); P4s.push_back(lep_i1);
                    if (Z2_lepindex[0] == j1) {
                        P4s.push_back(lep_j1); P4s.push_back(lep_j2);
                    } else {
                        P4s.push_back(lep_j2); P4s.push_back(lep_j1);
                    }
                } else if (Z1_lepindex[0] == j1) {
                    P4s.push_back(lep_j1); P4s.push_back(lep_j2);
                    if (Z2_lepindex[0] == i1) {
                        P4s.push_back(lep_i1); P4s.push_back(lep_i2);
                    } else {
                        P4s.push_back(lep_i2); P4s.push_back(lep_i1);
                    }
                } else if (Z1_lepindex[0] == j2) {
                    P4s.push_back(lep_j2); P4s.push_back(lep_j1);
                    if (Z2_lepindex[0] == i1) {
                        P4s.push_back(lep_i1); P4s.push_back(lep_i2);
                    } else {
                        P4s.push_back(lep_i2); P4s.push_back(lep_i1);
                    }
                }
                
                tmpIDs.push_back((*lep_id)[Z1_lepindex[0]]); tmpIDs.push_back((*lep_id)[Z1_lepindex[1]]);
                tmpIDs.push_back((*lep_id)[Z2_lepindex[0]]); tmpIDs.push_back((*lep_id)[Z2_lepindex[1]]);

                bool same4l=false;
                bool foundZ11=false; bool foundZ12=false; bool foundZ21=false; bool foundZ22=false;
                for(int l = 0; l < 4; l++){
                    if (lep_Hindex[l]==Z1_lepindex[0]) foundZ11 = true;
                    if (lep_Hindex[l]==Z1_lepindex[1]) foundZ12 = true;
                    if (lep_Hindex[l]==Z2_lepindex[0]) foundZ21 = true;
                    if (lep_Hindex[l]==Z2_lepindex[1]) foundZ22 = true;
                }
                same4l = (foundZ11 && foundZ12 && foundZ21 && foundZ22);
                
                if (signalRegion) { // Signal Region has priority
                    
                    if (!foundSRCandidate) same4l=false;
                    

                    if (  Z1DeltaM<=minZ1DeltaM_SR ) { 
                        
                        minZ1DeltaM_SR = Z1DeltaM;
                        
                        if (Z_Hindex[0]==Z1index && Z2SumPt<maxZ2SumPt_SR) continue;

                        Z_Hindex[0] = Z1index;
                        lep_Hindex[0] = Z1_lepindex[0];
                        lep_Hindex[1] = Z1_lepindex[1];
                    
                        maxZ2SumPt_SR = Z2SumPt;    
                        Z_Hindex[1] = Z2index;
                        lep_Hindex[2] = Z2_lepindex[0];
                        lep_Hindex[3] = Z2_lepindex[1];
                        
                        Z1Vec = Z1; Z2Vec = Z2; HVec = Z1+Z2;                   
                        massZ1 = Z1Vec.M(); massZ2 = Z2Vec.M(); mass4l = HVec.M();
                        
                        if (debug) cout<<" new best candidate SR: mass4l: "<<HVec.M()<<endl;
                        if ((HVec.M()>m4lLowCut)&&(HVec.M()<m4lHighCut))  {
                            foundHiggsCandidate=true;                    
                            foundSRCandidate=true;
                        }
                    }
                } else if (!foundSRCandidate) { // Control regions get second priority

                    if (  Z1DeltaM<=minZ1DeltaM_CR ) {                 

                        minZ1DeltaM_CR = Z1DeltaM;
                        
                        if (Z_Hindex[0]==Z1index && Z2SumPt<maxZ2SumPt_CR) continue;

                        Z_Hindex[0] = Z1index;
                        lep_Hindex[0] = Z1_lepindex[0];
                        lep_Hindex[1] = Z1_lepindex[1];
                        
                        maxZ2SumPt_CR = Z2SumPt;
                        Z_Hindex[1] = Z2index;
                        lep_Hindex[2] = Z2_lepindex[0];
                        lep_Hindex[3] = Z2_lepindex[1];

                        Z1Vec = Z1; Z2Vec = Z2; HVec = Z1+Z2;                   
                        massZ1 = Z1Vec.M(); massZ2 = Z2Vec.M(); mass4l = HVec.M();

                        if (debug) cout<<" new best candidate CR: mass4l: "<<HVec.M()<<endl;
                        if (HVec.M()>m4lLowCut&&HVec.M()<m4lHighCut) foundHiggsCandidate=true;                    
                    }
                } 

                if (debug) cout<<"Z_Hindex[0]: "<<Z_Hindex[0]<<" lep_Hindex[0]: "<<lep_Hindex[0]<<" lep_Hindex[1]: "<<lep_Hindex[1]
                               <<"Z_Hindex[1]: "<<Z_Hindex[1]<<" lep_Hindex[2]: "<<lep_Hindex[2]<<" lep_Hindex[3]: "<<lep_Hindex[3]<<endl;
                
            } // Zj
        } // Zi

        if (foundHiggsCandidate) {

            if (debug) cout<<" lep_Hindex[0]: "<<lep_Hindex[0]<<" lep_Hindex[1]: "<<lep_Hindex[1]
                           <<" lep_Hindex[2]: "<<lep_Hindex[2]<<" lep_Hindex[3]: "<<lep_Hindex[3]<<endl;
                    

            TLorentzVector Lep1, Lep2, Lep3, Lep4,  Jet1, Jet2, Jet1_2p5, Jet2_2p5;            
            TLorentzVector nullFourVector(0, 0, 0, 0);                 
            Lep1.SetPtEtaPhiM((*lep_pt)[lep_Hindex[0]],(*lep_eta)[lep_Hindex[0]],(*lep_phi)[lep_Hindex[0]],(*lep_mass)[lep_Hindex[0]]);
            Lep2.SetPtEtaPhiM((*lep_pt)[lep_Hindex[1]],(*lep_eta)[lep_Hindex[1]],(*lep_phi)[lep_Hindex[1]],(*lep_mass)[lep_Hindex[1]]);
            Lep3.SetPtEtaPhiM((*lep_pt)[lep_Hindex[2]],(*lep_eta)[lep_Hindex[2]],(*lep_phi)[lep_Hindex[2]],(*lep_mass)[lep_Hindex[2]]);
            Lep4.SetPtEtaPhiM((*lep_pt)[lep_Hindex[3]],(*lep_eta)[lep_Hindex[3]],(*lep_phi)[lep_Hindex[3]],(*lep_mass)[lep_Hindex[3]]);

            nZXCRFailedLeptons = 0;
            for(unsigned int i = 0; i <= 3; i++) {
                if ((!(abs((*lep_id)[lep_Hindex[i]])==11 && (*lep_RelIso)[lep_Hindex[i]]<isoCutEl)) &&
                    !(abs((*lep_id)[lep_Hindex[i]])==13 && (*lep_RelIso)[lep_Hindex[i]]<isoCutMu)){
                    nZXCRFailedLeptons++;
                }
            }
            if (debug) cout << nZXCRFailedLeptons<<" failing leptons in higgs candidate"<<endl;
            if (nZXCRFailedLeptons>0) { // at least one lepton has failed 
                passedZXCRSelection = true;
            } else { //  signal region candidate                    
                passedFullSelection = true;
            }

            idL1 = (*lep_id)[lep_Hindex[0]]; pTL1 = Lep1.Pt(); etaL1 = Lep1.Eta();
            idL2 = (*lep_id)[lep_Hindex[1]]; pTL2 = Lep2.Pt(); etaL2 = Lep2.Eta();
            idL3 = (*lep_id)[lep_Hindex[2]]; pTL3 = Lep3.Pt(); etaL3 = Lep3.Eta();       
            idL4 = (*lep_id)[lep_Hindex[3]]; pTL4 = Lep4.Pt(); etaL4 = Lep4.Eta();
            phiL1 = Lep1.Phi(); 
            phiL2 = Lep2.Phi();
            phiL3 = Lep3.Phi();
            phiL4 = Lep4.Phi();
            std::vector<TLorentzVector> P4s; vector<int> tmpIDs;             
            P4s.push_back(Lep1); P4s.push_back(Lep2);
            P4s.push_back(Lep3); P4s.push_back(Lep4);
            tmpIDs.push_back(idL1); tmpIDs.push_back(idL2);
            tmpIDs.push_back(idL3); tmpIDs.push_back(idL4);
            computeAngles(Lep1,idL1,Lep2,idL2,Lep3,idL3,Lep4,idL4,cosThetaStar,cosTheta1,cosTheta2,Phi,Phi1);

            TLorentzVector HVec = Lep1+Lep2+Lep3+Lep4;   
            TLorentzVector Z1Vec = Lep1+Lep2;   
            TLorentzVector Z2Vec = Lep3+Lep4;   
            
            massZ1 = (Lep1+Lep2).M(); massZ2 = (Lep3+Lep4).M(); 
            mass4l = HVec.M();
            pT4l = HVec.Pt(); eta4l = HVec.Eta(); rapidity4l = HVec.Rapidity(); phi4l = HVec.Phi();
            pTZ1 = Z1Vec.Pt(); pTZ2 = Z2Vec.Pt(); massZ1 = Z1Vec.M(); massZ2 = Z2Vec.M(); 

            if (abs(idL1)==11 && abs(idL3)==11) {mass4e=mass4l; mass4mu=-1.0; mass2e2mu=-1.0;}
            else if (abs(idL1)==13 && abs(idL3)==13) {mass4e=-1.0; mass4mu=mass4l; mass2e2mu=-1.0;}
            else if (abs(idL1)!=abs(idL3)) {mass4e=-1.0; mass4mu=-1.0; mass2e2mu=mass4l;}

                
            if(debug) cout<<"fill tree"<<endl;
            if(debug) cout<<endl;
            
            outTree->Fill();   
        }
        
        //outTree->Fill();

        lep_id->clear();
        lep_pt->clear(); 
        lep_eta->clear(); 
        lep_phi->clear();
        lep_mass->clear();
        lep_RelIso->clear(); 

    }
}

//------------------------------------------------------------------------------
void InitTree(TTree* outTree)
{

    outTree->Branch("lep_id",&lep_id);
    outTree->Branch("lep_pt",&lep_pt);
    outTree->Branch("lep_eta",&lep_eta);
    outTree->Branch("lep_phi",&lep_phi);
    outTree->Branch("lep_mass",&lep_mass);
    outTree->Branch("lep_RelIso",&lep_RelIso);

    outTree->Branch("mass4l",&mass4l);
    outTree->Branch("mass4e",&mass4e);
    outTree->Branch("mass4mu",&mass4mu);
    outTree->Branch("mass2e2mu",&mass2e2mu);
    outTree->Branch("massZ1",&massZ1);
    outTree->Branch("massZ2",&massZ2);
    outTree->Branch("cosTheta1",&cosTheta1);
    outTree->Branch("cosTheta2",&cosTheta2);
    outTree->Branch("cosThetaStar",&cosThetaStar);
    outTree->Branch("Phi",&Phi);
    outTree->Branch("Phi1",&Phi1);
    
    outTree->Branch("pT4l",&pT4l);
    outTree->Branch("eta4l",&eta4l);
    outTree->Branch("phi4l",&phi4l);
    outTree->Branch("pTZ1",&pTZ1);
    outTree->Branch("etaZ1",&etaZ1);
    outTree->Branch("phiZ1",&phiZ1);
    outTree->Branch("pTZ2",&pTZ2);
    outTree->Branch("etaZ2",&etaZ2);
    outTree->Branch("phiZ2",&phiZ2);

    outTree->Branch("pTL1",&pTL1);
    outTree->Branch("etaL1",&etaL1);
    outTree->Branch("phiL1",&phiL1);
    outTree->Branch("idL1",&idL1);
    outTree->Branch("pTL2",&pTL2);
    outTree->Branch("etaL2",&etaL2);
    outTree->Branch("phiL2",&phiL2);
    outTree->Branch("idL2",&idL2);
    outTree->Branch("pTL3",&pTL3);
    outTree->Branch("etaL3",&etaL3);
    outTree->Branch("phiL3",&phiL3);
    outTree->Branch("idL3",&idL3);
    outTree->Branch("pTL4",&pTL4);
    outTree->Branch("etaL4",&etaL4);
    outTree->Branch("phiL4",&phiL4);
    outTree->Branch("idL4",&idL4);
    outTree->Branch("passedFullSelection",&passedFullSelection);
    outTree->Branch("passedZXCRSelection",&passedZXCRSelection);
    outTree->Branch("nZXCRFailedLeptons",&nZXCRFailedLeptons);

}

//------------------------------------------------------------------------------

void LiteTreeProducer(const char *inputFile, const char *outputFile)
{
    gSystem->Load("libDelphes");
    
    TChain *chain = new TChain("Delphes");
    chain->Add(inputFile);
    
    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);

    TFile* outFile = new TFile(outputFile,"RECREATE");
    TTree* outTree = new TTree("LiteTree","LiteTree");

    InitTree(outTree);
    AnalyseEvents(treeReader,outTree);

    outFile->cd();
    outTree->Write("LiteTree",TObject::kOverwrite);
    outFile->Close();

    std::cout << "** Exiting..." << std::endl;
    
    delete treeReader;
    delete chain;
}

//------------------------------------------------------------------------------
