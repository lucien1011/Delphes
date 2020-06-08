/*
This macro is to pre-process the delphes tree to make muon fast simulation ntuples 
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
#include <iostream>
#else
class ExRootTreeReader;
#endif

//------------------------------------------------------------------------------

using namespace std;

std::vector<float> *Muon_Pt = 0; 
std::vector<float> *Muon_Eta = 0; 
std::vector<float> *Muon_Phi = 0;
std::vector<float> *Muon_IsolationVar = 0;
std::vector<float> *Muon_T = 0;
std::vector<float> *Muon_Charge = 0;
std::vector<float> *Muon_IsolationVarRhoCorr = 0;
std::vector<float> *Muon_SumPtCharged = 0;
std::vector<float> *Muon_SumPtNeutral = 0;   
std::vector<float> *Muon_SumPtChargedPU = 0;    
std::vector<float> *Muon_SumPt = 0;
std::vector<float> *Muon_GenPt = 0; 
std::vector<float> *Muon_GenEta = 0; 
std::vector<float> *Muon_GenPhi = 0;
std::vector<float> *Muon_GenMass = 0;
std::vector<float> *Muon_GenStatus = 0;
std::vector<float> *Muon_GenIsPU = 0;
std::vector<float> *Muon_GenCharge = 0;
std::vector<float> *Muon_GenD0 = 0;
std::vector<float> *Muon_GenDZ = 0;
std::vector<float> *Muon_GenVx = 0;
std::vector<float> *Muon_GenVy = 0;
std::vector<float> *Muon_GenVz = 0;
std::vector<float> *Muon_GenVt = 0;
std::vector<float> *Muon_TrackPt = 0; 
std::vector<float> *Muon_TrackEta = 0; 
std::vector<float> *Muon_TrackPhi = 0;
std::vector<float> *Muon_TrackCharge = 0;
std::vector<float> *Muon_TrackD0 = 0;
std::vector<float> *Muon_TrackDZ = 0;
std::vector<float> *Muon_TrackVx = 0;
std::vector<float> *Muon_TrackVy = 0;
std::vector<float> *Muon_TrackVz = 0;
std::vector<float> *Muon_TrackVt = 0;
std::vector<float> *Muon_TrackOuterx = 0;
std::vector<float> *Muon_TrackOutery = 0;
std::vector<float> *Muon_TrackOuterz = 0;
std::vector<float> *Muon_TrackOutert = 0;
std::vector<float> *Muon_TrackVxd = 0;
std::vector<float> *Muon_TrackVyd = 0;
std::vector<float> *Muon_TrackVzd = 0;
std::vector<float> *Muon_TrackL = 0;
std::vector<float> *Muon_TrackPtError = 0;
std::vector<float> *Muon_TrackPhiError = 0;
std::vector<float> *Muon_TrackTError = 0;
std::vector<float> *Muon_TrackD0Error = 0;
std::vector<float> *Muon_TrackDZError = 0;

//------------------------------------------------------------------------------

void AnalyseEvents(ExRootTreeReader* treeReader, TTree* outTree) {
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchTrack = treeReader->UseBranch("Track");
    TClonesArray *branchMuon = treeReader->UseBranch("Muon");
    
    Long64_t allEntries = treeReader->GetEntries();
    
    std::cout << "** Chain contains " << allEntries << " events" << std::endl;
    
    Long64_t entry;
    
    Int_t i, j, k, pdgCode;
    
    Track *track;
    Muon *muon;
    GenParticle *particle, *gentrack;

    for(entry = 0; entry < allEntries; ++entry) {
        treeReader->ReadEntry(entry);
        
        if(entry %1000 == 0) std::cout<<entry<< std::endl;

        for(i = 0; i < branchMuon->GetEntriesFast(); ++i) {
            muon = (Muon*) branchMuon->At(i);
            particle = (GenParticle*) muon->Particle.GetObject();
            for(k = 0; k < branchTrack->GetEntriesFast(); ++k) {
                track = (Track*) branchTrack->At(k);
                gentrack = (GenParticle*) track->Particle.GetObject();
                if(gentrack->GetUniqueID() == particle->GetUniqueID()) {
                    Muon_Pt->push_back(muon->PT);
                    Muon_Eta->push_back(muon->Eta);
                    Muon_Phi->push_back(muon->Phi);
                    Muon_IsolationVar->push_back(muon->IsolationVar);
                    Muon_T->push_back(muon->T);
                    Muon_Charge->push_back(muon->Charge);
                    Muon_IsolationVarRhoCorr->push_back(muon->IsolationVarRhoCorr);
                    Muon_SumPtCharged->push_back(muon->SumPtCharged);
                    Muon_SumPtNeutral->push_back(muon->SumPtNeutral);   
                    Muon_SumPtChargedPU->push_back(muon->SumPtChargedPU);
                    Muon_SumPt->push_back(muon->SumPt);
                    Muon_GenPt->push_back(particle->PT); 
                    Muon_GenEta->push_back(particle->Eta); 
                    Muon_GenPhi->push_back(particle->Phi);
                    Muon_GenMass->push_back(particle->Mass);
                    Muon_GenStatus->push_back(particle->Status);
                    Muon_GenIsPU->push_back(particle->IsPU);
                    Muon_GenCharge->push_back(particle->Charge);
                    Muon_GenD0->push_back(particle->D0);
                    Muon_GenDZ->push_back(particle->DZ);
                    Muon_GenVx->push_back(particle->X);
                    Muon_GenVy->push_back(particle->Y);
                    Muon_GenVz->push_back(particle->Z);
                    Muon_GenVt->push_back(particle->T);
                    Muon_TrackPt->push_back(track->PT); 
                    Muon_TrackEta->push_back(track->Eta); 
                    Muon_TrackPhi->push_back(track->Phi);
                    Muon_TrackCharge->push_back(track->Charge);
                    Muon_TrackD0->push_back(track->D0);
                    Muon_TrackDZ->push_back(track->DZ);
                    Muon_TrackVx->push_back(track->X);
                    Muon_TrackVy->push_back(track->Y);
                    Muon_TrackVz->push_back(track->Z);
                    Muon_TrackVt->push_back(track->T);
                    Muon_TrackOuterx->push_back(track->XOuter);
                    Muon_TrackOutery->push_back(track->YOuter);
                    Muon_TrackOuterz->push_back(track->ZOuter);
                    Muon_TrackOutert->push_back(track->TOuter);
                    Muon_TrackVxd->push_back(track->Xd);
                    Muon_TrackVyd->push_back(track->Yd);
                    Muon_TrackVzd->push_back(track->Zd);
                    Muon_TrackL->push_back(track->L);
                    Muon_TrackPtError->push_back(track->ErrorPT);
                    Muon_TrackPhiError->push_back(track->ErrorPhi);
                    Muon_TrackTError->push_back(track->ErrorT);
                    Muon_TrackD0Error->push_back(track->ErrorD0);
                    Muon_TrackDZError->push_back(track->ErrorDZ);
                }
            }
        }
        
        outTree->Fill();

        Muon_Pt->clear(); 
        Muon_Eta->clear(); 
        Muon_Phi->clear();
        Muon_IsolationVar->clear();
        Muon_T->clear();
        Muon_Charge->clear();
        Muon_IsolationVarRhoCorr->clear();
        Muon_SumPtCharged->clear();
        Muon_SumPtNeutral->clear();   
        Muon_SumPtChargedPU->clear();    
        Muon_SumPt->clear();
        Muon_GenPt->clear(); 
        Muon_GenEta->clear(); 
        Muon_GenPhi->clear();
        Muon_GenMass->clear();
        Muon_GenStatus->clear();
        Muon_GenIsPU->clear();
        Muon_GenCharge->clear();
        Muon_GenD0->clear();
        Muon_GenDZ->clear();
        Muon_GenVx->clear();
        Muon_GenVy->clear();
        Muon_GenVz->clear();
        Muon_GenVt->clear();
        Muon_TrackPt->clear(); 
        Muon_TrackEta->clear(); 
        Muon_TrackPhi->clear();
        Muon_TrackCharge->clear();
        Muon_TrackD0->clear();
        Muon_TrackDZ->clear();
        Muon_TrackVx->clear();
        Muon_TrackVy->clear();
        Muon_TrackVz->clear();
        Muon_TrackVt->clear();
        Muon_TrackOuterx->clear();
        Muon_TrackOutery->clear();
        Muon_TrackOuterz->clear();
        Muon_TrackOutert->clear();
        Muon_TrackVxd->clear();
        Muon_TrackVyd->clear();
        Muon_TrackVzd->clear();
        Muon_TrackL->clear();
        Muon_TrackPtError->clear();
        Muon_TrackPhiError->clear();
        Muon_TrackTError->clear();
        Muon_TrackD0Error->clear();
        Muon_TrackDZError->clear();
    }
}

//------------------------------------------------------------------------------
void InitTree(TTree* outTree)
{
    outTree->Branch("Muon_Pt",&Muon_Pt);
    outTree->Branch("Muon_Eta",&Muon_Eta);
    outTree->Branch("Muon_Phi",&Muon_Phi); 
    outTree->Branch("Muon_T",&Muon_T);
    outTree->Branch("Muon_Charge",&Muon_Charge);
    outTree->Branch("Muon_IsolationVar",&Muon_IsolationVar);
    outTree->Branch("Muon_IsolationVarRhoCorr",&Muon_IsolationVarRhoCorr);
    outTree->Branch("Muon_SumPtCharged",&Muon_SumPtCharged);
    outTree->Branch("Muon_SumPtNeutral",&Muon_SumPtNeutral);
    outTree->Branch("Muon_SumPtChargedPU",&Muon_SumPtChargedPU);
    outTree->Branch("Muon_SumPt",&Muon_SumPt);
    outTree->Branch("Muon_GenPt",&Muon_GenPt); 
    outTree->Branch("Muon_GenEta",&Muon_GenEta); 
    outTree->Branch("Muon_GenPhi",&Muon_GenPhi);
    outTree->Branch("Muon_GenMass",&Muon_GenMass);
    outTree->Branch("Muon_GenStatus",&Muon_GenStatus);
    outTree->Branch("Muon_GenIsPU",&Muon_GenIsPU);
    outTree->Branch("Muon_GenCharge",&Muon_GenCharge);
    outTree->Branch("Muon_GenD0",&Muon_GenD0);
    outTree->Branch("Muon_GenDZ",&Muon_GenDZ);
    outTree->Branch("Muon_GenVx",&Muon_GenVx);
    outTree->Branch("Muon_GenVy",&Muon_GenVy);
    outTree->Branch("Muon_GenVz",&Muon_GenVz);
    outTree->Branch("Muon_GenVt",&Muon_GenVt);
    outTree->Branch("Muon_TrackPt",&Muon_TrackPt); 
    outTree->Branch("Muon_TrackEta",&Muon_TrackEta); 
    outTree->Branch("Muon_TrackPhi",&Muon_TrackPhi);
    outTree->Branch("Muon_TrackCharge",&Muon_TrackCharge);
    outTree->Branch("Muon_TrackD0",&Muon_TrackD0);
    outTree->Branch("Muon_TrackDZ",&Muon_TrackDZ);
    outTree->Branch("Muon_TrackD0Error",&Muon_TrackD0Error);
    outTree->Branch("Muon_TrackDZError",&Muon_TrackDZError);
    outTree->Branch("Muon_TrackVx",&Muon_TrackVx);
    outTree->Branch("Muon_TrackVy",&Muon_TrackVy);
    outTree->Branch("Muon_TrackVz",&Muon_TrackVz);
    outTree->Branch("Muon_TrackVt",&Muon_TrackVt);
    outTree->Branch("Muon_TrackOuterx",&Muon_TrackOuterx);
    outTree->Branch("Muon_TrackOutery",&Muon_TrackOutery);
    outTree->Branch("Muon_TrackOuterz",&Muon_TrackOuterz);
    outTree->Branch("Muon_TrackOutert",&Muon_TrackOutert);
    outTree->Branch("Muon_TrackVxd",&Muon_TrackVxd);
    outTree->Branch("Muon_TrackVyd",&Muon_TrackVyd);
    outTree->Branch("Muon_TrackVzd",&Muon_TrackVzd);
    outTree->Branch("Muon_TrackL",&Muon_TrackL);
    outTree->Branch("Muon_TrackPtError",&Muon_TrackPtError); 
    outTree->Branch("Muon_TrackTError",&Muon_TrackTError); 
    outTree->Branch("Muon_TrackPhiError",&Muon_TrackPhiError);
}

//------------------------------------------------------------------------------

void MuonTreeProducer(const char *inputFile, const char *outputFile)
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
