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
#else
class ExRootTreeReader;
#endif

//------------------------------------------------------------------------------

std::vector<float> *Lep_Pt = 0; 
std::vector<float> *Lep_Eta = 0; 
std::vector<float> *Lep_Phi = 0;
std::vector<float> *Lep_IsolationVar = 0;
std::vector<float> *Lep_T = 0;
std::vector<float> *Lep_Charge = 0;
std::vector<float> *Lep_IsolationVarRhoCorr = 0;
std::vector<float> *Lep_SumPtCharged = 0;
std::vector<float> *Lep_SumPtNeutral = 0;   
std::vector<float> *Lep_SumPtChargedPU = 0;    
std::vector<float> *Lep_SumPt = 0;
std::vector<float> *Lep_GenPt = 0; 
std::vector<float> *Lep_GenEta = 0; 
std::vector<float> *Lep_GenPhi = 0;
std::vector<float> *Lep_GenMass = 0;
std::vector<float> *Lep_GenStatus = 0;
std::vector<float> *Lep_GenIsPU = 0;
std::vector<float> *Lep_GenCharge = 0;
std::vector<float> *Lep_GenD0 = 0;
std::vector<float> *Lep_GenDZ = 0;
std::vector<float> *Lep_GenVx = 0;
std::vector<float> *Lep_GenVy = 0;
std::vector<float> *Lep_GenVz = 0;
std::vector<float> *Lep_GenVt = 0;

//------------------------------------------------------------------------------

void AnalyseEvents(ExRootTreeReader* treeReader, TTree* outTree) {
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchTrack = treeReader->UseBranch("Track");
    TClonesArray *branchMuon = treeReader->UseBranch("Muon");
    
    Long64_t allEntries = treeReader->GetEntries();
    
    cout << "** Chain contains " << allEntries << " events" << endl;
    
    Long64_t entry;
    
    Int_t i, j, k, pdgCode;
    
    Track *track;
    Muon *muon;
    GenParticle *particle, *gentrack;

    for(entry = 0; entry < allEntries; ++entry) {
        treeReader->ReadEntry(entry);
        
        if(entry %1000 == 0) cout<<entry<< endl;

        for(i = 0; i < branchMuon->GetEntriesFast(); ++i) {
            muon = (Muon*) branchMuon->At(i);
            particle = (GenParticle*) muon->Particle.GetObject();
            for(k = 0; k < branchTrack->GetEntriesFast(); ++k) {
                track = (Track*) branchTrack->At(k);
                gentrack = (GenParticle*) track->Particle.GetObject();
                if(gentrack->GetUniqueID() == particle->GetUniqueID()) {
                    Lep_Pt->push_back(muon->PT);
                    Lep_Eta->push_back(muon->Eta);
                    Lep_Phi->push_back(muon->Phi);
                    Lep_IsolationVar->push_back(muon->IsolationVar);
                    Lep_T->push_back(muon->T);
                    Lep_Charge->push_back(muon->Charge);
                    Lep_IsolationVarRhoCorr->push_back(muon->IsolationVarRhoCorr);
                    Lep_SumPtCharged->push_back(muon->SumPtCharged);
                    Lep_SumPtNeutral->push_back(muon->SumPtNeutral);   
                    Lep_SumPtChargedPU->push_back(muon->SumPtChargedPU);
                    Lep_SumPt->push_back(muon->SumPt);
                    Lep_GenPt->push_back(particle->PT); 
                    Lep_GenEta->push_back(particle->Eta); 
                    Lep_GenPhi->push_back(particle->Phi);
                    Lep_GenMass->push_back(particle->Mass);
                    Lep_GenStatus->push_back(particle->Status);
                    Lep_GenIsPU->push_back(particle->IsPU);
                    Lep_GenCharge->push_back(particle->Charge);
                    Lep_GenD0->push_back(particle->D0);
                    Lep_GenDZ->push_back(particle->DZ);
                    Lep_GenVx->push_back(particle->X);
                    Lep_GenVy->push_back(particle->Y);
                    Lep_GenVz->push_back(particle->Z);
                    Lep_GenVt->push_back(particle->T);
                }
            }
        }
        
        outTree->Fill();

        Lep_Pt->clear(); 
        Lep_Eta->clear(); 
        Lep_Phi->clear();
        Lep_IsolationVar->clear();
        Lep_T->clear();
        Lep_Charge->clear();
        Lep_IsolationVarRhoCorr->clear();
        Lep_SumPtCharged->clear();
        Lep_SumPtNeutral->clear();   
        Lep_SumPtChargedPU->clear();    
        Lep_SumPt->clear();
        Lep_GenPt->clear(); 
        Lep_GenEta->clear(); 
        Lep_GenPhi->clear();
        Lep_GenMass->clear();
        Lep_GenStatus->clear();
        Lep_GenIsPU->clear();
        Lep_GenCharge->clear();
        Lep_GenD0->clear();
        Lep_GenDZ->clear();
        Lep_GenVx->clear();
        Lep_GenVy->clear();
        Lep_GenVz->clear();
        Lep_GenVt->clear();
    }
}

//------------------------------------------------------------------------------
void InitTree(TTree* outTree)
{
    outTree->Branch("Lep_Pt",&Lep_Pt);
    outTree->Branch("Lep_Eta",&Lep_Eta);
    outTree->Branch("Lep_Phi",&Lep_Phi); 
    outTree->Branch("Lep_T",&Lep_T);
    outTree->Branch("Lep_Charge",&Lep_Charge);
    outTree->Branch("Lep_IsolationVar",&Lep_IsolationVar);
    outTree->Branch("Lep_IsolationVarRhoCorr",&Lep_IsolationVarRhoCorr);
    outTree->Branch("Lep_SumPtCharged",&Lep_SumPtCharged);
    outTree->Branch("Lep_SumPtNeutral",&Lep_SumPtNeutral);
    outTree->Branch("Lep_SumPtChargedPU",&Lep_SumPtChargedPU);
    outTree->Branch("Lep_SumPt",&Lep_SumPt);
    outTree->Branch("Lep_GenPt",&Lep_GenPt); 
    outTree->Branch("Lep_GenEta",&Lep_GenEta); 
    outTree->Branch("Lep_GenPhi",&Lep_GenPhi);
    outTree->Branch("Lep_GenMass",&Lep_GenMass);
    outTree->Branch("Lep_GenStatus",&Lep_GenStatus);
    outTree->Branch("Lep_GenIsPU",&Lep_GenIsPU);
    outTree->Branch("Lep_GenCharge",&Lep_GenCharge);
    outTree->Branch("Lep_GenD0",&Lep_GenD0);
    outTree->Branch("Lep_GenDZ",&Lep_GenDZ);
    outTree->Branch("Lep_GenVx",&Lep_GenVx);
    outTree->Branch("Lep_GenVy",&Lep_GenVy);
    outTree->Branch("Lep_GenVz",&Lep_GenVz);
    outTree->Branch("Lep_GenVt",&Lep_GenVt);
}

//------------------------------------------------------------------------------

void MuonTreeProducer(const char *inputFile)
{
    gSystem->Load("libDelphes");
    
    TChain *chain = new TChain("Delphes");
    chain->Add(inputFile);
    
    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);

    TFile* outFile = new TFile("results.root","RECREATE");
    TTree* outTree = new TTree("LiteTree","LiteTree");

    InitTree(outTree);
    AnalyseEvents(treeReader,outTree);

    outFile->cd();
    outTree->Write("LiteTree",TObject::kOverwrite);
    outFile->Close();

    cout << "** Exiting..." << endl;
    
    delete treeReader;
    delete chain;
}

//------------------------------------------------------------------------------
