/*
This macro is to pre-process the delphes tree to make photon fast simulation ntuples 
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

std::vector<float> *Photon_Pt = 0; 
std::vector<float> *Photon_Eta = 0; 
std::vector<float> *Photon_Phi = 0;
std::vector<float> *Photon_E = 0;
std::vector<float> *Photon_T = 0;
std::vector<float> *Photon_IsolationVar = 0;
std::vector<float> *Photon_IsolationVarRhoCorr = 0;
std::vector<float> *Photon_SumPtCharged = 0;
std::vector<float> *Photon_SumPtNeutral = 0;   
std::vector<float> *Photon_SumPtChargedPU = 0;    
std::vector<float> *Photon_SumPt = 0;
std::vector<float> *Photon_EhadOverEem = 0;    
std::vector<float> *Photon_GenPt = 0; 
std::vector<float> *Photon_GenEta = 0; 
std::vector<float> *Photon_GenPhi = 0;
std::vector<float> *Photon_GenMass = 0;
std::vector<float> *Photon_GenStatus = 0;
std::vector<float> *Photon_GenIsPU = 0;
std::vector<float> *Photon_GenCharge = 0;
std::vector<float> *Photon_GenD0 = 0;
std::vector<float> *Photon_GenDZ = 0;
std::vector<float> *Photon_GenVx = 0;
std::vector<float> *Photon_GenVy = 0;
std::vector<float> *Photon_GenVz = 0;
std::vector<float> *Photon_GenVt = 0;

//------------------------------------------------------------------------------

void AnalyseEvents(ExRootTreeReader* treeReader, TTree* outTree) {
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
    
    Long64_t allEntries = treeReader->GetEntries();
    
    std::cout << "** Chain contains " << allEntries << " events" << std::endl;
    
    Long64_t entry;
    
    Int_t i, j, k, pdgCode;
    
    Photon *photon;
    GenParticle *particle;

    for(entry = 0; entry < allEntries; ++entry) {
        treeReader->ReadEntry(entry);
        
        if(entry %1000 == 0) std::cout<<entry<< std::endl;

        for(i = 0; i < branchPhoton->GetEntriesFast(); ++i) {
            photon = (Photon*) branchPhoton->At(i);
            particle = (GenParticle*) photon->Particle.GetObject();
            
            Photon_Pt->push_back(photon->PT);
            Photon_Eta->push_back(photon->Eta);
            Photon_Phi->push_back(photon->Phi);
            Photon_IsolationVar->push_back(photon->IsolationVar);
            Photon_T->push_back(photon->T);
            Photon_IsolationVarRhoCorr->push_back(photon->IsolationVarRhoCorr);
            Photon_SumPtCharged->push_back(photon->SumPtCharged);
            Photon_SumPtNeutral->push_back(photon->SumPtNeutral);   
            Photon_SumPtChargedPU->push_back(photon->SumPtChargedPU);
            Photon_SumPt->push_back(photon->SumPt);
            Photon_EhadOverEem->push_back(photon->EhadOverEem);
            
            Photon_GenPt->push_back(particle->PT); 
            Photon_GenEta->push_back(particle->Eta); 
            Photon_GenPhi->push_back(particle->Phi);
            Photon_GenMass->push_back(particle->Mass);
            Photon_GenStatus->push_back(particle->Status);
            Photon_GenIsPU->push_back(particle->IsPU);
            Photon_GenCharge->push_back(particle->Charge);
            Photon_GenD0->push_back(particle->D0);
            Photon_GenDZ->push_back(particle->DZ);
            Photon_GenVx->push_back(particle->X);
            Photon_GenVy->push_back(particle->Y);
            Photon_GenVz->push_back(particle->Z);
            Photon_GenVt->push_back(particle->T);
            }
        }
        
        outTree->Fill();

        Photon_Pt->clear(); 
        Photon_Eta->clear(); 
        Photon_Phi->clear();
        Photon_IsolationVar->clear();
        Photon_T->clear();
        Photon_Charge->clear();
        Photon_IsolationVarRhoCorr->clear();
        Photon_SumPtCharged->clear();
        Photon_SumPtNeutral->clear();   
        Photon_SumPtChargedPU->clear();    
        Photon_SumPt->clear();
        Photon_EhadOverEem->clear();

        Photon_GenPt->clear(); 
        Photon_GenEta->clear(); 
        Photon_GenPhi->clear();
        Photon_GenMass->clear();
        Photon_GenStatus->clear();
        Photon_GenIsPU->clear();
        Photon_GenCharge->clear();
        Photon_GenD0->clear();
        Photon_GenDZ->clear();
        Photon_GenVx->clear();
        Photon_GenVy->clear();
        Photon_GenVz->clear();
        Photon_GenVt->clear();
    }
}

//------------------------------------------------------------------------------
void InitTree(TTree* outTree)
{
    outTree->Branch("Photon_Pt",&Photon_Pt);
    outTree->Branch("Photon_Eta",&Photon_Eta);
    outTree->Branch("Photon_Phi",&Photon_Phi); 
    outTree->Branch("Photon_T",&Photon_T);
    outTree->Branch("Photon_Charge",&Photon_Charge);
    outTree->Branch("Photon_IsolationVar",&Photon_IsolationVar);
    outTree->Branch("Photon_IsolationVarRhoCorr",&Photon_IsolationVarRhoCorr);
    outTree->Branch("Photon_SumPtCharged",&Photon_SumPtCharged);
    outTree->Branch("Photon_SumPtNeutral",&Photon_SumPtNeutral);
    outTree->Branch("Photon_SumPtChargedPU",&Photon_SumPtChargedPU);
    outTree->Branch("Photon_SumPt",&Photon_SumPt);
    outTree->Branch("Photon_EhadOverEem",&Photon_EhadOverEem);

    outTree->Branch("Photon_GenPt",&Photon_GenPt); 
    outTree->Branch("Photon_GenEta",&Photon_GenEta); 
    outTree->Branch("Photon_GenPhi",&Photon_GenPhi);
    outTree->Branch("Photon_GenMass",&Photon_GenMass);
    outTree->Branch("Photon_GenStatus",&Photon_GenStatus);
    outTree->Branch("Photon_GenIsPU",&Photon_GenIsPU);
    outTree->Branch("Photon_GenCharge",&Photon_GenCharge);
    outTree->Branch("Photon_GenD0",&Photon_GenD0);
    outTree->Branch("Photon_GenDZ",&Photon_GenDZ);
    outTree->Branch("Photon_GenVx",&Photon_GenVx);
    outTree->Branch("Photon_GenVy",&Photon_GenVy);
    outTree->Branch("Photon_GenVz",&Photon_GenVz);
    outTree->Branch("Photon_GenVt",&Photon_GenVt);
}

//------------------------------------------------------------------------------

void PhotonTreeProducer(const char *inputFile, const char *outputFile)
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
