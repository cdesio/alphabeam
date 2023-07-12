//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************

#include "RunAction.hh"
#include "G4Run.hh"
#include "G4AnalysisManager.hh"
#include "globals.hh"
#include <map>
#include "CommandLineParser.hh"
#include "G4EventManager.hh"
#include "EventAction.hh"
#include "G4Event.hh"
#include "DetectorConstruction.hh"
#include "git_version.hh"
#include "G4SystemOfUnits.hh" 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
using namespace G4DNAPARSER;

RunAction::RunAction()
    : G4UserRunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::~RunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::BeginOfRunAction(const G4Run *)
{
    CommandLineParser *parser = CommandLineParser::GetParser();
    Command *command(0);
    if ((command = parser->GetCommandIfActive("-out")) == 0)
        return;

    // Open an output file
    G4String fileName{"output.root"};
    if (command->GetOption().empty() == false)
    {
        fileName = command->GetOption();
    }

    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("root");
    analysisManager->SetVerboseLevel(0);

    // open output file
    //
    G4bool fileOpen = analysisManager->OpenFile(fileName);
    if (!fileOpen)
    {
        G4cout << "\n---> HistoManager::book(): cannot open " << fileName << G4endl;
        return;
    }

    G4cout << "\n----> Histogram file is opened in " << fileName << G4endl;

    analysisManager->CreateNtuple("Info", "Info");
    analysisManager->CreateNtupleDColumn("NumPrimaries");
    analysisManager->CreateNtupleDColumn("Rmin");
    analysisManager->CreateNtupleDColumn("Rmax");
    analysisManager->CreateNtupleIColumn("Nrings");
    analysisManager->CreateNtupleIColumn("NperRing");
    analysisManager->CreateNtupleSColumn("GitHash");
    analysisManager->CreateNtupleDColumn("Rn220Desorption");
    analysisManager->CreateNtupleDColumn("Pb212Desorption");
    analysisManager->CreateNtupleDColumn("Pb212Leakage");
    analysisManager->FinishNtuple(0);

    analysisManager->CreateH1("0", "dose", 5000, 0, 5000);



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::EndOfRunAction(const G4Run *run)
{
    Write(run);
    auto fpEventAction = (EventAction *)G4EventManager::GetEventManager()->GetUserEventAction();

    G4int numPrimaries = run->GetNumberOfEvent();
    G4double RnDesorptionIN = fpEventAction->getRnDesorptionIN();
    G4cout << "Desorption of Rn220 from is " << (1 - RnDesorptionIN / numPrimaries) * 100 << "%, expect 40%." << G4endl;

    G4double PbDesorptionIN = fpEventAction->getPbDesorptionIN();
    G4cout << "Desorption of Pb212 from is " << (1 - PbDesorptionIN / numPrimaries) * 100 << "%, expect 55%." << G4endl;

    G4double PbLeakage = fpEventAction->getPbLeakage();
    G4double PbNoLeakage = fpEventAction->getPbNoLeakage();
    G4cout << "Leakage of Pb212 from is " << PbLeakage / (PbLeakage + PbNoLeakage) * 100 << "%, value depends on tumour size" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void RunAction::Write(const G4Run* run)
{
    CommandLineParser *parser = CommandLineParser::GetParser();
    Command *command(0);
    if ((command = parser->GetCommandIfActive("-out")) == 0)
        return;
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

    auto fpEventAction = (EventAction *)G4EventManager::GetEventManager()->GetUserEventAction();

    G4int numPrimaries = run->GetNumberOfEvent();
    G4double RnDesorptionIN = fpEventAction->getRnDesorptionIN();
    G4double PbDesorptionIN = fpEventAction->getPbDesorptionIN();
    G4double PbLeakage = fpEventAction->getPbLeakage();
    G4double PbNoLeakage = fpEventAction->getPbNoLeakage();

    analysisManager->FillNtupleDColumn(0,0, run->GetNumberOfEvent());
    analysisManager->FillNtupleDColumn(0, 1, Rmin/um);
    analysisManager->FillNtupleDColumn(0,2, Rmax/um);
    analysisManager->FillNtupleIColumn(0,3, Nrings);
    analysisManager->FillNtupleIColumn(0,4, NperRing);
    analysisManager->FillNtupleSColumn(0,5, kGitHash);
    analysisManager->FillNtupleDColumn(0,6, 1 - RnDesorptionIN / numPrimaries);
    analysisManager->FillNtupleDColumn(0,7, 1 - PbDesorptionIN / numPrimaries);
    analysisManager->FillNtupleDColumn(0,8, PbLeakage / (PbLeakage + PbNoLeakage));

    analysisManager->AddNtupleRow(0);

    analysisManager->Write();
    analysisManager->CloseFile();
    analysisManager->Clear();
    G4cout << "\n----> Histograms are saved" << G4endl;
}

