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

    // dose and strand breaks are calculated at 10 radial distances from the seed, in the nucleus volume

    analysisManager->CreateH1("DSB_0", "DSB_0", 336, 0, 336); // time (1hr)
    analysisManager->CreateH1("DSB_1", "DSB_1", 336, 0, 336); 
    analysisManager->CreateH1("DSB_2", "DSB_2", 336, 0, 336); 
    analysisManager->CreateH1("DSB_3", "DSB_3", 336, 0, 336); 
    analysisManager->CreateH1("DSB_4", "DSB_4", 336, 0, 336); 
    analysisManager->CreateH1("DSB_5", "DSB_5", 336, 0, 336); 
    analysisManager->CreateH1("DSB_6", "DSB_6", 336, 0, 336); 
    analysisManager->CreateH1("DSB_7", "DSB_7", 336, 0, 336); 
    analysisManager->CreateH1("DSB_8", "DSB_8", 336, 0, 336); 
    analysisManager->CreateH1("DSB_9", "DSB_9", 336, 0, 336);
    analysisManager->CreateH1("Dose_0", "Dose_0", 336, 0, 336);
    analysisManager->CreateH1("Dose_1", "Dose_1", 336, 0, 336);
    analysisManager->CreateH1("Dose_2", "Dose_2", 336, 0, 336);
    analysisManager->CreateH1("Dose_3", "Dose_3", 336, 0, 336);
    analysisManager->CreateH1("Dose_4", "Dose_4", 336, 0, 336);
    analysisManager->CreateH1("Dose_5", "Dose_5", 336, 0, 336);
    analysisManager->CreateH1("Dose_6", "Dose_6", 336, 0, 336);
    analysisManager->CreateH1("Dose_7", "Dose_7", 336, 0, 336);
    analysisManager->CreateH1("Dose_8", "Dose_8", 336, 0, 336);
    analysisManager->CreateH1("Dose_9", "Dose_9", 336, 0, 336);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::EndOfRunAction(const G4Run *run)
{
    Write(run);
    auto fpEventAction = (EventAction *)G4EventManager::GetEventManager()->GetUserEventAction();

    G4int numPrimaries = run->GetNumberOfEvent();
    G4double RnDeabsorptionIN = fpEventAction->getRnDeabsorptionIN();
    G4cout << "Deabsoption of Rn220 from is " << (1-RnDeabsorptionIN / numPrimaries) * 100 << "%, expect 40%." << G4endl;

    G4double PbDeabsorptionIN = fpEventAction->getPbDeabsorptionIN();
    G4cout << "Deabsoption of Pb212 from is " << (1-PbDeabsorptionIN / numPrimaries) * 100 << "%, expect 55%." << G4endl;

    G4double PbLeakage = fpEventAction->getPbLeakage();
    G4double PbNoLeakage = fpEventAction->getPbNoLeakage();
    G4cout << "Leakage of Pb212 from is " << PbLeakage / (PbLeakage + PbNoLeakage) * 100 << "%, value depends on tumour size" << G4endl;

    G4cout << "Activity of Radon 224 = " << numPrimaries/fpEventAction->getTotalRaDecayTime()/s << " s-1" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void RunAction::Write(const G4Run *)
{
    CommandLineParser *parser = CommandLineParser::GetParser();
    Command *command(0);
    if ((command = parser->GetCommandIfActive("-out")) == 0)
        return;
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

    analysisManager->Write();
    analysisManager->CloseFile();
    analysisManager->Clear();
    G4cout << "\n----> Histograms are saved" << G4endl;
}

void RunAction::saveDose(G4double inEdep, G4double time, G4int cp)
{
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

    G4double mass;
    mass = 997 * (4/3)*3.141593*5e-6*5e-6*5e-6; // sphere of water desnity water 997 kg/m3

    G4double Edep = ((inEdep) / joule);
    analysisManager->FillH1(8+cp, time/60/60, Edep / mass);
}
