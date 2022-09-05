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

    analysisManager->CreateH3("Energy3D", "Energy3D", 151, -7.5, 7.5, 151, -7.5, 7.5, 151, -7.5, 7.5);
    analysisManager->CreateH3("NumAlpha", "NumAlpha", 151, -7.5, 7.5, 151, -7.5, 7.5, 151, -7.5, 7.5);
    analysisManager->CreateH3("Dose", "Dose", 151, -7.5, 7.5, 151, -7.5, 7.5, 151, -7.5, 7.5);
    analysisManager->CreateH2("Energy2D", "Energy2D", 150, 0, 7.5, 100, 0, 10);
    analysisManager->CreateH2("Energy2D_Ra224", "Energy2D_Ra224", 150, 0, 7.5, 100, 0, 10);
    analysisManager->CreateH2("Energy2D_Rn220", "Energy2D_Rn220", 150, 0, 7.5, 100, 0, 10);
    analysisManager->CreateH2("Energy2D_Po216", "Energy2D_Po216", 150, 0, 7.5, 100, 0, 10);
    analysisManager->CreateH2("Energy2D_Pb212", "Energy2D_Pb212", 150, 0, 7.5, 100, 0, 10);
    analysisManager->CreateH2("Energy2D_Bi212", "Energy2D_Bi212", 150, 0, 7.5, 100, 0, 10);
    analysisManager->CreateH2("Energy2D_Po212", "Energy2D_Po212", 150, 0, 7.5, 100, 0, 10);
    analysisManager->CreateH2("Energy2D_Tl208", "Energy2D_Tl208", 150, 0, 7.5, 100, 0, 10);
    analysisManager->CreateH2("Energy2D_Pb208", "Energy2D_Pb208", 150, 0, 7.5, 100, 0, 10);
    analysisManager->CreateH1("Dose1D", "Dose1D", 150, 0, 7.5);
    analysisManager->CreateH1("Dose1D_Ra224", "Dose1D_Ra224", 150, 0, 7.5);
    analysisManager->CreateH1("Dose1D_Rn220", "Dose1D_Rn220", 150, 0, 7.5);
    analysisManager->CreateH1("Dose1D_Po216", "Dose1D_Po216", 150, 0, 7.5);
    analysisManager->CreateH1("Dose1D_Pb212", "Dose1D_Pb212", 150, 0, 7.5);
    analysisManager->CreateH1("Dose1D_Bi212", "Dose1D_Bi212", 150, 0, 7.5);
    analysisManager->CreateH1("Dose1D_Po212", "Dose1D_Po212", 150, 0, 7.5);
    analysisManager->CreateH1("Dose1D_Tl208", "Dose1D_Tl208", 150, 0, 7.5);
    analysisManager->CreateH1("Dose1D_Pb208", "Dose1D_Pb208", 150, 0, 7.5);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::EndOfRunAction(const G4Run *run)
{
    Write(run);
    auto fpEventAction = (EventAction *)G4EventManager::GetEventManager()->GetUserEventAction();

    G4double deabsorptionOUT = fpEventAction->getDeabsorptionOUT();
    G4double deabsorptionIN = fpEventAction->getDeabsorptionIN();
    G4cout << "Deabsoption of Radon from is " << deabsorptionOUT / (deabsorptionOUT + deabsorptionIN) * 100 << "%, expect 40%." << G4endl;

    G4double PbLeakage = fpEventAction->getPbLeakage();
    G4double PbNoLeakage = fpEventAction->getPbNoLeakage();
    G4cout << "Leakage of Pb212 from is " << PbLeakage / (PbLeakage + PbNoLeakage) * 100 << "%, value depends on tumour size" << G4endl;
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

void RunAction::saveDose(G4double inDose, G4double x, G4double y, G4double z, G4String parent)
{
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

    G4double mass;
    G4double dose;
    mass = 997 * (0.1e-3) * (0.1e-3) * (0.1e-3); // per cube of water, side length = 0.1mm, desnity water 997 kg/m3

    dose = ((inDose) / joule) / (mass / kg);
    analysisManager->FillH3(2, x, y, z, dose);

    // save dose distribution around source, assume radial symmetry and save radial distance from centre of source
    if ((-3 <= z) && (z <= 3)) // source is 6mm
    {
        analysisManager->FillH1(0, sqrt(x * x + y * y), dose);
        if (parent == "Ra224")
        {
            analysisManager->FillH1(1, sqrt(x * x + y * y), dose);
        }
        else if (parent == "Rn220")
        {
            analysisManager->FillH1(2, sqrt(x * x + y * y), dose);
        }
        else if (parent == "Po216")
        {
            analysisManager->FillH1(3, sqrt(x * x + y * y), dose);
        }
        else if (parent == "Pb212")
        {
            analysisManager->FillH1(4, sqrt(x * x + y * y), dose);
        }
        else if (parent == "Bi212")
        {
            analysisManager->FillH1(5, sqrt(x * x + y * y), dose);
        }
        else if (parent == "Po212")
        {
            analysisManager->FillH1(6, sqrt(x * x + y * y), dose);
        }
        else if (parent == "Tl208")
        {
            analysisManager->FillH1(7, sqrt(x * x + y * y), dose);
        }
        else if (parent == "Pb208")
        {
            analysisManager->FillH1(8, sqrt(x * x + y * y), dose);
        }
    }
}

void RunAction::saveKE(G4double inEnergy, G4double x, G4double y, G4double z, G4String parent)
{

    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

    analysisManager->FillH3(0, x, y, z, inEnergy);
    analysisManager->FillH3(1, x, y, z, 1);

    // save energy distribution around source, assume radial symmetry and save radial distance from centre of source
    if ((-3 <= z) && (z <= 3)) // source is 6mm
    {
        analysisManager->FillH2(0, sqrt(x * x + y * y), inEnergy, 1);
        if (parent == "Ra224")
        {
            analysisManager->FillH2(1, sqrt(x * x + y * y), inEnergy, 1);
        }
        else if (parent == "Rn220")
        {
            analysisManager->FillH2(2, sqrt(x * x + y * y), inEnergy, 1);
        }
        else if (parent == "Po216")
        {
            analysisManager->FillH2(3, sqrt(x * x + y * y), inEnergy, 1);
        }
        else if (parent == "Pb212")
        {
            analysisManager->FillH2(4, sqrt(x * x + y * y), inEnergy, 1);
        }
        else if (parent == "Bi212")
        {
            analysisManager->FillH2(5, sqrt(x * x + y * y), inEnergy, 1);
        }
        else if (parent == "Po212")
        {
            analysisManager->FillH2(6, sqrt(x * x + y * y), inEnergy, 1);
        }
        else if (parent == "Tl208")
        {
            analysisManager->FillH2(7, sqrt(x * x + y * y), inEnergy, 1);
        }
        else if (parent == "Pb208")
        {
            analysisManager->FillH2(8, sqrt(x * x + y * y), inEnergy, 1);
        }
    }
}
