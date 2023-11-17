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
//
//
#include "SteppingAction.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include <map>
#include "globals.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "CommandLineParser.hh"
#include "EventAction.hh"
#include "G4Ions.hh"

using namespace G4DNAPARSER;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::SteppingAction()
    : G4UserSteppingAction(), fpEventAction(0)
{
  fpEventAction = (EventAction *)G4EventManager::GetEventManager()->GetUserEventAction();
  fRunAction = (RunAction *)(G4RunManager::GetRunManager()->GetUserRunAction());

  CommandLineParser *parser = CommandLineParser::GetParser();
  Command *command(0);
  G4String fileName{"PSfile.bin"};
  if ((command = parser->GetCommandIfActive("-out")) != 0)

    if (command->GetOption().empty() == false)
    {
      fileName = command->GetOption() + ".bin";
    }

  PSfile.open(fileName, std::ios::out | std::ios::binary);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::~SteppingAction()
{
  PSfile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void SteppingAction::UserSteppingAction(const G4Step *step)
{
  // G4cout << "start of step" << G4endl;

  if (step->GetTrack()->GetParticleDefinition()->GetParticleName() == "anti_nu_e") // not anti neutrinos
    return;
  G4double dE = step->GetTotalEnergyDeposit();
  G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();

  // if (step->GetPreStepPoint() != nullptr)
  // {
    // G4cout << particleName << " " << step->GetTrack()->GetTrackID() << " " << step->GetPreStepPoint()->GetPhysicalVolume()->GetName() << " " << step->GetPostStepPoint()->GetPhysicalVolume()->GetName() << G4endl;
    // G4cout << particleName << " " << step->GetPreStepPoint()->GetKineticEnergy() << " " << step->GetPostStepPoint()->GetPhysicalVolume()->GetName() << " " << step->GetPostStepPoint()->GetPosition()<< G4endl;
  // }

  if (step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "water")
  {
    if (step->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "world")
    {
      step->GetTrack()->SetTrackStatus(fStopAndKill);

      return;
    }
  }

  CommandLineParser *parser = CommandLineParser::GetParser();
  Command *command(0);

  if ((command = parser->GetCommandIfActive("-out")) == 0)
    return;

  if (step->GetPreStepPoint() == nullptr) // is this right????
    return;

  G4String volumeNamePre = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();

  if (volumeNamePre == "water") // particle from water volume entering the cell - save details in PS file
  {
    if (step->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "voxel") // last step before entering cell
    {
      
      G4TouchableHandle theTouchable = step->GetPostStepPoint()->GetTouchableHandle();
      G4ThreeVector worldPos = step->GetPostStepPoint()->GetPosition();
      G4ThreeVector localPos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);

      G4ThreeVector worldMomentum = step->GetPostStepPoint()->GetMomentumDirection();
      //G4ThreeVector localMomentum = (*(theTouchable->GetHistory()->GetTopVolume()->GetRotation()))*worldMomentum; // rotate momentum direction by volume rotation

      G4double time = step->GetPostStepPoint()->GetGlobalTime();

      auto particleEnergy = step->GetPostStepPoint()->GetKineticEnergy();
      G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

      G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();
      G4int copyNo = step->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo();
      G4double steplength = step->GetStepLength();

      G4float particleID{0};
      if (G4StrUtil::contains(particleName,"e-"))
        particleID = 1;
      else if (G4StrUtil::contains(particleName,"gamma"))
        particleID = 2;
      else if (G4StrUtil::contains(particleName,"alpha"))
        particleID = 3;
      else
      {
        G4cout << particleName << " outside  not saved" << G4endl;
        return;
      }

      //G4cout << "particle " << particleID << " in layer " << copyNo << G4endl;
      float output[12];
      output[0] = localPos.x() / mm;
      output[1] = localPos.y() / mm;
      output[2] = localPos.z() / mm;
      output[3] = worldMomentum.x();
      output[4] = worldMomentum.y();
      output[5] = worldMomentum.z();
      output[6] = particleEnergy / MeV;
      output[7] = eventID;
      output[8] = particleID;
      output[9] = copyNo;
      output[10] = time / s;
      output[11] = ((const G4Ions*)( step->GetTrack()->GetParticleDefinition()))->GetExcitationEnergy();

      PSfile.write((char *)&output, sizeof(output));

      G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

      analysisManager->FillNtupleIColumn(1, 0, eventID);
      analysisManager->FillNtupleDColumn(1, 1, particleEnergy / MeV);
      analysisManager->FillNtupleDColumn(1, 2, dE / MeV);
      analysisManager->FillNtupleIColumn(1, 3, particleID);
      analysisManager->FillNtupleIColumn(1, 4, copyNo);
      analysisManager->FillNtupleDColumn(1, 5, worldPos.x() / nanometer);
      analysisManager->FillNtupleDColumn(1, 6, worldPos.y() / nanometer);
      analysisManager->FillNtupleDColumn(1, 7, worldPos.z() / nanometer);
      analysisManager->FillNtupleDColumn(1, 8, steplength);
      analysisManager->AddNtupleRow(1);

    
    }
  }

  if ((volumeNamePre == "voxel") && (step->IsFirstStepInVolume())) // save particles created in the cell or nucleus
  {

    if (step->GetPreStepPoint()->GetProcessDefinedStep() != nullptr)
      return; // if prestep process is nullptr this is the first step of particle created by interaction in the cell - only save those created by processes in cell

    if (step->GetTrack()->GetCreatorProcess()->GetProcessName() != "RadioactiveDecay")
      return; // only save products of radioactive decay other products are from parents which are saved on entering the cell and will be tracked in DNA simulation.


    G4TouchableHandle theTouchable = step->GetPreStepPoint()->GetTouchableHandle();
    G4ThreeVector worldPos = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector localPos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);

    G4ThreeVector worldMomentum = step->GetPreStepPoint()->GetMomentumDirection();
    G4ThreeVector localMomentum = (*(theTouchable->GetHistory()->GetTopVolume()->GetRotation()))*worldMomentum;

    G4double time = step->GetPreStepPoint()->GetGlobalTime();

    auto particleEnergy = step->GetPreStepPoint()->GetKineticEnergy();
    G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

    G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();
    G4int copyNo = step->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo();
    G4double steplength = step->GetStepLength();

    G4float particleID{0};
    if (G4StrUtil::contains(particleName,"e-"))
      particleID = 1;
    else if (G4StrUtil::contains(particleName,"gamma"))
      particleID = 2;
    else if (G4StrUtil::contains(particleName,"alpha"))
      particleID = 3;
    else
    {
      G4cout << particleName << " inside  not saved" << G4endl;
      return;
    }

    float output[12];
    output[0] = localPos.x() / mm;
    output[1] = localPos.y() / mm;
    output[2] = localPos.z() / mm;
    output[3] = localMomentum.x();
    output[4] = localMomentum.y();
    output[5] = localMomentum.z();
    output[6] = particleEnergy / MeV;
    output[7] = eventID;
    output[8] = particleID;
    output[9] = copyNo;
    output[10] = time / s;
    output[11] = ((const G4Ions*)( step->GetTrack()->GetParticleDefinition()))->GetExcitationEnergy();

    PSfile.write((char *)&output, sizeof(output));

    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

    analysisManager->FillNtupleDColumn(1, 1, particleEnergy / MeV);
    analysisManager->FillNtupleDColumn(1, 2, dE / MeV);
    analysisManager->FillNtupleIColumn(1, 3, particleID);
    analysisManager->FillNtupleIColumn(1, 4, copyNo);
    analysisManager->FillNtupleDColumn(1, 5, worldPos.x() / nanometer);
    analysisManager->FillNtupleDColumn(1, 6, worldPos.y() / nanometer);
    analysisManager->FillNtupleDColumn(1, 7, worldPos.z() / nanometer);
    analysisManager->FillNtupleDColumn(1, 8, steplength);
    analysisManager->AddNtupleRow(1);
  }

}

