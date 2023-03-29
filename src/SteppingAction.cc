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

  G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();

  if (step->GetTrack()->GetCreatorProcess() != nullptr)
  {
  G4cout << particleName << " Track ID = " << step->GetTrack()->GetTrackID() << " creator process = " << step->GetTrack()->GetCreatorProcess()->GetProcessName() << " KE = " << step->GetPreStepPoint()->GetKineticEnergy() << " parent = " << step->GetTrack()->GetParentID() << "position = " << step->GetPreStepPoint()->GetPosition() <<  " " << step->GetPreStepPoint()->GetPhysicalVolume()->GetName() <<"-"<<step->GetPostStepPoint()->GetPhysicalVolume()->GetName() <<G4endl;
  // G4cout << particleName << " " << step->GetPreStepPoint()->GetKineticEnergy() << " " << step->GetPostStepPoint()->GetPhysicalVolume()->GetName() << " " << step->GetPostStepPoint()->GetPosition()<< G4endl;
  }

  if (step->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "world")
  {
    step->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);

    return;
  }

  if ((particleName == "Ra224") && (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "RadioactiveDecay"))
  {
    fpEventAction->addDecayTimeRa(step->GetPostStepPoint()->GetGlobalTime());
  }

  if ((particleName == "Rn220") && (step->GetPreStepPoint()->GetKineticEnergy() == 0))
  {
    // deabsorption from the source through recoil only
    if (step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "seed")
    {
      fpEventAction->addRnDeabsorptionIN();
    }
  }
  if ((G4StrUtil::contains(particleName, "Pb212")) && (step->GetPreStepPoint()->GetKineticEnergy() == 0))
  {
    if (step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "seed")
    {
      fpEventAction->addPbDeabsorptionIN();
    }
  }

  G4double particleEnergy = step->GetPreStepPoint()->GetKineticEnergy();
  if ((particleName == "Pb212") && (particleEnergy == 0))
  {
    if (step->GetPreStepPoint()->GetPhysicalVolume()->GetName() != "seed")
    {
      if (G4UniformRand() >= 0.5)
      {
        step->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
        fpEventAction->addPbLeakage();
        return;
      }
      else
      {
        fpEventAction->addPbNoLeakage();
      }
    }
  }
  if ((particleName == "Pb208") && (particleEnergy == 0))
  {
    step->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
    return;
  }

  CommandLineParser *parser = CommandLineParser::GetParser();
  Command *command(0);

  if ((command = parser->GetCommandIfActive("-out")) == 0)
    return;

  if (step->GetPreStepPoint() == nullptr)
    return;

  // G4cout << step->GetTrack()->GetParticleDefinition()->GetParticleName() << " from r = " << std::pow(std::pow(step->GetPreStepPoint()->GetPosition().x(), 2) + std::pow(step->GetPreStepPoint()->GetPosition().y(), 2), 0.5) / mm << " to r = " << std::pow(std::pow(step->GetPostStepPoint()->GetPosition().x(), 2) + std::pow(step->GetPostStepPoint()->GetPosition().y(), 2), 0.5) / mm << G4endl;

  G4String volumeNamePre = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
  // if (volumeNamePre == "water") // particle entering cylindrical bands 
  // {
  //   if (step->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "cell")
  //   {
  //     if (step->GetPreStepPoint()->GetKineticEnergy() > 0)
  //     {
  //       savePoint(step, step->GetPostStepPoint());
  //     }
  //   }
  // }
  // if ((volumeNamePre == "cell") && (step->IsFirstStepInVolume())) // save particles created in the cell or nucleus
  // {
  //   if (step->GetPostStepPoint()->GetProcessDefinedStep() == nullptr)
  //     // if prestep process is nullptr this is the first step of particle created by interaction in the cell - only save those created by processes in cell

  //   if (step->GetTrack()->GetCreatorProcess()->GetProcessName() == "RadioactiveDecay")
  //     // only save products of radioactive decay other products are from parents which are saved on entering the cell and will be tracked in DNA simulation.
  //     if (step->GetPostStepPoint()->GetKineticEnergy() > 0)
  //     {
  //       savePoint(step, step->GetPreStepPoint());
  //     }
  //   }
  
  // save all steps in rings
  if ((volumeNamePre == "cell"))
{
      if (step->GetPreStepPoint()->GetKineticEnergy() > 0)
      {
        if ((step->IsFirstStepInVolume())&&(step->IsLastStepInVolume())) //only one step in volume
        {
            savePoint(step, step->GetPreStepPoint(), step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo());
            savePoint(step, step->GetPostStepPoint(), step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo());
        }
        else if (step->IsLastStepInVolume())
            savePoint(step, step->GetPostStepPoint(), step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo());
        else
            savePoint(step, step->GetPreStepPoint(),step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo());
      }
  }
}

void SteppingAction::savePoint(const G4Step *step, const G4StepPoint * point, const int copy)
{
  // G4TouchableHandle theTouchable = step->GetPostStepPoint()->GetTouchableHandle();
  G4ThreeVector worldPos = point->GetPosition();
  // G4ThreeVector localPos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);

  G4ThreeVector worldMomentum = point->GetMomentumDirection();
  // G4ThreeVector localMomentum = (*(theTouchable->GetHistory()->GetTopVolume()->GetRotation())) * worldMomentum; // rotate momentum direction by volume rotation

  G4double time = point->GetGlobalTime();

  auto particleEnergy = point->GetKineticEnergy();
  G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

  G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();

  G4float particleID{0};
  if (G4StrUtil::contains(particleName, "e-"))
    particleID = 1;
  else if (G4StrUtil::contains(particleName, "gamma"))
    particleID = 2;
  else if (G4StrUtil::contains(particleName, "alpha"))
    particleID = 3;
  else if (G4StrUtil::contains(particleName, "Rn220"))
    particleID = 4;
  else if (G4StrUtil::contains(particleName, "Po216"))
    particleID = 5;
  else if (G4StrUtil::contains(particleName, "Pb212"))
    particleID = 6;
  else if (G4StrUtil::contains(particleName, "Bi212"))
    particleID = 7;
  else if (G4StrUtil::contains(particleName, "Tl208"))
    particleID = 8;
  else if (G4StrUtil::contains(particleName, "Po212"))
    particleID = 9;
  else if (G4StrUtil::contains(particleName, "Pb208"))
    particleID = 10;
  else
  {
    G4cout << particleName << "  not saved" << G4endl;
    return;
  }
  double output[14];
  output[0] = worldPos.x() / mm;
  output[1] = worldPos.y() / mm;
  output[2] = worldPos.z() / mm;
  output[3] = worldMomentum.x();
  output[4] = worldMomentum.y();
  output[5] = worldMomentum.z();
  output[6] = particleEnergy / MeV;
  output[7] = eventID;
  output[8] = particleID;
  output[9] = copy;
  output[10] = time / s;
  output[11] = ((const G4Ions *)(step->GetTrack()->GetParticleDefinition()))->GetExcitationEnergy();
  output[12] = step->GetTrack()->GetTrackID();
  output[13] = step->GetTrack()->GetCreatorProcess()->GetProcessName() == "RadioactiveDecay" ? 1 : 0; // 1 if from radioactive decay

  PSfile.write((char *)&output, sizeof(output));

  G4cout << particleName << " saved at = " << std::pow(std::pow(worldPos.x(), 2) + std::pow(worldPos.y(), 2), 0.5) / mm << G4endl;
}
