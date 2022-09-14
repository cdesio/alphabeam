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

using namespace G4DNAPARSER;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::SteppingAction()
    : G4UserSteppingAction(), fpEventAction(0)
{
  fpEventAction = (EventAction *)G4EventManager::GetEventManager()->GetUserEventAction();
  fRunAction = (RunAction *)(G4RunManager::GetRunManager()->GetUserRunAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::~SteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void SteppingAction::UserSteppingAction(const G4Step *step)
{
  G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();

  if ((particleName == "Ra224") && (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="RadioactiveDecay"))
  {
  fpEventAction -> addDecayTimeRa(step->GetPostStepPoint()->GetGlobalTime()) ;
  }

  if ((particleName == "Rn220") && (step->GetPreStepPoint()->GetKineticEnergy()==0))
  {
    // deabsorption from the source through recoil only
    if (step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "world")
    {
      fpEventAction->addRnDeabsorptionOUT();
    }
  }
    if ((particleName.contains("Pb212"))&& (step->GetPreStepPoint()->GetKineticEnergy()==0))
  {
    if (step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "world")
    {
      fpEventAction->addPbDeabsorptionOUT();
    }
  }

  G4double particleEnergy = step->GetPreStepPoint()->GetKineticEnergy();
  if ((particleName == "Pb212") && (particleEnergy == 0))
  {
    if (step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "world")
    {
      if (G4UniformRand() >= 0.5)
      {
        step->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
        fpEventAction->addPbLeakage();
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
  }

  CommandLineParser *parser = CommandLineParser::GetParser();
  Command *command(0);

  if ((command = parser->GetCommandIfActive("-out")) == 0)
    return;

  if (particleName.contains("Ra224"))
  {
    fpEventAction->addToMap(step->GetTrack()->GetTrackID(), "Ra224");
  }
  else if (particleName.contains("Rn220"))
  {
    fpEventAction->addToMap(step->GetTrack()->GetTrackID(), "Rn220");
  }
  else if (particleName.contains("Po216"))
  {
    fpEventAction->addToMap(step->GetTrack()->GetTrackID(), "Po216");
  }
  else if (particleName.contains("Pb212"))
  {
    fpEventAction->addToMap(step->GetTrack()->GetTrackID(), "Pb212");
  }
  else if (particleName.contains("Bi212"))
  {
    fpEventAction->addToMap(step->GetTrack()->GetTrackID(), "Bi212");
  }
  else if (particleName.contains("Po212"))
  {
    fpEventAction->addToMap(step->GetTrack()->GetTrackID(), "Po212");
  }
  else if (particleName.contains("Tl208"))
  {
    fpEventAction->addToMap(step->GetTrack()->GetTrackID(), "Tl208");
  }
  else if (particleName.contains("Pb208"))
  {
    fpEventAction->addToMap(step->GetTrack()->GetTrackID(), "Pb208");
  }
  else if (particleName.contains("Bi213"))
  {
    fpEventAction->addToMap(step->GetTrack()->GetTrackID(), "Bi213");
  }
  else if (particleName.contains("At217"))
  {
    fpEventAction->addToMap(step->GetTrack()->GetTrackID(), "At217");
  }
  else
  {
    G4int parentID = step->GetTrack()->GetParentID();

    if (fpEventAction->getFromMap(parentID).contains("Ra224"))
    {
      fpEventAction->addToMap(step->GetTrack()->GetTrackID(), "Ra224");
    }
    else if (fpEventAction->getFromMap(parentID).contains("Rn220"))
    {
      fpEventAction->addToMap(step->GetTrack()->GetTrackID(), "Rn220");
    }
    else if (fpEventAction->getFromMap(parentID).contains("Po216"))
    {
      fpEventAction->addToMap(step->GetTrack()->GetTrackID(), "Po216");
    }
    else if (fpEventAction->getFromMap(parentID).contains("Pb212"))
    {
      fpEventAction->addToMap(step->GetTrack()->GetTrackID(), "Pb212");
    }
    else if (fpEventAction->getFromMap(parentID).contains("Bi212"))
    {
      fpEventAction->addToMap(step->GetTrack()->GetTrackID(), "Bi212");
    }
    else if (fpEventAction->getFromMap(parentID).contains("Po212"))
    {
      fpEventAction->addToMap(step->GetTrack()->GetTrackID(), "Po212");
    }
    else if (fpEventAction->getFromMap(parentID).contains("Tl208"))
    {
      fpEventAction->addToMap(step->GetTrack()->GetTrackID(), "Tl208");
    }
    else if (fpEventAction->getFromMap(parentID).contains("Pb208"))
    {
      fpEventAction->addToMap(step->GetTrack()->GetTrackID(), "Pb208");
    }
  }

  G4double edep = step->GetTotalEnergyDeposit();

  // Save Edep from all particles to calculate dose
  if ((edep > 0) && (step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "world"))
  {
    G4double positionX = (step->GetPreStepPoint()->GetPosition().x() + step->GetPostStepPoint()->GetPosition().x()) / 2 / mm;
    G4double positionY = (step->GetPreStepPoint()->GetPosition().y() + step->GetPostStepPoint()->GetPosition().y()) / 2 / mm;
    G4double positionZ = (step->GetPreStepPoint()->GetPosition().z() + step->GetPostStepPoint()->GetPosition().z()) / 2 / mm;

    fRunAction->saveEdep(edep, positionX, positionY, positionZ, fpEventAction->getFromMap(step->GetTrack()->GetTrackID()));
  }
  // Save KE of all alpha particles
  if ((particleName == "alpha") && (step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "world"))
  {
    G4double positionX = (step->GetPreStepPoint()->GetPosition().x() + step->GetPostStepPoint()->GetPosition().x()) / 2 / mm;
    G4double positionY = (step->GetPreStepPoint()->GetPosition().y() + step->GetPostStepPoint()->GetPosition().y()) / 2 / mm;
    G4double positionZ = (step->GetPreStepPoint()->GetPosition().z() + step->GetPostStepPoint()->GetPosition().z()) / 2 / mm;

    G4double time = (step->GetPreStepPoint()->GetGlobalTime() + step->GetPostStepPoint()->GetGlobalTime()) / 2 / s;

    G4double particleMeanEnergy = (step->GetPreStepPoint()->GetKineticEnergy() + step->GetPostStepPoint()->GetKineticEnergy()) / 2;
    fRunAction->saveKE(particleMeanEnergy, positionX, positionY, positionZ, fpEventAction->getFromMap(step->GetTrack()->GetTrackID()), time);
    // fRunAction->calculateDSB(particleMeanEnergy, step->GetStepLength(), positionX, positionY, positionZ, time);
  }
}
