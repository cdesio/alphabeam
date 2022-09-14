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
  if ((particleName.contains("Pb212")) && (step->GetPreStepPoint()->GetKineticEnergy() == 0))
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
  G4double edep = step->GetTotalEnergyDeposit();

  // Save Edep from all particles to calculate dose
  if ((edep > 0) && (step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "nucleus"))
  {
    G4double time = (step->GetPreStepPoint()->GetGlobalTime() + step->GetPostStepPoint()->GetGlobalTime()) / 2 / s;
    G4int cp = step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo();
    fRunAction->saveDose(edep, time, cp);
  }
  // Save SB of all alpha particles
  if ((particleName == "alpha") && (step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "nucleus"))
  {
    G4double time = (step->GetPreStepPoint()->GetGlobalTime() + step->GetPostStepPoint()->GetGlobalTime()) / 2 / s;

    G4double particleMeanEnergy = (step->GetPreStepPoint()->GetKineticEnergy() + step->GetPostStepPoint()->GetKineticEnergy()) / 2;
    G4int cp = step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo();
    G4cout << "here" << G4endl;

    calculateDSB(particleMeanEnergy, step->GetStepLength(), time, cp);
  }
}

void SteppingAction::calculateDSB(G4double KE, G4double stepLength, G4double time, G4int cp)
{
  // function to calculate the number of DSB per step
  G4double KEtoDSB[17][2] = {{0, 0},
                             {0.02666666666666667, 0.0011688818838326493},
                             {0.03, 0.0013335344030287663},
                             {0.10666666666666667, 0.002861315961221019},
                             {0.3, 0.004822935752598016},
                             {0.5300000000000001, 0.005741146233020532},
                             {0.7800000000000001, 0.006176385740656605},
                             {1.04, 0.006058728319567002},
                             {1.3083333333333333, 0.005461756820525401},
                             {1.84, 0.0045488435091638356},
                             {2.8799999999999994, 0.003275130975507074},
                             {3.9, 0.0025388106759176884},
                             {4.91, 0.002096712792670623},
                             {5.920000000000001, 0.0017191600802838665},
                             {6.93, 0.0014353052334767301},
                             {7.94, 0.0012337176891147356},
                             {8.94, 0.0011739286803916908}};

  G4int i = 0;
  while (KE > KEtoDSB[i][0])
  {
    ++i;
  }

  // index of points either side of E
  G4int minBound{i - 1};
  G4int maxBound{i};

  // linearly interpolate DSB to get value for KE
  G4double DSB = stepLength / nm * KEtoDSB[minBound][1] + ((KE - KEtoDSB[minBound][0]) / (KEtoDSB[maxBound][0] - KEtoDSB[minBound][0])) * abs(KEtoDSB[maxBound][1] - KEtoDSB[minBound][1]);

  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

  analysisManager->FillH1(cp, time / 60 / 60, DSB);
  return;
}