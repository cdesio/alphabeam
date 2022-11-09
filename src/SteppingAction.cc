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

    calculateDSB(particleMeanEnergy, step->GetStepLength(), time, cp);
  }
}

void SteppingAction::calculateDSB(G4double KE, G4double stepLength, G4double time, G4int cp)
{
  // function to calculate the number of DSB per step (IRT sync model)
    G4double KEtoSimpleDSB[16][2] = {{0, 0},
                                     {0.03333333333333333, 0.0012910796888576262},
                                     {0.10700000000000003, 0.0022926004077416817},
                                     {0.2993333333333333, 0.0032954991762141157},
                                     {0.5290000000000001, 0.003641997371113926},
                                     {0.7800000000000001, 0.003851114172218088},
                                     {1.0399999999999998, 0.0038328258902288687},
                                     {1.307666666666667, 0.003915121427206071},
                                     {1.836, 0.003502440088975107},
                                     {2.876, 0.0029836664610465645},
                                     {3.9, 0.0023459687884314732},
                                     {4.909999999999999, 0.0020190778664797944},
                                     {5.922000000000001, 0.0017789151617751196},
                                     {6.93, 0.0014588185233866448},
                                     {7.949333333333335, 0.0013897465179284274},
                                     {8.950000000000001, 0.00115589606962544}};

    G4double KEtoComplexDSB[16][2] = {{0, 0},
                                      {0.03333333333333333, 0.000235557063087997},
                                      {0.10700000000000003, 0.0008611800649338721},
                                      {0.2993333333333333, 0.0023226680399363164},
                                      {0.5290000000000001, 0.003405127911163989},
                                      {0.7800000000000001, 0.003406643154938256},
                                      {1.0399999999999998, 0.0030608295915904894},
                                      {1.307666666666667, 0.0027169840830291066},
                                      {1.836, 0.0019813627992904983},
                                      {2.876, 0.0012227892599067226},
                                      {3.9, 0.00065444523084271},
                                      {4.909999999999999, 0.0003764882183867844},
                                      {5.922000000000001, 0.00027534692121027477},
                                      {6.93, 0.00017964076149848385},
                                      {7.949333333333335, 0.00016706434762929373},
                                      {8.950000000000001, 0.0001225325097701195}};

  // KE values are the same for complex and simple, so can use the same index
  G4int i = 0;
  while (KE > KEtoSimpleDSB[i][0])
  {
    ++i;
  }

  // index of points either side of E
  G4int minBound{i - 1};
  G4int maxBound{i};

  // linearly interpolate DSB to get value for KE
  G4double simpleDSB = stepLength / nm * KEtoSimpleDSB[minBound][1] + ((KE - KEtoSimpleDSB[minBound][0]) / (KEtoSimpleDSB[maxBound][0] - KEtoSimpleDSB[minBound][0])) * abs(KEtoSimpleDSB[maxBound][1] - KEtoSimpleDSB[minBound][1]);

  G4double complexDSB = stepLength / nm * KEtoComplexDSB[minBound][1] + ((KE - KEtoComplexDSB[minBound][0]) / (KEtoComplexDSB[maxBound][0] - KEtoComplexDSB[minBound][0])) * abs(KEtoComplexDSB[maxBound][1] - KEtoComplexDSB[minBound][1]);

  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

  analysisManager->FillH1(cp, time / 60 / 60, simpleDSB);
  analysisManager->FillH1(cp+10, time / 60 / 60, complexDSB);
  return;
}