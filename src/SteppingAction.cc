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
  // function to calculate the number of DSB per step
  G4double KEtoSimpleDSB[17][2] = { {0, 0},
                                    {0.02666666666666667, 0.0010616628338083911},
                                    {0.03, 0.0011516155416274603},
                                    {0.10666666666666669, 0.0021639360306445173},
                                    {0.3, 0.0030253018618037614},
                                    {0.53, 0.003320716114160411},
                                    {0.7799999999999999, 0.0035208842090507153},
                                    {1.04, 0.0036304412586854224},
                                    {1.3083333333333333, 0.0034707972087794693},
                                    {1.84, 0.0031396687384708964},
                                    {2.8800000000000003, 0.0025469724755125107},
                                    {3.9, 0.002111859459586764},
                                    {4.91, 0.0017862752786142674},
                                    {5.919999999999999, 0.0015061988021849615},
                                    {6.93, 0.001287185641075164},
                                    {7.94, 0.001143362492762962},
                                    {8.94, 0.0010834915009852147} };

  G4double KEtoComplexDSB[17][2] = { {0, 0},
                                     {0.02666666666666667, 0.00010645727843204637},
                                     {0.03, 0.00018429546939030563},
                                     {0.10666666666666669, 0.000697391351508226},
                                     {0.3, 0.0017973074568268493},
                                     {0.53, 0.0024203714698786106},
                                     {0.7799999999999999, 0.002655930185781365},
                                     {1.04, 0.0024284360357074816},
                                     {1.3083333333333333, 0.0019913150796892188},
                                     {1.84, 0.0014109909236008616},
                                     {2.8800000000000003, 0.0007296379316051993},
                                     {3.9, 0.0004264609871747214},
                                     {4.91, 0.0003112096010395113},
                                     {5.919999999999999 ,0.00021332358389700178},
                                     {6.93, 0.00014876555806773325},
                                     {7.94, 9.058893448935311e-05},
                                     {8.94, 9.041053469685657e-05} };

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