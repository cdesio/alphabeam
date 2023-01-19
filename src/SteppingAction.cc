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

  if (step->GetPreStepPoint() == nullptr) // is this right????
    return;

  G4String volumeNamePre = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();

  if (volumeNamePre == "water") // particle from water volume entering the cell - save details in PS file
  {
    if (step->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "cell") // last step before entering cell
    {

      G4TouchableHandle theTouchable = step->GetPostStepPoint()->GetTouchableHandle();
      G4ThreeVector worldPos = step->GetPostStepPoint()->GetPosition();
      G4ThreeVector localPos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);

      G4ThreeVector worldMomentum = step->GetPostStepPoint()->GetMomentumDirection();
      G4ThreeVector localMomentum = (*(theTouchable->GetHistory()->GetTopVolume()->GetRotation()))*worldMomentum; // rotate momentum direction by volume rotation

      G4double time = step->GetPostStepPoint()->GetGlobalTime();

      auto particleEnergy = step->GetPostStepPoint()->GetKineticEnergy();
      G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

      G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();

      G4float particleID{0};
      if (particleName == "e-")
        particleID = 1;
      else if (particleName == "gamma")
        particleID = 2;
      else if (particleName == "alpha")
        particleID = 3;
      else
      {
        G4cout << particleName << " outside  not saved" << G4endl;
        return;
      }
      float output[11];
      output[0] = localPos.x() / mm;
      output[1] = localPos.y() / mm;
      output[2] = localPos.z() / mm;
      output[3] = localMomentum.x();
      output[4] = localMomentum.y();
      output[5] = localMomentum.z();
      output[6] = particleEnergy / MeV;
      output[7] = eventID;
      output[8] = particleID;
      output[9] = step->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo();
      output[10] = time / s;

      PSfile.write((char *)&output, sizeof(output));
    }
  }

  if ((volumeNamePre == "cell") && (step->IsFirstStepInVolume())) // save particles created in the cell or nucleus
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

    G4float particleID{0};
    if (particleName == "e-")
      particleID = 1;
    else if (particleName == "gamma")
      particleID = 2;
    else if (particleName == "alpha")
      particleID = 3;
    else
    {
      G4cout << particleName << " inside  not saved" << G4endl;
      return;
    }

    float output[11];
    output[0] = localPos.x() / mm;
    output[1] = localPos.y() / mm;
    output[2] = localPos.z() / mm;
    output[3] = localMomentum.x();
    output[4] = localMomentum.y();
    output[5] = localMomentum.z();
    output[6] = particleEnergy / MeV;
    output[7] = eventID;
    output[8] = particleID;
    output[9] = step->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo();
    output[10] = time / s;

    PSfile.write((char *)&output, sizeof(output));

  }

  // G4double edep = step->GetTotalEnergyDeposit();

  // // Save Edep from all particles to calculate dose
  // if ((edep > 0) && (step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "nucleus"))
  // {
  //   G4double time = (step->GetPreStepPoint()->GetGlobalTime() + step->GetPostStepPoint()->GetGlobalTime()) / 2 / s;
  //   G4int cp = step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo();
  //   fRunAction->saveDose(edep, time, cp);
  // }
  // // Save SB of all alpha particles
  // if ((particleName == "alpha") && (step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "nucleus"))
  // {
  //   G4double time = (step->GetPreStepPoint()->GetGlobalTime() + step->GetPostStepPoint()->GetGlobalTime()) / 2 / s;

  //   G4double particleMeanEnergy = (step->GetPreStepPoint()->GetKineticEnergy() + step->GetPostStepPoint()->GetKineticEnergy()) / 2;
  //   G4int cp = step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo();

  //   calculateDSB(particleMeanEnergy, step->GetStepLength(), time, cp);
  // }
  // G4cout << "end of step" << G4endl;
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
  analysisManager->FillH1(cp + 10, time / 60 / 60, complexDSB);
  return;
}