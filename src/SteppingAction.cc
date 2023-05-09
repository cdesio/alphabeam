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
#include <cmath>

using namespace G4DNAPARSER;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::SteppingAction(DetectorConstruction *pDetector)
    : G4UserSteppingAction(), fpEventAction(0), fDetector(pDetector)
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
  // G4cout << "start of step, volume = " << step->GetPreStepPoint()->GetPhysicalVolume()->GetName()<< " step length = " << step->GetStepLength()/nm << G4endl;

  if (step->GetTrack()->GetParticleDefinition()->GetParticleName() == "anti_nu_e") // not anti neutrinos
    return;

  G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();

  // if (step->GetTrack()->GetCreatorProcess() != nullptr)
  // {
  // G4cout << particleName << " Track ID = " << step->GetTrack()->GetTrackID() << " creator process = " << step->GetTrack()->GetCreatorProcess()->GetProcessName() << " KE = " << step->GetPreStepPoint()->GetKineticEnergy() << " parent = " << step->GetTrack()->GetParentID() << "position = " << step->GetPreStepPoint()->GetPosition() <<  " " << step->GetPreStepPoint()->GetPhysicalVolume()->GetName() <<"-"<<step->GetPostStepPoint()->GetPhysicalVolume()->GetName() <<G4endl;
  // // G4cout << particleName << " " << step->GetPreStepPoint()->GetKineticEnergy() << " " << step->GetPostStepPoint()->GetPhysicalVolume()->GetName() << " " << step->GetPostStepPoint()->GetPosition()<< G4endl;
  // }

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

  G4String volumeNamePre = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();

  // save all steps entering rings
  if ((volumeNamePre == "water") && (step->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "cell"))
  {
    if (step->GetPreStepPoint()->GetKineticEnergy() > 0)
    {
      // G4cout << "saving particle entering ring " << G4endl;

      G4ThreeVector worldPos = step->GetPostStepPoint()->GetPosition();
      G4double newX = (G4UniformRand() * .00017 * 2) - .00017;
      G4double newZ = (G4UniformRand() * .00017 * 2) - .00017;

      G4double radius = std::pow(worldPos.x() * worldPos.x() + worldPos.y() * worldPos.y(), 0.5);

      G4double newY = radius - fDetector->R[step->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo()];

      // pick position in box frame
      G4ThreeVector newPos = G4ThreeVector(newX, newY, newZ);
      G4ThreeVector newMomentum = transformDirection(step->GetPostStepPoint()->GetPosition(), step->GetPostStepPoint()->GetMomentumDirection());
      savePoint(step->GetTrack(), newPos, newMomentum, step->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo(), step->GetPostStepPoint()->GetKineticEnergy(), step->GetPostStepPoint()->GetGlobalTime());
    }
  }
  // save decay in box
  else if ((volumeNamePre == "cell") && (step->IsFirstStepInVolume())) // save particles created in the cell or nucleus if from radioactive decay as not simulated in RBE
  {
    if (step->GetPreStepPoint()->GetProcessDefinedStep() == nullptr)
    // if prestep process is nullptr this is the first step of particle created by interaction in the cell - only save those created by processes in cell
    {
      if (step->GetTrack()->GetCreatorProcess()->GetProcessName() == "RadioactiveDecay")
      {
        // only save products of radioactive decay other products are from parents which are saved on entering the cell and will be tracked in DNA simulation.
        G4cout << "saved Rdecay" << G4endl;

        G4ThreeVector worldPos = step->GetPreStepPoint()->GetPosition();
        G4double newX = (G4UniformRand() * .00017 * 2) - .00017;
        G4double newZ = (G4UniformRand() * .00017 * 2) - .00017;

        G4double radius = std::pow(worldPos.x() * worldPos.x() + worldPos.y() * worldPos.x(), 0.5);
        G4double newY = radius - fDetector->R[step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo()];

        // pick position in box frame
        G4ThreeVector newPos = G4ThreeVector(newX, newY, newZ);
        G4ThreeVector newMomentum = transformDirection(step->GetPreStepPoint()->GetPosition(), step->GetPreStepPoint()->GetMomentumDirection());

        savePoint(step->GetTrack(), newPos, newMomentum, step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo(), step->GetPreStepPoint()->GetKineticEnergy(), step->GetPreStepPoint()->GetGlobalTime());

        G4cout << step->GetTrack()->GetParticleDefinition()->GetParticleName() << " saved at " << newPos << G4endl;
      }
    }
  }
  else if (volumeNamePre == "cell") // is a step in the cell but not first, check if crosses virtual box boundary, if it does it is saved as if entering the box from the side faces
  {
    if (step->GetPreStepPoint()->GetKineticEnergy() > 0)
    {
      G4ThreeVector entryPosition = fpEventAction->particlePos[step->GetTrack()->GetTrackID()]; // look up position in box frame from last step

      G4ThreeVector deltaWorld = step->GetPostStepPoint()->GetPosition() - step->GetPreStepPoint()->GetPosition(); // change in position in world frame

      G4ThreeVector boxMomentum = transformDirection(step->GetPostStepPoint()->GetPosition(), step->GetPostStepPoint()->GetMomentumDirection()); // particle momentum in box frame

      G4ThreeVector delta = deltaWorld.mag() * boxMomentum; // change in position in box frame

      G4ThreeVector postStepBox = entryPosition + (fpEventAction->particleDist)[step->GetTrack()->GetTrackID()] + delta; // post step position in box frame

      if ((std::abs(postStepBox.x()) >= 0.00017) || (std::abs(postStepBox.z()) >= 0.00017)) // if >=0.00017 has crossed the boundary
      {
        // save particle, new position and distance saved

        G4ThreeVector preStepBox = entryPosition + (fpEventAction->particleDist)[step->GetTrack()->GetTrackID()]; // pre step point position in box frame

        // does step exit box in x and z?
        G4double tXneg = (-.00017 - preStepBox.x()) / boxMomentum.x();
        G4double tXpos = (.00017 - preStepBox.x()) / boxMomentum.x();

        G4double tYneg = (-.00017 - preStepBox.y()) / boxMomentum.y();
        G4double tYpos = (.00017 - preStepBox.y()) / boxMomentum.y();

        G4double tZneg = (-.00017 - preStepBox.z()) / boxMomentum.z();
        G4double tZpos = (.00017 - preStepBox.z()) / boxMomentum.z();

        tXneg = tXneg > 0 ? tXneg : 1e9;
        tXpos = tXpos > 0 ? tXpos : 1e9;
        tYneg = tYneg > 0 ? tYneg : 1e9;
        tYpos = tYpos > 0 ? tYpos : 1e9;
        tZneg = tZneg > 0 ? tZneg : 1e9;
        tZpos = tZpos > 0 ? tZpos : 1e9;

        G4double distanceToExit = std::min({tXneg, tXpos, tZneg, tZpos}); // shortest distance travelled to cross box surface

        if ((tYneg <= distanceToExit) || (tYpos <= distanceToExit)) // exit y
          return;

        G4double stepDistance = step->GetStepLength();

        G4ThreeVector newPos = preStepBox + (distanceToExit * boxMomentum);

        // check which side of the box was crossed and change sign as particle is entering adjacent box
        if (newPos.x() == 0.00017)
          newPos.setX(-0.00017);
        else if (newPos.x() == -0.00017)
          newPos.setX(+0.00017);

        if (newPos.z() == 0.00017)
          newPos.setZ(-0.00017);
        else if (newPos.z() == -0.00017)
          newPos.setZ(+0.00017);

        G4double percentageOfStep = distanceToExit / stepDistance;

        // calculate KE at point where crossing occurs
        G4double newKE = step->GetPreStepPoint()->GetKineticEnergy() - (step->GetPreStepPoint()->GetKineticEnergy() - step->GetPostStepPoint()->GetKineticEnergy()) * percentageOfStep;
        // calculate time at point where crossing occurs
        G4double newTime = step->GetPreStepPoint()->GetGlobalTime() + (step->GetDeltaTime() * percentageOfStep);

        savePoint(step->GetTrack(), newPos, step->GetPostStepPoint()->GetMomentumDirection(), step->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo(), newKE, newTime);
      }
      else
      {
        //  if not update distance travelled in box
        G4ThreeVector previousDelta = (fpEventAction->particleDist)[step->GetTrack()->GetTrackID()];

        // G4cout << "add delta " << delta << G4endl;
        fpEventAction->particleDist.erase(step->GetTrack()->GetTrackID());
        fpEventAction->particleDist.insert(std::pair<int, G4ThreeVector>(step->GetTrack()->GetTrackID(), previousDelta + delta));
      }
    }
  }
  // }
}

void SteppingAction::savePoint(const G4Track *track, G4ThreeVector newPos, G4ThreeVector boxMomentum, const int copy, G4double particleEnergy, G4double time)
{
  // save particle to phase space file in box reference frame

  fpEventAction->particlePos.erase(track->GetTrackID()); // erase current saved box entry position for this track

  fpEventAction->particlePos.insert(std::pair<int, G4ThreeVector>(track->GetTrackID(), newPos)); // add current box entry position for this track

  fpEventAction->particleDist.erase(track->GetTrackID());                                                  // erase distance travelled for this track in the box
  fpEventAction->particleDist.insert(std::pair<int, G4ThreeVector>(track->GetTrackID(), G4ThreeVector())); // made initial distance travelled zero

  G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

  G4String particleName = track->GetParticleDefinition()->GetParticleName();

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
  double output[12];
  output[0] = newPos.x() / mm;
  output[1] = newPos.y() / mm;
  output[2] = newPos.z() / mm;
  output[3] = boxMomentum.x();
  output[4] = boxMomentum.y();
  output[5] = boxMomentum.z();
  output[6] = particleEnergy / MeV;
  output[7] = eventID;
  output[8] = particleID;
  output[9] = copy;
  output[10] = time / s;
  output[11] = ((const G4Ions *)(track->GetParticleDefinition()))->GetExcitationEnergy();
  // output[12] = track->GetTrackID();
  // output[13] = track->GetCreatorProcess()->GetProcessName() == "RadioactiveDecay" ? 1 : 0; // 1 if from radioactive decay

  PSfile.write((char *)&output, sizeof(output));

  // G4cout << particleName << " saved at = " << newPos / mm << " with KE = " << particleEnergy << G4endl;
}

G4ThreeVector SteppingAction::transformDirection(G4ThreeVector position, G4ThreeVector worldMomentum)
{
  G4double theta = std::asin(position.x() / std::pow(position.y() * position.y() + position.x() * position.x(), 0.5));
  if ((position.x() > 0) && (position.y() > 0))
    // # positive-positive quadrant
    theta = theta;
  else if ((position.x() > 0) && (position.y() < 0))
    // # positive-negative quadrant
    theta = 3.14159 - theta;
  else if ((position.x() < 0) && (position.y() < 0))
    // # negative-negative quadrant
    theta = std::abs(theta) + 3.14159;
  else if ((position.x() < 0) && (position.y() > 0))
    // # negative-positive quadrant
    theta = theta;

  G4ThreeVector newMomentum = G4ThreeVector(worldMomentum.x() * std::cos(theta) - worldMomentum.y() * std::sin(theta), worldMomentum.x() * std::sin(theta) + worldMomentum.y() * std::cos(theta), worldMomentum.z());

  return newMomentum;
}