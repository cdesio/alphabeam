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
  // G4cout << particleName << " Track ID = " << step->GetTrack()->GetTrackID() << " creator process = " << step->GetTrack()->GetCreatorProcess()->GetProcessName() << " KE = " << step->GetPreStepPoint()->GetKineticEnergy() << " parent = " << step->GetTrack()->GetParentID() << "position = " << step->GetPreStepPoint()->GetPosition() << ", process = " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << " " << fpEventAction->parentParticle[step->GetTrack()->GetParentID()] << G4endl;
  // G4cout << particleName << " " << step->GetPreStepPoint()->GetKineticEnergy() << " " << step->GetPostStepPoint()->GetPhysicalVolume()->GetName() << " " << step->GetPostStepPoint()->GetPosition()<< G4endl;
  // }

  G4int TrackID = step->GetTrack()->GetTrackID();
  if ((fpEventAction->parentParticle.find(TrackID) == fpEventAction->parentParticle.end()) && (step->GetTrack()->GetCreatorProcess() != nullptr))
  {
    // track ID not found, save which to map, trackID: creator particle from decay (e-,alpha, gamma) for split of DNA damage by source
    if (step->GetTrack()->GetCreatorProcess()->GetProcessName() == "RadioactiveDecay")
    {
      if (G4StrUtil::contains(particleName, "[")) // to catch excited states
      {
        fpEventAction->parentParticle.insert(std::pair<G4int, G4int>(TrackID, particleOriginMap[particleName.substr(0, 5)]));
      }
      if ((particleName == "e-") || (particleName == "gamma") || (particleName == "alpha") || (particleName == "e+")) // save product and parent names
      {
        G4String parentName = reverseParticleOriginMap[fpEventAction->parentParticle[step->GetTrack()->GetParentID()]];
        G4String combinedName = particleName + parentName;

        fpEventAction->parentParticle.insert(std::pair<G4int, G4int>(TrackID, particleOriginMap[combinedName]));
      }
      else
      {
        fpEventAction->parentParticle.insert(std::pair<G4int, G4int>(TrackID, particleOriginMap[particleName]));
      }
    }
    else
    {
      // not radioactive decay so another process so parent ID should be in mapping
      G4int parentParticle = fpEventAction->parentParticle[step->GetTrack()->GetParentID()];
      // add current track with parent particle
      fpEventAction->parentParticle.insert(std::pair<G4int, G4int>(TrackID, parentParticle));
    }
  }

  if (step->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "world")
  {
    step->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);

    return;
  }

  if ((particleName == "Rn220") && (step->GetPreStepPoint()->GetKineticEnergy() == 0))
  {
    // Desorption from the source through recoil only
    if (step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "seed")
    {
      fpEventAction->addRnDesorptionIN();
    }
  }
  if ((G4StrUtil::contains(particleName, "Pb212")) && (step->GetPreStepPoint()->GetKineticEnergy() == 0))
  {
    if (step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "seed")
    {
      fpEventAction->addPbDesorptionIN();
    }
  }

  G4double particleEnergy = step->GetPreStepPoint()->GetKineticEnergy();
  if ((particleName == "Pb212") && (particleEnergy == 0))
  {
    if (step->GetPreStepPoint()->GetPhysicalVolume()->GetName() != "seed")
    {
      G4double alpha = 0.6931471806 / (10.64 * 60 * 60); // clearance rate for Pb due to vascular routes, equal to lambda pb resulting in a leakage of 50% Phys. Med. Biol. 65 (2020)
      G4double leakTime = step->GetPostStepPoint()->GetLocalTime() / s;
      G4double p = alpha * leakTime;
      p = std::pow(2.718, -1 * p);

      if (p <= G4UniformRand())
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
  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

  // Calculate dose for all rings
  if (volumeNamePre == "cell")
  {
    G4double radius = fDetector->R[step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo()] / um;
    G4double edep = step->GetTotalEnergyDeposit() / joule;

    G4double OutR = (radius + .150) * 1e-6; // m
    G4double InR = (radius - .150) * 1e-6;  // m
    G4double seedLength = 6.0e-3;           // mm

    G4double volumeCylinder = (3.14159 * seedLength * (OutR * OutR - InR * InR));
    G4double density = 1000; // water
    G4double massCylinder = density * volumeCylinder;

    analysisManager->FillH1(0, radius, edep / massCylinder);
  }

  // save all steps entering rings
  if ((volumeNamePre == "water") && (step->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "cell"))
  {
    if (step->GetPreStepPoint()->GetKineticEnergy() > 0)
    {

      G4ThreeVector worldPos = step->GetPostStepPoint()->GetPosition();
      G4double newX = (G4UniformRand() * .00015 * 2) - .00015;
      G4double newZ = (G4UniformRand() * .00015 * 2) - .00015;

      G4double radius = std::pow(worldPos.x() * worldPos.x() + worldPos.y() * worldPos.y(), 0.5);

      G4double newY = radius - fDetector->R[step->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo()];

      // pick position in box frame
      G4ThreeVector newPos = G4ThreeVector(newX, newY, newZ);
      G4ThreeVector newMomentum = transformDirection(step->GetPostStepPoint()->GetPosition(), step->GetPostStepPoint()->GetMomentumDirection());
      savePoint(step->GetTrack(), newPos, newMomentum, step->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo(), step->GetPostStepPoint()->GetKineticEnergy(), step->GetPostStepPoint()->GetGlobalTime(), fpEventAction->parentParticle[TrackID]);
    }
  }
  // save decay in box
  else if ((volumeNamePre == "cell") && (step->IsFirstStepInVolume()) && (particleName != "gamma")) // save particles created in the cell or nucleus if from radioactive decay as not simulated in RBE
  {
    if (step->GetPreStepPoint()->GetProcessDefinedStep() == nullptr)
    // if prestep process is nullptr this is the first step of particle created by interaction in the cell - only save those created by processes in cell not in other volumes
    {
      if ((step->GetTrack()->GetCreatorProcess()->GetProcessName() == "RadioactiveDecay") && (((const G4Ions *)(step->GetTrack()->GetParticleDefinition()))->GetExcitationEnergy() < 1e-10))
      {
        // only save products of radioactive decay other products are from parents which are saved on entering the cell and will be tracked in DNA simulation. Excited states are not saved as de-excitation is not simulated in RBE, products are saved to phase space file.
        // G4cout << "saved Rdecay" << G4endl;

        G4int parentID = step->GetTrack()->GetParentID();

        G4ThreeVector worldPos = step->GetPreStepPoint()->GetPosition();

        // decay products should start in the same place in box reference frame, check if is first product
        if (fpEventAction->decayPos.find(parentID) == fpEventAction->decayPos.end())
        {
          // parent ID not found, is first product, pick new position and save
          G4double newX = (G4UniformRand() * .00015 * 2) - .00015;
          G4double newZ = (G4UniformRand() * .00015 * 2) - .00015;

          G4double radius = std::pow(worldPos.x() * worldPos.x() + worldPos.y() * worldPos.y(), 0.5);
          G4double newY = radius - fDetector->R[step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo()];
          // pick position in box frame
          G4ThreeVector newPos = G4ThreeVector(newX, newY, newZ);

          // save
          fpEventAction->decayPos.insert(std::pair<int, G4ThreeVector>(parentID, newPos));

          G4ThreeVector newMomentum = transformDirection(worldPos, step->GetPreStepPoint()->GetMomentumDirection());

          savePoint(step->GetTrack(), newPos, newMomentum, step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo(), step->GetPreStepPoint()->GetKineticEnergy(), step->GetPreStepPoint()->GetGlobalTime(), fpEventAction->parentParticle[TrackID]);
        }
        else
        {
          // parent ID found, look up new position
          G4ThreeVector newPos = fpEventAction->decayPos[parentID];

          G4ThreeVector newMomentum = transformDirection(worldPos, step->GetPreStepPoint()->GetMomentumDirection());

          savePoint(step->GetTrack(), newPos, newMomentum, step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo(), step->GetPreStepPoint()->GetKineticEnergy(), step->GetPreStepPoint()->GetGlobalTime(), fpEventAction->parentParticle[TrackID]);
        }
      }
    }
  }
  else if ((volumeNamePre == "cell") && ((std::find(fpEventAction->tracks.begin(), fpEventAction->tracks.end(), TrackID) != fpEventAction->tracks.end())) && (particleName != "gamma")) // is a step in the cell but not first and it is a track which has previously been saved i.e. not a secondary created in the box, check if crosses virtual box boundary, if it does it is saved as if entering the box from the side faces.
  {
    if (step->GetPreStepPoint()->GetKineticEnergy() > 0)
    {

      G4ThreeVector entryPosition = fpEventAction->particlePos[step->GetTrack()->GetTrackID()]; // look up position in box frame from last step

      G4ThreeVector deltaWorld = step->GetPostStepPoint()->GetPosition() - step->GetPreStepPoint()->GetPosition(); // change in position in world frame

      G4ThreeVector boxMomentumPre = transformDirection(step->GetPreStepPoint()->GetPosition(), step->GetPreStepPoint()->GetMomentumDirection()); // particle momentum in box frame

      G4ThreeVector delta = deltaWorld.mag() * boxMomentumPre; // change in position in box frame

      G4ThreeVector postStepBox = entryPosition + (fpEventAction->particleDist)[step->GetTrack()->GetTrackID()] + delta; // post step position in box frame

      if ((std::abs(postStepBox.x()) >= 0.00015) || (std::abs(postStepBox.y()) >= 0.00015) || (std::abs(postStepBox.z()) >= 0.00015)) // if >=0.00015 has crossed the boundary
      {
        // save particle, new position and distance saved
        G4ThreeVector preStepBox = entryPosition + (fpEventAction->particleDist)[step->GetTrack()->GetTrackID()]; // pre step point position in box frame

        G4double distanceToExit = calculateDistanceToExitBox(preStepBox, boxMomentumPre);
        // if (std::abs(preStepBox.y()) >= 0.00015)
        // {
        //   return;
        // }

        if (distanceToExit == DBL_MAX)
        // exit y
        // update start position to y exit point and zero distance travelled, in case scattering changes direction
        {
          G4double tYneg = (-.00015 - preStepBox.y()) / boxMomentumPre.y();
          G4double tYpos = (.00015 - preStepBox.y()) / boxMomentumPre.y();
          tYneg = tYneg > 1e-10 ? tYneg : DBL_MAX;
          tYpos = tYpos > 1e-10 ? tYpos : DBL_MAX;

          G4double distanceToExit = std::min({tYpos, tYneg}); // shortest distance travelled to cross y box surface

          G4ThreeVector newPos = preStepBox + (distanceToExit * boxMomentumPre);

          fpEventAction->particleDist.erase(step->GetTrack()->GetTrackID());
          fpEventAction->particleDist.insert(std::pair<int, G4ThreeVector>()); // travelled (0,0,0) from the new starting position

          fpEventAction->particlePos.erase(step->GetTrack()->GetTrackID()); // erase current saved box entry position for this track
          fpEventAction->particlePos.insert(std::pair<int, G4ThreeVector>(step->GetTrack()->GetTrackID(), newPos));

          return;
        }
        else
        {
          G4double stepDistance = step->GetStepLength();

          G4ThreeVector newPos = preStepBox + (distanceToExit * boxMomentumPre);

          // check which side of the box was crossed and change sign as particle is entering adjacent box
          if ((newPos.x() > 0) && (std::abs(newPos.x() - 0.00015) < 1e-10))
            newPos.setX(-0.00015);
          else if ((newPos.x() < 0) && (std::abs(newPos.x() + 0.00015) < 1e-10))
            newPos.setX(+0.00015);

          if ((newPos.z() > 0) && (std::abs(newPos.z() - 0.00015) < 1e-10))
            newPos.setZ(-0.00015);
          else if ((newPos.z() < 0) && (std::abs(newPos.z() + 0.00015) < 1e-10))
            newPos.setZ(+0.00015);

          G4double percentageOfStep = distanceToExit / stepDistance;

          G4double percentageAccountedFor = percentageOfStep;

          // calculate KE at point where crossing occurs
          G4double newKE = step->GetPreStepPoint()->GetKineticEnergy() - (step->GetPreStepPoint()->GetKineticEnergy() - step->GetPostStepPoint()->GetKineticEnergy()) * percentageOfStep;

          // calculate time at point where crossing occurs
          G4double newTime = step->GetPreStepPoint()->GetGlobalTime() + (step->GetDeltaTime() * percentageOfStep);

          savePoint(step->GetTrack(), newPos, boxMomentumPre, step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo(), newKE, newTime, fpEventAction->parentParticle[TrackID]);

          while (percentageAccountedFor < 1)
          {
            G4ThreeVector restOfStep = newPos + (stepDistance - distanceToExit) * boxMomentumPre;

            if ((std::abs(restOfStep.x()) < 0.00015) && (std::abs(restOfStep.y()) < 0.00015) && (std::abs(restOfStep.z()) < 0.00015))
            {
              // remainder of step is contained in the adjacent box
              // save remainder of track travel to next box
              fpEventAction->particleDist.erase(step->GetTrack()->GetTrackID());
              fpEventAction->particleDist.insert(std::pair<int, G4ThreeVector>(step->GetTrack()->GetTrackID(), G4ThreeVector()));

              fpEventAction->particlePos.erase(step->GetTrack()->GetTrackID()); // erase current saved box entry position for this track

              fpEventAction->particlePos.insert(std::pair<int, G4ThreeVector>(step->GetTrack()->GetTrackID(), newPos + (stepDistance - distanceToExit) * boxMomentumPre)); // add current box entry position for this track
              percentageAccountedFor = 1;
            }
            else
            {
              // crosses another boundary
              // find which boundary
              G4double distanceToExitRemainder = calculateDistanceToExitBox(newPos, boxMomentumPre);

              if (distanceToExitRemainder == DBL_MAX) // exit y
              // update start position to y exit point and zero distance travelled, in case scattering changes direction
              {
                G4double tYneg = (-.00015 - preStepBox.y()) / boxMomentumPre.y();
                G4double tYpos = (.00015 - preStepBox.y()) / boxMomentumPre.y();
                tYneg = tYneg > 1e-10 ? tYneg : DBL_MAX;
                tYpos = tYpos > 1e-10 ? tYpos : DBL_MAX;

                G4double distanceToExit = std::min({tYpos, tYneg}); // shortest distance travelled to cross y box surface

                G4ThreeVector newPos = preStepBox + (distanceToExit * boxMomentumPre);

                fpEventAction->particleDist.erase(step->GetTrack()->GetTrackID());
                fpEventAction->particleDist.insert(std::pair<int, G4ThreeVector>()); // travelled (0,0,0) from the new starting position

                fpEventAction->particlePos.erase(step->GetTrack()->GetTrackID()); // erase current saved box entry position for this track
                fpEventAction->particlePos.insert(std::pair<int, G4ThreeVector>(step->GetTrack()->GetTrackID(), newPos));

                percentageAccountedFor = 1;

                return;
              }

              newPos += (distanceToExitRemainder * boxMomentumPre);

              // check which side of the box was crossed and change sign as particle is entering adjacent box
              if ((newPos.x() > 0) && (std::abs(newPos.x() - 0.00015) < 1e-10))
                newPos.setX(-0.00015);
              else if ((newPos.x() < 0) && (std::abs(newPos.x() + 0.00015) < 1e-10))
                newPos.setX(+0.00015);

              if ((newPos.z() > 0) && (std::abs(newPos.z() - 0.00015) < 1e-10))
                newPos.setZ(-0.00015);
              else if ((newPos.z() < 0) && (std::abs(newPos.z() + 0.00015) < 1e-10))
                newPos.setZ(+0.00015);

              percentageOfStep = distanceToExitRemainder / stepDistance;

              newKE = newKE - (newKE - step->GetPostStepPoint()->GetKineticEnergy()) * percentageOfStep;

              newTime += (step->GetDeltaTime() * percentageOfStep);
              distanceToExit += distanceToExitRemainder;

              savePoint(step->GetTrack(), newPos, boxMomentumPre, step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo(), newKE, newTime, fpEventAction->parentParticle[TrackID]);

              percentageAccountedFor += percentageOfStep;
            }
          }
        }
      }

      else
      {
        //  if position hasn't crossed the box bounday update distance travelled in box
        G4ThreeVector previousDelta = (fpEventAction->particleDist)[step->GetTrack()->GetTrackID()];

        // G4cout << "add delta " << delta << G4endl;
        fpEventAction->particleDist.erase(step->GetTrack()->GetTrackID());
        fpEventAction->particleDist.insert(std::pair<int, G4ThreeVector>(step->GetTrack()->GetTrackID(), previousDelta + delta));
      }
    }
  }
  else if ((volumeNamePre == "cell") && ((std::find(fpEventAction->tracks.begin(), fpEventAction->tracks.end(), TrackID) != fpEventAction->tracks.end())) && (particleName == "gamma")) // is a step in the cell but not first and it is a track which has previously been saved i.e. not a secondary created in the box, check if crosses virtual box boundary, if it does it is saved as if entering the box from the side faces.

  // gamma steps are not limited by the strp limiter process, so the crossing of boxes does not occur. It can be assumed that between steps the path is a straight line and there is no energy loss. The entry point is saved on entrance to the volume.
  {

    G4ThreeVector preStepBox = fpEventAction->particlePos[step->GetTrack()->GetTrackID()]; // look up position in box frame from last step

    G4ThreeVector deltaWorld = step->GetPostStepPoint()->GetPosition() - step->GetPreStepPoint()->GetPosition(); // change in position in world frame

    G4ThreeVector boxMomentumPre = transformDirection(step->GetPreStepPoint()->GetPosition(), step->GetPreStepPoint()->GetMomentumDirection()); // particle momentum in box frame

    G4ThreeVector delta = deltaWorld.mag() * boxMomentumPre; // change in position in box frame

    G4double distanceToExit = calculateDistanceToExitBox(preStepBox, boxMomentumPre);

    // G4cout << "distanceToExit" << distanceToExit << G4endl;
    if (distanceToExit == DBL_MAX)
    // exit y
    // update start position to y exit point and zero distance travelled, in case scattering changes direction
    {
      G4double tYneg = (-.00015 - preStepBox.y()) / boxMomentumPre.y();
      G4double tYpos = (.00015 - preStepBox.y()) / boxMomentumPre.y();
      tYneg = tYneg > 1e-10 ? tYneg : DBL_MAX;
      tYpos = tYpos > 1e-10 ? tYpos : DBL_MAX;

      G4double distanceToExit = std::min({tYpos, tYneg}); // shortest distance travelled to cross y box surface

      G4ThreeVector newPos = preStepBox + (distanceToExit * boxMomentumPre);

      fpEventAction->particleDist.erase(step->GetTrack()->GetTrackID());
      fpEventAction->particleDist.insert(std::pair<int, G4ThreeVector>()); // travelled (0,0,0) from the new starting position

      fpEventAction->particlePos.erase(step->GetTrack()->GetTrackID()); // erase current saved box entry position for this track
      fpEventAction->particlePos.insert(std::pair<int, G4ThreeVector>(step->GetTrack()->GetTrackID(), newPos));
    }
    else if (distanceToExit > deltaWorld.mag())
    {
      G4ThreeVector boxMomentumPost = transformDirection(step->GetPostStepPoint()->GetPosition(), step->GetPostStepPoint()->GetMomentumDirection()); // particle momentum in box frame for next step

      // step is contained within box save end point as start of next step
      savePoint(step->GetTrack(), preStepBox + delta, boxMomentumPost, step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo(), step->GetPreStepPoint()->GetKineticEnergy(), step->GetPreStepPoint()->GetGlobalTime(), fpEventAction->parentParticle[TrackID]);
    }
    else
    {
      // crosses x or z boundary
      G4double stepDistance = step->GetStepLength();

      G4ThreeVector newPos = preStepBox + (distanceToExit * boxMomentumPre);

      // check which side of the box was crossed and change sign as particle is entering adjacent box
      if ((newPos.x() > 0) && (std::abs(newPos.x() - 0.00015) < 1e-10))
        newPos.setX(-0.00015);
      else if ((newPos.x() < 0) && (std::abs(newPos.x() + 0.00015) < 1e-10))
        newPos.setX(+0.00015);

      if ((newPos.z() > 0) && (std::abs(newPos.z() - 0.00015) < 1e-10))
        newPos.setZ(-0.00015);
      else if ((newPos.z() < 0) && (std::abs(newPos.z() + 0.00015) < 1e-10))
        newPos.setZ(+0.00015);

      G4double percentageOfStep = distanceToExit / stepDistance;

      G4double percentageAccountedFor = percentageOfStep;

      // calculate KE at point where crossing occurs
      G4double newKE = step->GetPreStepPoint()->GetKineticEnergy() - (step->GetPreStepPoint()->GetKineticEnergy() - step->GetPostStepPoint()->GetKineticEnergy()) * percentageOfStep;

      // calculate time at point where crossing occurs
      G4double newTime = step->GetPreStepPoint()->GetGlobalTime() + (step->GetDeltaTime() * percentageOfStep);

      savePoint(step->GetTrack(), newPos, boxMomentumPre, step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo(), newKE, newTime, fpEventAction->parentParticle[TrackID]);

      while (percentageAccountedFor < 1)
      {
        G4ThreeVector restOfStep = newPos + (stepDistance - distanceToExit) * boxMomentumPre;

        if ((std::abs(restOfStep.x()) < 0.00015) && (std::abs(restOfStep.y()) < 0.00015) && (std::abs(restOfStep.z()) < 0.00015))
        {
          // remainder of step is contained in the adjacent box
          // save remainder of track travel to next box
          fpEventAction->particleDist.erase(step->GetTrack()->GetTrackID());
          fpEventAction->particleDist.insert(std::pair<int, G4ThreeVector>(step->GetTrack()->GetTrackID(), G4ThreeVector()));

          fpEventAction->particlePos.erase(step->GetTrack()->GetTrackID()); // erase current saved box entry position for this track

          fpEventAction->particlePos.insert(std::pair<int, G4ThreeVector>(step->GetTrack()->GetTrackID(), newPos + (stepDistance - distanceToExit) * boxMomentumPre)); // add current box entry position for this track
          percentageAccountedFor = 1;
        }
        else
        {
          // crosses another boundary
          // find which boundary
          G4double distanceToExitRemainder = calculateDistanceToExitBox(newPos, boxMomentumPre);

          if (distanceToExitRemainder == DBL_MAX) // exit y
          // update start position to y exit point and zero distance travelled, in case scattering changes direction
          {
            G4double tYneg = (-.00015 - preStepBox.y()) / boxMomentumPre.y();
            G4double tYpos = (.00015 - preStepBox.y()) / boxMomentumPre.y();
            tYneg = tYneg > 1e-10 ? tYneg : DBL_MAX;
            tYpos = tYpos > 1e-10 ? tYpos : DBL_MAX;

            G4double distanceToExit = std::min({tYpos, tYneg}); // shortest distance travelled to cross y box surface

            newPos = preStepBox + (distanceToExit * boxMomentumPre);

            fpEventAction->particleDist.erase(step->GetTrack()->GetTrackID());
            fpEventAction->particleDist.insert(std::pair<int, G4ThreeVector>()); // travelled (0,0,0) from the new starting position

            fpEventAction->particlePos.erase(step->GetTrack()->GetTrackID()); // erase current saved box entry position for this track
            fpEventAction->particlePos.insert(std::pair<int, G4ThreeVector>(step->GetTrack()->GetTrackID(), newPos));

            percentageAccountedFor = 1;

            return;
          }

          newPos += (distanceToExitRemainder * boxMomentumPre);

          // check which side of the box was crossed and change sign as particle is entering adjacent box

          // G4cout << newPos.x() << G4endl;
          if ((newPos.x() > 0) && (std::abs(newPos.x() - 0.00015) < 1e-10))
            newPos.setX(-0.00015);
          else if ((newPos.x() < 0) && (std::abs(newPos.x() + 0.00015) < 1e-10))
            newPos.setX(+0.00015);

          if ((newPos.z() > 0) && (std::abs(newPos.z() - 0.00015) < 1e-10))
            newPos.setZ(-0.00015);
          else if ((newPos.z() < 0) && (std::abs(newPos.z() + 0.00015) < 1e-10))
            newPos.setZ(+0.00015);

          percentageOfStep = distanceToExitRemainder / stepDistance;

          newKE = newKE - (newKE - step->GetPostStepPoint()->GetKineticEnergy()) * percentageOfStep;

          newTime += (step->GetDeltaTime() * percentageOfStep);
          distanceToExit += distanceToExitRemainder;

          savePoint(step->GetTrack(), newPos, boxMomentumPre, step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo(), newKE, newTime, fpEventAction->parentParticle[TrackID]);

          percentageAccountedFor += percentageOfStep;
        }
      }
    }
  }
}

void SteppingAction::savePoint(const G4Track *track, G4ThreeVector newPos, G4ThreeVector boxMomentum, const int copy, G4double particleEnergy, G4double time, G4int originParticle)
{
  // save particle to phase space file in box reference frame

  fpEventAction->particlePos.erase(track->GetTrackID()); // erase current saved box entry position for this track

  fpEventAction->particlePos.insert(std::pair<int, G4ThreeVector>(track->GetTrackID(), newPos)); // add current box entry position for this track

  fpEventAction->particleDist.erase(track->GetTrackID());                                                  // erase distance travelled for this track in the box
  fpEventAction->particleDist.insert(std::pair<int, G4ThreeVector>(track->GetTrackID(), G4ThreeVector())); // made initial distance travelled zero

  G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

  G4String particleName = track->GetParticleDefinition()->GetParticleName();

  G4float particleID = particleMap[particleName];
  if (particleID == 0)
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
  output[11] = originParticle;

  PSfile.write((char *)&output, sizeof(output));
  fpEventAction->tracks.push_back(track->GetTrackID());

  // G4cout << particleName << " saved at = " << newPos / mm << " with KE = " << particleEnergy << " with momentum " << boxMomentum << " TracKID = " << track->GetTrackID() << " originParticle " << originParticle << " copy " << copy << G4endl;
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

G4double SteppingAction::calculateDistanceToExitBox(G4ThreeVector preStepPosition, G4ThreeVector preStepMomentumDirection)
{
  // does step exit box in x and z?
  G4double tXneg = (-.00015 - preStepPosition.x()) / preStepMomentumDirection.x();
  G4double tXpos = (.00015 - preStepPosition.x()) / preStepMomentumDirection.x();

  G4double tYneg = (-.00015 - preStepPosition.y()) / preStepMomentumDirection.y();
  G4double tYpos = (.00015 - preStepPosition.y()) / preStepMomentumDirection.y();

  G4double tZneg = (-.00015 - preStepPosition.z()) / preStepMomentumDirection.z();
  G4double tZpos = (.00015 - preStepPosition.z()) / preStepMomentumDirection.z();

  // G4cout << tXneg << " " << tXpos << " " << tYneg << " " << tYpos << " " << tZneg << " " << tZpos << " " << G4endl;

  tXneg = tXneg > 1e-10 ? tXneg : DBL_MAX;
  tXpos = tXpos > 1e-10 ? tXpos : DBL_MAX;
  tYneg = tYneg > 1e-10 ? tYneg : DBL_MAX;
  tYpos = tYpos > 1e-10 ? tYpos : DBL_MAX;
  tZneg = tZneg > 1e-10 ? tZneg : DBL_MAX;
  tZpos = tZpos > 1e-10 ? tZpos : DBL_MAX;

  G4double distanceToExit = std::min({tXneg, tXpos, tZneg, tZpos}); // shortest distance travelled to cross box surface

  // G4cout << tXneg << " " << tXpos << " " << tYneg << " " << tYpos << " " << tZneg << " " << tZpos << " " << G4endl;
  if ((tYneg <= distanceToExit) || (tYpos <= distanceToExit)) // exit y
  {
    return DBL_MAX;
  }

  return distanceToExit;
}