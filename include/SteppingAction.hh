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
#pragma once
#include "G4UserSteppingAction.hh"
#include "G4String.hh"
#include <fstream>
#include <iostream>
#include "RunAction.hh"
#include "G4Track.hh"
class EventAction;
class RunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction(DetectorConstruction *fpDet);
  ~SteppingAction() override;

  void UserSteppingAction(const G4Step *step) override;

private:
  EventAction *fpEventAction;
  RunAction *fRunAction;
  std::ofstream PSfile;
  void savePoint(const G4Track *track, G4ThreeVector worldPos, G4ThreeVector worldMomentum, const int copy, G4double KE, G4double time, G4int orignParticle);
  G4ThreeVector transformDirection(G4ThreeVector position, G4ThreeVector worldMomentum);
  DetectorConstruction *fDetector;
  G4double calculateDistanceToExitBox(G4ThreeVector preStepPosition, G4ThreeVector preStepMomentumDirection);

  std::map<G4String, G4int> particleMap{
      {"e-", 1},
      {"gamma", 2},
      {"alpha", 3},
      {"Rn220", 4},
      {"Po216", 5},
      {"Pb212", 6},
      {"Bi212", 7},
      {"Tl208", 8},
      {"Po212", 9},
      {"Pb208", 10},
      {"e+", 11}
      };

  std::map<G4String, G4int> particleOriginMap{
      {"Ra224", 0},
      {"Rn220", 1},
      {"Po216", 2},
      {"Pb212", 3},
      {"Bi212", 4},
      {"Tl208", 5},
      {"Po212", 6},
      {"Pb208", 7},
      {"alphaRa224", 8},
      {"alphaRn220", 9},
      {"alphaPo216", 10},
      {"alphaBi212", 11},
      {"alphaPo212", 12},
      {"e-Rn220", 13},
      {"e-Po216", 14},
      {"e-Pb212", 15},
      {"e-Bi212", 16},
      {"e-Tl208", 17},
      {"e-Po212", 18},
      {"e-Pb208", 19},
      {"gammaRn220", 20},
      {"gammaPo216", 21},
      {"gammaPb212", 22},
      {"gammaBi212", 23},
      {"gammaTl208", 24},
      {"gammaPo212", 25},
      {"gammaPb208", 26},
      {"e+", 27}};

  std::map<G4int, G4String> reverseParticleOriginMap{
      {0, "Ra224"},
      {1, "Rn220"},
      {2, "Po216"},
      {3, "Pb212"},
      {4, "Bi212"},
      {5, "Tl208"},
      {6, "Po212"},
      {7, "Pb208"},
  };
};
