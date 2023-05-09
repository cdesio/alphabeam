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
//
/// \file EventAction.hh
/// \brief Definition of the EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include <map>
class EventAction : public G4UserEventAction
{
public:
  EventAction();
  ~EventAction();

public:
  virtual void BeginOfEventAction(const G4Event *);
  virtual void EndOfEventAction(const G4Event *);

  void addRnDeabsorptionIN() { RnDeabsorptionIN++; }
  G4int getRnDeabsorptionIN() { return RnDeabsorptionIN; }

  void addPbDeabsorptionIN() { PbDeabsorptionIN++; }
  G4int getPbDeabsorptionIN() { return PbDeabsorptionIN; }

  void addPbLeakage() { PbLeakage++; }
  void addPbNoLeakage() { PbNoLeakage++; }
  G4int getPbLeakage() { return PbLeakage; }
  G4int getPbNoLeakage() { return PbNoLeakage; }

  void addDecayTimeRa(G4double val) { totalRaDecayTime += val; }
  G4double getTotalRaDecayTime() { return totalRaDecayTime; }

  std::map<G4int, G4ThreeVector> particlePos;
  std::map<G4int, G4ThreeVector> particleDist;

private:
  G4int RnDeabsorptionIN{0};
  G4int PbDeabsorptionIN{0};
  G4int PbLeakage{0};
  G4int PbNoLeakage{0};
  G4double totalRaDecayTime{0};
};

#endif