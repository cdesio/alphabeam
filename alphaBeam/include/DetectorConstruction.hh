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

#pragma once
#include "G4VUserDetectorConstruction.hh"
#include <memory>
#include "G4RotationMatrix.hh"
#include "DetectorMessenger.hh"

class G4VPhysicalVolume;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class DetectorConstruction
    : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    ~DetectorConstruction() override;
    G4VPhysicalVolume *Construct() override;
    void set_spacing (G4double);
    void set_startZ(G4double);
    void set_ndiv_X (G4int);
    void set_ndiv_Y(G4int);
    void set_ndiv_Z(G4int);
    DetectorMessenger* fDetectorMessenger;
    G4double get_spacing() { return spacing; }
    G4double get_start_Z() { return start_Z; }
private:

    G4double spacing;
    G4double start_Z;

    G4int ndiv_X;
    G4int ndiv_Y;
    G4int ndiv_Z;
};
