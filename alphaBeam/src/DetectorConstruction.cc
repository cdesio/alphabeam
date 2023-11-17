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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class


#include "DetectorConstruction.hh"

#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4VisAttributes.hh"
#include "G4NistManager.hh"
#include <G4SystemOfUnits.hh>
#include "G4UserLimits.hh"
#include "G4UnitsTable.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PVPlacement.hh"

#include "CommandLineParser.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"

#include "G4ProductionCuts.hh"
#include "G4RunManager.hh"
#include "DetectorMessenger.hh"

#include "RunAction.hh"

using namespace G4DNAPARSER;

#define countof(x) (sizeof(x) / sizeof(x[0]))

using namespace std;
using CLHEP::angstrom;
using CLHEP::degree;
using CLHEP::micrometer;
using CLHEP::mm;
using CLHEP::nanometer;

static G4VisAttributes visGrey(true, G4Colour(0.839216, 0.839216, 0.839216));
static G4VisAttributes invisGrey(false, G4Colour(0.839216, 0.839216, 0.839216));
static G4VisAttributes visRed(false, G4Colour(1, 0, 0));
static G4VisAttributes visBlue(true, G4Colour(0.0, 0.0, 1.0));


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction()
{
  // R = {155 * micrometer, 175 * micrometer, 195 * micrometer, 215 * micrometer, 235 * micrometer, 255 * micrometer, 275 * micrometer, 295 * micrometer, 315 * micrometer, 335 * micrometer};

  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::Construct()
{

  /***************************************************************************/
  //                               World
  /***************************************************************************/

  G4NistManager *man = G4NistManager::Instance();
  G4Material *waterMaterial = man->FindOrBuildMaterial("G4_WATER");
  G4Material *S_Steel = G4NistManager::Instance()->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  G4Material *air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");


  G4double nucleusSize = 300*nm;
  G4double margin = 20*nm;

  G4Box *solidWorld = new G4Box("world", 100 * mm, 100 * mm, 100 * mm);

  G4Box *solidWater = new G4Box("water", 10 * mm, 10 * mm, 10 * mm);

  G4Box *solidVoxel = new G4Box("voxel", nucleusSize/2+margin, nucleusSize/2 + margin, nucleusSize/2 + margin);
  
  G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld,
                                                    air,
                                                    "world");

  G4PVPlacement *physiWorld = new G4PVPlacement(0,
                                                G4ThreeVector(),
                                                "world",
                                                logicWorld,
                                                0,
                                                false,
                                                0);

  G4LogicalVolume *logicWater = new G4LogicalVolume(solidWater,
                                                    waterMaterial,
                                                    "water");
  G4PVPlacement *physiWater = new G4PVPlacement(0,
                                                G4ThreeVector(),
                                                logicWater,
                                                "water",
                                                logicWorld,
                                                0,
                                                false,
                                                0);


  G4LogicalVolume *logicVoxel = new G4LogicalVolume(solidVoxel,
                                                         waterMaterial,
                                                         "voxel");
  G4UserLimits* userLimits = new G4UserLimits();
  userLimits->SetMaxAllowedStep(3 * nm);
  G4int noVoxels =0;
  G4double spacing = 0.5;
  for (G4int i=-10; i <10; i++){
      for (G4int j = -5; j<5; j++){
        for (G4int k=0; k<60; k++){
            

            // G4RotationMatrix* rot = new G4RotationMatrix(-1*theta, 
            //                                   0,
            //                                   0) ;

            G4PVPlacement *physiCell = new G4PVPlacement(0,
                                                     G4ThreeVector(0.5*um + i*spacing*um, 0.5*um + j*spacing* um, 5* um + k*spacing * um),
                                                     logicVoxel,
                                                     "voxel",
                                                     logicWater,
                                                     0,
                                                     k,
                                                     0);
        noVoxels++; 
        }


      }
    
 }
  G4cout << "placed " << noVoxels << " voxels. " << G4endl;
  logicVoxel->SetVisAttributes(&visBlue);
  logicWorld->SetVisAttributes(&invisGrey);
  logicWater->SetVisAttributes(&invisGrey);
  
  return physiWorld;
}
