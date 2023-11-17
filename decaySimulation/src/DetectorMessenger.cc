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
/// \file DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction *Det)
    : G4UImessenger(), fDetector(Det)//, fPosMin(0), fPosMax(0)
{
  // fPosMin = new G4UIcmdWithADoubleAndUnit("/det/cellPos_min", this);
  // fPosMin->SetGuidance("Set min cell positions");
  // fPosMin->SetParameterName("min", false);
  // fPosMin->SetRange("min>0.");
  // fPosMin->SetUnitCategory("Length");
  // fPosMin->AvailableForStates(G4State_PreInit, G4State_Idle);
  // fPosMin->SetToBeBroadcasted(false);

  // fPosMax = new G4UIcmdWithADoubleAndUnit("/det/cellPos_max", this);
  // fPosMax->SetGuidance("Set max cell positions");
  // fPosMax->SetParameterName("max", false);
  // fPosMax->SetRange("max>0.");
  // fPosMax->SetUnitCategory("Length");
  // fPosMax->AvailableForStates(G4State_PreInit, G4State_Idle);
  // fPosMax->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{

  // delete fPosMin;
  // delete fPosMax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void DetectorMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
// {
//   if(command == fPosMin ) 
//   {
//     fDetector->SetMin(fPosMin->GetNewDoubleValue(newValue));
//      }
//     if(command == fPosMax ) 
//   {
//     fDetector->SetMax(fPosMax->GetNewDoubleValue(newValue));
//      }
// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......