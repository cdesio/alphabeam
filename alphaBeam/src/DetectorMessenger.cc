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
#include "G4UIcmdWithAnInteger.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction *Det)
    : G4UImessenger(), fDetector(Det), spacing(0), start_Z(0), ndiv_X(0), ndiv_Y(0), ndiv_Z(0)
{
  start_Z = new G4UIcmdWithADoubleAndUnit("/det/set_startZ",this);
  start_Z->SetGuidance("Set starting Z coords of voxels");
  start_Z->SetParameterName("start_Z",false);
  start_Z->SetDefaultValue(5);
  start_Z->SetDefaultUnit("micrometer");
  start_Z->AvailableForStates(G4State_PreInit,G4State_Idle);
  start_Z->SetToBeBroadcasted(false);

  spacing = new G4UIcmdWithADoubleAndUnit("/det/set_spacing",this);
  spacing->SetGuidance("Set radial spacing of boxes");
  spacing->SetParameterName("spacing",false);
  spacing->SetRange("spacing>0.");
  spacing->SetDefaultValue(0.5);
  spacing->SetDefaultUnit("micrometer");
  spacing->AvailableForStates(G4State_PreInit,G4State_Idle);
  spacing->SetToBeBroadcasted(false);

  ndiv_X = new G4UIcmdWithAnInteger("/det/set_ndiv_X",this);
  ndiv_X->SetGuidance("Set no. divisions in X");
  ndiv_X->SetParameterName("ndiv_X",false);
  ndiv_X->SetDefaultValue(10);
  ndiv_X->AvailableForStates(G4State_PreInit,G4State_Idle);
  ndiv_X->SetToBeBroadcasted(false);
  
  ndiv_Z = new G4UIcmdWithAnInteger("/det/set_ndiv_Z",this);
  ndiv_Z->SetGuidance("Set no. divisions in z");
  ndiv_Z->SetParameterName("ndiv_Z",false);
  ndiv_Z->SetDefaultValue(100);
  ndiv_Z->AvailableForStates(G4State_PreInit,G4State_Idle);
  ndiv_Z->SetToBeBroadcasted(false);
  
  ndiv_Y = new G4UIcmdWithAnInteger("/det/set_ndiv_Y",this);
  ndiv_Y->SetGuidance("Set no. divisions in Y");
  ndiv_Y->SetParameterName("ndiv_theta",false);
  ndiv_Y->SetDefaultValue(10);
  ndiv_Y->AvailableForStates(G4State_PreInit,G4State_Idle);
  ndiv_Y->SetToBeBroadcasted(false);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{


delete ndiv_X;
delete spacing;
delete start_Z;
delete ndiv_Y;
delete ndiv_Z;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
  if( command == spacing )
  {
     fDetector->set_spacing(spacing->GetNewDoubleValue(newValue));
  }
  if( command == ndiv_Z )
  {
     fDetector->set_ndiv_Z(ndiv_Z->GetNewIntValue(newValue));
  }
  if( command == ndiv_Y )
  {
     fDetector->set_ndiv_Y(ndiv_Y->GetNewIntValue(newValue));
  }
  if( command == ndiv_X )
  {
     fDetector->set_ndiv_X(ndiv_X->GetNewIntValue(newValue));
  }
  if (command == start_Z)
  {
     fDetector->set_startZ(start_Z->GetNewDoubleValue(newValue));
  }}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......