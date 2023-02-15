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
  if ((G4StrUtil::contains(particleName,"Pb212")) && (step->GetPreStepPoint()->GetKineticEnergy() == 0))
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


  G4double edep = step->GetTotalEnergyDeposit();

  // Save Edep from all particles to calculate dose
  if ((edep > 0) && (step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "cell"))
  {
    G4TouchableHandle theTouchable = step->GetPreStepPoint()->GetTouchableHandle();
    G4ThreeVector worldPos = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector localPos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);

    if ((abs(localPos.x())<=150*nm) && (abs(localPos.y())<=150*nm) && (abs(localPos.z())<=150*nm))// check within cell not margin
    {
    G4double time = (step->GetPreStepPoint()->GetGlobalTime() + step->GetPostStepPoint()->GetGlobalTime()) / 2 / s;
    G4int cp = step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo();

    G4double mass = 997 * 300e-9 * 300e-9 * 300e-9; // cube of water desnity water 997 kg/m3
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

    analysisManager->FillH1(20 + cp, time / 60 / 60, edep / joule  / mass);
 
    }
  }
  // Save SB of all alpha particles
  if ((particleName == "alpha") && (step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "cell")) // save for alpha steps in the cell
  {
    G4TouchableHandle theTouchable = step->GetPreStepPoint()->GetTouchableHandle();
    G4ThreeVector worldPos = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector localPos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);

    if ((abs(localPos.x())<=150*nm) && (abs(localPos.y())<=150*nm) && (abs(localPos.z())<=150*nm))// check within cell not margin
    {
    G4double time = (step->GetPreStepPoint()->GetGlobalTime() + step->GetPostStepPoint()->GetGlobalTime()) / 2 / s;

    G4double particleMeanEnergy = (step->GetPreStepPoint()->GetKineticEnergy() + step->GetPostStepPoint()->GetKineticEnergy()) / 2;
    G4int cp = step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo();


    calculateDSB(particleMeanEnergy, step->GetStepLength(), time, cp);
    }
  }


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
      if (G4StrUtil::contains(particleName,"e-"))
        particleID = 1;
      else if (G4StrUtil::contains(particleName,"gamma"))
        particleID = 2;
      else if (G4StrUtil::contains(particleName,"alpha"))
        particleID = 3;
      else if (G4StrUtil::contains(particleName,"Rn220"))
        particleID = 4;
      else if (G4StrUtil::contains(particleName,"Po216"))
        particleID = 5;
      else if (G4StrUtil::contains(particleName,"Pb212"))
        particleID = 6;
      else if (G4StrUtil::contains(particleName,"Bi212"))
        particleID = 7;
      else if (G4StrUtil::contains(particleName,"Tl208"))
        particleID = 8;
      else if (G4StrUtil::contains(particleName,"Po212"))
        particleID = 9;
      else if (G4StrUtil::contains(particleName,"Pb208"))
        particleID = 10;
      else
      {
        G4cout << particleName << " outside  not saved" << G4endl;
        return;
      }
      float output[12];
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
      output[11] = ((const G4Ions*)( step->GetTrack()->GetParticleDefinition()))->GetExcitationEnergy();

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
    if (G4StrUtil::contains(particleName,"e-"))
      particleID = 1;
    else if (G4StrUtil::contains(particleName,"gamma"))
      particleID = 2;
    else if (G4StrUtil::contains(particleName,"alpha"))
      particleID = 3;
    else if (G4StrUtil::contains(particleName,"Rn220"))
      particleID = 4;
    else if (G4StrUtil::contains(particleName,"Po216"))
      particleID = 5;
    else if (G4StrUtil::contains(particleName,"Pb212"))
      particleID = 6;
    else if (G4StrUtil::contains(particleName,"Bi212"))
      particleID = 7;
    else if (G4StrUtil::contains(particleName,"Tl208"))
      particleID = 8;
    else if (G4StrUtil::contains(particleName,"Po212"))
      particleID = 9;
    else if (G4StrUtil::contains(particleName,"Pb208"))
      particleID = 10;
    else
    {
      G4cout << particleName << " inside  not saved" << G4endl;
      return;
    }

    float output[12];
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
    output[11] = ((const G4Ions*)( step->GetTrack()->GetParticleDefinition()))->GetExcitationEnergy();

    PSfile.write((char *)&output, sizeof(output));

  }

}

void SteppingAction::calculateDSB(G4double KE, G4double stepLength, G4double time, G4int cp)
{
  // function to calculate the number of DSB per step (IRT model)
    G4double KEtoSimpleDSB[19][2] = {{0, 0},
                                    { 0.0323994848612319 , 0.000822124951113371 },
                                    { 0.07331541275594013 , 0.001260892261700282 },
                                    { 0.10526875591635568 , 0.0015894685658994524 },
                                    { 0.21507798335636447 , 0.0019816137839281486 },
                                    { 0.2977516237592551 , 0.0021355788894433366 },
                                    { 0.4322209190842615 , 0.002404038864169155 },
                                    { 0.5272446136259942 , 0.0025020526253341103 },
                                    { 0.7781060949171443 , 0.0027010955277722306 },
                                    { 1.0405580968411194 , 0.002826203719055142 },
                                    { 1.305946937526972 , 0.0027134382480948184 },
                                    { 1.8361135787388494 , 0.0024749385171951746 },
                                    { 2.875233873666144 , 0.0021252485845039645 },
                                    { 3.897210833430267 , 0.0017181087287708374 },
                                    { 4.912156939921248 , 0.001474603986659365 },
                                    { 5.923284731371911 , 0.0012570639073881825 },
                                    { 6.931295280324434 , 0.0011539019920032583 },
                                    { 7.937943873452893 , 0.0009839364297663405 },
                                    { 8.943067543805933 , 0.0009202625486618789 }};

    G4double KEtoComplexDSB[19][2] = {{0, 0},
                                      { 0.0323994848612319 , 0.0004661538727662737 },
                                      { 0.07331541275594013 , 0.001043588338656949 },
                                      { 0.10526875591635568 , 0.0014392680179136298 },
                                      { 0.21507798335636447 , 0.0023925828245363865 },
                                      { 0.2977516237592551 , 0.0028391964878529694 },
                                      { 0.4322209190842615 , 0.0033766871078409525 },
                                      { 0.5272446136259942 , 0.0035028713592227282 },
                                      { 0.7781060949171443 , 0.0035090382944870774 },
                                      { 1.0405580968411194 , 0.003279310696866271 },
                                      { 1.305946937526972 , 0.0029336215302320128 },
                                      { 1.8361135787388494 , 0.0022277537138450667 },
                                      { 2.875233873666144 , 0.0013256274508608485 },
                                      { 3.897210833430267 , 0.0009556822035093806 },
                                      { 4.912156939921248 , 0.0007197188643735522 },
                                      { 5.923284731371911 , 0.0005418645525326428 },
                                      { 6.931295280324434 , 0.000447865508992609 },
                                      { 7.937943873452893 , 0.0003379691041638497 },
                                      { 8.943067543805933 , 0.0003270240995646832 }};

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