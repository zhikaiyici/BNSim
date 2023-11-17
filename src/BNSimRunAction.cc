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
// $Id: BNSimRunAction.cc 99560 2016-09-27 07:03:29Z gcosmo $
//
/// \file BNSimRunAction.cc
/// \brief Implementation of the BNSimRunAction class

#include "BNSimRunAction.hh"
#include "BNSimPrimaryGeneratorAction.hh"
#include "BNSimDetectorConstruction.hh"
// #include "B1Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
//#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "BNSimSteppingAction.hh"
#include "BNSimEventAction.hh"
#include "UserDataInput.hh"

#include <fstream>
#include <iomanip>
#include <cmath>
#include <sstream>
/*#include <direct.h>
#include <io.h>*/

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BNSimRunAction::BNSimRunAction()
: G4UserRunAction(),
  neutroncount(0),
  effofneucap(0.),
  EnergyDeposit(0),
  EnergyDepositOf2nd(0)
  /*fEdep(0.),
  fEdep2(0.)*/
{ 
  /*
  // add new units for dose
  // 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray); 

  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fEdep);
  accumulableManager->RegisterAccumulable(fEdep2);
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BNSimRunAction::~BNSimRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BNSimRunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  //// reset accumulables to their initial values
  //G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  //accumulableManager->Reset();

//  G4int nofEvents = run->GetNumberOfEvent();
//  if (nofEvents == 0) return;
//
// // Run conditions
////  note: There is no primary generator action object for "master"
////        run manager for multi-threaded mode.
//  const BNSimPrimaryGeneratorAction* generatorAction
//      = static_cast<const BNSimPrimaryGeneratorAction*>
//      (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
//  G4String runCondition;
//  if (generatorAction)
//  {
//      const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
//      runCondition += particleGun->GetParticleDefinition()->GetParticleName();
//      runCondition += " of ";
//      G4double particleEnergy = particleGun->GetParticleEnergy();
//      runCondition += G4BestUnit(particleEnergy, "Energy");
//  }
//
//  G4cout
//      << G4endl
//      << " The run consists of " << nofEvents << " " << runCondition
//      << G4endl;
//  getchar();

  neutroncount = 0;
  effofneucap = 0.;
  EnergyDeposit.clear();
  EnergyDepositOf2nd.clear();
  BNSimSteppingAction::ResetNeutornCount();

  G4int numofevents = UserDataInput::GetNumberOfEvents();
  G4double dtctrx = UserDataInput::GetDectorDimensionX() / mm;
  G4double dtctry = UserDataInput::GetDectorDimensionY() / mm;
  G4double dtctrz = UserDataInput::GetDectorDimensionZ() / um;

  G4long seconds = time(NULL); // 格林威治时间
  seconds = seconds + 8 * 3600; // 北京时间
  G4int secondnow = seconds % 60;
  G4int minutes = (seconds - secondnow) / 60;
  G4int minutenow = minutes % 60;
  G4int hours = (minutes - minutenow) / 60;
  G4int hournow = hours % 24;
  G4cout << G4endl << " Initialization completed. " << numofevents << " event(s) will be simulated."
         << G4endl
         << " The dimension of detector is " << dtctrx << " mm × " << dtctry << " mm × " << dtctrz << " um" 
         << G4endl
         << G4endl
         << " Simulating..." 
         << G4endl
         << G4endl
         << " Time now: " << setw(2) << hournow << ":" << setw(2) << minutenow << ":" << setw(2) << secondnow << "."
         << G4endl;
  //getchar();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BNSimRunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  neutroncount = BNSimSteppingAction::GetNeutronCount();
  effofneucap = 100. * neutroncount / nofEvents;

  /*// 判断输出目标文件夹是否存在，不存在则创建
  if (_access("output", 4) != 0 )
  {
      _mkdir("output");
      //G4cout << "已创建output" << G4endl;
  }*/

  ofstream EDEPFof2nd, EDEPF, EFFOFNEUF;

  G4double dtctrz = UserDataInput::GetDectorDimensionZ() / um;
  G4double numofevent = nofEvents;
  ostringstream thickness, eventnumber;
  thickness << dtctrz;
  eventnumber << setprecision(1) << numofevent;
  G4String Thickness = thickness.str();
  G4String EventNumber = eventnumber.str();
  G4String EDEPFof2ndname = "output/edepof2nd_" + Thickness + "_um_" + EventNumber + ".txt";
  EDEPFof2nd.open(EDEPFof2ndname, ios_base::out);
  for (list<G4double>::iterator i = EnergyDepositOf2nd.begin(); i != EnergyDepositOf2nd.end(); i++)
  {
      if (*i != 0)
      {
          EDEPFof2nd << *i << G4endl;
      }
  }
  EDEPFof2nd.close();

  G4String EDEPFname = "output/edep_" + Thickness + "_um_" + EventNumber + ".txt";
  EDEPF.open(EDEPFname, ios::out);
  for (list<G4double>::iterator i = EnergyDeposit.begin(); i != EnergyDeposit.end(); i++)
  {
      if (*i != 0)
      {
          EDEPF << *i << G4endl;
      }
  }
  EDEPF.close();

  G4String EFFOFNEUFname = "output/effofneucap_" + EventNumber + ".txt";
  EFFOFNEUF.open(EFFOFNEUFname, ios_base::app);
  EFFOFNEUF << dtctrz << setw(10) << effofneucap << "%" << G4endl;
  EFFOFNEUF.close();

  /*
  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Compute dose = total energy deposit in a run and its variance
  //
  G4double edep  = fEdep.GetValue();
  G4double edep2 = fEdep2.GetValue();
  
  G4double rms = edep2 - edep*edep/nofEvents;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;  

  const BNSimDetectorConstruction* detectorConstruction
   = static_cast<const BNSimDetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double mass = detectorConstruction->GetScoringVolume()->GetMass();
  G4double dose = edep/mass;
  G4double rmsDose = rms/mass;
  */

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const BNSimPrimaryGeneratorAction* generatorAction
   = static_cast<const BNSimPrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }
        
  // Print
  //  
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }
  
  G4cout
     << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     << " The count of captured neutron is " << neutroncount
     << G4endl
     << " The efficiency of neutron capture is " << effofneucap << "%"
     << G4endl
     //<< " Cumulated dose per run, in scoring volume : " 
     //<< G4BestUnit(dose,"Dose") << " rms = " << G4BestUnit(rmsDose,"Dose")
     //<< G4endl
     << "------------------------------------------------------------"
     << G4endl
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//void BNSimRunAction::AddEdep(G4double edep)
//{
//  fEdep  += edep;
//  fEdep2 += edep*edep;
//}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

