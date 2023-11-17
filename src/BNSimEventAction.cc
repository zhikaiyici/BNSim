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
// $Id: BNSimEventAction.cc 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file BNSimEventAction.cc
/// \brief Implementation of the BNSimEventAction class

#include "BNSimEventAction.hh"
#include "BNSimRunAction.hh"
#include "BNSimSteppingAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

#include "BNSimSteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

////vector<vector<G4double>> BNSimEventAction::EnergyDeposit;// = *(new vector<vector<G4double>>);
//list<vector<G4double>> BNSimEventAction::EnergyDeposit = *(new list<vector<G4double>>);

BNSimEventAction::BNSimEventAction(BNSimRunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdep(0.),
  EnergyDepositOf2nd(0)
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BNSimEventAction::~BNSimEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BNSimEventAction::BeginOfEventAction(const G4Event* event)
{    
  fEdep = 0.;
  BNSimSteppingAction::Reset();
  //G4cout << "----------------------BeginOfEventAction------------------------------" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BNSimEventAction::EndOfEventAction(const G4Event* event)
{   
  // accumulate statistics in run action
  if (fEdep != 0)
  {
      fRunAction->PushBackEnergyDepositOfEvent(fEdep);
  }

  EnergyDepositOf2nd = BNSimSteppingAction::GetEnergyDepositOf2nd();
  vector<G4double>::iterator itr = EnergyDepositOf2nd.begin();
  //while (itr != EnergyDepositOf2nd.end())
  //{
  //    if (*itr == 0)
  //    {
  //        itr = EnergyDepositOf2nd.erase(itr);
  //    }
  //    else
  //        itr++;
  //}
  //itr = EnergyDepositOf2nd.begin();
  while (itr != EnergyDepositOf2nd.end())
  {
      if (*itr != 0)
      {
          fRunAction->PushBackEnergyDepositOf2nd(*itr);
      }
      itr++;
  }

  G4int eventID = event->GetEventID();
  G4int numberofevents = UserDataInput::GetNumberOfEvents();
  if ((eventID + 1) % (numberofevents / 10) == 0)
  {
      G4double per = (1. * eventID + 1) / (numberofevents * 0.01);
      //G4cout << " numberofevents: "<< numberofevents << " eventID: "<< eventID <<G4endl;
      G4long seconds = time(NULL); // 格林威治时间
      seconds = seconds + 8 * 3600; // 北京时间
      G4int secondnow = seconds % 60;
      G4int minutes = (seconds - secondnow) / 60;
      G4int minutenow = minutes % 60;
      G4int hours = (minutes - minutenow) / 60;
      G4int hournow = hours % 24;
      G4cout << " Time now: " << setw(2) << hournow << ":" << setw(2) << minutenow << ":" << setw(2) << secondnow
             << ". " << setw(3) << per << "% of simulation completed."
             << G4endl;
      //getchar();
  }
  //G4cout << "-----------------------EndOfEventAction------------------------------" << G4endl << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
