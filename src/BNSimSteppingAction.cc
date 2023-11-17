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
// $Id: BNSimSteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file BNSimSteppingAction.cc
/// \brief Implementation of the BNSimSteppingAction class

#include "BNSimSteppingAction.hh"
#include "BNSimEventAction.hh"
#include "BNSimDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

#include "G4SystemOfUnits.hh"

#include "UserDataInput.hh"

#include <algorithm>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

vector<G4double> BNSimSteppingAction::EnergyDepositOf2nd;// = *(new vector<G4double>);
vector<G4int> BNSimSteppingAction::trackIDs;// = *(new vector<G4int>);
G4int BNSimSteppingAction::neutroncount = 0;

BNSimSteppingAction::BNSimSteppingAction(BNSimEventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0),
  numof2ndprtclofneutron(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BNSimSteppingAction::~BNSimSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BNSimSteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringVolume) { 
    const BNSimDetectorConstruction* detectorConstruction
      = static_cast<const BNSimDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();   
  }

  // push_back trackID和parentID应放在check scoring volume之前，以防止有次级粒子在scoring volume外而漏记
  // 导致后续无法找到三级粒子的母粒子
  G4Track* theTrack = step->GetTrack();
  G4int trackID = theTrack->GetTrackID();
  G4int parentID = theTrack->GetParentID();
  if (parentID > 0)
  {
      trackIDs.push_back(trackID);
      trackIDs.push_back(parentID);
  }

  // get volume of the current step
  G4TouchableHandle touch = step->GetPreStepPoint()->GetTouchableHandle();
  G4LogicalVolume* volume = touch->GetVolume()->GetLogicalVolume();

  // check if we are in scoring volume
  if (volume != fScoringVolume) return;

  G4TrackStatus trackStatus = theTrack->GetTrackStatus();

  //G4cout << "TrackStatus: " << trackStatus << G4endl;
  //G4cout << "trackID: " << trackID << G4endl;
  //G4cout << "parentID: " << parentID << G4endl;
  
  //G4String modelname = theTrack->GetCreatorModelName();
  //G4cout << "modelname: " << modelname << G4endl;

  // 判断是否为中子，准备好vector EnergyDepositOf2nd用于储存次级粒子能量沉积
  if (parentID == 0)
  {
      //numof2ndprtclofneutron = step->GetNumberOfSecondariesInCurrentStep();
      //EnergyDepositOf2nd.resize(numof2ndprtclofneutron); //获取次级粒子数量并以此resize vector EnergyDepositOf2nd
      if (trackStatus == fStopAndKill) //若中子在探测器内终止则中子俘获计数+1
          neutroncount += 1;
  }

  //G4cout << "numof2ndprtclofneutron: " << numof2ndprtclofneutron << G4endl;
  //G4cout << "EnergyDepositOf2nd大小：" << EnergyDepositOf2nd.size() << G4endl;
  //G4cout << "neutroncount: " << neutroncount << G4endl;
  //G4String strPrtclName = theTrack->GetDefinition()->GetParticleName();
  //G4cout << "particalname: " << strPrtclName << G4endl;

  // 有时中子碰撞多次，每步都产生了次级粒子，而GetNumberOfSecondariesInCurrentStep只能获取当前步的次级粒子数量
  // 此时vector容量小于需记录的粒子数，根据trackID再次resize EnergyDepositOf2nd
  if (parentID == 1 && trackID - 1 > EnergyDepositOf2nd.size())
  {
      EnergyDepositOf2nd.resize(trackID - 1);
      //G4cout << "EnergyDepositOf2nd.size(): " << EnergyDepositOf2nd.size() << G4endl;
      //getchar();
  }
  numof2ndprtclofneutron = EnergyDepositOf2nd.size();
  /* if外添加 
        numof2ndprtclofneutron = EnergyDepositOf2nd.size();
     后,此语句可完全替代
        numof2ndprtclofneutron = step->GetNumberOfSecondariesInCurrentStep();
        EnergyDepositOf2nd.resize(numof2ndprtclofneutron);

     Geant4 10.3之前版本Step没有成员函数GetNumberOfSecondariesInCurrentStep()
*/

  //G4cout << "trackIDs: ";
  //for (int i = 0; i < trackIDs.size(); i++)
  //{
  //    G4cout << trackIDs[i] << "  ";
  //}
  //G4cout << G4endl;
  //G4cout << "trackIDs大小：" << trackIDs.size() << G4endl;
  //G4cout << "EnergyDepositOf2nd大小：" << EnergyDepositOf2nd.size() << G4endl;

  //G4float stepLength = step->GetStepLength();		//cm
  //G4StepPoint* prestepPoint = step->GetPreStepPoint();
  //G4StepPoint* poststepPoint = step->GetPostStepPoint();
  //G4StepStatus poststepStatus = poststepPoint->GetStepStatus();
  //G4float PreKineticEnergy = prestepPoint->GetKineticEnergy(); 	//MeV
  //G4float PostKineticEnergy = poststepPoint->GetKineticEnergy(); 	//MeV
  //G4cout << "PreKineticEnergy: " << PreKineticEnergy / keV << G4endl;
  //G4cout << "PostKineticEnergy: " << PostKineticEnergy / keV << G4endl;
  //G4cout << "stepLength: " << stepLength /um << G4endl;
  //G4ThreeVector prepoint = prestepPoint->GetPosition();
  //G4ThreeVector postpoint = poststepPoint->GetPosition();
  //G4cout << "prepoint: " << prepoint / um << G4endl;
  //G4cout << "postpoint: " << postpoint / um << G4endl;

  ////G4double dtctrx = UserDataInput::GetDectorDimensionX();
  ////G4double dtctry = UserDataInput::GetDectorDimensionY();
  ////G4double dtctrz = UserDataInput::GetDectorDimensionZ();
  ////if (abs(postpoint.getZ()) >= dtctrz / 2. || abs(postpoint.getY()) >= dtctry / 2. || abs(postpoint.getX()) >= dtctrx / 2.)
  ////    G4cout << "粒子射出探测器" << G4endl;
  ////else
  ////{
  ////    G4cout << "粒子被吸收" << G4endl;
  ////    getchar();
  ////}

  //if (poststepStatus == fGeomBoundary)
  //{
  //    G4cout << "粒子射出探测器" << G4endl;
  //}
  //else if (poststepStatus != fGeomBoundary && trackStatus == fAlive)
  //{
  //    G4cout << "粒子继续碰撞" << G4endl;
  //}
  //else
  //{
  //    G4cout << "粒子被吸收" << G4endl;
  //    //getchar();
  //}

  // if (strPrtclType ==)
  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit() / keV;

  //
  // 统计次级粒子的能量沉积
  if (parentID == 1)
      EnergyDepositOf2nd[trackID - 2] += edepStep;
  else if (parentID > 1 && parentID <= numof2ndprtclofneutron + 1)
      EnergyDepositOf2nd[parentID - 2] += edepStep;
  else if (parentID > numof2ndprtclofneutron + 1)
  {
      G4int newparentID = FindParentID(trackID, trackIDs);
      while (newparentID > numof2ndprtclofneutron + 1)
      {
          G4int parentID_temp = newparentID;
          //G4cout << "parentID_temp: " << parentID_temp << G4endl;
          newparentID = FindParentID(newparentID, trackIDs);
          if (newparentID == 1)
          {
              //G4cout << "parentID_temp: " << parentID_temp << G4endl;
              //getchar();
              newparentID = parentID_temp;
              break;
          }
      }
      EnergyDepositOf2nd[newparentID - 2] += edepStep;
  }

  //G4cout << "EnergyDepositOf2nd: ";
  //for (int i = 0; i < EnergyDepositOf2nd.size(); i++)
  //{
  //    G4cout << EnergyDepositOf2nd[i] << "  "; 
  //}
  //G4cout << G4endl;

  //统计每个event的能量沉积，即中子反应产生的所有次级粒子的能量沉积的总和
  fEventAction->AddEdep(edepStep);

  //G4cout << "edepStep: " << edepStep << G4endl << G4endl;
  //getchar();

}

void BNSimSteppingAction::Reset()
{
    EnergyDepositOf2nd.clear();
    trackIDs.clear();
}

G4int BNSimSteppingAction::FindParentID(G4int trackID, vector<G4int> trackIDs)
{
    //G4cout << "trackIDs: ";
    //for (int i = 0; i < trackIDs.size(); i++)
    //{
    //    G4cout << trackIDs[i] << "  ";
    //}
    //G4cout << G4endl;
    //G4cout << "find parentID trackID: " << trackID << G4endl;
    //getchar();
    vector <int>::iterator iElement = find(trackIDs.begin(), trackIDs.end(), trackID);
    G4int nPosition = 0;
    if (iElement != trackIDs.end())
    {
        nPosition = distance(trackIDs.begin(), iElement);
        return trackIDs[nPosition + 1];
    }
    else
    {
        G4cout << "Unable to find trackID: " << trackID << G4endl;
        G4cout << "trackIDs: ";
        vector<int>::iterator itr = trackIDs.begin();
        while (itr != trackIDs.end())
        {
            G4cout << *itr << "  ";
            itr++;
        }
        G4cout << G4endl;
        //getchar();
        return EnergyDepositOf2nd.size() + 1;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

