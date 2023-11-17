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
// $Id: BNSimDetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file BNSimDetectorConstruction.cc
/// \brief Implementation of the BNSimDetectorConstruction class

#include "BNSimDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4VisAttributes.hh"

#include "UserDataInput.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BNSimDetectorConstruction::BNSimDetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BNSimDetectorConstruction::~BNSimDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* BNSimDetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  UserDataInput userdata;

  //     
  // World
  //
  G4double world_sizeX = 1.1 * userdata.GetDectorDimensionX();
  G4double world_sizeY = 1.1 * userdata.GetDectorDimensionY();
  G4double world_sizeZ = 2. * userdata.GetDistanceOfSD() + userdata.GetDectorDimensionZ();
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking 

  //
  // 探测器几何
  //

  // 定义探测器材料hBN
  G4double z, a, density;
  G4String name, symbol;
  G4int numberofelements, numberofatoms;

  G4Element* elB = new G4Element(name = "Boron",
	  symbol = "B",
	  z = 5.0, a = 10.811 * g / mole);
  /*G4Element* elB = new G4Element(name = "Boron",
	  symbol = "B",
	  z = 7.0, a = 14.0067 * g / mole);*/
  G4Element* elN = new G4Element(name = "Nitrogen",
	  symbol = "N",
	  z = 7.0, a = 14.0067 * g / mole);

  //G4int iz, in, numberofisotopes;
  //
  //// 数据来源
  //// a: http://www.nuclear.csdb.cn/texing.html
  //// 丰度：https://www.nndc.bnl.gov/nudat2/
  //
  //G4Isotope* isoB10 = new G4Isotope(name = "B10",
  //    iz = 5, in = 10, a = 10.012937068613 * g / mole);
  //G4Isotope* isoB11 = new G4Isotope(name = "B11",
  //    iz = 5, in = 11, a = 11.0093054826847 * g / mole);
  //G4Element* elB = new G4Element(name = "Boron",
  //    symbol = "B",
  //    numberofisotopes = 2);
  //elB->AddIsotope(isoB10, 19.9 * perCent);
  //elB->AddIsotope(isoB11, 80.1 * perCent);
  //
  //G4Isotope* isoN14 = new G4Isotope(name = "N14",
  //    iz = 7, in = 14, a = 14.0030739869773 * g / mole);
  //G4Isotope* isoN15 = new G4Isotope(name = "N15",
  //    iz = 7, in = 15, a = 15.0001088574001 * g / mole);
  //G4Element* elN = new G4Element(name = "Nitrogen",
  //    symbol = "N",
  //    numberofisotopes = 2);
  //elN->AddIsotope(isoN14, 99.636 * perCent);
  //elN->AddIsotope(isoN15, 0.364 * perCent);

  //G4Element* elB = nist->FindOrBuildElement(z = 5);
  //G4Element* elN = nist->FindOrBuildElement(z = 7);

  G4Material* hBN = new G4Material(name = "h-BN",
	  density = 2.28 * g / cm3,
	  numberofelements = 2);
  hBN->AddElement(elB, numberofatoms = 1);
  hBN->AddElement(elN, numberofatoms = 1);

  G4double dtctrx = userdata.GetDectorDimensionX();
  G4double dtctry = userdata.GetDectorDimensionY();
  G4double dtctrz = userdata.GetDectorDimensionZ();

  G4Box* solidDetector = new G4Box("Detector", 0.5 * dtctrx, 0.5 * dtctry, 0.5 * dtctrz);

  G4LogicalVolume* logicDetector = new G4LogicalVolume(solidDetector, hBN, "Detector"); 

  logicDetector->SetVisAttributes(new G4VisAttributes(G4Color(300. / 255.0, 200. / 255.0, 100. / 255.0)));

  G4VPhysicalVolume* physDetector =
	  new G4PVPlacement(0, G4ThreeVector(), logicDetector, "Detector", logicWorld, false, 0, checkOverlaps);

  // Set Detector as scoring volume 
  fScoringVolume = logicDetector;
  
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
