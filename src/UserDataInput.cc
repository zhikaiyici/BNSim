
#include "UserDataInput.hh"

#include "G4SystemOfUnits.hh"

#include <fstream>
#include <iomanip>

using namespace std;

G4int UserDataInput::numberofevents = 1000;

G4String UserDataInput::sourse = "Point";

G4double UserDataInput::stod = 20.;
G4double UserDataInput::dtctrx = 2.;
G4double UserDataInput::dtctry = 2.;
G4double UserDataInput::dtctrz = 50.;
//G4double UserDataInput::StartThickness = 5.;
//G4double UserDataInput::ThicknessIncrement = 5.;
//G4double UserDataInput::EndThickness = 50.;

UserDataInput::UserDataInput()
{
}

UserDataInput::~UserDataInput()
{
}

void UserDataInput::ReadInputData()
{
	G4String prmtrfl = "InputSet.txt";
	ifstream PrameterFile(prmtrfl, ios_base::in);		//输入文件
	if (PrameterFile.good() != 1)
		G4cout << "FAILED to open file " << prmtrfl << ". ERROR!" << G4endl;
	else
	{
		G4cout << "File " << prmtrfl << " has been opened SUCCESSFULLY." << G4endl;

		G4String input;

		while (PrameterFile >> input)			//读取参数
		{
			if (input == "Number_Of_Events:")
			{
				PrameterFile >> numberofevents;
				G4cout <<  numberofevents << " events will be simulated." << G4endl;
			}
			else if (input == "Sourse_Type:")
			{
				PrameterFile >> sourse;
				G4cout << "The sourse tpye is " << sourse << " sourse." << G4endl;
			}
			else if (input == "Distance_of_SD:")
			{
				PrameterFile >> stod;
				G4cout << "The distance of SD is " << stod << " cm." << G4endl;
			}
			else if (input == "XYZ:")
			{
				PrameterFile >> dtctrx >> dtctry >> dtctrz;
				G4cout << "The dimension of detector is " << dtctrx << " mm × " << dtctry << " mm × " << dtctrz << " um" << G4endl;
			}
			//else if (input == "Start_Thickness:")
			//{
			//	PrameterFile >> StartThickness;
			//	G4cout << "" << StartThickness << " um" << G4endl;
			//}
			//else if (input == "Increment_Of_Thickness:")
			//{
			//	PrameterFile >> ThicknessIncrement;
			//	G4cout << "" << ThicknessIncrement << " um" << G4endl;
			//}
			//else if (input == "End_Thickness:")
			//{
			//	PrameterFile >> EndThickness;
			//	G4cout << "" << EndThickness << " um" << G4endl;
			//}
		}
		PrameterFile.close();
	}
}
