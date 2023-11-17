#ifndef UserDataInput_h
#define UserDataInput_h 1

#include "globals.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

class UserDataInput
{
public:
	UserDataInput();
	~UserDataInput();

	void static ReadInputData(); //读取外部输入文件

	inline static G4int GetNumberOfEvents() { return numberofevents; }

	inline static G4String GetSourseType() { return sourse; }

	inline static G4double GetDistanceOfSD() { return stod * cm; }
	inline static G4double GetDectorDimensionX() { return dtctrx * mm; }
	inline static G4double GetDectorDimensionY() { return dtctry * mm; }
	inline static G4double GetDectorDimensionZ() { return dtctrz * um; }

private:

	static G4int numberofevents; //要模拟的离子数

	static G4String sourse; //源的类型
	
	static G4double stod; //源到探测器距离
	static G4double dtctrx, dtctry, dtctrz; //探测器的长宽高
	//static G4double StartThickness, ThicknessIncrement, EndThickness;
};

#endif // !UserDataInput_h