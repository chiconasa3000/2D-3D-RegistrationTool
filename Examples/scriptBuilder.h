#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
#include <stdlib.h>
#include "utils.h"

using namespace std;

class ScriptBuilder{

private:
	string origin_volume = "";
	string target_volume = "";
	string tipoScript = "";
	string inputTipoProy = "";
	string comman = "";
	int numTests = 1;
	int indexTest = 1;
	float rotVol[3] = {0.0, 0.0, 0.0};
	float trasVol[3] = {0.0, 0.0, 0.0};
	float generalScale=1.0;
    bool compareVols = false;
    bool writestatistics = false;
	bool randMode = false;
	int genthreshold = 0;
public:
	void setThreshold(int umbral);
	void setIndexTest(int indexTest);
	void setNumTests(int numTests);
	void setTipoProy(string tipoProy);
	void asignarScript(string nombreScript);
	void buildScript();
	ScriptBuilder();
	void setRotation(float rx, float ry, float rz);
	void setTranslation(float tx, float ty, float tz);
	void setScale(float scale);
    void setCompareVols(bool flagCompVols);
    void setWriteStatistics(bool writestatistics);
	void setRandMode(bool randmode);
	void setTargetVolume(string inputVol);	
	void setOriginVolume(string inputVol);	

	void GetStdoutFromCommand(string command);

};
