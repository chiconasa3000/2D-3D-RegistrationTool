#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
using namespace std;

class ScriptBuilder{

private:
	string tipoScript = "";
	string inputTipoProy = "";
	string comman = "";
	int numTests = 1;
	int indexTest = 1;
	 
public:
	void setIndexTest(int indexTest);
	void setNumTests(int numTests);
	void setTipoProy(string tipoProy);
	void asignarScript(string nombreScript);
	void buildScript();
	ScriptBuilder();

	string GetStdoutFromCommand(string command);

};
