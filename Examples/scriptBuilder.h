#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>

using namespace std;

class ScriptBuilder{

private:
	string tipoScript = "";
	string comman = "";
	 
public:
	void asignarScript(string nombreScript);
	void buildScript();
	ScriptBuilder();

	string GetStdoutFromCommand(string command);

};
