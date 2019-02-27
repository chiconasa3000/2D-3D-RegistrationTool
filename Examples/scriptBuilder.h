#include <iostream>
#include <string>
class ScriptBuilder{

private:
	std::string tipoScript = "";
	std::string comman = "";
	 
public:
	void asignarScript(string nombreScript);
	void buildScript();
	ScriptBuilder();

	string GetStdoutFromCommand(string command);

};
