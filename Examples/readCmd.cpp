#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <sstream>

using namespace std;

string GetStdoutFromCommand(string cmd){

	string data;
	FILE * stream;
	const int max_buffer = 256;
	char buffer[max_buffer];
	cmd.append(" 2>&1");

	stream = popen(cmd.c_str(),"r");
	if(stream){
		while(!feof(stream))
			if(fgets(buffer, max_buffer, stream) != NULL){
				std::cout<<buffer<<endl;
				data.append(buffer);
			}
			pclose(stream);
	}
	return data;

}

int main(){
	//vector of step tolerances
    std::vector<float> stepSize={2.0};

    for(int i=0;i <stepSize.size(); i++){

		char comman[200];
        string command = "./MultiImageRegistration ";
		strcpy(comman, command.c_str());
        string movingImage = "../Cpu-Midas-Journal-800/bestData/pelvisSegmIntensity.mha ";
		strcat(comman,movingImage.c_str());
		string numImages = "2 ";
		strcat(comman, numImages.c_str());
        string fixed1Image = "../Cpu-Midas-Journal-800/bestData/pelvisDRRG0LspCenter.mha ";
		strcat(comman, fixed1Image.c_str());
		string focal1Point = "0 1990 0 ";
		strcat(comman, focal1Point.c_str());
        string fixed2Image = "../Cpu-Midas-Journal-800/bestData/pelvisDRRG90PsrCenter.mha ";
		strcat(comman, fixed2Image.c_str());
		string focal2Point = "-1990 0 0 ";
		strcat(comman, focal2Point.c_str());
        string stepTolerance = "0.01 ";
        strcat(comman, stepTolerance.c_str());

		ostringstream ostr;
        ostr << stepSize[i];
        std::string strSizeValue = ostr.str()+" ";
        //string stepSize = ostr.str()+" ";
        strcat(comman, strSizeValue.c_str());

		string schedule = "4 6 4 2 1 ";
		strcat(comman, schedule.c_str());
        string outputDir ="../Cpu-Midas-Journal-800/bestData/outDirNewUmbral";
		strcat(comman, outputDir.c_str());

		string outputTextRegistration = GetStdoutFromCommand(comman);
		
		//el archivo del log tendra los parametros que fueron usados
		char nameLogRegistro[100];
		string cabezera = "LogRegisterIteration_";
		strcpy(nameLogRegistro, cabezera.c_str());

		//numero de imagenes
		replace(numImages.begin(), numImages.end(), ' ', '_');
		strcat(nameLogRegistro, numImages.c_str());

		//step tolerance
		replace(stepTolerance.begin(), stepTolerance.end(), ' ', '_');
		strcat(nameLogRegistro, stepTolerance.c_str());

		//size step
        replace(strSizeValue.begin(), strSizeValue.end(), ' ', '_');
        strcat(nameLogRegistro, strSizeValue.c_str());

		//schedule
		replace(schedule.begin(), schedule.end(), ' ', '_');
		strcat(nameLogRegistro, schedule.c_str());
		strcat(nameLogRegistro, ".txt");
		
		ofstream out(nameLogRegistro);
		out << outputTextRegistration << endl;
		out.close();
	}
	return 0;
}
