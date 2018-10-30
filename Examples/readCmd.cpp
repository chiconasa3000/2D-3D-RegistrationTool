#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <sstream>

#include "subprocess.hpp"
using namespace std;
namespace sp = subprocess;



void GetStdoutFromCommand(){
    auto p = sp::Popen({"./../build-Cpu-midas-journal-800-MyCmake-Default/MultiImageRegistration", "../Cpu-midas-journal-800/bestData/settingRai/pelvisSegmIntensityRai.mha", "2", "../Cpu-midas-journal-800/bestData/settingRai/brokenAp_r0t0_iso225145110_res0.mha", "0","-1990", "0", "../Cpu-midas-journal-800/bestData/settingRai/brokenLt_r0t0_iso225145110_res0.mha", "-2067.5", "0", "0", "0.01", "2.0", "4", "6", "4", "2", "1", "../Cpu-midas-journal-800/bestData/outDirNewUmbral"});
    //auto p = sp::Popen({"./MultiImageRegistration", "../bestData/settingRai/pelvisSegmIntensityRai.mha", "2", "../bestData/settingRai/brokenAp_r0t0_iso225145110_res0.mha", "0","-1990", "0", "../bestData/settingRai/brokenLt_r0t0_iso225145110_res0.mha", "-2067.5", "0", "0", "0.01", "2.0", "4", "6", "4", "2", "1", "../bestData/outDirNewUmbral"});
    auto obuf = p.communicate().first;
    std::cout << "Data : " << obuf.buf.data() << std::endl;
    //return "";
    
    //auto obuf = sp::check_output({"./../build-Cpu-midas-journal-800-MyCmake-Default/MultiImageRegistration", "../Cpu-midas-journal-800/bestData/settingRai/pelvisSegmIntensityRai.mha", "2", "../Cpu-midas-journal-800/bestData/settingRai/brokenAp_r0t0_iso225145110_res0.mha", "0","-1990", "0", "../Cpu-midas-journal-800/bestData/settingRai/brokenLt_r0t0_iso225145110_res0.mha", "-2067.5", "0", "0", "0.01", "2.0", "4", "6", "4", "2", "1", "../Cpu-midas-journal-800/bestData/outDirNewUmbral"});
    //std::cout << "Data : " << obuf.buf.data() << std::endl;
    //std::cout << "Data len: " << obuf.length << std::endl;
    //return "";
}

int main(){
    // GetStdoutFromCommand();
    //vector of step tolerances
    
    std::vector<float> stepTolerances={2.0};
    std::vector<float> focalPoint1={0, 1000, 0};
    std::vector<float> focalPoint2={-1000,0,0};
    std::vector<float> schedules={6, 4, 2, 1};
    for(int i=0;i <stepTolerances.size(); i++){
        
        char comman[200];
        
        string registrationCommand = "./../build-Cpu-midas-journal-800-MyCmake-Default/MultiImageRegistration";
        
        string movingImage = "../Cpu-midas-journal-800/bestData/settingRai/volumenAsTorax/newBrokenTwistedRAI.mha";
        
        string numImages = std::to_string(2);
        
        string fixed1Image = "../Cpu-midas-journal-800/bestData/settingRai/volumenAsTorax/newbrokenAP_scd1000_iso226_-253_112.mha";
        
        string fixed2Image = "../Cpu-midas-journal-800/bestData/settingRai/volumenAsTorax/newbrokenLT_scd1000_iso452_126_4.mha";
        
        string stepSize = std::to_string(0.01);
        
        std::string stepTolerance = std::to_string(stepTolerances[i]);
        
        std::vector<string> focalPoint1Str;
        std::vector<string> focalPoint2Str;
        
        for(int i=0; i<3; i++){
            focalPoint1Str.push_back(std::to_string(focalPoint1[i]));
            focalPoint2Str.push_back(std::to_string(focalPoint2[i]));
        }
        
        std::string numberLevels = std::to_string(schedules.size());
        
        std::vector<string> schedulesStr;
        
        for(int i=0; i<schedules.size(); i++){
            schedulesStr.push_back(std::to_string(schedules[i]));
        }
        
        string outputDir ="../Cpu-midas-journal-800/bestData/outDirNewUmbral";
        
        auto p = sp::Popen({registrationCommand.c_str(), movingImage.c_str(), numImages.c_str(), fixed1Image.c_str(), focalPoint1Str[0].c_str(),focalPoint1Str[1].c_str(),focalPoint1Str[2].c_str(),fixed2Image.c_str(), focalPoint2Str[0].c_str(), focalPoint2Str[1].c_str(), focalPoint2Str[2].c_str(), stepSize.c_str(), stepTolerance.c_str(), numberLevels.c_str(), schedulesStr[0].c_str(), schedulesStr[1].c_str(), schedulesStr[2].c_str(), schedulesStr[3].c_str(), outputDir.c_str()});
        //auto p = sp::Popen({"./MultiImageRegistration", "../bestData/settingRai/pelvisSegmIntensityRai.mha", "2", "../bestData/settingRai/brokenAp_r0t0_iso225145110_res0.mha", "0","-1990", "0", "../bestData/settingRai/brokenLt_r0t0_iso225145110_res0.mha", "-2067.5", "0", "0", "0.01", "2.0", "4", "6", "4", "2", "1", "../bestData/outDirNewUmbral"});
        auto obuf = p.communicate().first;
        //std::cout << "Data : " << obuf.buf.data() << std::endl;
        
        //el archivo del log tendra los parametros que fueron usados
        char nameLogRegistro[100];
        string cabezera = "LogRegisterIteration_";
        strcpy(nameLogRegistro, cabezera.c_str());
        
        //numero de imagenes
        strcat(nameLogRegistro, numImages.c_str());
        
        //step tolerance
        strcat(nameLogRegistro, stepTolerance.c_str());
        
        //size step
        strcat(nameLogRegistro, stepSize.c_str());
        
        //schedule
        strcat(nameLogRegistro, numberLevels.c_str());
        strcat(nameLogRegistro, ".txt");
        
        ofstream out(nameLogRegistro);
        out << obuf.buf.data() << endl;
        out.close();
    }
    return 0;
}

