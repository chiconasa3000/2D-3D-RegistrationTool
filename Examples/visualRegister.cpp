#include <fstream>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include "utils.h"


//Visualizacion de Registro

int main(int argc, char *argv[]){

	//Lectura del archivo
	int indRegistro = 0;
	int numLevels = 0;
	char *movingImage = NULL;

	bool ok = false;	
	while(argc > 1){
		ok = false;
		
		if((ok==false) && (strcmp(argv[1], "-movingImage")==0)){
			argc--; argv++;
			ok = true;
			movingImage = argv[1];
			argc--; argv++;
		}

		if((ok==false) && (strcmp(argv[1], "-indRegistro")==0)){
			argc--; argv++;
			ok = true;
			indRegistro = atoi(argv[1]);
			argc--; argv++;
		}

		if((ok==false) && (strcmp(argv[1], "-numLevels")==0)){
			argc--; argv++;
			ok = true;
			numLevels = atoi(argv[1]);
			argc--; argv++;
		}

	}
	
	//Indice de registro en texto
	std::string strIndReg = std::to_string(indRegistro);

	//Armando la ruta y el nombre del log
	std::string logfilename("");
	
	//Nombre de la ruta de resultados de registro
	std::string dirResultados = "../outputData/resultsReg_" + strIndReg + "/";

	//Ruta de la imagen movible
	//const unsigned int Dimensions = 3;
	
	//typedef itk::Image<short int,Dimensions> FixedImageType;
	//typedef itk::Image<short int, Dimensions> MovingImageType;

	//Lectura de la imagen Movible
	//typedef itk::ImageFileReader<MovingImageType> MovingImageReaderType;  	
	//MovingImageReaderType::Pointer movingReader = MovingImageReaderType::New();
	//movingReader->SetFileName(movingImage);
		
	

	//Creacion del directorio donde se almacenaran los resultados
	std::string outDir = "../outputData/visImages";
	itksys::SystemTools::MakeDirectory(outDir);
	Utilitarios * util = new Utilitarios();
	
	//Recorriendo el nro de niveles
	for(int i=0; i < numLevels; i++){

		//Armando la ruta y el nombre del log
		logfilename = "level" + std::to_string(i) + ".txt";
		std::string currentFileLog(dirResultados + logfilename);

		//Limpiando la primera columna de cada archivo
		std::string cmdClearFirstColLog("awk '{$1=$2=\"\"; print $0}' " + currentFileLog +" > " +dirResultados + "n"+ logfilename);
		std::system(cmdClearFirstColLog.c_str());

		std::ifstream fileParametros(dirResultados + "n" + logfilename);


		float vx,vy,vz,tx,ty,tz,sg,rx,ry,rz;
		//Recorre cada linea del archivo
		std::string lineParametros;
		std::string cmdGenFixedImages="";
		int cont = 0;

		while(std::getline(fileParametros,lineParametros)){
			
			vx = 0; vy=0; vz=0; tx=0; ty=0; tz=0; sg=0; rx=0; ry=0; rz=0;
			
			std::string strcont = std::to_string(cont);
			
			std::istringstream currentParamTransf(lineParametros);
			
			
			if(!(currentParamTransf >> vx >> vy >> vz >> tx >> ty >> tz >> sg))
			{
				std::cout << "Dentro de While"<< std::endl;		

				break;
			}

			//Procesamos los valores actuales de transformacion
			std::cout << vx << " " << vy << " " << vz << " " << tx << " " << ty << " " << tz << " " << sg << std::endl;
			
			//Conversion de versor a grados	
			util->convertVersorToEuler(vx,vy,vz,rx,ry,rz);
	
			//Conversion a str de un flotante
			std::string str_rx = std::to_string(rx);
			std::string str_ry = std::to_string(ry);
			std::string str_rz = std::to_string(rz);
			std::string str_tx = std::to_string(tx);
			std::string str_ty = std::to_string(ty);
			std::string str_tz = std::to_string(tz);
			std::string str_sg = std::to_string(sg);
		
			cmdGenFixedImages = "";
			//Generacion de la primera imagen fija AP
			cmdGenFixedImages = cmdGenFixedImages +  "./genVirtualImage -v -p AP -dc 90 0 0 -foc 0 -1000 0" + " -rx " 
			+ str_rx + " -ry " + str_ry + " -rz " + str_rz + " -t " + str_tx + " " + str_ty + " " + str_tz + " " 
			+ " -sg " + str_sg + 
			" -o pelvisHealthy_ap_vis_" + strcont + ".mha " + "-scd -200 -threshold 0 -inputVol " + movingImage + 
			+ " -logFileName ../outputData/" + "visImages/log_vis_ap"+ strcont + " ";	
			
			std::cout << cmdGenFixedImages << std::endl;
			//std::system(cmdGenFixedImages.c_str());

			//Generacion de la segunda imagen fija ML				
				

			cont++;
		}
	}
}

