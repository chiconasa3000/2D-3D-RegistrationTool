#include "scriptBuilder.h"

ScriptBuilder::ScriptBuilder(){
	//TODO: To do something :p
}

void ScriptBuilder::setThreshold(int number){
	this->genthreshold = number;
}

void ScriptBuilder::setOriginVolume(string input){
	this->origin_volume = input;
}


void ScriptBuilder::setTargetVolume(string input){
	this->target_volume = input;
}



void ScriptBuilder::setIndexTest(int indexTest){
	this->indexTest = indexTest;
}


void ScriptBuilder::setNumTests(int numTests){
	this->numTests = numTests;
}

void ScriptBuilder::setTipoProy(string tipoProy){
	this->inputTipoProy = tipoProy;
}

void ScriptBuilder::setRotation(float rx, float ry, float rz){
	rotVol[0] = rx; 
	rotVol[1] = ry;
	rotVol[2] = rz;
}

void ScriptBuilder::setScale(float scale){
	this->generalScale = scale;
}

void ScriptBuilder::setTranslation(float tx, float ty, float tz){
	trasVol[0] = tx;
	trasVol[1] = ty;
	trasVol[2] = tz;
}

void ScriptBuilder::setCompareVols(bool flagCompVols){
	this->compareVols = flagCompVols;
}

void ScriptBuilder::setWriteStatistics(bool flagWriteStast){
    this->writestatistics = flagWriteStast;
}

void ScriptBuilder::setRandMode(bool randmode){
	this->randMode = randmode;
}

void ScriptBuilder::asignarScript(string nombreScript){

    if(nombreScript.compare("MultiImageRegistration")==0){
		comman = "./"+nombreScript+" ";
		tipoScript = "MultiImageRegistration";
		cout<<"MultiImageRegistrationActivado"<<endl;
	}else if(nombreScript.compare("CreateImageSetSimilarity")==0){
		comman = "./"+nombreScript+" ";
		tipoScript = "CreateImageSetSimilarity";

	}else if(nombreScript.compare("genVirtualImage")==0){
		comman = "./"+nombreScript+" ";
		tipoScript = "genVirtualImage";

	}else
		std::cerr << "Comando no encontrado" << std::endl;

}

void ScriptBuilder::buildScript(){
	//Creando Clase de ayuda para la comparacion de volumenes
	Utilitarios *util = new Utilitarios();


	if(tipoScript.compare("MultiImageRegistration")==0){
				
		//Modelo 3D a registrar
		string movingImage = "-movingImage " +origin_volume + " ";
		comman += movingImage;

		//Nro de imagenes a registrar
		string numImages = "-numFixedImages 2 ";
		comman += numImages;

		//1er Imagen 2D Fija (DRR o imagen virtual)
		string fixed1Image = "-f0 ../outputData/virtualImages/pelvisHealthy_ap_"+to_string(indexTest)+".mha ";
		comman += fixed1Image;

		//Punto Focal de la 1era Imagen 2D
		string focal1Point = "0 -1000 0 ";
		comman += focal1Point;

		//2da Imagen 2D Fija (DRR o imagen virtual)
		string fixed2Image = "-f1 ../outputData/virtualImages/pelvisHealthy_ml_"+to_string(indexTest)+".mha ";
		comman += fixed2Image;

		//Punto Focal de la 2da Imagen 2D
		string focal2Point = "-1000 0 0 ";
		comman += focal2Point;

		//Tolerancia de la metrica para terminar la optimización
		string stepTolerance = "-steptolerance 0.02 ";
		comman += stepTolerance;

		//Tamanio de Paso
		string stepSize = "-stepsize 6.0 ";
		comman += stepSize;

		//Nro de Niveles de Resolucion y 
		//sus respectivos factores de escala en cada nivel de resolución
		int nroLevels = 4;
		string numLevels = "-numLevels "+std::to_string(nroLevels)+" ";
		comman += numLevels;

		string schedule = "-schedule 4 3 2 1 ";
		comman += schedule;

		//TODO: Create Directory for every test
		//Directorio de Salida de los resultados del registro
		string strDir = "../outputData/resultsReg_"+to_string(indexTest);
		string outputDir ="-outputDirectory " + strDir + " ";
		comman += outputDir;

		//Construimos un archivo que almacena todo el stream del comando ejecutado
		string nameLogRegistro = "LogMultiImageRegistration_"+to_string(indexTest);

		//Nombre del Archivo de Registro
		string logFileName = "-logfilename " + strDir +"/"+ nameLogRegistro + " ";
		comman += logFileName;	

		string activeNewVol = "-writeFinalVol ";
		comman += activeNewVol;

		//string activeStatistics = "-writeStatistics ";
		//comman += activeStatistics;


		//Para conseguir el stream cuando ejecutamos el comando que hemos construido
		GetStdoutFromCommand(comman);
		string fileNuevoVolumen = "../outputData/resultsReg_";
		string fileDeforVolumen = "../outputData/ImagesDefs/Images/";
		string logFileNameTest = strDir + "/" + nameLogRegistro;
		//El volumen reconstruido, la distancia de housdorff y las estadisticas
		//seran activadas ya que requiere una carga adicional para el registro
		if(compareVols){
			util->compareVols(logFileNameTest, fileNuevoVolumen, fileDeforVolumen, indexTest);
		}
		if(writestatistics){
			util->createStats(nroLevels, logFileNameTest, fileNuevoVolumen, indexTest);
		}
		
	}else if(tipoScript.compare("CreateImageSetSimilarity")==0){
		//Activar modo Verbose
		comman += "-v ";

		if(randMode){
			//Mode Random
			comman += "-rnd ";
			comman += "-rnd_sem ";
		}
		//Directorio de Salida
		comman += "-folderName ";
		//TODO: Create Directory for every test
		//Directorio de Salida de los resultados del registro
		string outputDir ="../outputData/ImagesDefs ";
		comman += outputDir;
	
		//Numero de Imagenes a generar
		string numImagesToGenerate = to_string(numTests);
		comman += "-numImages " + numImagesToGenerate + " ";
	
		//Transformacion del modelo
		comman += "-rx " + to_string(rotVol[0]) + " -ry " + to_string(rotVol[1]) +  " -rz " + to_string(rotVol[2])+" ";
		comman += "-t "+ to_string(trasVol[0]) +" "+to_string(trasVol[1])+" "+to_string(trasVol[2])+" ";	
		comman += "-sg "+ to_string(generalScale)+" ";
		
		//Modelo 3D a registrar
        	string movingImage = "-inputVol " + target_volume + " ";
		comman += movingImage;

		//Construimos un archivo que almacena todo el stream del comando ejecutado
		string nameLogRegistro;
		string cabezera = "LogCreateDefImageWithSimilarity_";
		nameLogRegistro += cabezera;
		
		replace(outputDir.begin(), outputDir.end(), ' ', '/');

		//nameLogRegistro += ".txt";
		
		//Asignando el nombre del archivo log
		string logfilename = "-logFileName " + nameLogRegistro;
		comman += logfilename;		
		//Escribiendo el archivo con el stream del comando ejecutado
		//ofstream out(outputDir + nameLogRegistro);
		//out << outputTextRegistration << endl;
		//out.close();
	
		//Ejecución de POPEN
		//Para conseguir el stream cuando ejecutamos el comando que hemos construido
		GetStdoutFromCommand(comman);

		


	}else if(tipoScript.compare("genVirtualImage")==0){
		//Volumen de Entrada
		string volEntrada = "../outputData/ImagesDefs/Images/imagenDef_"+to_string(indexTest)+".mha ";

		//Conseguir Tamanio y Resolucion de la imagen Deformada
		//util->getSizeAndSpacingFromImage(targeVolume);

		//Modo Verbose activado
		string verboseOn = "-v ";
		comman += verboseOn;
	
		//Tipo de Proyeccion
		string tipoProy = "-p "+inputTipoProy+" ";
		comman += tipoProy;

		//Cosenos Directores
		string directCosine = "";

		//Posicion de Punto Focal
		string puntoFocal = "";

		//Tamanio de la Image Virtual
		string tamanio = "";
		//Resolucion de la Imagen Virtual
		string resolucionImagen = "";

		//Nombre de la Imagen Virtual
		string nameVirtualImage = "";


		//Separacion por el tipo de proyeccion
		if(inputTipoProy.compare("AP")==0){
			//Cosenos Directores
			directCosine = "-dc 90 0 0 ";
			comman += directCosine;

			//Posicion de Punto Focal
			puntoFocal = "-foc 0 -1000 0 ";
			comman += puntoFocal;

			//Tamanio de la Image Virtual
           	 	//tamanio = "-size 334 214 ";
			//comman += tamanio;
			
			//Resolucion de la Imagen Virtual
			//resolucionImagen = "-res 1 1 ";
			//comman += resolucionImagen;

			//Nombre de la Imagen Virtual
			nameVirtualImage = "-o pelvisHealthy_ap_"+to_string(indexTest)+" ";
			comman += nameVirtualImage;
			
			//Distancia de Fuente a Isocentro
			string sourceToIsocenterDistance = "-scd -200 ";
			comman += sourceToIsocenterDistance;



		}else if(inputTipoProy.compare("ML")==0){
			//Cosenos Directores
			directCosine = "-dc 90 90 0 ";
			comman += directCosine;

			//Posicion de Punto Focal
			puntoFocal = "-foc -1000 0 0 ";
			comman += puntoFocal;

			//Tamanio de la Image Virtual
            		//tamanio = "-size 179 214 ";
			//comman += tamanio;
			
			//Resolucion de la Imagen Virtual
			//resolucionImagen = "-res 1 1 ";
			//comman += resolucionImagen;

			//Nombre de la Imagen Virtual
			nameVirtualImage = "-o pelvisHealthy_ml_"+to_string(indexTest)+" ";
			comman += nameVirtualImage;

			//Distancia de Fuente a Isocentro
			string sourceToIsocenterDistance = "-scd 200 ";
			comman += sourceToIsocenterDistance;

		}

		//Umbral de Hounsfield
		string threshold = "-threshold " + std::to_string(this->genthreshold) + " ";
		comman += threshold;

		//Volumen de Entrada
		string volumenToProject = "-inputVol " + volEntrada;
		comman += volumenToProject;


		//TODO: Create Directory for every test
		//Directorio de Salida de los resultados del registro
		string outputDir ="../outputData/virtualImages ";
	
		
		//Construimos un archivo que almacena todo el stream del comando ejecutado
        	string nameLogRegistro = "LogVirtualImages_"+inputTipoProy+"_"+to_string(indexTest);
        	nameLogRegistro += ".txt";

		replace(outputDir.begin(), outputDir.end(), ' ', '/');
	
		string logfilename = "-logFileName " + outputDir + nameLogRegistro;	
		comman += logfilename;
		//Escribiendo el archivo con el stream del comando ejecutado
		//ofstream out(outputDir + nameLogRegistro);
		//out << outputTextRegistration << endl;
		//out.close();
	
		//Ejecución del comando
		//Para conseguir el stream cuando ejecutamos el comando que hemos construido
		GetStdoutFromCommand(comman);

		
	}else
		std::cerr << "Comando no encontrado" << std::endl;
}

void ScriptBuilder::GetStdoutFromCommand(string command){

	/*string data;
	FILE * stream;
	const int max_buffer = 300;
	char buffer[max_buffer];
	command.append(" 2>&1");

	//popen es el encargado de ejecutar el comando y devolver el stream
	//de la salida del comandoen este caso es en modo lectura
	stream = popen(command.c_str(),"r");
	if(stream){
		while(!feof(stream))
			if(fgets(buffer, max_buffer, stream) != NULL){
				//guardamos cada linea que imprime la salida del comando (solo para visualizar)
				std::cout<<buffer<<endl;
				//salvamos el stream actual (al final tendremos el stream completo de la salida del comando)
				data.append(buffer);
			}
		pclose(stream);
	}*/

	std::system(command.c_str());	
	//return data;

}
