#include "scriptBuilder.h"

ScriptBuilder::ScriptBuilder(){
	//TODO: To do something :p
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


	if(tipoScript.compare("MultiImageRegistration")==0){
	
		//Modelo 3D a registrar
		string movingImage = origin_volume + " ";
		comman += movingImage;

		//Nro de imagenes a registrar
		string numImages = "2 ";
		comman += numImages;

		//1er Imagen 2D Fija (DRR o imagen virtual)
		string fixed1Image = "../outputData/virtualImages/pelvisHealthy_ap_"+to_string(indexTest)+".mha ";
		comman += fixed1Image;

		//Punto Focal de la 1era Imagen 2D
		string focal1Point = "0 -1000 0 ";
		comman += focal1Point;

		//2da Imagen 2D Fija (DRR o imagen virtual)
		string fixed2Image = "../outputData/virtualImages/pelvisHealthy_ml_"+to_string(indexTest)+".mha ";
		comman += fixed2Image;

		//Punto Focal de la 2da Imagen 2D
		string focal2Point = "-1000 0 0 ";
		comman += focal2Point;

		//Tolerancia de la metrica para terminar la optimización
		string stepTolerance = "0.02 ";
		comman += stepTolerance;

		//Tamanio de Paso
		string stepSize = "4.0 ";
		comman += stepSize;

		//Nro de Niveles de Resolucion y 
		//sus respectivos factores de escala en cada nivel de resolución
		//string schedule = "4 6 4 2 1 ";
		string schedule = "4 6 3 2 1 ";
		comman += schedule;

		//TODO: Create Directory for every test
		//Directorio de Salida de los resultados del registro
		string outputDir ="../outputData/resultsReg_"+to_string(indexTest) + " ";
		comman += outputDir;

		//Ejecución de POPEN
		//Para conseguir el stream cuando ejecutamos el comando que hemos construido
		GetStdoutFromCommand(comman);


		//Construimos un archivo que almacena todo el stream del comando ejecutado
		string nameLogRegistro;
		string cabezera = "LogMultiImageRegistration_"+to_string(indexTest);
		nameLogRegistro += cabezera;

		//datos adicionales en el nombre del archivo a LogRegisterIteration

		//LogRegisterIteration + numero de imagenes
		replace(numImages.begin(), numImages.end(), ' ', '_');
		nameLogRegistro += numImages;

		//LogRegisterIteration + numero de imagenes + step tolerance
		replace(stepTolerance.begin(), stepTolerance.end(), ' ', '_');
		nameLogRegistro += stepTolerance;

		//LogRegisterIteration + numero de imagenes + step tolerance + size step
		replace(stepSize.begin(), stepSize.end(), ' ', '_');
		nameLogRegistro += stepSize;

		//LogRegisterIteration + numero de imagenes + step tolerance + size step + schedule
		replace(schedule.begin(), schedule.end(), ' ', '_');
		nameLogRegistro +=  schedule;
		nameLogRegistro += ".txt";

		//Nombre del Archivo de Registro
		string logFileName = outputDir + "/" + nameLogRegistro + " ";
		comman += logFileName;	
		
		//Escribiendo el archivo con el stream del comando ejecutado
		//ofstream out(outputDir+"/" + nameLogRegistro);
		//out << outputTextRegistration << endl;
		//out.close();


	}else if(tipoScript.compare("CreateImageSetSimilarity")==0){
		//Activar modo Verbose
		comman += "-v ";

		//Mode Random
		comman += "-rnd ";
		comman += "-rnd_sem ";

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
		string movingImage = "-inputVol" + target_volume + " ";
		comman += movingImage;

		//Construimos un archivo que almacena todo el stream del comando ejecutado
		string nameLogRegistro;
		string cabezera = "LogCreateDefImageWithSimilarity_" + to_string(indexTest);
		nameLogRegistro += cabezera;
		
		replace(outputDir.begin(), outputDir.end(), ' ', '/');

		nameLogRegistro += ".txt";
		
		//Asignando el nombre del archivo log
		string logfilename = "-logFileName" + outputDir + nameLogRegistro;
		comman += logfilename;		
		//Escribiendo el archivo con el stream del comando ejecutado
		//ofstream out(outputDir + nameLogRegistro);
		//out << outputTextRegistration << endl;
		//out.close();
		
		//Ejecución de POPEN
		//Para conseguir el stream cuando ejecutamos el comando que hemos construido
		GetStdoutFromCommand(comman);




	}else if(tipoScript.compare("genVirtualImage")==0){
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
			tamanio = "-size 452 225 ";
			comman += tamanio;
			//Resolucion de la Imagen Virtual
			resolucionImagen = "-res 0.6641 1 ";
			comman += resolucionImagen;

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
			tamanio = "-size 253 225 ";
			comman += tamanio;
			//Resolucion de la Imagen Virtual
			resolucionImagen = "-res 0.6641 1 ";
			comman += resolucionImagen;

			//Nombre de la Imagen Virtual
			nameVirtualImage = "-o pelvisHealthy_ml_"+to_string(indexTest)+" ";
			comman += nameVirtualImage;

			//Distancia de Fuente a Isocentro
			string sourceToIsocenterDistance = "-scd 200 ";
			comman += sourceToIsocenterDistance;

		}

		//Umbral de Hounsfield
		string threshold = "-threshold 0 ";
		comman += threshold;

		//Volumen de Entrada
		string volumenToProject = "-inputVol ../outputData/ImagesDefs/Images/imagenDef_"+to_string(indexTest)+".mha ";
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
