#include "scriptBuilder.h"

ScriptBuilder::ScriptBuilder(){
	//TODO: To do something :p
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
		string movingImage = "../inputData/pelvisSegmIntensityLPI.mha ";
		comman += movingImage;

		//Nro de imagenes a registrar
		string numImages = "2 ";
		comman += numImages;

		//1er Imagen 2D Fija (DRR o imagen virtual)
		string fixed1Image = "../outputData/virtualImages/pelvisHealthy_ap.mha ";
		comman += fixed1Image;

		//Punto Focal de la 1era Imagen 2D
		string focal1Point = "0 1000 0 ";
		comman += focal1Point;

		//2da Imagen 2D Fija (DRR o imagen virtual)
		string fixed2Image = 	"../outputData/virtualImages/pelvisHealthy_ml.mha ";
		comman += fixed2Image;

		//Punto Focal de la 2da Imagen 2D
		string focal2Point = "1000 0 0 ";
		comman += focal2Point;

		//Tolerancia de la metrica para terminar la optimización
		string stepTolerance = "0.02 ";
		comman += stepTolerance;

		//Tamanio de Paso
		string stepSize = "4.0 ";
		comman += stepSize;

		//Nro de Niveles de Resolucion y 
		//sus respectivos factores de escala en cada nivel de resolución
		string schedule = "4 6 4 2 1 ";
		comman += schedule;

		//TODO: Create Directory for every test
		//Directorio de Salida de los resultados del registro
		string outputDir ="../outputData/resultsReg";
		comman += outputDir;


		//Ejecución de POPEN
		//Para conseguir el stream cuando ejecutamos el comando que hemos construido
		string outputTextRegistration = GetStdoutFromCommand(comman);


		//Construimos un archivo que almacena todo el stream del comando ejecutado
		string nameLogRegistro;
		string cabezera = "LogRegisterIteration_";
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

		//Escribiendo el archivo con el stream del comando ejecutado
		char * c_nameLogRegistro = &nameLogRegistro[0];
		ofstream out(outputDir + nameLogRegistro);
		out << outputTextRegistration << endl;
		out.close();


	}else if(tipoScript.compare("CreateImageSetSimilarity")==0){



	}else if(tipoScript.compare("genVirtualImage")==0){

	}else
		std::cerr << "Comando no encontrado" << std::endl;
}

string ScriptBuilder::GetStdoutFromCommand(string command){
	
	string data;
	FILE * stream;
	const int max_buffer = 256;
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
	}
	return data;

}
