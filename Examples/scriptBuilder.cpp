#include "scriptBuilder.h"

ScriptBuilder::ScriptBuilder(){
	//TODO: To do something :p
}

void ScriptBuilder::asignarScript(string nombreScript){

	if(strcmp(nombreScript, "MultiImageRegistration")==0){
		comman = "./"+nombreScript+" ";
		tipoScript = "MultiImageRegistration";

	}else if(strcmp(nombreScript, "CreateImageSetSimilarity")==0){
		comman = "./"+nombreScript+" ";
		tipoScript = "CreateImageSetSimilarity";

	}else if(strcmp(nombreScript, "genVirtualImage")==0){
		comman = "./"+nombreScript+" ";
		tipoScript = "genVirtualImage";

	}else
		std::cerr << "Comando no encontrado" << std::endl;

}

void ScriptBuilder::buildScript(){


	if(strcmp(tipoScript, "MultiImageRegistration")==0){

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
		string focal2Point = "-1000 0 0 ";
		comman += focal2Point;

		//Tolerancia de la metrica para terminar la optimización
		string stepTolerance = "0.01 ";
		comman += stepTolerance;

		//Nro de Niveles de Resolucion y 
		//sus respectivos factores de escala en cada nivel de resolución
		string schedule = "4 6 4 2 1 ";
		comman += schedule;

		//TODO: Create Directory for every test
		//Directorio de Salida de los resultados del registro
		string outputDir ="../outputData/resultReg";
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
		nameLogRegistr += stepTolerance;

		//LogRegisterIteration + numero de imagenes + step tolerance + size step
		replace(strSizeValue.begin(), strSizeValue.end(), ' ', '_');
		nameLogRegistro += strSizeValue;

		//LogRegisterIteration + numero de imagenes + step tolerance + size step + schedule
		replace(schedule.begin(), schedule.end(), ' ', '_');
		nameLogRegistro +=  schedule;
		nameLogRegistro += ".txt";

		//Escribiendo el archivo con el stream del comando ejecutado
		ofstream out(nameLogRegistro);
		out << outputTextRegistration << endl;
		out.close();


	}else if(strcmp(tipoScript, "CreateImageSetSimilarity")==0){



	}else if(strcmp(tipoScript, "genVirtualImage")==0){

	}else
		std::cerr << "Comando no encontrado" << std::endl;
}
