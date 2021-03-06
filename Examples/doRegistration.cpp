#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <sstream>

using namespace std;


/////////////////////////////////////////////////////////////////////////////
//
// GENERAR UN LOG DE LOS RESULTADOS DEL REGISTRO
//
/////////////////////////////////////////////////////////////////////////////


//GetStdoutFromCommand perimite conseguir todo los impreso
//en consola al ejecutarse el comando dado por cmd
string GetStdoutFromCommand(string cmd){

	string data;
	FILE * stream;
	const int max_buffer = 256;
	char buffer[max_buffer];
	cmd.append(" 2>&1");

	//popen es el encargado de ejecutar el comando y devolver el stream
	//de la salida del comandoen este caso es en modo lectura
	stream = popen(cmd.c_str(),"r");
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

int main(){
	//vector of tamaños de paso en la optimizacion
	std::vector<float> stepSize={2.0};

	for(int i=0;i <stepSize.size(); i++){
		
		//Contruimos el comando a ejecutarse guardado en comman
		char comman[200];

		//Comando principal
		string command = "./MultiImageRegistration ";
		strcpy(comman, command.c_str());
		
		//Modelo 3D a registrar
		string movingImage = "../../Cpu-midas-journal-800/bestData/settingRai/volumenAsTorax/newBrokenTwistedRAI.mha ";
		strcat(comman,movingImage.c_str());
		
		//Nro de imagenes a registrar
		string numImages = "2 ";
		strcat(comman, numImages.c_str());

		//1er Imagen 2D Fija (DRR o imagen virtual)
		string fixed1Image = "../../Cpu-midas-journal-800/bestData/settingRai/volumenAsTorax/newbrokenAP_scd1000_iso226_-253_112.mha ";
		strcat(comman, fixed1Image.c_str());

		//Punto Focal de la 1era Imagen 2D
		string focal1Point = "0 1000 0 ";
		strcat(comman, focal1Point.c_str());

		//2da Imagen 2D Fija (DRR o imagen virtual)
		string fixed2Image = "../../Cpu-midas-journal-800/bestData/settingRai/volumenAsTorax/newbrokenLT_scd1000_iso452_126_4.mha ";
		strcat(comman, fixed2Image.c_str());

		//Punto Focal de la 2da Imagen 2D
		string focal2Point = "-1000 0 0 ";
		strcat(comman, focal2Point.c_str());

		//Tolerancia de la metrica para terminar la optimización
		string stepTolerance = "0.01 ";
		strcat(comman, stepTolerance.c_str());

		//Leer el vector de tamaños de paso
		//Insertar el actual tamaño de paso leído
		//(convertir de texto a valor numérico)
		ostringstream ostr;
		ostr << stepSize[i];
		std::string strSizeValue = ostr.str()+" ";
		strcat(comman, strSizeValue.c_str());
		
		//Nro de Niveles de Resolucion y 
		//sus respectivos factores de escala en cada nivel de resolución
		string schedule = "4 6 4 2 1 ";
		strcat(comman, schedule.c_str());

		//Directorio de Salida de los resultados del registro
		string outputDir ="../../Cpu-midas-journal-800/bestData/outDirNewUmbral";
		strcat(comman, outputDir.c_str());


		//Ejecución de POPEN
		//Para conseguir el stream cuando ejecutamos el comando que hemos construido
		string outputTextRegistration = GetStdoutFromCommand(comman);


		//Construimos un archivo que almacena todo el stream del comando ejecutado
		char nameLogRegistro[100];
		string cabezera = "LogRegisterIteration_";
		strcpy(nameLogRegistro, cabezera.c_str());

		//datos adicionales en el nombre del archivo a LogRegisterIteration
		
		//LogRegisterIteration + numero de imagenes
		replace(numImages.begin(), numImages.end(), ' ', '_');
		strcat(nameLogRegistro, numImages.c_str());

		//LogRegisterIteration + numero de imagenes + step tolerance
		replace(stepTolerance.begin(), stepTolerance.end(), ' ', '_');
		strcat(nameLogRegistro, stepTolerance.c_str());

		//LogRegisterIteration + numero de imagenes + step tolerance + size step
		replace(strSizeValue.begin(), strSizeValue.end(), ' ', '_');
		strcat(nameLogRegistro, strSizeValue.c_str());

		//LogRegisterIteration + numero de imagenes + step tolerance + size step + schedule
		replace(schedule.begin(), schedule.end(), ' ', '_');
		strcat(nameLogRegistro, schedule.c_str());
		strcat(nameLogRegistro, ".txt");

		//Escribiendo el archivo con el stream del comando ejecutado
		ofstream out(nameLogRegistro);
		out << outputTextRegistration << endl;
		out.close();
	}
	return 0;
}
