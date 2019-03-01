#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include "scriptBuilder.h"
using namespace std;


/////////////////////////////////////////////////////////////////////////////
//
// GENERAR UN LOG DE LOS RESULTADOS DEL REGISTRO
//
/////////////////////////////////////////////////////////////////////////////


int main(int argc, char *argv[]){	

	//Lectura del numero de imagenes para las pruebas (Imagenes Deformadas)
	int numImagenes = 2;

	//Imagen de Ingreso : InputData
	char *input_volume = NULL;

	//Bandera de lectura satisfactoria
	bool ok;


	//Initialization de variables
	while(argc > 1){
		ok = false;
		if((ok == false) && (strcmp(argv[1], "-numImag") == 0))
		{
			argc--; argv++;
			ok = true;
			numImagenes = atoi(argv[1]);
			argc--; argv++;
		}
		if(ok == false)
		{
			if(input_volume == NULL){
				input_volume = argv[1];
				argc--;
				argv++;
				std::cerr << "Imagen Leida"<< std::endl;
			}else
				std::cerr << "Error: No se puede leer la imagen de entrada" << argv[1] << std::endl;
		}

	}

	ScriptBuilder *scriptbuilder = new ScriptBuilder();
	scriptbuilder->setNumTests(numImagenes);


	scriptbuilder->asignarScript("CreateImageSetSimilarity");	
	scriptbuilder->buildScript();	


	//Recorrer el numero de pruebas
	for(int currentIndexTest = 0; currentIndexTest < numImagenes; currentIndexTest++){

		//Actualizando el actual index para seleccionar una especifica
		//imagen deformada y generar especificas imagenes virtuales	
		scriptbuilder->setIndexTest(currentIndexTest);

		//Generacion de Imagen AP
		scriptbuilder->asignarScript("genVirtualImage");
		scriptbuilder->setTipoProy("AP");
		scriptbuilder->buildScript();

		//Generacion de Imagen LT
		scriptbuilder->asignarScript("genVirtualImage");
		scriptbuilder->setTipoProy("ML");
		scriptbuilder->buildScript();
		
		//Generar el registro de imagenes
		scriptbuilder->asignarScript("MultiImageRegistration");
		scriptbuilder->buildScript();

		//Leer Transformacion de Salida
		
		//No haremos comprobacion de Tipo de Transformacion
		//ya que sabemos que es de Similaridad

		//No haremos comprobacion de multiplicidad de transformaciones
		//ya que sabemos que una sola transformacion es guardada en el archivo
		itk::TransformFileReader::Pointer transformReader = itk::TransformFileReader::New();
		transformReader->SetFileName(argv[0]);
		
		std::cout<<"Lectura de Archivo de Transformacion del "<<indexTest<<" registro"<<std::endl;
		try{
			transformReader->Update();
		}catch(itk::ExceptionObject &e)
		{
			std::cout << e.GesDescription() << std:endl;
		}
		itk::TransformFileReader::TransformListType* transformList = transformReader->GetTransformList();	
		itk::TransformFileReaer::TransformPointer baseTransform = transformList->front();
	}

	return 0;
}
