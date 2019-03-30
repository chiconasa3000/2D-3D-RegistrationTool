#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include "scriptBuilder.h"
#include <itkTransformFileReader.h>
#include "utils.h"
using namespace std;


/////////////////////////////////////////////////////////////////////////////
//
// GENERAR UN LOG DE LOS RESULTADOS DEL REGISTRO
//
/////////////////////////////////////////////////////////////////////////////


int main(int argc, char *argv[]){	

	//Valores de entrada para la transformacion inicial
	float rx = 0.0;
	float ry = 0.0;
	float rz = 0.0;

	float tx = 0.0;
	float ty = 0.0;
	float tz = 0.0;

	float sg = 1.0;

	//Error Acumulado Al cuadrado
	float rx_error = 0.0;
	float ry_error = 0.0;
	float rz_error = 0.0;

	float tx_error = 0.0;
	float ty_error = 0.0;
	float tz_error = 0.0;

	float sg_error = 0.0;
				
	//Lectura del numero de imagenes para las pruebas (Imagenes Deformadas)
	int numImagenes = 2;

	//Imagen de Ingreso : InputData
	char *input_volume = NULL;

	//Bandera de lectura satisfactoria
	bool ok;


	//Initialization de variables
	while(argc > 1){
		ok = false;
		if((ok == false) && (strcmp(argv[1], "-rx") == 0)){
			argc--; argv++;
			ok = true;
			rx = atof(argv[1]);
			argc--; argv++;			
		}

		if((ok == false) && (strcmp(argv[1], "-ry") == 0)){
			argc--; argv++;
			ok = true;
			ry = atof(argv[1]);
			argc--; argv++;			
		}

		if((ok == false) && (strcmp(argv[1], "-rz") == 0)){
			argc--; argv++;
			ok = true;
			rz = atof(argv[1]);
			argc--; argv++;			
		}

		if ((ok == false) && (strcmp(argv[1], "-t") == 0))
		{
			argc--; argv++;
			ok = true;
			tx=atof(argv[1]);
			argc--; argv++;
			ty=atof(argv[1]);
			argc--; argv++;
			tz=atof(argv[1]);
			argc--; argv++;
		}


		if((ok == false) && (strcmp(argv[1], "-sg") == 0)){
			argc--; argv++;
			ok = true;
			sg = atof(argv[1]);
			argc--; argv++;			
		}



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
	scriptbuilder->setInputVolume(input_volume);
	scriptbuilder->setRotation(rx,ry,rz);
	scriptbuilder->setTranslation(tx,ty,tz);
	scriptbuilder->setScale(sg);
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
		
		//Formando el nombre del archivo de transformacion
		string currentTransformFile = "../outputData/resultsReg_"+to_string(currentIndexTest)+"/outTransform.txt";
				
		itk::TransformFileReader::Pointer transformReader = itk::TransformFileReader::New();
		transformReader->SetFileName(currentTransformFile);

		std::cout<<"Lectura de Archivo de Transformacion del "<<currentIndexTest<<" registro"<<std::endl;
		try{
			transformReader->Update();
		}catch(itk::ExceptionObject &e)
		{
			std::cout << e.GetDescription() << std::endl;
		}

		itk::TransformFileReader::TransformListType* transformList = transformReader->GetTransformList();	
		itk::TransformFileReader::TransformPointer baseTransform = transformList->front();
		
		std::cout<<"Current Parameters Transform"<<std::endl;
		std::cout<< baseTransform->GetParameters() << std::endl;
		
		//adaptar los euler a su forma versor
		Utilitarios * util = new Utilitarios();
		
		double vx,vy,vz,newangle;
		util->unirVectorWithAngle(rx,ry,rz,vx,vy,vz,newangle);
		
		rx_error += vx - baseTransform->GetParameters()[0];
		ry_error += vy - baseTransform->GetParameters()[1];
		rz_error += vz - baseTransform->GetParameters()[2];

		tx_error += tx - baseTransform->GetParameters()[3];
		ty_error += ty - baseTransform->GetParameters()[4];
		tz_error += tz - baseTransform->GetParameters()[5];

		sg_error += sg - baseTransform->GetParameters()[6];
		std::cout << "newangle versor: " << newangle << std::endl;	
		std::cout << "rx_error: " << rx_error << std::endl;
		std::cout << "ry_error: " << ry_error << std::endl;
		std::cout << "rz_error: " << rz_error << std::endl;

		std::cout << "tx_error: " << tx_error << std::endl;
		std::cout << "ty_error: " << ty_error << std::endl;
		std::cout << "tz_error: " << tz_error << std::endl;
 		
		std::cout << "sg_error: " << sg_error << std::endl;
		
		

	}

	return 0;
}
