#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include "scriptBuilder.h"
#include <itkTransformFileReader.h>
#include "utils.h"
#include <fstream>
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


	ofstream myfile;
	myfile.open("RMSE_Registro.txt");	

	//Recorrer el numero de pruebas
	for(int currentIndexTest = 0; currentIndexTest < numImagenes; currentIndexTest++){
		
		string cindex = to_string(currentIndexTest);

		//Lectura de Archivos de el volumen Transformado Aleatoriamente
		string deformRandomTransform = "../outputData/ImagesDefs/TransformFiles/transfSim_"+ cindex +".txt";

		itk::TransformFileReader::Pointer transformReader_2 = itk::TransformFileReader::New();
		transformReader_2->SetFileName(deformRandomTransform);

		std::cout<<"Lectura de Transformacion de la Deformacion "<<currentIndexTest<<" registro"<<std::endl;
		try{
			transformReader_2->Update();
		}catch(itk::ExceptionObject &e)
		{
			std::cout << e.GetDescription() << std::endl;
		}

		itk::TransformFileReader::TransformListType* transformList_2 = transformReader_2->GetTransformList();	
		itk::TransformFileReader::TransformPointer baseTransform_2 = transformList_2->front();

		std::cout<<"Parameters Deformed Volume Transform"<<std::endl;
		std::cout<< baseTransform_2->GetParameters() << std::endl;


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
		string currentTransformFile = "../outputData/resultsReg_"+ cindex +"/outTransform.txt";

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
		//Utilitarios * util = new Utilitarios();

		//double vx,vy,vz,newangle;
		//util->unirVectorWithAngle(rx,ry,rz,vx,vy,vz,newangle);

		rx_error += pow(baseTransform_2->GetParameters()[0] - baseTransform->GetParameters()[0],2.0);
		ry_error += pow(baseTransform_2->GetParameters()[1] - baseTransform->GetParameters()[1],2.0);
		rz_error += pow(baseTransform_2->GetParameters()[2] - baseTransform->GetParameters()[2],2.0);

		tx_error += pow(baseTransform_2->GetParameters()[3] - baseTransform->GetParameters()[3],2.0);
		ty_error += pow(baseTransform_2->GetParameters()[4] - baseTransform->GetParameters()[4],2.0);
		tz_error += pow(baseTransform_2->GetParameters()[5] - baseTransform->GetParameters()[5],2.0);

		sg_error += pow(baseTransform_2->GetParameters()[6]- baseTransform->GetParameters()[6],2.0);
		
		//std::cout << "newangle versor: " << newangle << std::endl;	
		/*
		std::cout << "rx_error: " << rx_error << std::endl;
		std::cout << "ry_error: " << ry_error << std::endl;
		std::cout << "rz_error: " << rz_error << std::endl;

		std::cout << "tx_error: " << tx_error << std::endl;
		std::cout << "ty_error: " << ty_error << std::endl;
		std::cout << "tz_error: " << tz_error << std::endl;

		std::cout << "sg_error: " << sg_error << std::endl;
		*/
		myfile << "rx_error: " << rx_error << std::endl;
		myfile << "ry_error: " << ry_error << std::endl;
		myfile << "rz_error: " << rz_error << std::endl;

		myfile << "tx_error: " << tx_error << std::endl;
		myfile << "ty_error: " << ty_error << std::endl;
		myfile << "tz_error: " << tz_error << std::endl;

		myfile << "sg_error: " << sg_error << std::endl;




	}
	/*
	std::cout << "Rx_final_Error: " << sqrt(rx_error/numImagenes) << std::endl;
	std::cout << "Ry_final_Error: " << sqrt(ry_error/numImagenes) << std::endl;
	std::cout << "Rz_final_Error: " << sqrt(rz_error/numImagenes) << std::endl;
	std::cout << "Tx_final_Error: " << sqrt(tx_error/numImagenes) << std::endl;
	std::cout << "Ty_final_Error: " << sqrt(ty_error/numImagenes) << std::endl;
	std::cout << "Tz_final_Error: " << sqrt(tz_error/numImagenes) << std::endl;
	std::cout << "Sg_final_Error: " << sqrt(sg_error/numImagenes) << std::endl;
	*/
	myfile << "Rx_final_Error: " << sqrt(rx_error/numImagenes) << std::endl;
	myfile << "Ry_final_Error: " << sqrt(ry_error/numImagenes) << std::endl;
	myfile << "Rz_final_Error: " << sqrt(rz_error/numImagenes) << std::endl;
	myfile << "Tx_final_Error: " << sqrt(tx_error/numImagenes) << std::endl;
	myfile << "Ty_final_Error: " << sqrt(ty_error/numImagenes) << std::endl;
	myfile << "Tz_final_Error: " << sqrt(tz_error/numImagenes) << std::endl;
	myfile << "Sg_final_Error: " << sqrt(sg_error/numImagenes) << std::endl;


	

	return 0;
}
