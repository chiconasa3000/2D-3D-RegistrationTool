#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include "scriptBuilder.h"
#include <itkTransformFileReader.h>
#include <fstream>
#include "itkTimeProbe.h"

using namespace std;


/////////////////////////////////////////////////////////////////////////////
//
// GENERAR UN LOG DE LOS RESULTADOS DEL REGISTRO
//
/////////////////////////////////////////////////////////////////////////////


int main(int argc, char *argv[]){	
	
	//Acumulador de distancia de hausdorff
	float hausdorffDistanceValue = 0.0;
	float hausdorffDistanceAcum = 0.0;
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

	//Volumen Referencia : InputData
	char *origin_volume = NULL;
	
	//Volumen Movible:
	char *target_volume = NULL;
	
	//Bandera de Comparacion de Volumenes
	bool comparevolumes = false;	 
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
		if((ok == false) && (strcmp(argv[1], "-targetVol") == 0))
		{
			argc--; argv++;
			ok = true;
			target_volume = argv[1];
			argc--; argv++;
		}

		if((ok == false) && (strcmp(argv[1], "-originVol") == 0))
		{
			argc--; argv++;
			ok = true;
			origin_volume = argv[1];
			argc--; argv++;
		}
		if((ok == false) && (strcmp(argv[1], "-compareVols") == 0))
		{
			argc--; argv++;
			ok = true;
			comparevolumes = true;
		}

	}

	ScriptBuilder *scriptbuilder = new ScriptBuilder();
	scriptbuilder->setNumTests(numImagenes);
	scriptbuilder->setOriginVolume(origin_volume);
	scriptbuilder->setTargetVolume(target_volume);
	scriptbuilder->setRotation(rx,ry,rz);
	scriptbuilder->setTranslation(tx,ty,tz);
	scriptbuilder->setScale(sg);
	scriptbuilder->asignarScript("CreateImageSetSimilarity");	
	scriptbuilder->buildScript();	
	

	//Archivo log para los valores de distancia de hausdorff
	ofstream hausdorffDistances;
	hausdorffDistances.open("HausdorffDistances.txt");
	//Archivo log para los errores de parametros de transformacion	
	ofstream myfile;
	myfile.open("RMSE_Registro.txt");	

	itk::TimeProbe cputimer;
	cputimer.Start();
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
		//Comparacion de Volumenes solo en caso que este activado
		if(comparevolumes){
			scriptbuilder->setCompareVols(true);
		}
		scriptbuilder->buildScript();
		
		//Obtener las lineas correspondientes a la distancia de hausdorff
		std::string nameLogRegistro = "LogMultiImageRegistration_"+to_string(currentIndexTest);
		std::string logfilename = "../outputData/resultsReg_"+to_string(currentIndexTest) + "/" + nameLogRegistro; 
		
		FILE * stream; //salida en el cmd
		const int max_buffer = 6; //nro de cifras del valor de hausdorff
		char buffer[max_buffer]; //string que guarda el valor de hausdorff
		std::string cmdDistHauss = "cat " + logfilename + " | tail -n 1 | awk '{print $2}'"; 
		cmdDistHauss.append(" 2>&1");
		//Ejecucion y lectura del comando: Save Hausdorff Distance Value
		stream = popen(cmdDistHauss.c_str(),"r");
		//La lectura es una sola vez
		if(stream && fgets(buffer, max_buffer, stream) != NULL){
			//agregamos el valor al actual acumulador de hausdorff distance
 			hausdorffDistanceValue = atof(buffer);
			hausdorffDistanceAcum += hausdorffDistanceValue;
		}

		//Almacenamos el actual valor de hausdorff en un archivo
		//std::string cmdDistHauss("cat " + logfilename + " | tail -n 1 | awk '{print $2}' >> valuesHaussdorffDistance.txt");
		//std::system(cmdDistHauss.c_str());
		hausdorffDistances << "HaussDist "<< currentIndexTest << " : " << hausdorffDistanceValue << std::endl;	



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
		//double vx,vy,vz,newangle;
		//util->unirVectorWithAngle(rx,ry,rz,vx,vy,vz,newangle);
				
		double t_rx, t_ry, t_rz, t_tx, t_ty, t_tz, t_sg; 		

		t_rx = pow(baseTransform_2->GetParameters()[0] - baseTransform->GetParameters()[0],2.0);
		t_ry = pow(baseTransform_2->GetParameters()[1] - baseTransform->GetParameters()[1],2.0);
		t_rz = pow(baseTransform_2->GetParameters()[2] - baseTransform->GetParameters()[2],2.0);

		t_tx = pow(baseTransform_2->GetParameters()[3] - baseTransform->GetParameters()[3],2.0);
		t_ty = pow(baseTransform_2->GetParameters()[4] - baseTransform->GetParameters()[4],2.0);
		t_tz = pow(baseTransform_2->GetParameters()[5] - baseTransform->GetParameters()[5],2.0);

		t_sg = pow(baseTransform_2->GetParameters()[6]- baseTransform->GetParameters()[6],2.0);
		
		rx_error += t_rx;
		ry_error += t_ry;
		rz_error += t_rz;
		tx_error += t_tx;
		ty_error += t_ty;
		tz_error += t_tz;
		sg_error += t_sg;
			
		myfile << "Registro num " << currentIndexTest << std::endl;	
		myfile << "rx_error: " << t_rx << std::endl;
		myfile << "ry_error: " << t_ry << std::endl;
		myfile << "rz_error: " << t_rz << std::endl;

		myfile << "tx_error: " << t_tx << std::endl;
		myfile << "ty_error: " << t_ty << std::endl;
		myfile << "tz_error: " << t_tz << std::endl;

		myfile << "sg_error: " << t_sg << std::endl;
		myfile << std::endl;
		
	}
	
	hausdorffDistanceAcum /= numImagenes;
	hausdorffDistances << "AverageTestHausdorffDistance: "<< hausdorffDistanceAcum;
	myfile << "Rx_final_Error: " << sqrt(rx_error/numImagenes) << std::endl;
	myfile << "Ry_final_Error: " << sqrt(ry_error/numImagenes) << std::endl;
	myfile << "Rz_final_Error: " << sqrt(rz_error/numImagenes) << std::endl;
	myfile << "Tx_final_Error: " << sqrt(tx_error/numImagenes) << std::endl;
	myfile << "Ty_final_Error: " << sqrt(ty_error/numImagenes) << std::endl;
	myfile << "Tz_final_Error: " << sqrt(tz_error/numImagenes) << std::endl;
	myfile << "Sg_final_Error: " << sqrt(sg_error/numImagenes) << std::endl;

	cputimer.Stop();
	
	myfile << "Registrations took " << cputimer.GetMean() << " seconds.\n" << std::endl;
	hausdorffDistances.close();
	myfile.close();
	return 0;
}
