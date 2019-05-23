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
#include <iomanip>
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
	
	int threshold = 0;
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

	//Bandera de Escritura de Estadisticas
	bool writestatistics = false;
	
	//Bandera de modo aleatorio de parametros iniciales
	bool randmode = false;

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

		if((ok == false) && (strcmp(argv[1], "-threshold") == 0))
		{
			argc--; argv++;
			ok = true;
			threshold = atoi(argv[1]); 
			argc--; argv++;
		}


		if((ok == false) && (strcmp(argv[1], "-writeStatistics") == 0))
		{
			argc--; argv++;
			ok = true;
			writestatistics = true;
		}
		if((ok == false) && (strcmp(argv[1], "-randMode") == 0))
		{
			argc--; argv++;
			ok = true;
			randmode = true;
		}

	}

	ScriptBuilder *scriptbuilder = new ScriptBuilder();
	scriptbuilder->setNumTests(numImagenes);
	scriptbuilder->setOriginVolume(origin_volume);
	scriptbuilder->setTargetVolume(target_volume);
	scriptbuilder->setRotation(rx,ry,rz);
	scriptbuilder->setTranslation(tx,ty,tz);
	scriptbuilder->setScale(sg);
	scriptbuilder->setThreshold(threshold);
	//Asignacion de banderas para comparacion de volumenes y estadisticas
	if(comparevolumes){
		scriptbuilder->setCompareVols(true);
	}
	if(writestatistics){
		scriptbuilder->setWriteStatistics(true);
	}
	if(randmode){
		scriptbuilder->setRandMode(true);
	}
	//Creacion de Imagenes Deformadas
	scriptbuilder->asignarScript("CreateImageSetSimilarity");	
	scriptbuilder->buildScript();	


	//Archivo log para los valores de distancia de hausdorff
	ofstream hausdorffDistances;
	hausdorffDistances.open("../outputData/HausdorffDistances.txt");
	//Archivo log para los errores de parametros de transformacion	
	ofstream myfile;
	myfile.open("../outputData/RMSE_Registro.txt");	
	
	//utilitarios para la creacion de estadisticas	
	Utilitarios *utils = new Utilitarios();

	//Buffer para acumular etiquetas de errores
	std::string srx="", sry="", srz="", stx="", sty="", stz="", ssg="",csrx="";	
	std::string etq="";
	//Cabezeras iniciales en el log de errores
	etq += "#Level\t#Test";	
	srx += "Rotx\t1";
	sry += "Roty\t2";
	srz += "Rotz\t3";
	stx += "Trax\t4";
	sty += "Tray\t5";
	stz += "Traz\t6";
	ssg += "Sca\t7";
	csrx += "Rotx\t1";


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
		scriptbuilder->buildScript();
		
		//Lectura de las Distancias de Hausdorff de cada registro dentro del archivo Log
		
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
		
		float gt_rx, gt_ry, gt_rz, gt_tx, gt_ty, gt_tz, gt_sg, rg_rx, rg_ry, rg_rz, rg_tx, rg_ty, rg_tz, rg_sg;	
		float ngt_rx,ngt_ry,ngt_rz, nrg_rx, nrg_ry, nrg_rz;	
		//Captura de valores de transformacion
		gt_rx = baseTransform_2->GetParameters()[0];
		gt_ry = baseTransform_2->GetParameters()[1];
		gt_rz = baseTransform_2->GetParameters()[2];
		gt_tx = baseTransform_2->GetParameters()[3];
		gt_ty = baseTransform_2->GetParameters()[4];
		gt_tz = baseTransform_2->GetParameters()[5];
		gt_sg = baseTransform_2->GetParameters()[6];

		//cast from versor to euler	
		utils->convertVersorToEuler(gt_rx, gt_ry, gt_rz, ngt_rx, ngt_ry, ngt_rz);		
		rg_rx = baseTransform->GetParameters()[0];
		rg_ry = baseTransform->GetParameters()[1];
		rg_rz = baseTransform->GetParameters()[2];
		rg_tx = baseTransform->GetParameters()[3];
		rg_ty = baseTransform->GetParameters()[4];
		rg_tz = baseTransform->GetParameters()[5];
		rg_sg = baseTransform->GetParameters()[6];
	
		utils->convertVersorToEuler(rg_rx, rg_ry, rg_rz, nrg_rx, nrg_ry, nrg_rz);
	
		//Archivo log de cada test para almacenar los valores de los parametros de transformacion el antes y el despues
		string filenameTransValue = "../outputData/resultsReg_" + cindex + "/valueTransf" + cindex + ".txt";
		ofstream fileTransValue;
		fileTransValue.open(filenameTransValue);
		
		//Escritura de los valores de los parametros de transformacion
		fileTransValue << "#\tGroundTruth\tRegistration"<<std::endl;
		fileTransValue << "Rx\t" + to_string(ngt_rx) + "\t" + to_string(nrg_rx) << std::endl;
		fileTransValue << "Ry\t" + to_string(ngt_ry) + "\t" + to_string(nrg_ry) << std::endl;
		fileTransValue << "Rz\t" + to_string(ngt_rz) + "\t" + to_string(nrg_rz) << std::endl;
		fileTransValue << "Tx\t" + to_string(gt_tx) + "\t" + to_string(rg_tx) << std::endl;
		fileTransValue << "Ty\t" + to_string(gt_ty) + "\t" + to_string(rg_ty) << std::endl;
		fileTransValue << "Tz\t" + to_string(gt_tz) + "\t" + to_string(rg_tz) << std::endl;
		fileTransValue << "Sg\t" + to_string(gt_sg) + "\t" + to_string(rg_sg) << std::endl;
		fileTransValue.close();

		float t_rx, t_ry, t_rz, t_tx, t_ty, t_tz, t_sg; 		

		t_rx = abs(ngt_rx - nrg_rx);
		t_ry = abs(ngt_ry - nrg_ry);
		t_rz = abs(ngt_rz - nrg_rz);

		t_tx = abs(gt_tx - rg_tx);
		t_ty = abs(gt_ty - rg_ty);
		t_tz = abs(gt_tz - rg_tz);

		t_sg = abs(gt_sg - rg_sg);
		
		rx_error += pow(t_rx,2.0);
		ry_error += pow(t_ry,2.0);
		rz_error += pow(t_rz,2.0);
		tx_error += pow(t_tx,2.0);
		ty_error += pow(t_ty,2.0);
		tz_error += pow(t_tz,2.0);
		sg_error += pow(t_sg,2.0);
		
		etq  = etq + "\t" + "#Error" + cindex;		
		srx += "\t" + std::to_string(t_rx);	
		sry += "\t" + std::to_string(t_ry);		
		srz += "\t" + std::to_string(t_rz);
		stx += "\t" + std::to_string(t_tx);
		sty += "\t" + std::to_string(t_ty);
		stz += "\t" + std::to_string(t_tz);
		ssg += "\t" + std::to_string(t_sg);		
		csrx += "\t" + std::to_string(t_rx);	
		
		//Generar las graficas de valores de transformacion por cada prueba
		std::string dirres = "../outputData/resultsReg_"+to_string(currentIndexTest) + "/"; 
		utils->createStatsOfTransValues(dirres,filenameTransValue,currentIndexTest);
	
	}

	hausdorffDistanceAcum /= numImagenes;
	hausdorffDistances << "AverageTestHausdorffDistance: "<< hausdorffDistanceAcum;
	
	etq = etq + "\t" + "#FinalError";		
	srx += "\t" + std::to_string(sqrt(rx_error/numImagenes));	
	sry += "\t" + std::to_string(sqrt(ry_error/numImagenes));		
	srz += "\t" + std::to_string(sqrt(rz_error/numImagenes));
	stx += "\t" + std::to_string(sqrt(tx_error/numImagenes));
	sty += "\t" + std::to_string(sqrt(ty_error/numImagenes));
	stz += "\t" + std::to_string(sqrt(tz_error/numImagenes));
	ssg += "\t" + std::to_string(sqrt(sg_error/numImagenes));		
	csrx += "\t" + std::to_string(sqrt(rx_error/numImagenes));

	myfile << etq << std::endl;
	myfile << srx << std::endl;
	myfile << sry << std::endl;
	myfile << srz << std::endl;
	myfile << stx << std::endl;
	myfile << sty << std::endl;
	myfile << stz << std::endl;
	myfile << ssg << std::endl;
	myfile << csrx << std::endl;

	cputimer.Stop();

	//myfile << "Registrations took " << cputimer.GetMean() << " seconds.\n" << std::endl;
	hausdorffDistances.close();
	myfile.close();
	//Una vez creado el archivo de errores Podemos llamar al plot
	utils->createStatsBarHausdorff();
	utils->createStatsBoxPlotsTypeTransParams();
	utils->createStatsOfErrors(numImagenes);
	return 0;
}
