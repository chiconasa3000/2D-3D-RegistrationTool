#include <fstream>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include "GenDRRs/genVirtGen.h"
#include <itkSubtractImageFilter.h>
#include <itkFlipImageFilter.h>
//Visualizacion de Registro

int main(int argc, char *argv[]){

	//Lectura del archivo
	int indRegistro = 0;
	int numLevels = 0;
	char *movingImage = NULL;
	char *dirResReg = NULL;


	char *output_name = NULL;		//virtual image
	char *volFixedImage = NULL;
	char *type_projection = (char*)"OR";		//tipo de proyeccion [AP(anteroposterior) , ML(mediolateral)]

	// Translation parameter of the isocenter in mm
	float tx = 0.;
	float ty = 0.;
	float tz = 0.;

	//Direction Cosines (rotation for orientation volume)
	float dcx = 0.;
	float dcy = 0.;
	float dcz = 0.;


	// CT volume rotation around isocenter along x,y,z axis in degrees
	float rx = 0.;
	float ry = 0.;
	float rz = 0.;

	//Scale
	float sg = 1.;

	int dxx = 0;				//pixels number virtual image in x
	int dyy = 0;				//pixels number virtual image in y

	float im_sx = 0;			//virtual image spacing in x
	float im_sy = 0;			//virtual image spacing in y

	float focalPointx = 0.;			//focalPoint in x
	float focalPointy = 1000.0;		//focalPoint in y
	float focalPointz = 0.;			//focalPoint in z

	float threshold = 0.;			//virtual image threshold

	bool customized_2DCX = false;		//flag for central of 2d image
	float o2Dx;				//virtual image origin in x
	float o2Dy;				//virtual image origin in y

	float scd = 1000.0;			//distance from source to isocenter



	bool ok = false;	
	while(argc > 1){
		ok = false;
		if((ok ==false) && (strcmp(argv[1], "-dirResultsReg") == 0)){
			argc--; argv++;
			ok == true;
			dirResReg = argv[1];
			argc--; argv++;
		}

		if((ok==false) && (strcmp(argv[1], "-movingImage")==0)){
			argc--; argv++;
			ok = true;
			movingImage = argv[1];
			argc--; argv++;
		}

		if((ok==false) && (strcmp(argv[1], "-indRegistro")==0)){
			argc--; argv++;
			ok = true;
			indRegistro = atoi(argv[1]);
			argc--; argv++;
		}

		if((ok==false) && (strcmp(argv[1], "-numLevels")==0)){
			argc--; argv++;
			ok = true;
			numLevels = atoi(argv[1]);
			argc--; argv++;
		}
		
		if ((ok == false) && (strcmp(argv[1], "-p") == 0))
		{
			argc--; argv++;
			ok = true;
			type_projection = argv[1];
			argc--; argv++;
		}

		if((ok==false) && strcmp(argv[1], "-volFixedImage")==0)
		{
			argc--; argv++;
			ok = true;
			volFixedImage = argv[1];
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

		if ((ok == false) && (strcmp(argv[1], "-dc") == 0))
		{
			argc--; argv++;
			ok = true;
			dcx=atof(argv[1]);
			argc--; argv++;
			dcy=atof(argv[1]);
			argc--; argv++;
			dcz=atof(argv[1]);
			argc--; argv++;
		}

		if ((ok == false) && (strcmp(argv[1], "-rx") == 0))
		{
			argc--; argv++;
			ok = true;
			rx=atof(argv[1]);
			argc--; argv++;
		}

		if ((ok == false) && (strcmp(argv[1], "-ry") == 0))
		{
			argc--; argv++;
			ok = true;
			ry=atof(argv[1]);
			argc--; argv++;
		}

		if ((ok == false) && (strcmp(argv[1], "-rz") == 0))
		{
			argc--; argv++;
			ok = true;
			rz=atof(argv[1]);
			argc--; argv++;
		}

		if ((ok == false) && (strcmp(argv[1], "-sg") == 0))
		{
			argc--; argv++;
			ok = true;
			sg=atof(argv[1]);
			argc--; argv++;
		}

		if ((ok == false) && (strcmp(argv[1], "-foc") == 0))
		{
			argc--; argv++;
			ok = true;
			focalPointx = atof(argv[1]);
			argc--; argv++;
			focalPointy = atof(argv[1]);
			argc--; argv++;
			focalPointz = atof(argv[1]);
			argc--; argv++;
		}

		if ((ok == false) && (strcmp(argv[1], "-threshold") == 0))
		{
			argc--; argv++;
			ok = true;
			threshold=atof(argv[1]);
			argc--; argv++;
		}




		if ((ok == false) && (strcmp(argv[1], "-res") == 0))
		{
			argc--; argv++;
			ok = true;
			im_sx=atof(argv[1]);
			argc--; argv++;
			im_sy=atof(argv[1]);
			argc--; argv++;
		}

		if ((ok == false) && (strcmp(argv[1], "-size") == 0))
		{
			argc--; argv++;
			ok = true;
			dxx = atoi(argv[1]);
			argc--; argv++;
			dyy = atoi(argv[1]);
			argc--; argv++;
		}
		if((ok == false) && (strcmp(argv[1], "-2dcx")==0))
		{
			argc--; argv++;
			ok = true;
			o2Dx = atof(argv[1]);
			argc--; argv++;
			o2Dy = atof(argv[1]);
			argc--; argv++;
			customized_2DCX = true;		
		}
		if ((ok == false) && (strcmp(argv[1], "-scd") == 0))
		{
			argc--; argv++;
			ok = true;
			scd = atof(argv[1]);
			argc--; argv++;
		}
		if ((ok == false) && (strcmp(argv[1], "-o") == 0))
		{
			argc--; argv++;
			ok = true;
			output_name = argv[1];
			argc--; argv++;
		}

	}
	
	//Indice de registro en texto
	std::string strIndReg = std::to_string(indRegistro);

	//Armando la ruta y el nombre del log
	std::string logfilename("");
	
	//Nombre de la ruta de resultados de registro
	std::string dirResRegStr = dirResReg;
	std::string dirResultados = dirResRegStr + "resultsReg_" + strIndReg + "/";

	Utilitarios * util = new Utilitarios();

	const unsigned int Dimension = 3;
	//tipo de pixel by default para las imagenes
	typedef short int InputPixelType;
	typedef unsigned char OutputPixelType;
	typedef itk::Image<InputPixelType, Dimension> FixedImageType;
	typedef itk::Image<InputPixelType, Dimension> MovingImageType;
	typedef itk::Image<OutputPixelType, 3> OutputImageType;
	MovingImageType::Pointer imageVolTemp;
	MovingImageType::Pointer imageVolRef;
	
	//Tipo para la substraccion
	typedef itk::SubtractImageFilter<MovingImageType,MovingImageType, MovingImageType>  SubtracterType;
	//Tipo para el rescalado de contraste
	typedef itk::RescaleIntensityImageFilter< MovingImageType, OutputImageType > RescaleFilterType;
	//Tipo para escritura de diferencia
	typedef itk::ImageFileWriter<OutputImageType> WriterType;
	
    	//el directorio donde estaran los resultados de las diferencias
    	itksys::SystemTools::MakeDirectory("../outputData/difimages");
	itksys::SystemTools::MakeDirectory("../outputData/logimages");


	GenVirtGen *genVolRef = new GenVirtGen(type_projection,"../outputData/logImages/logVolRef.txt");
	genVolRef->readMovingImage(volFixedImage);	
	genVolRef->initTransform(tx,ty,tz,dcx,dcy,dcz,rx,ry,rz,sg);
	genVolRef->initInterpolator(focalPointx,focalPointy,focalPointz,100);
	genVolRef->initResampleFilter(dxx, dyy, im_sx,im_sy,customized_2DCX,o2Dx,o2Dy,scd);
	
	std::string outVolRef = "pelvisHealthyVolRef";

	imageVolRef = genVolRef->returnResultImage((char*)outVolRef.c_str());	
	//genVolRef->writeResultImage((char*)outVolRef.c_str());

	//Recorriendo el nro de niveles
	for(int i=0; i < numLevels; i++){
		//Numero de niveles actual
		std::string currentLvl = std::to_string(i);

		//Armando la ruta y el nombre del log
		logfilename = "level" + currentLvl + ".txt";
		std::string currentFileLog(dirResultados + logfilename);

		//Limpiando la primera columna de cada archivo
		std::string cmdClearFirstColLog("awk '{$1=$2=\"\"; print $0}' " + currentFileLog +" > " +dirResultados + "n"+ logfilename);
		std::system(cmdClearFirstColLog.c_str());

		std::ifstream fileParametros(dirResultados + "n" + logfilename);


		float vx,vy,vz; //,tx,ty,tz,sg,rx,ry,rz;
		//Recorre cada linea del archivo
		std::string lineParametros;
		std::string cmdGenFixedImages="";
		int cont = 0;

		while(std::getline(fileParametros,lineParametros)){
			vx = 0; vy=0; vz=0; tx=0; ty=0; tz=0; sg=0; rx=0; ry=0; rz=0;
			
			std::string strcont = std::to_string(cont);
			
			std::istringstream currentParamTransf(lineParametros);
			
				
			if(!(currentParamTransf >> vx >> vy >> vz >> tx >> ty >> tz >> sg))
			{
				std::cout << "Dentro de While"<< std::endl;		

				break;
			}

			//Procesamos los valores actuales de transformacion
			std::cout << vx << " " << vy << " " << vz << " " << tx << " " << ty << " " << tz << " " << sg << std::endl;
			
			//Conversion de versor a grados	
			util->convertVersorToEuler(vx,vy,vz,rx,ry,rz);
	
			cmdGenFixedImages = "";
			
			std::string logdifimage = "../outputData/logimages/log" + currentLvl + strcont + ".txt";

			GenVirtGen *genImagVirtual = new GenVirtGen(type_projection,(char*)logdifimage.c_str());
			genImagVirtual->readMovingImage(movingImage);	
			genImagVirtual->initTransform(tx,ty,tz,dcx,dcy,dcz,rx,ry,rz,sg);
			genImagVirtual->initInterpolator(focalPointx,focalPointy,focalPointz,threshold);
			genImagVirtual->initResampleFilter(dxx, dyy, im_sx,im_sy,customized_2DCX,o2Dx,o2Dy,scd);
			std::string outputFile = "pelvisHealthy_ap_vis_" + currentLvl + strcont;
		        //char *ofile = outputFile.c_str();	
			imageVolTemp = genImagVirtual->returnResultImage((char*)outputFile.c_str());
			//genImagVirtual->writeResultImage((char*)outputFile.c_str());
			
			//Realizamos la substraccion de la imagen de referencia y la imagen fija trasladada
			

			//the difference is between the projection and the current fixed image in this case the last level
			SubtracterType::Pointer subtracter = SubtracterType::New();
			subtracter->SetInput1( imageVolTemp );
			subtracter->SetInput2( imageVolRef); 
			
			//Rescalando la imagen a un umbral de 0 y 255			

			typedef itk::RescaleIntensityImageFilter< MovingImageType, OutputImageType > RescaleFilterType;
			//threshold the image with the correct intensity values
			RescaleFilterType::Pointer rescaler1 = RescaleFilterType::New();
			//rescaler1->SetOutputMinimum(   0 );
			//rescaler1->SetOutputMaximum( 255 );
			rescaler1->SetInput( subtracter->GetOutput());
			rescaler1->Update(); //Imagen movible proyeccion con  0-255 y threshold 0 por el interpolador

			//Girando la imagen en eje horizontal
			using FlipImageFilterType = itk::FlipImageFilter<OutputImageType>;

			FlipImageFilterType::Pointer flipFilter = FlipImageFilterType::New();
			flipFilter->SetInput(rescaler1->GetOutput());
			FlipImageFilterType::FlipAxesArrayType flipAxes;
			flipAxes[0] = true;
			flipAxes[1] = true;
			flipFilter->SetFlipAxes(flipAxes);

			std::string outDifFile = "../outputData/difimages/diference_" + currentLvl + strcont + ".png";
			WriterType::Pointer subtractionWriter = WriterType::New();
			subtractionWriter->SetFileName((char*) outDifFile.c_str() );
			subtractionWriter->SetInput( flipFilter->GetOutput() );
			std::cout << "Writing subtraction file " << subtractionWriter->GetFileName() << std::endl;
			try
			{
				subtractionWriter->Update();
			}
			catch( itk::ExceptionObject & e )
			{
				std::cerr << e.GetDescription() << std::endl;
			}

			cont++;


		}
	}

}

