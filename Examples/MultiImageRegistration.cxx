#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include <itkCommand.h>
//#include <itkEuler3DTransform.h>
#include <itkSimilarity3DTransform.h>

#include <itkFRPROptimizer.h>
//#include "itkVersorRigid3DTransformOptimizer.h"
//#include "itkPowellOptimizer.h"

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMultiResolutionPyramidImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include "itkTimeProbesCollectorBase.h"
#include "utils.h"
#include "itkNormalizedGradientCorrelationMultiImageToImageMetric.h"

#include "itkMultiResolutionMultiImageToImageRegistrationMethod.h"
#include "itkPatchedRayCastInterpolateImageFunction.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <itksys/SystemTools.hxx>
#include <vector>
#include "itkNormalizeImageFilter.h"
/**
 * MultiImageRegistration
 *
 * This example is a basic application which takes a moving 3D image and
 * registers it to a given number of fixed 2D projections. Registration is
 * done maximising the Normalized Gradient Correlation metric with a rigid
 * transform. Three resolution levels are used, using downsampling factors of
 * 4 and 2 for the first and second levels and the original resolutions for
 * the last one.
 *
 */


//----------------------------------------------------------------------------
// Registration observer class
//----------------------------------------------------------------------------
template <class TFixedImage, class TMovingImage>
class RegistrationObserver : public itk::Command
{
	public:
		typedef  RegistrationObserver<TFixedImage,TMovingImage> Self;
		typedef  itk::Command             Superclass;
		typedef  itk::SmartPointer<Self>  Pointer;
		itkNewMacro( Self );
		std::ofstream *logRegistration;
		std::string logoptimizer = "";
	protected:
		RegistrationObserver() {};

	public:
		typedef itk::MultiResolutionMultiImageToImageRegistrationMethod<TFixedImage,TMovingImage> RegistrationType;
		typedef RegistrationType*   RegistrationPointer;
		typedef itk::FRPROptimizer  OptimizerType;

		void Execute(itk::Object *caller, const itk::EventObject & event)
		{
			RegistrationPointer registration = dynamic_cast<RegistrationPointer>( caller );
			if( itk::IterationEvent().CheckEvent( & event ) )
			{	
				logoptimizer = "Resolution level " + std::to_string(registration->GetCurrentLevel()) + "\n";
				logRegistration->write(logoptimizer.c_str(),logoptimizer.size());			
				std::cout << std::endl << "Resolution level " << registration->GetCurrentLevel() << std::endl;
				//Reducir a la mitad la tolerancia de parametros y el tamaño de paso en la busqueda de valores
				//los dos a la mitad en cada nivel de resolucion, en la ultima ambos deben ser de muy bajo valor
				OptimizerType* optimizer = dynamic_cast<OptimizerType*>(registration->GetOptimizer() );

				optimizer->SetStepLength( 0.5 * optimizer->GetStepLength() );
				//std::cout << "StepLength set to " << optimizer->GetStepLength() << std::endl;
				optimizer->SetStepTolerance( 0.5 * optimizer->GetStepTolerance() );
				//std::cout << "StepTolerance set to " << optimizer->GetStepTolerance() << std::endl;
			}
		}

		//again call the execution method
		void Execute(const itk::Object* caller, const itk::EventObject & event)
		{}
		void buildOfstream(std::ofstream & log){
			logRegistration = &log;
		}

};


//----------------------------------------------------------------------------
// Optimizer observer class
//----------------------------------------------------------------------------
class OptimizerObserver : public itk::Command
{
	public:
		typedef  OptimizerObserver     Self;
		typedef  itk::Command             Superclass;
		typedef  itk::SmartPointer<Self>  Pointer;
		itkNewMacro( Self );
		std::ofstream *logOptimizer;	
		std::string logoptimizer = "";
	protected:
		OptimizerObserver() {};


	public:
		//typedef itk::VersorRigid3DTransformOptimizer  OptimizerType;
		typedef itk::FRPROptimizer  OptimizerType;

		typedef const OptimizerType*      OptimizerPointer;

		void Execute(itk::Object *caller, const itk::EventObject & event)
		{
			Execute( (const itk::Object *) caller, event);
		}

		//again call the execution method
		void Execute(const itk::Object* caller, const itk::EventObject & event)
		{
			OptimizerPointer optimizer = dynamic_cast<OptimizerPointer>( caller );
			if( itk::IterationEvent().CheckEvent( & event ) )
			{
				//std::cout << optimizer;
				/*std::cout << "Iteration " << optimizer->GetCurrentIteration()
				  << "/" << optimizer->GetMaximumIteration() << " Position: " <<
				  optimizer->GetCurrentPosition() << " Value: " <<
				  optimizer->GetCurrentCost() << std::endl;*/

				logoptimizer = "Iteration: " + std::to_string(optimizer->GetCurrentIteration()) + " "
					+ "Similarity: " + std::to_string(optimizer->GetValue()) + " "
					+ "Position: " + std::to_string(optimizer->GetCurrentPosition()[0]) + " " +
					std::to_string(optimizer->GetCurrentPosition()[1]) + " " +
					std::to_string(optimizer->GetCurrentPosition()[2]) + " " +
					std::to_string(optimizer->GetCurrentPosition()[3]) + " " +
					std::to_string(optimizer->GetCurrentPosition()[4]) + " " +
					std::to_string(optimizer->GetCurrentPosition()[5]) + " " +
					std::to_string(optimizer->GetCurrentPosition()[6]) + "\n";
				logOptimizer->write(logoptimizer.c_str(),logoptimizer.size());		
				std::cout << " Similarity: " << optimizer->GetValue() << std::endl;
			}
			else if( itk::StartEvent().CheckEvent( & event ) )
			{
				//std::cout << "Optimization started ..." << std::endl;
			}
			else if( itk::EndEvent().CheckEvent( & event ) )
			{
				//std::cout << "Optimization ended." << std::endl;
			}

		}

		void buildOfstream(std::ofstream & log){
			logOptimizer = &log;
		}

};


//----------------------------------------------------------------------------
// Main function
//----------------------------------------------------------------------------
int main(int argc, char* argv[] )
{
	char * movingImage = NULL;
	int numberFixedImages = 0;

	//Declaracion Dinamica del Entorno de FixedImages
	//
	std::vector < char * > fixedImages;	
	std::vector< std::vector <int>> fp;


	float steptolerance = 0.;
	float steplength = 0.;

	int numLevels = 0;

	std::vector< int > resolutions;

	char * outputDirectory = NULL;
	char * inputTransform = NULL;
	char * logFileName = NULL;

	bool writeFinalVol = false;
	bool writeStatistics = false;
	bool ok;	
	while(argc > 1){
		ok = false;

		if ((ok == false) && (strcmp(argv[1], "-movingImage") == 0))
		{
			argc--; argv++;
			ok = true;
			movingImage = argv[1];
			argc--; argv++;
		}
        if ((ok == false) && (strcmp(argv[1], "-numFixedImages") == 0))
		{
			argc--; argv++;
			ok = true;
			numberFixedImages = atoi(argv[1]);
			argc--; argv++;

			//Obviamente debemos confiar del anterior parametro
			if(numberFixedImages != 0){
				int contFixedImages = 0;

				//Iterar sobre el nro de Fixed Images
				while(contFixedImages < numberFixedImages){

					std::string str_fixed = "-f"+std::to_string(contFixedImages); 	
                    if ((ok == true) && (strcmp(argv[1], str_fixed.c_str()) == 0)){
						argc--; argv++;
						ok = true;
						fixedImages.push_back(argv[1]);
						//Lectura de su respectivo punto focal
                        argc--; argv++;

                        std::vector<int> fp_temp;
                        for(int j=0; j<3; j++){
                            fp_temp.push_back(atoi(argv[1])); argc--; argv++;
                        }
                        fp.push_back(fp_temp);
                        contFixedImages++;
					}
				}
			}

		}

		if ((ok == false) && (strcmp(argv[1], "-steptolerance") == 0))
		{
			argc--; argv++;
			ok = true;
			steptolerance = atof(argv[1]);
			argc--; argv++;
		}

		if ((ok == false) && (strcmp(argv[1], "-stepsize") == 0))
		{
			argc--; argv++;
			ok = true;
			steplength = atof(argv[1]);
			argc--; argv++;
		}

		if ((ok == false) && (strcmp(argv[1], "-numLevels") == 0))
		{
			argc--; argv++;
			ok = true;
			numLevels = atoi(argv[1]);
			argc--; argv++;
			//Confiamos en el que anterior parametro existe
			if(numLevels != 0 && (ok == true) && (strcmp(argv[1],"-schedule")==0)){
				int contLevels = 0;
				while(contLevels < numLevels){
					argc--; argv++;
					ok = true;
					resolutions.push_back(atoi(argv[1]));
					contLevels++;
				}
				argc--; argv++;
			}else{
				std::cerr<<"Falta Nro de Niveles "<<std::endl; 
			}

		}

		if ((ok == false) && (strcmp(argv[1], "-outputDirectory") == 0))
		{
			argc--; argv++;
			ok = true;
			outputDirectory = argv[1];
			argc--; argv++;
		}

		if ((ok == false) && (strcmp(argv[1], "-inpuTransform") == 0))
		{
			argc--; argv++;
			ok = true;
			inputTransform = argv[1];
			argc--; argv++;
		}
		if ((ok == false) && (strcmp(argv[1], "-logfilename") == 0))
		{
			argc--; argv++;
			ok = true;
			logFileName = argv[1];
			argc--; argv++;
		}
		if ((ok == false) && (strcmp(argv[1], "-writeFinalVol") == 0)){
			argc--; argv++;
			ok = true;
			writeFinalVol = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-writeStatistics") == 0)){
			argc--; argv++;
			ok = true;
			writeStatistics = true;
		}
	}	
	/*
	   if( argc < 13 )
	   {
	   std::cerr << "Usage: " << argv[0] << " [movingImage] [N] [fixedImage1] "
	   << "[focalPoint1_x] [focalPoint1_y] [focalPoint1_z] ... "
	   << "[focalPointN_x] [focalPointN_y] [focalPointN_z] "
	   << "[tolerance] [stepsize]"
	   << "[levels_n] [schedule_1,schedule_2...,schedule_n]"
	   << "[outputDirectory] <inputTransform> [logFileName]" << std::endl;
	   return EXIT_FAILURE;
	   }
	   */	
	//---------------------------------------------------------------------------
	//Ofstream General para todo el proceso de registro	
	//---------------------------------------------------------------------------

	//Collector de Tiempos
	itk::TimeProbesCollectorBase timer;

	const unsigned int Dimensions = 3;

	typedef itk::Image<short int,Dimensions>  FixedImageType;
	typedef itk::Image<short int,Dimensions>  MovingImageType;

	//create the type of the registration
	typedef itk::MultiResolutionMultiImageToImageRegistrationMethod< FixedImageType,MovingImageType >    RegistrationType;

	RegistrationType::Pointer registration = RegistrationType::New();

	typedef RegistrationObserver< FixedImageType, MovingImageType> RegistrationObserverType;
	RegistrationObserverType::Pointer registrationObserver = RegistrationObserverType::New();
	//the registration add the registraton observer in order to show some datas
	registration->AddObserver( itk::IterationEvent(), registrationObserver );
	
	//Variable para normalizacion de imagenes con mismo contraste
	/*using FixedNormalizeFilterType = itk::NormalizeImageFilter< FixedImageType, FixedImageType >;
	using MovingNormalizeFilterType = itk::NormalizeImageFilter< MovingImageType, MovingImageType>;
	FixedNormalizeFilterType::Pointer fixedNormalizer = FixedNormalizeFilterType::New();
	MovingNormalizeFilterType::Pointer movingNormalizer = MovingNormalizeFilterType::New();	
	*/
	//----------------------------------------------------------------------------
	// Create the transform
	//----------------------------------------------------------------------------
	typedef itk::Similarity3DTransform< double> TransformType;

	TransformType::Pointer transform = TransformType::New();
	//transform->SetComputeZYX(true);
	transform->SetIdentity();
	registration->SetTransform( transform );


	//----------------------------------------------------------------------------
	// Load the moving image
	//----------------------------------------------------------------------------
	typedef itk::ImageFileReader<MovingImageType> MovingImageReaderType;
	//pyramid
	typedef itk::MultiResolutionPyramidImageFilter< MovingImageType, MovingImageType  > MovingImagePyramidType;

	MovingImageReaderType::Pointer movingReader = MovingImageReaderType::New();
	movingReader->SetFileName( movingImage );

	try
	{
		movingReader->Update();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cout << e.GetDescription() << std::endl;
		return EXIT_FAILURE;
	}

	//movingNormalizer->SetInput(movingReader->GetOutput());
	registration->SetMovingImage( movingReader->GetOutput() ); //Imagen Movible consu informacion de intensidad completa 

	MovingImagePyramidType::Pointer movingPyramidFilter = MovingImagePyramidType::New();
	//it is only the allocate movin pyramid and given to the registration
	registration->SetMovingImagePyramid( movingPyramidFilter );

	//----------------------------------------------------------------------------
	// Load the fixed images
	//----------------------------------------------------------------------------
	typedef FixedImageType::ConstPointer          FixedImageConstPointer;
	typedef itk::ImageFileReader<FixedImageType>  FixedImageReaderType;

	typedef itk::PatchedRayCastInterpolateImageFunction<MovingImageType, double> InterpolatorType;
	typedef InterpolatorType::InputPointType FocalPointType;

	typedef itk::MultiResolutionPyramidImageFilter< FixedImageType, FixedImageType  >  FixedImagePyramidType;

	const unsigned int FImgTotal = numberFixedImages;
	argc--;
	argv++;

	float vfocalPoint[FImgTotal][3];

	//for every fixed image
	for( unsigned int f=0; f<FImgTotal; f++ )
	{
		FixedImageReaderType::Pointer fixedReader = FixedImageReaderType::New();
		fixedReader->SetFileName( fixedImages[f] );

		try
		{
			fixedReader->Update();
		}
		catch( itk::ExceptionObject & e )
		{
			std::cout << e.GetDescription() << std::endl;
			return EXIT_FAILURE;
		}

		FixedImageConstPointer fixedImage = fixedReader->GetOutput();
		
		//normalizar la imagen fija
		//fixedNormalizer->SetInput(fixedImage);
	
		registration->AddFixedImage( fixedImage ); //Imagen Fija leida con threhosld 0 con 0 255 seteado en el generador

		registration->AddFixedImageRegion( fixedImage->GetBufferedRegion() );

		InterpolatorType::Pointer interpolator = InterpolatorType::New();
		FocalPointType focalPoint;

		vfocalPoint[f][0] = focalPoint[0] = fp[f][0];
		vfocalPoint[f][1] = focalPoint[1] = fp[f][1];
		vfocalPoint[f][2] = focalPoint[2] = fp[f][2];

		interpolator->SetFocalPoint( focalPoint );
		interpolator->SetTransform( transform );
		interpolator->SetThreshold( 0.0 );


		/* Fijando al interpolador la nueva direccion para el volumen acorde a la direccion de la imagen fija*/	
		//typedef MovingImageType::DirectionType directionType;
		//directionType::Pointer directionPointer = directionType::New();
		//directionPointer = fixedImage->GetDirection();

		interpolator->SetDirectionFixed( fixedImage->GetDirection() );


		registration->AddInterpolator( interpolator );

		FixedImagePyramidType::Pointer fixedPyramidFilter = FixedImagePyramidType::New();
		registration->AddFixedImagePyramid( fixedPyramidFilter );

	}


	//----------------------------------------------------------------------------
	// Create the multi metric
	//----------------------------------------------------------------------------
	typedef itk::NormalizedGradientCorrelationMultiImageToImageMetric<FixedImageType, MovingImageType> MultiMetricType;
	MultiMetricType::Pointer multiMetric = MultiMetricType::New();

	//it only allocate the metric
	registration->SetMultiMetric( multiMetric );


	//----------------------------------------------------------------------------
	// Create the optimizer
	//----------------------------------------------------------------------------
	//typedef itk::VersorRigid3DTransformOptimizer OptimizerType;
	typedef itk::FRPROptimizer OptimizerType;


	OptimizerType::Pointer optimizer = OptimizerType::New();

	unsigned int ParTotal = transform->GetNumberOfParameters();

	//La escala inicial que tendran los parametros en este caso
	//sera la sexta parte del total 100/6 = 16.6 o cuarta parte 100/4 = 25.0
	OptimizerType::ScalesType scales( ParTotal );
	scales.Fill(itk::NumericTraits<OptimizerType::ScalesType::ValueType>::One);
	scales[0] = 1.0*10.0;
	scales[1] = 1.0*10.0;
	scales[2] = 1.0*10.0;
	scales[3] = 1.0/256.3;
	scales[4] = 1.0/256.3;
	scales[5] = 1.0/256.3;
	scales[6] = 1.0;

	optimizer->SetScales( scales );

	optimizer->SetMaximize( true ); //true
	//optimizer->SetRelaxationFactor(0.5);
	//optimizer->SetMaximumStepLength(0.2000);
	//optimizer->SetMinimumStepLength(0.0001);
	//optimizer->SetNumberOfIterations(2500);

	optimizer->SetMaximumIteration( 1000 ); //100
	optimizer->SetMaximumLineIteration(10); //10 - 4
	optimizer->SetValueTolerance( 1e-3 );
	optimizer->SetUseUnitLengthGradient( true );
	optimizer->SetToPolakRibiere();
	optimizer->SetCatchGetValueException( true );
	optimizer->SetMetricWorstPossibleValue(-itk::NumericTraits<MultiMetricType::MeasureType>::infinity() );

	//These parameters will be halved at the begginng of each resolution level
	optimizer->SetStepTolerance(steptolerance);
	optimizer->SetStepLength( steplength );

	//print the information in the optimizer.
	OptimizerObserver::Pointer optimizerObserver = OptimizerObserver::New();
	optimizer->AddObserver( itk::IterationEvent(), optimizerObserver );
	optimizer->AddObserver( itk::StartEvent(), optimizerObserver );
	optimizer->AddObserver( itk::EndEvent(), optimizerObserver );	

	//set the optimizer in the registration process
	registration->SetOptimizer( optimizer );

	//----------------------------------------------------------------------------
	// Set the moving and fixed images' schedules
	//----------------------------------------------------------------------------
	const unsigned int ResolutionLevels = numLevels; 

	RegistrationType::ScheduleType movingSchedule( ResolutionLevels,Dimensions);
	RegistrationType::ScheduleType fixedSchedule( ResolutionLevels,Dimensions );

	//set the scales according to the vector of the scales
	for(int i=0; i<ResolutionLevels; i++){
		for(int j=0; j<Dimensions; j++){
			if(j!=2){
				fixedSchedule[i][j]=resolutions[i];
				movingSchedule[i][j]=resolutions[i];
			}else{
				fixedSchedule[i][j]=1;
				movingSchedule[i][j]=resolutions[i];
			}
		}
	}

	//set both schedules from fixed and moving images in the registration process
	registration->SetSchedules( fixedSchedule, movingSchedule );

	//----------------------------------------------------------------------------
	// Get the output directory
	//----------------------------------------------------------------------------
	std::string outDir = outputDirectory;
	//el directorio donde estaran los resultados del registro
	std::string fname;
	itksys::SystemTools::MakeDirectory((fname + outDir));

	//----------------------------------------------------------------------------
	//Nombre del archivo de registro 
	//----------------------------------------------------------------------------
	std::ofstream logregistro( logFileName );

	//Set the observer in the ofstream register
	registrationObserver->buildOfstream(logregistro);
	optimizerObserver->buildOfstream(logregistro);


	//----------------------------------------------------------------------------
	// Save projected images
	//----------------------------------------------------------------------------
	typedef itk::ResampleImageFilter<MovingImageType, FixedImageType>  ResamplerType;
	typedef itk::SubtractImageFilter<FixedImageType,FixedImageType, FixedImageType>  SubtracterType;
	typedef itk::ImageFileWriter<FixedImageType> WriterType;
	typedef itk::RescaleIntensityImageFilter< MovingImageType, MovingImageType > RescaleFilterType;


	//traverse every fixed image
	for( unsigned int f=0; f<FImgTotal; f++ )
	{
		//the moving image is updating with the last transformation
		ResamplerType::Pointer resampler = ResamplerType::New();
		resampler->SetInput( registration->GetMovingImage() ); //Input Imagen movible con toda su informacion de intensidad
		resampler->SetInterpolator( registration->GetMultiInterpolator()[f] );
		resampler->SetTransform( transform );
		resampler->SetUseReferenceImage( true );
		resampler->SetReferenceImage( registration->GetFixedMultiImage()[f] );


		//threshold the image with the correct intensity values
		RescaleFilterType::Pointer rescaler1 = RescaleFilterType::New();
		rescaler1->SetOutputMinimum(   0 );
		rescaler1->SetOutputMaximum( 255 );
		rescaler1->SetInput( resampler->GetOutput());
		rescaler1->Update(); //Imagen movible proyeccion con  0-255 y threshold 0 por el interpolador


		std::stringstream strStream;
		strStream << outDir << "/initial_projection";

		//impresion del punto focal
		strStream << "_"<< vfocalPoint[f][0] <<"_"<<vfocalPoint[f][1]<<"_"<<vfocalPoint[f][2];
		strStream.width(3);
		strStream.fill('FF');
		strStream << f;
		strStream.width(0);
		strStream << ".mha";

		//the projection is condition
		WriterType::Pointer projectionWriter = WriterType::New();
		projectionWriter->SetFileName( strStream.str().c_str() );
		projectionWriter->SetInput( rescaler1->GetOutput() );

		logregistro << "Writing projection file " << projectionWriter->GetFileName() << std::endl;
		try
		{
			projectionWriter->Update();
		}
		catch( itk::ExceptionObject & e )
		{
			std::cerr << e.GetDescription() << std::endl;
		}

		//the difference is between the projection and the current fixed image in this case the last level
		SubtracterType::Pointer subtracter = SubtracterType::New();
		subtracter->SetInput1( registration->GetFixedMultiImage()[f] ); //imagen fija generada con threhsold 0 y ya rescalado previamente de 0 255 en el generador 
		subtracter->SetInput2( rescaler1->GetOutput() ); //Imagen movible proyeccion con rescala de  0-255 y threshold 0 por el interpolador


		std::stringstream strStream2;
		strStream2 << outDir << "/initial_subtraction";
		strStream2 << "_"<< vfocalPoint[f][0] <<"_"<<vfocalPoint[f][1]<<"_"<<vfocalPoint[f][2];
		strStream2.width(3);
		strStream2.fill('SS');
		strStream2 << f;
		strStream2.width(0);
		strStream2 << ".mha";

		WriterType::Pointer subtractionWriter = WriterType::New();
		subtractionWriter->SetFileName( strStream2.str().c_str() );
		subtractionWriter->SetInput( subtracter->GetOutput() );
		logregistro << "Writing subtraction file " << subtractionWriter->GetFileName() << std::endl;
		try
		{
			subtractionWriter->Update();
		}
		catch( itk::ExceptionObject & e )
		{
			std::cerr << e.GetDescription() << std::endl;
		}

		std::stringstream strStream3;
		strStream3 << outDir << "/fixedImage";
		strStream3 << f;
		strStream3.width(0);
		strStream3 << ".mha";

		WriterType::Pointer fixedWriter = WriterType::New();
		fixedWriter->SetFileName( strStream3.str().c_str() );
		fixedWriter->SetInput( registration->GetFixedMultiImage()[f] );  //imagen fija generada con threhsold 100 (por ser virtual viene de un volumen  0 255) 

		logregistro << "Writing subtraction file " << fixedWriter->GetFileName() << std::endl; 
		try
		{
			fixedWriter->Update();
		}
		catch( itk::ExceptionObject & e )
		{
			std::cerr << e.GetDescription() << std::endl;
		}


	}


	//----------------------------------------------------------------------------
	// Read the transform file, if given
	//----------------------------------------------------------------------------

	if( inputTransform != NULL ) // if there are still arguments in the command line
	{
		itk::TransformFileReader::Pointer transformReader = itk::TransformFileReader::New();
		transformReader->SetFileName( inputTransform );

		logregistro << "Reading transform file " << transformReader->GetFileName() << std::endl;
		try
		{
			transformReader->Update();
		}
		catch( itk::ExceptionObject & e )
		{
			std::cout << e.GetDescription() << std::endl;
			return EXIT_FAILURE;
		}

		itk::TransformFileReader::TransformListType* transformList = transformReader->GetTransformList();

		if( transformList->size() > 1 )
		{ 
			std::cout << "Only one transform expected from file. "
				<< "Using first available transform" << std::endl;
		}

		itk::TransformFileReader::TransformPointer baseTransform = transformList->front();

		//the transformation file has a specify format
		std::string className = baseTransform->GetNameOfClass();
		if( className.compare("Similarity3DTransform") != 0 )
		{
			std::cerr << "Transform class must be Similarity3DTransform." << std::endl;
			std::cerr << "Found " << className << " instead."
				<< std::endl;
			return EXIT_FAILURE;
		}
		//get separated paremters as fixed paremeters as normal parameters
		transform->SetFixedParameters( baseTransform->GetFixedParameters() );
		transform->SetParameters( baseTransform->GetParameters() );
	}

	//this will be the initial parameters
	registration->SetInitialTransformParameters( transform->GetParameters() );


	//Print parameters using in the registration

	logregistro << "StepTolerance "<< optimizer->GetStepTolerance()<<std::endl; 
	logregistro <<  "StepLength "<< optimizer->GetStepLength()<<std::endl; 
	logregistro <<  "Schedules"<< resolutions[0]<<" "<<resolutions[1]<<" "<<resolutions[2]<<std::endl; 

	//Create stream for important parameters
	std::ostringstream strStepTol, strStepLen;
	strStepTol << steptolerance;
	strStepLen << steplength;

	//----------------------------------------------------------------------------
	// Save the initial matrix file
	//----------------------------------------------------------------------------
	std::string inputMatrixPath = outDir;
	inputMatrixPath.append( "/inMatrix" );
	inputMatrixPath.append("_"+strStepTol.str()+"_"+strStepLen.str());
	inputMatrixPath.append(".txt");

	std::ofstream inputMatrixFile;
	logregistro <<  "Saving input matrix file " << inputMatrixPath << std::endl;
	inputMatrixFile.open( inputMatrixPath.c_str() );

	if(  inputMatrixFile.is_open() )
	{
		//print matrix input information
		TransformType::MatrixType m = transform->GetMatrix();
		TransformType::OffsetType off = transform->GetOffset();
		inputMatrixFile << m[0][0] << " " << m[0][1] << " " << m[0][2] << " " << off[0] << std::endl;
		inputMatrixFile << m[1][0] << " " << m[1][1] << " " << m[1][2] << " " << off[1] << std::endl;
		inputMatrixFile << m[2][0] << " " << m[2][1] << " " << m[2][2] << " " << off[2] << std::endl;
		inputMatrixFile << "0.0 0.0 0.0 1.0" << std::endl;
		inputMatrixFile.close();
	}
	else
	{
		std::cerr << "Could not open input matrix file " << inputMatrixPath << std::endl;
	}


	//----------------------------------------------------------------------------
	// Start the registration
	//----------------------------------------------------------------------------

	try
	{
		//measure the total time
		timer.Start("Tiempo de Registro");
		registration->Update();
		timer.Stop("Tiempo de Registro");
		//logregistro << 	"Total Registration time " << cputimer.GetMean() << " mean.\n"<< std::endl;
		//logregistro << "TotalTime Process Object: "<< registration->GetMTime() << std::endl;
		//logregistro <<  "TimeTransform: " << transform->GetMTime() << std::endl;
		//logregistro << "TimeMetric: "<< multiMetric->GetMTime() << std::endl;	
		//logregistro << "TimeOptmizer: " << optimizer->GetMTime() << std::endl;
		//std::cout << "CPU Registration took " << cputimer.GetMean() << " mean.\n"<< std::endl;
		//std::cout << "TotalTime Process Object: "<< registration->GetMTime() << std::endl;

	}
	catch( itk::ExceptionObject & e )
	{
		std::cout << e.GetDescription() << std::endl;
		return EXIT_FAILURE;
	}


	//---------------------------------------------------------------------------
	// General Information after Registration Process
	//---------------------------------------------------------------------------

	typedef RegistrationType::ParametersType ParametersType;
	ParametersType finalParameters = registration->GetLastTransformParameters();
	//Parametro para convertir de angulos sexagesimales a radianes
	const double dtr = ( atan(1.0) * 4.0 ) / 180.0;
	const double rtd = 180 / ( atan(1.0) * 4.0);
	//La salida en Euler es en radianes asi que debemos convertir a grados
	//para observar el angulo en grados sexagesimales

	const double RotationAlongX = finalParameters[0]; // Convert radian to degree
	const double RotationAlongY = finalParameters[1];
	const double RotationAlongZ = finalParameters[2];
	const double TranslationAlongX = finalParameters[3];
	const double TranslationAlongY = finalParameters[4];
	const double TranslationAlongZ = finalParameters[5];
	const double EscalaXYZ = finalParameters[6];
	const int numberOfIterations = optimizer->GetCurrentIteration();

	const double bestValue = optimizer->GetValue();

	logregistro << "Result = " << std::endl;
	logregistro << " Rotation Along X = " << RotationAlongX  << " deg" << std::endl;
	logregistro << " Rotation Along Y = " << RotationAlongY  << " deg" << std::endl; 
	logregistro << " Rotation Along Z = " << RotationAlongZ  << " deg" << std::endl;
	logregistro << " Translation X = " << TranslationAlongX  << " mm" << std::endl;
	logregistro << " Translation Y = " << TranslationAlongY  << " mm" << std::endl;
	logregistro << " Translation Z = " << TranslationAlongZ  << " mm" << std::endl;
	logregistro << " Escala value = " << EscalaXYZ << std::endl;
	logregistro << " Number Of Iterations = " << numberOfIterations << std::endl;
	logregistro << " Metric value  = " << bestValue          << std::endl;

	//----------------------------------------------------------------------------
	// Save the transform files
	//----------------------------------------------------------------------------

	//output transformation
	std::string outputTransformPath = outDir;
	outputTransformPath.append( "/outTransform" );
	outputTransformPath.append(".txt");

	itk::TransformFileWriter::Pointer transformWriter = itk::TransformFileWriter::New();
	transformWriter->SetInput( transform );
	transformWriter->SetFileName( outputTransformPath.c_str() );

	logregistro << "Saving output transform file " << transformWriter->GetFileName() << std::endl;
	try
	{
		transformWriter->Update();
	}
	catch( itk::ExceptionObject &e )
	{
		std::cout << e.GetDescription() << std::endl;
	}

	std::string outputMatrixPath = outDir;
	outputMatrixPath.append( "/outMatrix" );
	outputMatrixPath.append("_"+strStepTol.str()+"_"+strStepLen.str());
	outputMatrixPath.append(".txt");

	std::ofstream outputMatrixFile;
	logregistro << "Saving output matrix file " << outputMatrixPath << std::endl;
	outputMatrixFile.open( outputMatrixPath.c_str() );

	if(  outputMatrixFile.is_open() )
	{
		//print ouput matrix information
		TransformType::MatrixType m = transform->GetMatrix();
		TransformType::OffsetType off = transform->GetOffset();
		outputMatrixFile << m[0][0] << " " << m[0][1] << " " << m[0][2] << " " << off[0] << std::endl;
		outputMatrixFile << m[1][0] << " " << m[1][1] << " " << m[1][2] << " " << off[1] << std::endl;
		outputMatrixFile << m[2][0] << " " << m[2][1] << " " << m[2][2] << " " << off[2] << std::endl;
		outputMatrixFile << "0.0 0.0 0.0 1.0" << std::endl;
		outputMatrixFile.close();
	}
	else
	{
		std::cerr << "Could not open output matrix file " << outputMatrixPath << std::endl;
	}


	//----------------------------------------------------------------------------
	// Save projected images
	//----------------------------------------------------------------------------

	//traverse every fixed image
	for( unsigned int f=0; f<FImgTotal; f++ )
	{
		//the moving image is updating with the last transformation
		ResamplerType::Pointer resampler = ResamplerType::New();
		resampler->SetInput( registration->GetMovingImage() ); //Imagen Movible entra con toda su informacion de intensidad
		resampler->SetInterpolator( registration->GetMultiInterpolator()[f] );
		resampler->SetTransform( transform );
		resampler->SetUseReferenceImage( true );
		resampler->SetReferenceImage( registration->GetFixedMultiImage()[f] );

		RescaleFilterType::Pointer rescaler1 = RescaleFilterType::New();
		rescaler1->SetOutputMinimum(   0 );
		rescaler1->SetOutputMaximum( 255 );
		rescaler1->SetInput( resampler->GetOutput());
		rescaler1->Update(); //Imagen movible proyeccion sin  0-255 ya que esta por defecto en la imagen y threshold 100 por el interpolador


		std::stringstream strStream;
		strStream << outDir << "/projection";

		//impresion del punto focal
		strStream << "_"<< vfocalPoint[f][0] <<"_"<<vfocalPoint[f][1]<<"_"<<vfocalPoint[f][2];

		//impresion Resolutions Levels, schedules, optimizer->GetStepLength optimizer->GetStepTolerance

		strStream<<"_"<<strStepTol.str();
		strStream<<"_"<<strStepLen.str();
		strStream.width(3);
		strStream.fill('FF');
		strStream << f;
		strStream.width(0);
		strStream << ".mha";

		//the projection is condition
		WriterType::Pointer projectionWriter = WriterType::New();
		projectionWriter->SetFileName( strStream.str().c_str() );
		//projectionWriter->SetInput( resampler->GetOutput() );
		projectionWriter->SetInput( rescaler1->GetOutput() ); //Imagen Movible proyeccion de 0 a 255 (por deformada) y threshold 100 por ser rescalada
		logregistro << "Writing Final projection file " << projectionWriter->GetFileName() << std::endl;

		try
		{
			projectionWriter->Update();
		}
		catch( itk::ExceptionObject & e )
		{
			std::cerr << e.GetDescription() << std::endl;
		}

		//the difference is between the projection and the current fixed image in this case the last level
		SubtracterType::Pointer subtracter = SubtracterType::New();
		subtracter->SetInput1( registration->GetFixedMultiImage()[f] ); //Imagen Fija despues del registro sigue siendo la misma de entrada 
		//imagen fija generada con rescale de 0 a 255 y threshold 0 

		subtracter->SetInput2( rescaler1->GetOutput() ); //Imagen movible proyeccion con rescale de  0-255 y threshold 0 por el interpolador


		std::stringstream strStream2;
		strStream2 << outDir << "/subtraction";
		strStream2 << "_"<< vfocalPoint[f][0] <<"_"<<vfocalPoint[f][1]<<"_"<<vfocalPoint[f][2];
		strStream2 << "_" << strStepTol.str();
		strStream2 << "_" << strStepLen.str();
		strStream2.width(3);
		strStream2.fill('SS');
		strStream2 << f;
		strStream2.width(0);
		strStream2 << ".mha";

		WriterType::Pointer subtractionWriter = WriterType::New();
		subtractionWriter->SetFileName( strStream2.str().c_str() );
		subtractionWriter->SetInput( subtracter->GetOutput() );
		logregistro << "Writing Final subtraction file " <<subtractionWriter->GetFileName() << std::endl;
		try
		{
			subtractionWriter->Update();
		}
		catch( itk::ExceptionObject & e )
		{
			std::cerr << e.GetDescription() << std::endl;
		}

	}

	//Aplicando Final Transformacion al volume con referencia al volume fijo
	//La referencia lo unico que se puede cambiar es espaciado y size lo cual 
	//no cambiaria el volumen pero si su resolucion. Asi que todo se dejara por defecto
	//con el mismo volumen de entrada (size, origen, espaciado, direcccion)


	//El tipo de pixel y el nro de dimensiones de salida seran las mismas al volumen de entrada
	if(writeFinalVol){	
		timer.Start("Tiempo de Transformacion"); 	
		using  WriterTypeVol = itk::ImageFileWriter<FixedImageType>;
		using ResampleFilterType = itk::ResampleImageFilter< MovingImageType, MovingImageType>;
		ResampleFilterType::Pointer resamplerVol = ResampleFilterType::New();
		WriterTypeVol::Pointer writer = WriterTypeVol::New();

		resamplerVol->SetTransform(transform);
		resamplerVol->SetInput(registration->GetMovingImage());
		resamplerVol->SetSize(registration->GetMovingImage()->GetBufferedRegion().GetSize());
		resamplerVol->SetUseReferenceImage( true );
		resamplerVol->SetReferenceImage( registration->GetMovingImage() );

		writer->SetFileName(outDir + "/newVolumen.mha");
		writer->SetInput(resamplerVol->GetOutput());
		logregistro<<"Escribiendo Volumen de Salida "<<writer->GetFileName() << std::endl;
		try{
			writer->Update();
		}
		catch(itk::ExceptionObject & e){
			std::cerr << e.GetDescription() << std::endl;
		}
		timer.Stop("Tiempo de Transformacion");

		timer.ExpandedReport(logregistro);
	
	}
    if(writeStatistics){
	Utilitarios *util = new Utilitarios();
	//util->createStats(4,"../outputData/log.txt","../outputData",0);
    }
	
	logregistro.close();

	//----------------------------------------------------------------------------
	// End of the example
	//----------------------------------------------------------------------------
	return EXIT_SUCCESS;
}

