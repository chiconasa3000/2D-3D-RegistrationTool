#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNormalizedGradientCorrelationImageToImageMetric.h"
#include "itkTranslationTransform.h"
#include "itkPatchedRayCastInterpolateImageFunction.h"
#include <itkTransformFileReader.h>
#include "utils.h"
int main( int argc, char * argv[] )
{

	char * movingImageName = NULL;
	//char * outputSpaceParams = NULL;
	char * transfMovingImage = NULL;
	char * typePlot = NULL;	
	int numberFixedImages = 0;

	//Declaracion Dinamica del Entorno de FixedImages
	//
	std::vector < char * > fixedImages;	
	std::vector< std::vector <int> > fp;

	//Rotacion
	//float rx = 0.;
	//float ry = 0.;
	//float rz = 0.;  
	//Translation
	float tx = 0.;
	float ty = 0.;
	float tz = 0.;
	//Scale
	float sg = 1.;

	//Rangos
	float rangex_max = 0;
	float rangex_min = 0;
	float rangey_max = 0;
	float rangey_min = 0;
	float stepx = 0;
	float stepy = 0;	
	
	//Indice de imagen a evaluar
	int indImagenEval = 0;

	//Pueden ser 2 unicas posibilidades para la visualizacion o 2D con un solo parametro
	//o en 3D alterando dos parametros

	//Indice Parametro 1
	int indParam1 = 0;
	//Indice Parametro 2
	int indParam2 = 0;
		
	
	//combineImages
	bool combineImages = false;		

	//Debemos considerar la lectura desde el mismo archivo y manualmente
	bool ok;	
	while(argc > 1){
		ok = false;

		if ((ok == false) && (strcmp(argv[1], "-movingImage") == 0))
		{
			argc--; argv++;
			ok = true;
			movingImageName = argv[1];
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

		/*if ((ok == false) && (strcmp(argv[1], "-outputSpaceParams") == 0))
		{
			argc--; argv++;
			ok = true;
			outputSpaceParams = argv[1];
			argc--; argv++;
		}*/
		if ((ok == false) && (strcmp(argv[1], "-transfMovingImage") == 0))
		{
			argc--; argv++;
			ok = true;
			transfMovingImage = argv[1];
			argc--; argv++;
		}

		if((ok == false) && (strcmp(argv[1], "-indImagenEval") == 0)){
			argc--; argv++;
			ok = true;
			indImagenEval = atoi(argv[1]);
			argc--; argv++;
		}

		if((ok == false) && (strcmp(argv[1], "-typePlot") == 0)){
			argc--; argv++;
			ok = true;
			typePlot = argv[1];
			argc--; argv++;
		}

		if((ok == false) && (strcmp(argv[1],"-indParam")==0)){	
			argc--; argv++;
			ok = true;
			indParam1 = atoi(argv[1]);
			argc--; argv++;
			indParam2 = atoi(argv[1]);
			argc--; argv++;

		}

		if((ok == false) && (strcmp(argv[1], "-rangex") == 0)){
			argc--; argv++;
			ok = true;
			rangex_min = atof(argv[1]);
			argc--; argv++;
			rangex_max = atof(argv[1]);
			argc--; argv++;	
		}
		if((ok == false) && (strcmp(argv[1], "-step") == 0)){
			argc--; argv++;
			ok = true;
			stepx = atof(argv[1]);
			argc--; argv++;
			stepy = atof(argv[1]);
			argc--; argv++;

		}

	
		if((ok == false) && (strcmp(argv[1], "-rangey") == 0)){
			argc--; argv++;
			ok = true;
			rangey_min = atof(argv[1]);
			argc--; argv++;
			rangey_max = atof(argv[1]);
			argc--; argv++;	
		}
		if((ok == false) && (strcmp(argv[1], "-combineImages") == 0)){
			argc--; argv++;
			ok = true;
			combineImages = true;
		}
	
		


	}	

	
	const     unsigned int   Dimensions = 3;
	typedef   short int          PixelType;
	typedef itk::Image<PixelType,Dimensions>  FixedImageType;
	typedef itk::Image<PixelType,Dimensions>  MovingImageType;

	typedef itk::ImageFileReader<MovingImageType> MovingImageReaderType;
	typedef FixedImageType::ConstPointer          FixedImageConstPointer;

	//Vector de imagenes fijas
	typedef std::vector<FixedImageConstPointer> FixedMultiImageType;

	FixedMultiImageType m_FixedMultiImage;
	
	typedef itk::ImageFileReader<FixedImageType>  FixedImageReaderType;
	
	typedef itk::PatchedRayCastInterpolateImageFunction<MovingImageType, double> InterpolatorType;
	typedef typename InterpolatorType::Pointer   InterpolatorPointer;
	typedef std::vector<InterpolatorPointer>    MultiInterpolatorType;

	/** Add, set or get the Interpolators. */
	MultiInterpolatorType       m_MultiInterpolator;

	typedef InterpolatorType::InputPointType FocalPointType;
    	typedef typename FixedImageType::RegionType     FixedImageRegionType;
    	typedef std::vector<FixedImageRegionType>       FixedMultiImageRegionType;

	FixedMultiImageRegionType m_FixedMultiImageRegion;

	//Vector de Metricas
	typedef itk::NormalizedGradientCorrelationImageToImageMetric< FixedImageType, MovingImageType >  MetricType;
	typedef typename MetricType::Pointer MetricPointer;
	typedef std::vector<MetricPointer> MultiMetricType;	
	
	MultiMetricType m_MultiMetric;

	//Transformacion
	typedef itk::Similarity3DTransform< double> TransformType;
	TransformType::Pointer transform = TransformType::New();
	transform->SetIdentity();

	typedef typename TransformType::ParametersType  TransformParametersType;

	const unsigned int FImgTotal = numberFixedImages;

	MovingImageReaderType::Pointer movingReader = MovingImageReaderType::New();
	
	//Lectura de la imagen Movible
	movingReader->SetFileName( movingImageName );
	try
	{
		movingReader->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception catched !" << std::endl;
		std::cerr << excep << std::endl;
	}
	
	MovingImageType::ConstPointer movingImage = movingReader->GetOutput();

  	float vfocalPoint[FImgTotal][3];

	//Lectura de las Imagenes Fijas
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
		
		m_FixedMultiImage.push_back(fixedReader->GetOutput());


		InterpolatorType::Pointer interpolator = InterpolatorType::New();
		FocalPointType focalPoint;

		vfocalPoint[f][0] = focalPoint[0] = fp[f][0];
		vfocalPoint[f][1] = focalPoint[1] = fp[f][1];
		vfocalPoint[f][2] = focalPoint[2] = fp[f][2];

		interpolator->SetFocalPoint( focalPoint );
		interpolator->SetTransform( transform );
		interpolator->SetThreshold( 100.0 );
		interpolator->SetDirectionFixed( m_FixedMultiImage[f]->GetDirection() );
		interpolator->SetInputImage(movingImage);
		m_MultiInterpolator.push_back(interpolator);
		m_FixedMultiImageRegion.push_back( m_FixedMultiImage[f]->GetBufferedRegion() );
		
		//Cada imagen fija tendra su propia metrica
		MetricType::Pointer metricTmp = MetricType::New();
		
		metricTmp->SetTransform( transform );
		metricTmp->SetInterpolator( m_MultiInterpolator[f] );
		metricTmp->SetFixedImage(  m_FixedMultiImage[f]  );
		metricTmp->SetFixedImageRegion(m_FixedMultiImageRegion[f]);
		metricTmp->SetMovingImage( movingImage );

		try
		{
			metricTmp->Initialize();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception catched !" << std::endl;
			std::cerr << excep << std::endl;
			return EXIT_FAILURE;
		}

		m_MultiMetric.push_back(metricTmp);	
	
		
		
	}
		
	//Lectura del archivo de transformacion de la plantilla
	std::string transfPlant = transfMovingImage;

	itk::TransformFileReader::Pointer transformReader = itk::TransformFileReader::New();
	transformReader->SetFileName(transfPlant);

	std::cout<<"Lectura de Transformacion de la Deformacion " << std::endl;
	try{
		transformReader->Update();	
	}
	catch(itk::ExceptionObject &e){
		std::cout << e.GetDescription() << std::endl;
	}

	itk::TransformFileReader::TransformListType* transformList = transformReader->GetTransformList();
	itk::TransformFileReader::TransformPointer baseTransform = transformList->front();
	
	std::cout<<"Parameters  de Transformacion de la Plantilla" << std::endl;
	std::cout<<baseTransform->GetParameters()<<std::endl;


	TransformType::ParametersType similarityParameters;
	similarityParameters = transform->GetParameters();
	
	float *vectParams = new float[7];	
	for(int i=0; i<7; i++){	
		similarityParameters[i] = baseTransform->GetParameters()[i];
		vectParams[i] = similarityParameters[i];
	}
	std::ofstream outputTransfByOneParameter("../CostFunction/MeanSquaresMetricOutput.txt");	
	double valueMetric = 0.0;

	//En caso de que los indices de parametros correspondan a rotaciones
	//Seran dados en grados

	Utilitarios *util = new Utilitarios();
	double vx,vy,vz,newangle;		
	float rx,ry,rz;

	//Precalculate the rotations fo the outTransform in Euler form
	util->convertVersorToEuler(vectParams[0], vectParams[1], vectParams[2],rx,ry,rz);
	std::cout << "Rx: "<<rx<<"Ry: "<< ry << "Rz: " << rz << std::endl; 
	
	util->unirVectorWithAngle(rx,ry,rz,vx,vy,vz,newangle);
	std::cout << "Vx: " << vx << "Vy: " << vy << "Vz: " << vz << std::endl; 
	

	MetricType::Pointer metric = m_MultiMetric[indImagenEval];

	if(strcmp(typePlot, "2D") == 0){
		//Alterar por un solo parametro 
		for( float dx = rangex_min; dx >= rangex_max; dx+=stepx )
		{
			//Limpiamos los versores para una nueva conversion
			vx = 0; vy = 0; vz = 0; newangle=0;
			if(indParam1>=0 && indParam1<=2){
				//Convertir el actual versor a rotacion
				if(indParam1 == 0){
					util->unirVectorWithAngle(dx,ry,rz,vx,vy,vz,newangle);
					similarityParameters[indParam1] = vx;
				}else if(indParam1 == 1){
					util->unirVectorWithAngle(rx,dx,rz,vx,vy,vz,newangle);
					similarityParameters[indParam1] = vy;
				}else if(indParam1 == 2){
					util->unirVectorWithAngle(rx,ry,dx,vx,vy,vz,newangle);
					similarityParameters[indParam1] = vz;

				}
			}else{

				similarityParameters[indParam1] = (dx >= -1 && dx <= 0.0) ? -1 : dx;
				similarityParameters[indParam1] = (dx <=  1 && dx >= 0.0) ? +1 : dx;
			}
			transform->SetParameters(similarityParameters);
			
			//Siempre resetera el valor de la metrica a cero
			valueMetric = 0.0;
			
			if(!combineImages){
				valueMetric = metric->GetValue(similarityParameters);
			}else{
				for(unsigned int itMet=0; itMet < FImgTotal; itMet++){
					valueMetric += m_MultiMetric[itMet]->GetValue(similarityParameters);
				}

			}
			std::cout << dx << "   " << valueMetric << std::endl;
			outputTransfByOneParameter << dx << "   " << valueMetric << std::endl;

		}

		outputTransfByOneParameter.close();
		std::string cmdPlot2d("gnuplot ../CostFunction/MeanSquares2D.gnup");
		std::system(cmdPlot2d.c_str());

	}else if(strcmp(typePlot,"3D")==0){
		//Alterar por un solo parametro 
		for( float dx = rangex_min; dx >= rangex_max; dx+=stepx )
		{
			for( float dy = rangey_min; dy >= rangey_max; dy+=stepy )
			{
				similarityParameters[indParam1] = (dx > -0.1 && dx < 0) ? -0.1 : dx;
				similarityParameters[indParam1] = (dx <  0.1 && dx > 0) ? +0.1 : dx;
				similarityParameters[indParam2] = (dy > -0.1 && dy < 0) ? -0.1 : dy;
				similarityParameters[indParam2] = (dy <  0.1 && dy > 0) ? +0.1 : dy;	
				transform->SetParameters(similarityParameters);
				//Siempre resetera el valor de la metrica a cero
				valueMetric = 0.0;
	
				if(!combineImages){
					valueMetric = metric->GetValue(similarityParameters);
				}else{
					for(unsigned int itMet=0; itMet < FImgTotal; itMet++){
						valueMetric += m_MultiMetric[itMet]->GetValue(similarityParameters);
					}

				}

				std::cout << dx << "   " << dy << "   " << valueMetric << std::endl;
				outputTransfByOneParameter << dx << "   " << dy << "   " << valueMetric << std::endl;
			}
		}
		outputTransfByOneParameter.close();
		std::string cmdPlot3D("gnuplot ../CostFunction/MeanSquaresMetricFigures.gnup");
		std::system(cmdPlot3D.c_str());

	}

	return EXIT_SUCCESS;
}
