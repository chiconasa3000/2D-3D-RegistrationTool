/* Generador de imagenes virtuales usando interpolador Patched */

#include <iostream>
#include "itkImage.h"
#include "itkTimeProbesCollectorBase.h"
#include "itkImageFileReader.h"
#include "itkResampleImageFilter.h"
#include "itkSimilarity3DTransform.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkImageFileWriter.h"
#include "../itkPatchedRayCastInterpolateImageFunction.h"
#include "itkMatrix.h"
#include "gestorMatrixRot.h"
#include "itkVersor.h"
#include <string>
#include <itksys/SystemTools.hxx>
// funcion de menu principal 

void menu(){
	std::cout<<"Generador de imagenes virtuales PATCHED GENERATOR"<< std::endl;
}

int main(int argc, char *argv[]){

	//variable definition
	char *input_name = NULL;		//input volume
	char *output_name = NULL;		//virtual image
	char *type_projection = NULL;		//tipo de proyeccion [AP(anteroposterior) , ML(mediolateral)]

	float scd = 1000.0;			//distance from source to isocenter

	// CT volume rotation around isocenter along x,y,z axis in degrees
	float rx = 0.;
	float ry = 0.;
	float rz = 0.;

	//Direction Cosines (rotation for orientation volume)
	float dcx = 0.;
	float dcy = 0.;
	float dcz = 0.;

	// Translation parameter of the isocenter in mm
	float tx = 0.;
	float ty = 0.;
	float tz = 0.;

	bool customized_iso  = false; 		//flag for iso given by user
	bool customized_2DCX = false;		//flag for central of 2d image
	float cx = 0.;				//virtual image isocenter in x
	float cy = 0.;				//virtual image isocenter in y
	float cz = 0.;				//virtual image isocenter in z

	float threshold = 0.;			//virtual image threshold

	int dxx = 452;				//pixels number virtual image in x
	int dyy = 225;				//pixels number virtual image in y

	float im_sx = 1;			//virtual image spacing in x
	float im_sy = 1;			//virtual image spacing in y

	float o2Dx;				//virtual image origin in x
	float o2Dy;				//virtual image origin in y

	float focalPointx = 0.;			//focalPoint in x
	float focalPointy = 1000.0;		//focalPoint in y
	float focalPointz = 0.;			//focalPoint in z


	bool ok;				//sanity of parameters
	float rprojection = 0.; 		//projection angle: AP view by default
	bool verbose = false;			//information flag

	itk::TimeProbesCollectorBase timer;	//Time Record


	//initialization variables
	//In order to evaluate every parameter the condition flag is always false
	//so whatever quantity of parameters the last parameter put on true the flag
	//condition
	//
	//So, the first argument argv[0] is the executable command
	while(argc > 1){
		ok = false;

		if ((ok == false) && (strcmp(argv[1], "-v") == 0))
		{
			argc--; argv++;
			ok = true;
			verbose = true;
		}

		if ((ok == false) && (strcmp(argv[1], "-p") == 0))
		{
			argc--; argv++;
			ok = true;
			type_projection = argv[1];
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


		if ((ok == false) && (strcmp(argv[1], "-rp") == 0))
		{
			argc--; argv++;
			ok = true;
			rprojection=atof(argv[1]);
			argc--; argv++;
		}

		if ((ok == false) && (strcmp(argv[1], "-iso") == 0))
		{
			argc--; argv++;
			ok = true;
			cx=atof(argv[1]);
			argc--; argv++;
			cy=atof(argv[1]);
			argc--; argv++;
			cz=atof(argv[1]);
			argc--; argv++;
			customized_iso = true;
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

		if ((ok == false) && (strcmp(argv[1], "-scd") == 0))
		{
			argc--; argv++;
			ok = true;
			scd = atof(argv[1]);
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
		if ((ok == false) && (strcmp(argv[1], "-o") == 0))
		{
			argc--; argv++;
			ok = true;
			output_name = argv[1];
			argc--; argv++;
		}

		if (ok == false) 
		{

			if (input_name == NULL) 
			{
				input_name = argv[1];
				argc--;
				argv++;
			}
			else 
				std::cerr << "ERROR: Can not parse the image " << argv[1] << std::endl;

		}


	}
	if (verbose) 
	{
		if (input_name)  std::cout << "Input image: "  << input_name  << std::endl;
		if (output_name) std::cout << "Output image: " << output_name << std::endl;
	}


	//function con las terceras parte
	menu();	
	const unsigned int Dimension = 3;

	//Parameters of fixed image


	//tipo de pixel by default para las imagenes
	typedef short int InputPixelType;
	typedef unsigned char OutputPixelType;

	//typedef itk::Image<short int, Dimensions> FixedImageType;
	typedef itk::Image<InputPixelType, Dimension> MovingImageType;
	typedef itk::Image<OutputPixelType, Dimension> OutputImageType;

	MovingImageType::Pointer image;

	// Reading the Moving Image
	if (input_name) 
	{
		timer.Start("Loading Input Image");
		typedef itk::ImageFileReader< MovingImageType >  ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( input_name );

		try { 
			reader->Update();
		} 

		catch( itk::ExceptionObject & err ) 
		{ 
			std::cerr << "ERROR: ExceptionObject caught !" << std::endl; 
			std::cerr << err << std::endl; 
			return EXIT_FAILURE;
		} 

		image = reader->GetOutput();
		timer.Stop("Loading Input Image");
	}
	else 
	{
		std::cerr << "Input image file missing !" << std::endl;
		return EXIT_FAILURE;
	}


	//force origin to cero position
	/*	MovingImageType::PointType ctOrigin;
		ctOrigin[0] = 0.0;
		ctOrigin[1] = 0.0;
		ctOrigin[2] = 0.0;
		image->SetOrigin(ctOrigin);
		*/
	if (verbose) 
	{
		unsigned int i;
		const itk::Vector<double, 3> spacing = image->GetSpacing();  
		std::cout << std::endl << "Input ";

		MovingImageType::RegionType region = image->GetBufferedRegion();
		region.Print(std::cout);

		std::cout << "  Resolution: [";
		for (i=0; i<Dimension; i++) 
		{
			std::cout << spacing[i];
			if (i < Dimension-1) std::cout << ", ";
		}
		std::cout << "]" << std::endl;

		const itk::Point<double, 3> origin = image->GetOrigin();
		std::cout << "  Origin: [";
		for (i=0; i<Dimension; i++) 
		{
			std::cout << origin[i];
			if (i < Dimension-1) std::cout << ", ";
		}
		std::cout << "]" << std::endl<< std::endl;
		//std::cout<< "Isocenter: " <<isocenter[0] <<", "<< isocenter[1] << ", "<< isocenter[2] << std::endl;
		//std::cout << "Transform: " << transform << std::endl;
	}

	//The resample filter enables coordinates for each of the pixels in DRR image.
	//these coordinates are used by interpolator to determine the equatio of each

	//ray which trough the volume

	typedef itk::ResampleImageFilter<MovingImageType, MovingImageType> FilterType;

	FilterType::Pointer filter = FilterType::New();
	filter->SetInput(image);
	filter->SetDefaultPixelValue(0);

	typedef itk::Similarity3DTransform< double > TransformType;
	TransformType::Pointer transform = TransformType::New();
	transform->SetIdentity();	
	//constant for casting degrees into radians format of rotation projection
	const double dtr = ( atan(1.0) * 4.0 ) / 180.0;

	//set general transformations transform
	TransformType::OutputVectorType translation;

	translation[0] = tx;
	translation[1] = ty;
	translation[2] = tz;

	transform->SetTranslation( translation );


	HelperRot helperRot; 
	HelperRot helperRotInit;

	helperRot.initRotX(dtr*(dcx));
	helperRot.initRotY(dtr*(dcy));
	helperRot.initRotZ(dtr*(dcz));

	helperRotInit.initRotX(dtr*(rx));
	helperRotInit.initRotY(dtr*(ry));
	helperRotInit.initRotZ(dtr*(rz));
	
	helperRot.composeMatrixRot();
	helperRotInit.composeMatrixRot();

	std::cout<<"Matriz de Rotacion para la imagen Movible"<<std::endl;
	helperRotInit.printRotg();
	std::cout<<std::endl;

	//itk::Versor<double> vectRot;
	//vectRot.SetRotationAroundX(dtr*(-90.0));	
	transform->SetMatrix(helperRotInit.getRotg());


	//Read image properties in order to build our isocenter
	MovingImageType::PointType imOrigin = image->GetOrigin();
	MovingImageType::SpacingType imRes = image->GetSpacing();

	typedef MovingImageType::RegionType InputImageRegionType;
	typedef MovingImageType::SizeType InputImageSizeType;

	InputImageRegionType imRegion = image->GetBufferedRegion();
	InputImageSizeType imSize = imRegion.GetSize();

	TransformType::InputPointType isocenter;

	//Consider the center of the volume

	if (customized_iso)
	{
		// Isocenter location given by the user.
		isocenter[0] = imOrigin[0] + imRes[0] * cx; 
		isocenter[1] = imOrigin[1] + imRes[1] * cy; 
		isocenter[2] = imOrigin[2] + imRes[2] * cz;
	}
	else
	{
		// Set the center of the image as the isocenter.
		isocenter[0] = imOrigin[0] + imRes[0] * static_cast<double>( imSize[0] ) / 2.0; 
		isocenter[1] = imOrigin[1] + imRes[1] * static_cast<double>( imSize[1] ) / 2.0; 
		isocenter[2] = imOrigin[2] + imRes[2] * static_cast<double>( imSize[2] ) / 2.0;
	}

	//transform->SetCenter(isocenter);

	//Instance of the interpolator
	typedef itk::PatchedRayCastInterpolateImageFunction<MovingImageType, double> InterpolatorType;
	typedef InterpolatorType::InputPointType FocalPointType;

	InterpolatorType::Pointer interpolator = InterpolatorType::New();

	//REMEMBER: It is according to Interpolator Class
	//take care about this
	interpolator->SetThreshold(threshold);

	//FocalPoint
	FocalPointType focalPoint;
	focalPoint[0] = focalPointx;
	focalPoint[1] = focalPointy;
	focalPoint[2] = focalPointz;
	//Set the focal point in interpolator
	interpolator->SetFocalPoint(focalPoint);

	//TODO: class need this function
	//interpolator->SetProjectionAngle( dtr * rprojection );

	interpolator->SetTransform(transform);
	//TODO: class need this function
	//interpolator->Initialize();

	//insert the interpolator into the filter
	filter->SetInterpolator(interpolator);

	//Setting properties of fixed image

	MovingImageType::SizeType size;
	double spacing[Dimension];
	double origin[Dimension];

	//Number of pixels of Fixed Image
	size[0] = dxx;
	size[1] = dyy;
	size[2] = 1;

	//Resolution of Fixed Image
	spacing[0] = im_sx;
	spacing[1] = im_sy;
	spacing[2] = 1.0;

	if(!customized_2DCX)
	{
		//center virtual image by default
		o2Dx = ((double) dxx - 1.)/2.;
		o2Dy = ((double) dyy - 1.)/2.;
	}

	//TODO: Need to classify orientation types (sign and scd)
	//Hint: consider the origin like a other input parameter
	if(strcmp(type_projection,"AP")==0){
		// Compute the origin (inmm) of the 2D Image
		origin[0] = im_sx * o2Dx;
		origin[1] = -scd;
		origin[2] = im_sy * o2Dy;
	}else if(strcmp(type_projection,"ML")==0){
		// Compute the origin (inmm) of the 2D Image
		origin[0] = -scd;
		origin[1] = -im_sx * o2Dx;
		origin[2] = im_sy * o2Dy;
	}

	//set identity in direction cosine
	itk::Versor< double > rotationDirect;
	const double angleRadians = (-180.0) * vnl_math::pi/180.0;
	rotationDirect.SetRotationAroundX(angleRadians);
	const OutputImageType::DirectionType direction = image->GetDirection();
	const OutputImageType::DirectionType newDirection = direction * helperRot.getRotg();
	//OutputImageType::DirectionType direction;
	//direction.SetIdentity();

	filter->SetOutputDirection(newDirection);
	//set properties of the virtual image
	filter->SetSize(size);
	filter->SetOutputSpacing(spacing);
	filter->SetOutputOrigin(origin);

	//transform in the filter
	filter->SetTransform(transform);

	//Virtual Image Properties information
	if(verbose)
	{
		std::cout << "Rotation: " 
			<< rx << ", " 
			<< ry << ", " 
			<< rz << std::endl;

		std::cout << "Traslation: " 
			<< tx << ", " 
			<< ty << ", " 
			<< tz << std::endl;


		std::cout << "Output image size: " 
			<< size[0] << ", " 
			<< size[1] << ", " 
			<< size[2] << std::endl;

		std::cout << "Output image spacing: " 
			<< spacing[0] << ", " 
			<< spacing[1] << ", " 
			<< spacing[2] << std::endl;

		std::cout << "Output image origin: "
			<< origin[0] << ", " 
			<< origin[1] << ", " 
			<< origin[2] << std::endl;
		std::cout << "Focal point image: "
			<< focalPoint[0] << ", " 
			<< focalPoint[1] << ", " 
			<< focalPoint[2] << std::endl;

	}

	//applying the resample filter (generation)
	timer.Start("DRR generation");
	filter->Update();
	timer.Stop("DRR generation");

	if (output_name) 
	{
		timer.Start("DRR post-processing");

		// The output of the filter can then be passed to a writer to
		// save the DRR image to a file.

		typedef itk::RescaleIntensityImageFilter< 
		MovingImageType, OutputImageType > RescaleFilterType;
		RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
		rescaler->SetOutputMinimum(   0 );
		rescaler->SetOutputMaximum( 255 );
		rescaler->SetInput( filter->GetOutput() );
		rescaler->Update();

		
		// Out of some reason, the computed projection is upsided-down.
		// Here we use a FilpImageFilter to flip the images in y direction.
		/*typedef itk::FlipImageFilter< OutputImageType > FlipFilterType;
		FlipFilterType::Pointer flipFilter = FlipFilterType::New();

		typedef FlipFilterType::FlipAxesArrayType FlipAxesArrayType;
		FlipAxesArrayType flipArray;
		flipArray[0] = 0;
		flipArray[1] = 0;

		flipFilter->SetFlipAxes( flipArray );
		flipFilter->SetInput( rescaler->GetOutput() );
		flipFilter->Update();
		*/
		timer.Stop("DRR post-processing");
		

		/*typedef itk::ImageFileWriter< OutputImageType >  WriterType;
		  WriterType::Pointer writer = WriterType::New();

		// Now we are ready to write the projection image.
		writer->SetFileName( output_name );
		writer->SetInput( rescaler->GetOutput() );

		try 
		{ 
		std::cout << "Writing image: " << output_name << std::endl;
		writer->Update();
		} 
		catch( itk::ExceptionObject & err ) 
		{ 
		std::cerr << "ERROR: ExceptionObject caught !" << std::endl; 
		std::cerr << err << std::endl; 
		} */
		std::string fname;
	        typedef itk::ImageFileWriter< OutputImageType >  WriterType;
		  WriterType::Pointer writer = WriterType::New();

		itksys::SystemTools::MakeDirectory( (fname + "../outputData"+"/virtualImages/").c_str() );

		fname = fname + "../outputData" + "/virtualImages/" + output_name;
		if(Dimension == 2)
		{
			fname += ".png";
		}
		else
		{
			fname += ".mha";
		}
		writer->SetFileName( fname.c_str() );
		//writer->SetInput(dynamic_cast<const MovingImageType *>(filter->GetOutput()));
		writer->SetInput(rescaler->GetOutput());

		try{
			std::cout << "Writing Virtual Image" << fname.c_str() << std::endl;
			writer->Update();
		}catch( itk::ExceptionObject & err ) 
		{ 
			std::cerr << "ERROR: ExceptionObject caught !" << std::endl; 
			std::cerr << err << std::endl; 
		} 
	}

	timer.Report();
	return 0;
}

