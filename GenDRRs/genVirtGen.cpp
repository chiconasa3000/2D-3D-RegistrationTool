#include "genVirtGen.h"

GenVirtGen::GenVirtGen(char *typeProjection,char *logFileName){
	type_projection = typeProjection;
	setLogFile(logFileName);
}

int GenVirtGen::readMovingImage(char *input_name){
	//comprobar existencia de imagen movible para proyeccion
	if(input_name){
	
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
		//std::cout<<"Imagen Movible a proyectar" << image<<std::endl;
		timer.Stop("Loading Input Image");

		unsigned int i;
		const itk::Vector<double, 3> spacing = image->GetSpacing();  
		MovingImageType::RegionType region = image->GetBufferedRegion();

		logregistro << "  Resolution: [";	
		for (i=0; i<Dimension; i++) 
		{
			logregistro << spacing[i];
			if (i < Dimension-1){ logregistro << ", ";}
		}
		logregistro << "]" << std::endl;

		const itk::Point<double, 3> origin = image->GetOrigin();
		logregistro << "  Origin: ["; 
		for (i=0; i<Dimension; i++) 
		{
			logregistro << origin[i];	
			if (i < Dimension-1) 
				logregistro << ", ";
		}
		logregistro << "]" << std::endl<< std::endl;

		return EXIT_SUCCESS;

	}else{
		std::cerr << "Input image file missing !" << std::endl;
		return EXIT_FAILURE;

	}	
}


void GenVirtGen::setLogFile(char *logfile){
	this->logregistro = std::ofstream(logfile);	
}

void GenVirtGen::initResampleFilter(int &dxx, int &dyy, float &im_sx, float &im_sy, bool customized_2DCX, float o2Dx, float o2Dy,float scd){
	filter = FilterType::New();
	filter->SetInput(image);
	filter->SetDefaultPixelValue(0);
	
	//insert the interpolator into the filter
	filter->SetInterpolator(interpolator);
	
	//Seteando dimension y espaciado de acuerdo a los parametros ingresados
	copyPropertiesInputImage(dxx,dyy,im_sx,im_sy);	


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
		origin[0] = +im_sx * o2Dx;
		origin[1] = +scd;
		origin[2] = -im_sy * o2Dy;
	}else if(strcmp(type_projection,"ML")==0){
		// Compute the origin (inmm) of the 2D Image
		origin[0] = -scd;
		origin[1] = +im_sx * o2Dx;
		origin[2] = -im_sy * o2Dy;
	}else{
		//byDefault
		origin[0] = im_sx * o2Dx;
		origin[1] = im_sy * o2Dy;
		origin[2] = scd;

	}

	//set identity in direction cosine
	const MovingImageType::DirectionType direction = image->GetDirection();
	const MovingImageType::DirectionType newDirection = direction * helperRot.getRotg();
	filter->SetOutputDirection(newDirection);
	
	//Set properties of the virtual image
	filter->SetSize(size);
	filter->SetOutputSpacing(spacing);
	filter->SetOutputOrigin(origin);

	//transform in the filter
	filter->SetTransform(transform);

	//applying the resample filter (generation)
	timer.Start("DRR generation");
	filter->Update();
	timer.Stop("DRR generation");

	logregistro << "Output image size: " << size[0] << ", " << size[1] << ", " << size[2] << std::endl;
	logregistro << "Output image spacing: " << spacing[0] << ", " << spacing[1] << ", " << spacing[2] << std::endl;
	logregistro << "Output image origin: "<< origin[0] << ", " << origin[1] << ", " << origin[2] << std::endl;

}

void GenVirtGen::initTransform(float tx, float ty, float tz, float dcx, float dcy, float dcz,float rx, float ry, float rz,float sg){
	transform = TransformType::New();
	transform->SetIdentity();	
	
	//Translate vector in Similarity Transform
	TransformType::OutputVectorType translation;

	translation[0] = tx;
	translation[1] = ty;
	translation[2] = tz;

	transform->SetTranslation( translation );
	//std::cout<<"Transform: "<<transform << std::endl;
	
	//Matriz de Rotacion de Dirección Euleriana a pesar de que es de similaridad (Ojo but Works)
	helperRot.initRotX(dtr*(dcx));
	helperRot.initRotY(dtr*(dcy));
	helperRot.initRotZ(dtr*(dcz));
	helperRot.composeMatrixRot();
	//helperRot.printRotg();

	//Rotacion Normal para el volumen en Versor
	typedef TransformType::VersorType VersorType;
	typedef VersorType::VectorType VectorType;
	VersorType rotation;
	VectorType axis;

	Utilitarios *util = new Utilitarios();
	
	util->convertEulerToVersor(rx,ry,rz,ax,ay,az,angle);		
	axis[0] = ax; axis[1] = ay; axis[2] = az;
	rotation.Set(axis, angle);
	transform->SetRotation(rotation);
	
	//Escala en Similarity Transform
	TransformType::ScaleType scale;
	scale = sg;
	transform->SetScale(scale);	
	similarityParameters = transform->GetParameters();

	logregistro << "Rotation: " << rx << ", " << ry << ", " << rz << std::endl;
	logregistro << "Versor: "<< similarityParameters[0] << ", " << similarityParameters[1] << ", "<< similarityParameters[2] << std::endl;
	logregistro << "Traslation: " << similarityParameters[3] << ", " << similarityParameters[4] << ", " << similarityParameters[5] << std::endl;
	logregistro << "Scale " << similarityParameters[6] << std::endl;
}

void GenVirtGen::copyPropertiesInputImage(int &dxx, int &dyy, float &im_sx, float &im_sy){
	//Read image properties in case of user don't insert size, spacing origin or isocenter
	MovingImageType::PointType imOrigin = image->GetOrigin();
	MovingImageType::SpacingType imRes = image->GetSpacing();

	typedef MovingImageType::RegionType InputImageRegionType;
	typedef MovingImageType::SizeType InputImageSizeType;

	InputImageRegionType imRegion = image->GetBufferedRegion();
	InputImageSizeType imSize = imRegion.GetSize();

	if(im_sx == 0 && im_sy == 0 && dxx == 0 && dyy == 0){
		//capturar del propio volumen
		if(strcmp(type_projection,"AP")==0){
			dxx = imSize[0];
			dyy = imSize[2];
			im_sx = imRes[0];
			im_sy = imRes[2];
		}else
		if(strcmp(type_projection,"ML")==0){
			dxx = imSize[1];
			dyy = imSize[2];
			im_sx = imRes[1];
			im_sy = imRes[2];
		}else{
			dxx = 500;
			dyy = 500;
			im_sx = 1;
			im_sy = 1;

		}	
	}


}


void GenVirtGen::initInterpolator(float focalPointx, float focalPointy, float focalPointz, float threshold){
	interpolator = InterpolatorType::New();

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
	interpolator->SetTransform(transform);
	logregistro << "Focal point image: "<< focalPoint[0] << ", " << focalPoint[1] << ", " << focalPoint[2] << std::endl;


}

void GenVirtGen::writeResultImage(char *output_name){

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
		
		timer.Stop("DRR post-processing");
		
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
	timer.ExpandedReport(logregistro);	
	logregistro.close();

}




itk::Image<short int,3>::Pointer GenVirtGen::returnResultImage(char *output_name){

	if (output_name) 
	{
		
		return filter->GetOutput();	
	}	
}
