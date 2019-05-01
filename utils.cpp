

#include "utils.h"


Utilitarios::Utilitarios(){
}

void Utilitarios::convertEulerToVersor(float &rx, float &ry, float &rz, double &ax, double &ay, double &az, double &angle){
	
	const double dtr = (atan(1.0) * 4.0)/180.0;

	//Conversion de Euler a Versor
	double c1 = cos(rx*dtr/2);
	double s1 = sin(rx*dtr/2);

	double c2 = cos(ry*dtr/2);
	double s2 = sin(ry*dtr/2);

	double c3 = cos(rz*dtr/2);
	double s3 = sin(rz*dtr/2);

	double aw = c1*c2*c3 - s1*s2*s3;
	ax = c1*c2*s3 + s1*s2*c3;
	ay = s1*c2*c3 + c1*s2*s3;
	az = c1*s2*c3 - s1*c2*s3;

	angle = 2*acos(aw);
	double norm = ax*ax+ay*ay+az*az;

	if(norm < 0.001){
		az = 1;
		ay = ax = 0;
	}else{
		norm = sqrt(norm);
		ax /= norm;
		ay /= norm;
		az /= norm;
	}
}

void Utilitarios::unirVectorWithAngle(float &rx, float &ry, float &rz, double &vx, double &vy, double &vz, double &newangulo){

	double ax,ay,az,angle;
	//Usar la funcion de conversion de Euler a Versor
	convertEulerToVersor(rx,ry,rz,ax,ay,az,angle);
	double axisnorm = sqrt(ax*ax + ay*ay + az*az);
	if(axisnorm == 0.0){
		std::cerr<< "Error norma de vector del versor es cercano a  cero"<<std::endl;
	}else{

	double cosangle = cos(angle/2.0);
	double sinangle = sin(angle/2.0);
	
	double factor = sinangle/ axisnorm;

	vx = ax * factor;
	vy = ay * factor;
	vz = az * factor;
	
	newangulo = cosangle;	}
}

void Utilitarios::convertVersorToEuler(float &x, float &y, float &z, float &w, float &rx, float &ry, float &rz){
	float t0 = +2.0 * (w * x + y * z);
        float t1 = +1.0 - 2.0 * (x * x + y * y);
        rx = atan(t0/ t1);

        float t2 = +2.0 * (w * y - z * x);
 	t2 = (t2 > +1.0) ? +1.0 : t2;
	t2 = (t2 < -1.0) ? -1.0 : t2;
        ry = asin(t2);

        float t3 = +2.0 * (w * z + x * y);
        float t4 = +1.0 - 2.0 * (y * y + z * z);
        rz = atan(t3/ t4);
}

void Utilitarios::compareVols(std::string  newvolume, std::string  imagendef, int numTest){
	const unsigned int Dimensions = 3;
	typedef itk::Image<short int, Dimensions> imagenDefType;	
	typedef itk::Image<short int, Dimensions> newvolumeType;

	//Leyendo las imagenes
	typedef itk::ImageFileReader<imagenDefType> defReaderType;
	defReaderType::Pointer defreader = defReaderType::New();
	defreader->SetFileName(imagendef.c_str());

	typedef itk::ImageFileReader<newvolumeType> newvolumeReaderType;
	newvolumeReaderType::Pointer newvolumereader = newvolumeReaderType::New();
	newvolumereader->SetFileName(newvolume.c_str());

	//Variable para las imagenes

	typedef imagenDefType::ConstPointer imageDefConstPointer;
	typedef newvolumeType::ConstPointer newvolConstPointer;


	try{
		defreader->Update();
	}catch(itk::ExceptionObject & e){
		std::cout << e.GetDescription() << std::endl;
		return;
	}

	try{
		newvolumereader->Update();
	}catch(itk::ExceptionObject & e){
		std::cout << e.GetDescription() << std::endl;
		return;
	}

	//Capturando en variables
	imageDefConstPointer imageDeformable = defreader->GetOutput();
	newvolConstPointer newVolume = newvolumereader->GetOutput();	


	typedef itk::LinearInterpolateImageFunction<newvolumeType,double> InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();

	typedef itk::Euler3DTransform <double> TransformType;
	TransformType::Pointer etransform = TransformType::New();
	etransform->SetComputeZYX(true);

	TransformType::OutputVectorType translation;
	translation[0] = 0.0;
	translation[1] = 0.0;
	translation[2] = 0.0;		
	etransform->SetTranslation(translation);

	typedef itk::ResampleImageFilter<imagenDefType, newvolumeType > ResampleImageFilterType;
	ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();
	resample->SetTransform(etransform);
	resample->SetInterpolator(interpolator);

	resample->SetSize(newVolume->GetBufferedRegion().GetSize());
	resample->SetOutputOrigin(newvolumereader->GetOutput()->GetOrigin());
	resample->SetOutputSpacing(newvolumereader->GetOutput()->GetSpacing());
	resample->SetOutputDirection(newvolumereader->GetOutput()->GetDirection());
	resample->SetDefaultPixelValue(0);

	resample->SetInput(imageDeformable);


	newvolumeType::Pointer output = resample->GetOutput();

	std::cout << "Output size: " << newVolume->GetBufferedRegion().GetSize() << std::endl;
	std::cout << "Spacing: " << newVolume->GetSpacing()<<std::endl;
	std::cout << "Origin: " << newVolume->GetOrigin()<<std::endl;

	std::cout << "Output size: " << output->GetBufferedRegion().GetSize() << std::endl;
	std::cout << "Spacing: " << output->GetSpacing()<<std::endl;
	std::cout << "Origin: " << output->GetOrigin()<<std::endl; 

	typedef itk::ImageFileWriter<newvolumeType> WriterType;	
	std::cout << "Writing output... " << std::endl;
	WriterType::Pointer outputWriter = WriterType::New();
	outputWriter->SetFileName("output.mha");
	outputWriter->SetInput(output);
	outputWriter->Update();




}
