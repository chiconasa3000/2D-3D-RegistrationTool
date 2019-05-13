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


void Utilitarios::getBaseElementsVersor(double &vx, double &vy, double &vz, double &rx, double &ry, double &rz){
	//Conseguir el factor K
	double k = sqrt(pow(vx,2) + pow(vy,2) + pow(vz,2));
	//Consiguiendo el angulo en su forma inicial
	float angle = asin(k)*2;
	
	//consiguiendo los axis ejes del versor
	float ax = vx/k;
	float ay = vy/k;
	float az = vz/k;
	
	//usar el convert Versor to Euler con los valores dados
	convertVersorToEuler(ax,ay,az, angle, rx, ry, rz);
		
}

void Utilitarios::convertVersorToEuler(float &x, float &y, float &z, float &w, double &rx, double &ry, double &rz){
	float t0 = +2.0 * (w * x + y * z);
	float t1 = +1.0 - 2.0 * (x * x + y * y);
	rz = atan(t0/ t1);

	float t2 = +2.0 * (w * y - z * x);
	t2 = (t2 > +1.0) ? +1.0 : t2;
	t2 = (t2 < -1.0) ? -1.0 : t2;
	rx = asin(t2);

	float t3 = +2.0 * (w * z + x * y);
	float t4 = +1.0 - 2.0 * (y * y + z * z);
	ry = atan(t3/ t4);
}

void Utilitarios::compareVols(std::string logfilename, std::string  dirnewvolume, std::string  dirimagedef, int indexTest){
	const unsigned int Dimensions = 3;
	typedef itk::Image<short int, Dimensions> imagenDefType;	
	typedef itk::Image<short int, Dimensions> newvolumeType;

	//Leyendo las imagenes
	typedef itk::ImageFileReader<imagenDefType> defReaderType;
	defReaderType::Pointer defreader = defReaderType::New();
	std::string filenameimagedef = dirimagedef + "imagenDef_" + std::to_string(indexTest) + ".mha";
	defreader->SetFileName(filenameimagedef.c_str());
	typedef itk::ImageFileReader<newvolumeType> newvolumeReaderType;
	newvolumeReaderType::Pointer newvolumereader = newvolumeReaderType::New();
	std::string filenamenewvolume = dirnewvolume + std::to_string(indexTest) + "/newVolumen.mha";
	newvolumereader->SetFileName(filenamenewvolume.c_str());

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


	newvolumeType::ConstPointer imagedefnewsize = resample->GetOutput();
	/*
	   std::cout << "Output size: " << newVolume->GetBufferedRegion().GetSize() << std::endl;
	   std::cout << "Spacing: " << newVolume->GetSpacing()<<std::endl;
	   std::cout << "Origin: " << newVolume->GetOrigin()<<std::endl;

	   std::cout << "Output size: " << output->GetBufferedRegion().GetSize() << std::endl;
	   std::cout << "Spacing: " << output->GetSpacing()<<std::endl;
	   std::cout << "Origin: " << output->GetOrigin()<<std::endl; 
	   */
	typedef itk::ImageFileWriter<newvolumeType> WriterType;	
	std::cout << "Writing output... " << std::endl;
	WriterType::Pointer outputWriter = WriterType::New();
	std::string dirsalida = dirnewvolume + std::to_string(indexTest) + "/";
	outputWriter->SetFileName(dirsalida + "imagedefnewsize.mha");
	outputWriter->SetInput(imagedefnewsize);
	outputWriter->Update();

	//aplicar Hausdorff distance entre ambos volumenes
	computeHausdorffDistance(logfilename, imagedefnewsize, newVolume);	


}

void Utilitarios::computeHausdorffDistance(std::string logfilename, typename itk::Image<short int,3>::ConstPointer a, typename itk::Image<short int, 3>::ConstPointer b){
	//Open the logfilename and append this data
	std::ofstream logfile;
	logfile.open(logfilename, std::ios_base::app);	
	
	const unsigned int Dimensions = 3;
	typedef itk::Image<short int, Dimensions> imageType;

	typedef itk::HausdorffDistanceImageFilter < imageType, imageType > HausdorffDistanceImageFilterType;
	HausdorffDistanceImageFilterType::Pointer hausdorffFilter = HausdorffDistanceImageFilterType::New();
	hausdorffFilter->SetInput1(a.GetPointer());
	hausdorffFilter->SetInput2(b.GetPointer());
	/*
	   typedef itk::StatisticsImageFilter < FloatImageType > StatisticsFilterType;
	   StatisticsFilterType ::Pointer statistics1 = StatisticsFilterType::New();
	   statistics1->SetInput(reader->GetOutput());

	   StatisticsFilterType ::Pointer statistics2 = StatisticsFilterType::New();
	   statistics2->SetInput(normalizeFilter->GetOutput());
	   */
	hausdorffFilter->Update();

	logfile << "HausdorffDistance: " << hausdorffFilter->GetHausdorffDistance() << std::endl;
	logfile << "HD_Average: " << hausdorffFilter->GetAverageHausdorffDistance() << std::endl;
	logfile.close();

}
void Utilitarios::createStats(int numLevels, std::string logfilename, std::string dir_resultados, int indexTest){

    //Eliminar datos innecesarios del log de resultados
    std::string cmdPickResults("cat " + logfilename + " | sed -n -e :a -e '1,36!{P;N;D;};N;ba' | sed 1,10d > newlog.txt");
    int result = std::system(cmdPickResults.c_str());

    //Separar el log de resultados por cada nivel de resolucion
    std::string directorioResultados = dir_resultados + std::to_string(indexTest) +"/";
    std::string cmdDivideResults("awk '/^Resolution/{close(file);file = \"" + directorioResultados + "\" $2 $NF \".txt\"; next}  /./{print >> file}' newlog.txt");
    std::system(cmdDivideResults.c_str());

    //Recorrer cada archivo de registro
    for(int i=0;i<numLevels;i++){

        //Solo seleccionar ciertas columnas para el plot
        std::string nameloginLevel(directorioResultados + "level" + std::to_string(i));
        std::string cmdChooseCols("gawk -i inplace '{print $2 \"\\t\" $4 \"\\t\" $6 \"\\t\" $8 \"\\t\" $10}' " + nameloginLevel+".txt");
        std::system(cmdChooseCols.c_str());

        //writing the plot with respective logfile
        std::string namePlotFile(directorioResultados + "plot_temp" + std::to_string(i) + ".gnup");

        std::ofstream plot_temp(namePlotFile);
        plot_temp << "set terminal svg size 600,400 dynamic  enhanced font 'arial,10' mousing name \"metricevolution\" butt dashlength 1.0"  << "\n";
        plot_temp << "set xlabel \"Iteration No.\"" << "\n";
        plot_temp << "set ylabel \"NormalizedGradientCorrelation\" " << "\n";
	plot_temp << "set title \"Metric Evolution Level " << std::to_string(i) << "\""<<"\n";
        plot_temp << "set output \"" + directorioResultados + "ImageRegistrationProgressMetric"+std::to_string(i)+".svg\""<<"\n";
        plot_temp << "plot \"" + nameloginLevel +
            ".txt\" using 1:2 notitle with lines lt 1, \"" + nameloginLevel +
            ".txt\" using 1:2 notitle with points lt 0 pt 12 ps 1" << "\n";

        plot_temp.close();

        std::string nameCmdPlot("gnuplot " + namePlotFile);
        std::system(nameCmdPlot.c_str());

    }
}

void Utilitarios::createStatsOfTransValues(std::string dirRes, std::string logfilename,int numTest){

	//writing the plot with respective logfile
        std::string namePlotFile( dirRes + "plotValueTrans" + std::to_string(numTest) + ".gnup");

	std::ofstream plot_temp(namePlotFile);
	plot_temp << "set terminal svg size 600,400 dynamic  enhanced font 'arial,10' mousing name \"BoxPlotTransformation\" butt dashlength 1.0" << "\n";
	plot_temp << "set output \"" + dirRes + "ValueTraslationScale"+std::to_string(numTest)+".svg\""<<"\n";
	plot_temp << "green = \"#80bfaa\"; skyblue = \"#55a0d5\""<<"\n";
	plot_temp << "set yrange [-20:20]"<<"\n";	
	plot_temp << "set style data histogram"<<"\n";
	plot_temp << "set style histogram cluster gap 1"<<"\n";
	plot_temp << "set style fill solid"<<"\n";
	plot_temp << "set boxwidth 0.9"<<"\n";
	plot_temp << "set xtics format \"\""<<"\n";
	plot_temp << "set grid ytics"<<"\n";
	plot_temp << "set ylabel \"Valor de Parametros de Transformacion (mm)\""<<"\n";
	plot_temp << "set xlabel \"Parametros de Traslacion y Escala\""<<"\n";
	plot_temp << "set title \"Diferencia entre los valores del GroundTruth y el de Registro de los Parametros de Transformacion\""<<"\n";
	plot_temp << "plot \"../outputData/valueTransf"<< std::to_string(numTest) <<".txt" <<"\" using 2:xtic(1) with histogram title \"GroundTruth\" linecolor rgb green, \\"<<"\n";
	plot_temp << "'' using 3 with histogram title \"Registration\" linecolor rgb skyblue, \\"<<"\n";
	plot_temp << "'' u 0:2:2 with labels font \"Helvetica, 10\" offset char -3,0.5 title \"\", \\"<<"\n";
	plot_temp << "'' u 0:3:3 with labels font \"Helvetica, 10\" offset char 3,0.5 title \"\" \\"<<"\n";
	//aqui es lo diferente para cada archivo las 3 ultimas lineas
	//lo demas se repite auque seria bueno solo modificar esto.
	plot_temp.close();

	std::string nameCmdPlot("gnuplot " + namePlotFile);
        std::system(nameCmdPlot.c_str());
	
}


void Utilitarios::createStatsOfErrors(int numImags){
	
	//obtener al mayor error de RMSE_registro para considerar el mayor rango para el spyder
	
	std::string command("sort -t, -nk3 ../outputData/RMSE_Registro.txt | tail -1 | awk '{$1=\"\"; $2=\"\"; print $0}' | tr -s ' '  '\\n' | sort -n | tail -1");
	std::string maxRange = "";
	FILE * stream;
	const int max_buffer = 6;
	char buffer[max_buffer];
	command.append(" 2>&1");
	stream = popen(command.c_str(),"r");
	if(stream && fgets(buffer,max_buffer,stream) != NULL){
		maxRange.append(buffer);
	}

	double nmaxRange = atof(maxRange.c_str()) + 1;		


	std::string defineRange("");	
    defineRange += "sed -i -e '/a1_max =/c\\a1_max = "+ std::to_string(nmaxRange) + "'"+
        " -e '/a2_max =/c\\a2_max = "+ std::to_string(nmaxRange) + "'" +
        " -e '/a3_max =/c\\a3_max = "+ std::to_string(nmaxRange) + "'" +
        " -e '/a4_max =/c\\a4_max = "+ std::to_string(nmaxRange) + "'" +
        " -e '/a5_max =/c\\a5_max = "+ std::to_string(nmaxRange) + "'" +
        " -e '/a6_max =/c\\a6_max = "+ std::to_string(nmaxRange) + "'" +
        " -e '/a7_max =/c\\a7_max = "+ std::to_string(nmaxRange) + "'" + " ../cmdPlots/spyder4.gnup";
	std::system(defineRange.c_str());

	int numColsPlot = numImags + 3;
	std::string fixNumColsCmd("");
	fixNumColsCmd += "sed -i '/do/c\\do for [COL=3:"+std::to_string(numColsPlot)+"] {' ../cmdPlots/spyder4.gnup"; 
	std::system(fixNumColsCmd.c_str());

	std::string fileoutput("");
	fileoutput += "sed -i '/set output/c\\set output sprintf(\"../outputData/spyder\%d.svg\",tag)' ../cmdPlots/spyder4.gnup";
	std::system(fileoutput.c_str());

	//writing the plot with respective logfile
	std::string namePlotFile("../cmdPlots/spyder4.gnup");
	std::string nameCmdPlot("gnuplot " + namePlotFile);
	std::system(nameCmdPlot.c_str());

}

void Utilitarios::createStatsBarHausdorff(){
    std::string cmdDelLastLineHausDist("sed -n -e :a -e '1,1!{P;N;D;};N;ba' ../outputData/HausdorffDistances.txt > ../outputData/hdgeneral.txt");
	std::system(cmdDelLastLineHausDist.c_str());

	std::string namePlotFile("../cmdPlots/barhdist.gnup");
	std::string nameCmdPlot("gnuplot "+namePlotFile);
	std::system(nameCmdPlot.c_str());

}

void Utilitarios::createStatsBoxPlotsTypeTransParams(){
	std::string cmdTransposeRMSE("");
	cmdTransposeRMSE += "cat ../outputData/RMSE_Registro.txt |sed -n -e :a -e '1,1!{P;N;D;};N;ba' | sed 1,1d | awk '{$1=$2=\"\"; print $0}' | awk -f ../cmdPlots/transposeTable.awk > ../outputData/newRMSE_Registro.txt";
	std::system(cmdTransposeRMSE.c_str());
	
	std::string namePlotFile("../cmdPlots/boxplotTranfParams.gnup");
	std::string nameCmdPlot("gnuplot "+namePlotFile);
	std::system(nameCmdPlot.c_str());

}
