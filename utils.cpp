#include "utils.h"

Utilitarios::Utilitarios(){
}


/*void Utilitarios::getSizeAndSpacingFromImage(std::string imageFileName){

	

}*/

void Utilitarios::round2Decimals(float &number){
	number = ceil(number*100)/100;
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

	double aw = c1*c2*c3 + s1*s2*s3;//cambio de signo
	ax = c1*c2*s3 - s1*s2*c3;//cambio de signo
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
	//double k = sqrt(pow(vx,2) + pow(vy,2) + pow(vz,2));
	//Consiguiendo el angulo en su forma inicial
	//float angle = asin(k)*2;
	//consiguiendo los axis ejes del versor
	//float ax = vx/k;
	//float ay = vy/k;
	//float az = vz/k;	
	//double newangle = cos(angle/2.0);
	//usar el convert Versor to Euler con los valores dados
	//convertVersorToEuler(vx,vy,vz, newangle, rx, ry, rz);
	//Convertir de radianes a grados


}


void Utilitarios::convertVersorToEuler(float &vx, float &vy, float &vz, float &rx, float &ry, float &rz){

	//Obteniendo el angulo nuevo con respecto a Vx, Vy, Vz es decir el super W :P
	double k = sqrt(pow(vx,2) + pow(vy,2) + pow(vz,2));
	float angle = asin(k)*2;
	double w = cos(angle/2.0); //super W

	double t0 = +2.0 * (w * vx + vy * vz);
	double t1 = +1.0 - 2.0 * (vx * vx + vy * vy);
	rz = atan2(t0, t1);

	double t2 = +2.0 * (w * vy - vz * vx);
	t2 = (t2 > +1.0) ? +1.0 : t2;
	t2 = (t2 < -1.0) ? -1.0 : t2;
	rx = asin(t2);

	double t3 = +2.0 * (w * vz + vx * vy);
	double t4 = +1.0 - 2.0 * (vy * vy + vz * vz);
	ry = atan2(t3,t4);

	//convert radians to degrees
	rx = rx * rtd;
	ry = ry * rtd;
	rz = rz * rtd;	

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

	std::ofstream plot_temp;

	//Recorrer cada archivo de registro
	for(int i=0;i<numLevels;i++){
		//Solo seleccionar ciertas columnas para el plot
		std::string nameloginLevel(directorioResultados + "level" + std::to_string(i));
		std::string cmdChooseCols("gawk -i inplace '{$1=\"\"; $3=\"\"; $5=\"\"; print $0}' " + nameloginLevel + ".txt");
		//std::string cmdChooseCols("gawk -i inplace '{print $2 \"\\t\" $4 \"\\t\" $6 \"\\t\" $8 \"\\t\" $10}' " + nameloginLevel+".txt");
		std::system(cmdChooseCols.c_str());

		//writing el plot para Metrica vs Iteracion
		std::string namePlotFile(directorioResultados + "plot_temp" + std::to_string(i) + ".gnup");

		plot_temp.open(namePlotFile, std::ios_base::app);
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

		//Plot para la Trans vs Iteración (Traslación Z)

		std::string plotTraceTransf(directorioResultados + "traceTrans" + std::to_string(i) + ".gnup");

		plot_temp.open(plotTraceTransf, std::ios_base::app);
		plot_temp << "set terminal svg size 600,400 dynamic  enhanced font 'arial,10' mousing name \"traceevolution\" butt dashlength 1.0"  << "\n";
		plot_temp << "set xlabel \"Rotacion en X.\"" << "\n";
		plot_temp << "set ylabel \"Rotacion en Y\" " << "\n";
		plot_temp << "set title \"Param Transf Evolution Level " << std::to_string(i) << "\""<<"\n";
		plot_temp << "set parametric"<<"\n";
		plot_temp << "set size square"<<"\n";
		plot_temp << "set output \"" + directorioResultados + "TraceTransf"+std::to_string(i)+".svg\""<<"\n";
		plot_temp << "plot \"" + nameloginLevel +
			".txt\" using 3:4 notitle with lines lt 1, \"" + nameloginLevel +
			".txt\" using 3:4 notitle with points lt 0 pt 12 ps 1" << "\n";
		plot_temp.close();
	
		nameCmdPlot = "gnuplot " + plotTraceTransf;
		std::system(nameCmdPlot.c_str());



	}
}

void Utilitarios::createStatsOfTransValues(std::string dirRes, std::string logfilename,int numTest){

	//writing the plot with respective logfile
	//dirRes = ../outputData/resultsReg0
	std::string nameOutputImagePlot("ValueTraslationScale" + std::to_string(numTest) + ".svg");	
	std::string newPlotParamsValuesDir("set output \"" + dirRes + nameOutputImagePlot);
	std::string namePlotFile(dirRes + "plotValueTrans"  + std::to_string(numTest) + ".gnup");

	//Replace the name of the plot and the outputname graph 
	std::string replaceNewNamePlot("");
	std::string replaceNameOutputPlotCmd("gsed '/set output/c\\" + newPlotParamsValuesDir + "\"' ../cmdPlots/plotValueTransGen.gnup > " + namePlotFile);
	std::system(replaceNameOutputPlotCmd.c_str());

	//Change the name of inputFile respectivily
	std::string nameNewReadingFilePlot("plot \"../outputData/resultsReg_"+std::to_string(numTest) + "/valueTransf" + std::to_string(numTest) + ".txt\"\\\\");
	std::string replaceInputFileNameCmd("gsed -i '/plot/c\\" + nameNewReadingFilePlot + "' "+namePlotFile);
	std::system(replaceInputFileNameCmd.c_str());

	std::string nameCmdPlot("gnuplot " + namePlotFile);
	std::system(nameCmdPlot.c_str());

}


void Utilitarios::createStatsOfErrors(int numImags){

	//obtener al mayor error de RMSE_registro para considerar el mayor rango para el spyder

	//std::string command("sort -nk3 ../outputData/RMSE_Registro.txt | tail -1 | awk '{$1=\"\"; $2=\"\"; print $0}' | tr -s ' '  '\\n' | sort -n | tail -1");
	std::string command(   "cat ../outputData/RMSE_Registro.txt | sed 1,1d | awk '{$1=$2=\"\"; print $0}' | tr -s ' '  '\\n' | sort -n | tail -1");
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
	defineRange += "gsed -i 's/a1_max =.*/a1_max = "+ std::to_string(nmaxRange) + "/g"+
		"; s/a2_max =.*/a2_max = "+ std::to_string(nmaxRange) + "/g" +
		"; s/a3_max =.*/a3_max = "+ std::to_string(nmaxRange) + "/g" +
		"; s/a4_max =.*/a4_max = "+ std::to_string(nmaxRange) + "/g" +
		"; s/a5_max =.*/a5_max = "+ std::to_string(nmaxRange) + "/g" +
		"; s/a6_max =.*/a6_max = "+ std::to_string(nmaxRange) + "/g" +
		"; s/a7_max =.*/a7_max = "+ std::to_string(nmaxRange) + "/g'" + " ../cmdPlots/spyder4.gnup";
	std::system(defineRange.c_str());

	int numColsPlot = numImags + 3;
	std::string fixNumColsCmd("");
	fixNumColsCmd += "gsed -i '/do/c\\do for [COL=3:"+std::to_string(numColsPlot)+"] {' ../cmdPlots/spyder4.gnup";
	std::system(fixNumColsCmd.c_str());

	std::string fileoutput("");
	fileoutput += "gsed -i '/set output/c\\set output sprintf(\"../outputData/spyder\%d.svg\",tag)' ../cmdPlots/spyder4.gnup";
	std::system(fileoutput.c_str());

	//writing the plot with respective logfile
	std::string namePlotFile("../cmdPlots/spyder4.gnup");
	std::string nameCmdPlot("gnuplot " + namePlotFile);
	std::system(nameCmdPlot.c_str());

}

void Utilitarios::createStatsBarHausdorff(){
	std::string cmdDelLastLineHausDist("gsed -n -e :a -e '1,1!{P;N;D;};N;ba' ../outputData/HausdorffDistances.txt > ../outputData/hdgeneral.txt");
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
