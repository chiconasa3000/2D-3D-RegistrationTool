#include "genVirtGen.h"


int main(int argc,char *argv[]){
	char *imageToProject = NULL;
	char *logFileName = NULL;
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


	bool ok;
	while(argc > 1){
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-p") == 0))
		{
			argc--; argv++;
			ok = true;
			type_projection = argv[1];
			argc--; argv++;
		}


		if((ok==false) && strcmp(argv[1], "-imageToProject")==0)
		{
			argc--; argv++;
			ok = true;
			imageToProject = argv[1];
			argc--; argv++;
		}

		if((ok==false) && strcmp(argv[1], "-logFileName")==0)
		{
			argc--; argv++;
			ok = true;
			logFileName = argv[1];
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

	}
	
	GenVirtGen *genImagVirtual = new GenVirtGen();
	genImagVirtual->readMovingImage(imageToProject);
		
	//el reporte log debe ser al ultimo despues de haber iniciado todas las variables a imprimir
	//setlogfile deberia estar dentro de printself por seguridad	
	genImagVirtual->setLogFile(logFileName);
	genImagVirtual->printSelf();
	genImagVirtual->initResampleFilter();
	genImagVirtual->initTransform(tx,ty,tz, dcx,dcy,dcz,rx,ry,rz,sg);
	genImagVirtual->initInterpolator(focalPointx,focalPointy,focalPointz,threshold);
	return 1;
}


