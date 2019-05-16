
#include <math.h>
#include <iostream>
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkResampleImageFilter.h"
#include "itkEuler3DTransform.h"
#include "itkHausdorffDistanceImageFilter.h"
#include <string>

class Utilitarios{

private:
	const double dtr = (atan(1.0)*4.0)/180.0;
	const double rtd = 180.0/(atan(1.0)*4.0);
public:
	Utilitarios();
	void convertEulerToVersor(float &rx, float &ry, float &rz, double &ax, double &ay, double &az, double &angle);
    void convertVersorToEuler(double &x, double &y, double &z, double &rx, double &ry, double &rz);
	void unirVectorWithAngle(float &rx, float &ry, float&rz, double &vx, double &vy, double &vz, double &newangle);
	void compareVols(std::string logfilename, std::string newvolume, std::string imagendef, int numTest);
	void computeHausdorffDistance(std::string logfilename, typename itk::Image<short int, 3>::ConstPointer a, itk::Image<short int,3>::ConstPointer b);
    	void createStats(int numLevels, std::string logfilename, std::string dir_resultados, int indexTest);
	void createStatsOfTransValues(std::string dirres,std::string logfilename,int numTest);
	void getBaseElementsVersor(double &vx, double &vy, double &vz, double &rx, double &ry, double &rz);

	void createStatsOfErrors(int numTest);
	void createStatsBarHausdorff();
	void createStatsBoxPlotsTypeTransParams();


};
