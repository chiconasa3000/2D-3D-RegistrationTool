
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

public:
	Utilitarios();
	void convertEulerToVersor(float &rx, float &ry, float &rz, double &ax, double &ay, double &az, double &angle);
	void convertVersorToEuler(float &x, float &y, float &z, float &w, float &rx, float &ry, float &rz);
	void unirVectorWithAngle(float &rx, float &ry, float&rz, double &vx, double &vy, double &vz, double &newangle);
	void compareVols(std::string logfilename, std::string newvolume, std::string imagendef, int numTest);
	void computeHausdorffDistance(std::string logfilename, typename itk::Image<short int, 3>::ConstPointer a, itk::Image<short int,3>::ConstPointer b);
    void createStats(int numLevels, std::string logfilename, std::string dir_resultados, int indexTest);
};
