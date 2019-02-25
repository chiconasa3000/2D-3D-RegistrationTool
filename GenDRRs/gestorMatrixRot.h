#include <iostream>
#include <math.h>
#include "itkMatrix.h"

//funciones basicas de control de rotaciones

class HelperRot{

private:
	typedef itk::Matrix<double,3,3> MatrixType;
	MatrixType matrixRotG;	
	MatrixType matrixRotX;
	MatrixType matrixRotY;
	MatrixType matrixRotZ;

public:
	HelperRot();
	void composeMatrixRot();
	void initRotX(float angleRadians);
	void initRotY(float angleRadians);
	void initRotZ(float angleRadians);
	void printRotx();
	void printRoty();
	void printRotz();
	void printRotg();
	MatrixType getRotx(){return matrixRotX;};
	MatrixType getRoty(){return matrixRotY;};
	MatrixType getRotz(){return matrixRotZ;};
	MatrixType getRotg(){return matrixRotG;};

};

