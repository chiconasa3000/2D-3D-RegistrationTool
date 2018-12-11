#include "gestorMatrixRot.h"

HelperRot::HelperRot(){
	matrixRotX.SetIdentity();
	matrixRotY.SetIdentity();
	matrixRotZ.SetIdentity();

void HelperRot::initRotX(float angle){

	matrixRotX(0,0) = 1.0;
	matrixRotX(0,1) = 0.0;
	matrixRotX(0,2) = 0.0;
	matrixRotX(1,0) = 0.0;
	matrixRotX(1,1) = cos(angle);
	matrixRotX(1,2) = -sin(angle);
	matrixRotX(2,0) = 0.0;
	matrixRotX(2,1) = sin(angle);
	matrixRotX(2,2) = cos(angle);

}

void HelperRot::initRotY(float angle){

	matrixRotX(0,0) = cos(angle);
	matrixRotX(0,1) = 0.0;
	matrixRotX(0,2) = sin(angle);
	matrixRotX(1,0) = 0.0;
	matrixRotX(1,1) = 1.0;
	matrixRotX(1,2) = 0.0;
	matrixRotX(2,0) = -sin(angle);
	matrixRotX(2,1) = 0.0;
	matrixRotX(2,2) = cos(angle);
}

void HelperRot::initRotZ(float angle){

	matrixRotX(0,0) = cos(angle);
	matrixRotX(0,1) = -sin(angle);
	matrixRotX(0,2) = 0.0;
	matrixRotX(1,0) = sin(angle);
	matrixRotX(1,1) = cos(angle);
	matrixRotX(1,2) = 0.0;
	matrixRotX(2,0) = 0.0;
	matrixRotX(2,1) = 0.0;
	matrixRotX(2,2) = 1.0;

}

void HelperRot::printRotx(){
	cout<<getRotx();
}
void HelperRot::printRoty(){
	cout<<getRoty();
}
void HelperRot::printRotz(){
	cout<<getRotz();
}
