#include "gestorMatrixRot.h"

HelperRot::HelperRot(){
	matrixRotX.SetIdentity();
	matrixRotY.SetIdentity();
	matrixRotZ.SetIdentity();
	matrixRotG.SetIdentity();
}
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

	matrixRotY(0,0) = cos(angle);
	matrixRotY(0,1) = 0.0;
	matrixRotY(0,2) = sin(angle);
	matrixRotY(1,0) = 0.0;
	matrixRotY(1,1) = 1.0;
	matrixRotY(1,2) = 0.0;
	matrixRotY(2,0) = -sin(angle);
	matrixRotY(2,1) = 0.0;
	matrixRotY(2,2) = cos(angle);
}

void HelperRot::initRotZ(float angle){

	matrixRotZ(0,0) = cos(angle);
	matrixRotZ(0,1) = -sin(angle);
	matrixRotZ(0,2) = 0.0;
	matrixRotZ(1,0) = sin(angle);
	matrixRotZ(1,1) = cos(angle);
	matrixRotZ(1,2) = 0.0;
	matrixRotZ(2,0) = 0.0;
	matrixRotZ(2,1) = 0.0;
	matrixRotZ(2,2) = 1.0;

}

void HelperRot::printRotg(){
	std::cout<<getRotg();
}

void HelperRot::printRotx(){
	std::cout<<getRotx();
}
void HelperRot::printRoty(){
	std::cout<<getRoty();
}
void HelperRot::printRotz(){
	std::cout<<getRotz();
}
void HelperRot::composeMatrixRot(){
	matrixRotG = getRotx() * getRoty() * getRotz();
}

