/* Generador de imagenes virtuales usando interpolador Patched */

#include <iostream>
#include <itkImage.h>
#include "../itkPatchedRayCastInterpolateImageFunction.h"
// funcion de menu principal 


void menu(){
	std::cout<<"Generador de imagenes virtuales : Parametros necesitados"<< std::endl;
}

int main(int argc, char *argv[]){
	//function con las terceras parte
	menu();	
	const unsigned int Dimensions = 3;

	typedef itk::Image<short int,Dimensions> FixedImageType;
	typedef itk::Image<short int,Dimensions> MovingImageType;
	typedef itk::PatchedRayCastInterpolateImageFunction<MovingImageType, double> InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	return 0;
}

