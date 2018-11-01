
/* Generador de imagenes virtuales usando interpolador Patched */


#include <iostream>


// funcion de menu principal 

namespace{
void menu(){
	std::cout<<"Generador de imagenes virtuales : Parametros necesitados"<< std::endl;
}	



int main(int argc, char *argv[]){
	//function con las terceras parte
	
	const unsigned int Dimensions = 3;

	typedef itk::Image<short int, Dimensions> FixedImageType;
	typedef itk::Image<short int,Dimensions> MovingImageType;
	typedef itk::PatchedRayCastInterpolateImageFunction<MovingimageType, double> InterpolateType;
	IntepolatorType::Pointer interpolator = InterpolatorType::New();
	return 0;
}
}
