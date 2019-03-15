

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

void Utilitarios::convertVersorToEuler(float &x, float &y, float &z, float &w, float &rx, float &ry, float &rz){
	float t0 = +2.0 * (w * x + y * z);
        float t1 = +1.0 - 2.0 * (x * x + y * y);
        rx = atan(t0/ t1);

        float t2 = +2.0 * (w * y - z * x);
 	t2 = (t2 > +1.0) ? +1.0 : t2;
	t2 = (t2 < -1.0) ? -1.0 : t2;
        ry = asin(t2);

        float t3 = +2.0 * (w * z + x * y);
        float t4 = +1.0 - 2.0 * (y * y + z * z);
        rz = atan(t3/ t4);
}
