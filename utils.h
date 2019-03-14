
#include <math.h>


class Utilitarios{

private:
	const double dtr = (atan(1.0)*4.0)/180.0;

public:
	Utilitarios();
	void convertEulerToVersor(float &rx, float &ry, float &rz, double &ax, double &ay, double &az, double &angle);
	void convertVersorToEuler(float &x, float &y, float &z, float &w, float &rx, float &ry, float &rz);


};
