#ifndef GLOBAL_H
#define GLOBAL_H
#include"lunamatrix.h"
#include"neuropoint.h"
#include"neurogeometry.h"
#define PI 3.1415926
//class NeuroGeometry;
//extern NeuroGeometry* neurogeo;


double angle2radian(double angle);
LMatrix4 calRotateMatrix(const vec3& axis1,float cosA,float sinA);
std::ostream& operator<<(std::ostream& out,vec3 temp);
string vec3tostring(const vec3& temp);
string double2string(double temp);
string int2string(int temp);
double caldistance(const vec3&a,const vec3& b);
vec3 string2vec3(const string& s1,const string& s2,const string&s3);
double string2double(const string & str);
int string2int(const string& str );
std::ostream& operator<<( ostream& out,const NeuroPoint& temp);
#endif // GLOBAL_H
