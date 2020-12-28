#include"vec3.h"

double vec3::operator[](int i)const{
    return this->data[i];
}
double& vec3::operator[](int i){
    return this->data[i];
}

vec3 vec3::operator-(const vec3& temp)const{
       return vec3(data[0]-temp[0],data[1]-temp[1],data[2]-temp[2]);
}
vec3 vec3::operator/(float k)const{
    return vec3(data[0]/k,data[1]/k,data[2]/k);
}
bool vec3::operator==(const vec3&temp){
    for(int i=0;i<3;i++){
        if(data[i]!=temp[i])return false;
    }
    return true;
}

vec3 vec3::normal()const{
    double len=this->length();
    return vec3(data[0]/len,data[1]/len,data[2]/len);
}
double vec3::length()const{
    return std::sqrt(this->length_squared());
}
double vec3::length_squared()const{
    return data[0]*data[0]+data[1]*data[1]+data[2]*data[2];
}
vec3 vec3::cross(const vec3&v)const{
//    return vec3(
//                data[1]*temp[2]-temp[1]*data[2],
//                -(data[0]*temp[2]-temp[0]*data[2]),
//                data[0]*temp[1]-temp[0]*data[1]
//                );
    return vec3(data[1] * v[2] - data[2] * v[1],
                data[2] * v[0] - data[0] * v[2],
                data[0] * v[1] - data[1] * v[0]);
}
float vec3::dot(const vec3& temp)const{
    return data[0]*temp[0]+data[1]*temp[1]+data[2]*temp[2];
}

