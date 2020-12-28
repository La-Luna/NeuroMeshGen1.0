#ifndef MY_VEC3_H
#define MY_VEC3_H
#include<iostream>
#include"math.h"
#include<cmath>


class vec3{
public:
    vec3(){}
    vec3(double a,double b,double c){
        data[0]=a;
        data[1]=b;
        data[2]=c;
    }
    vec3(const vec3& temp){
        for(int i=0;i<3;i++){
            data[i]=temp[i];
        };
    }

public:
    double operator[](int i)const;
    double& operator[](int i);
    vec3 operator-(const vec3& temp)const;
    vec3 operator/(float k)const;
    vec3 normal()const;
    double length()const;
    double length_squared()const;
    vec3 cross(const vec3&temp)const;
    bool operator==(const vec3&temp);
    float dot(const vec3&temp)const;
 public:   
    double  data[3];
};



class vec4{
public:
    vec4(){}
    vec4(double a,double b,double c,double d){
        data[0]=a;
        data[1]=b;
        data[2]=c;
        data[3]=d;
    }
    vec4(const vec4& temp){
        for(int i=0;i<4;i++){
            data[i]=temp[i];
        };
    }

public:
    double operator[](int i)const{return data[i];}
    double& operator [](int i){return data[i];}


public:
    double data[4];
};

#endif
