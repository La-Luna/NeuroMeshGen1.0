#ifndef NEUROPOINT_H
#define NEUROPOINT_H
#include"vec3.h"
#include<iostream>
using namespace std;


class NeuroPoint{
public:
    NeuroPoint(){}
    NeuroPoint(int _n,int t, vec3 p ,double r,int pn):n(_n),type(t),pos(p),radius(r),parentn(pn){}
    NeuroPoint(const NeuroPoint& temp){
        n=temp.n;
        type=temp.type;
        pos=temp.pos;
        radius=temp.radius;
        parentn=temp.parentn;
    }



public:
    int n;
    int type;
    vec3 pos;
    double radius;
    int parentn;

};


#endif
