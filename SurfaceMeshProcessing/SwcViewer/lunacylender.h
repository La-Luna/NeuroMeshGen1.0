#ifndef LUNACYLENDER_H
#define LUNACYLENDER_H
#include"object.h"
#include"vec3.h"
#include<vector>
#include"lunamatrix.h"
class LunaCylender:public Object
{
public:
    LunaCylender();
    LunaCylender(const vec3&p1,const vec3&p2,float r1,float r2,int level);

    void move(const LMatrix4& rotmat,const LMatrix4& transmat);
    static void addCylenderNum();

public:
    std::vector<vec3> circle1;
    std::vector<vec3 > circle2;
    std::vector<int> renderindex;
    static int cylendersnum;

};

#endif // LUNACYLENDER_H
