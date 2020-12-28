#include "lunacylender.h"
#include"lunamatrix.h"
#include"global.h"
 int LunaCylender::cylendersnum=0;
LunaCylender::LunaCylender()
{

}
LunaCylender::LunaCylender(const vec3& p1,const vec3& p2,float r1,float r2,int level){
    //构建圆柱体的表面vertex
    float dtheta=360.0/level;
    vec3 p1_p2=p1-p2;
    float segmentlength=p1_p2.length();
    //float segmentlength=1;
    float tx,ty,tz;
    ty=segmentlength/2.0;


    for(int i=0;i<=level;i++){
        float thetatemp=angle2radian(dtheta*i);
        tx=r1*std::cos( thetatemp);
        tz=r1*std::sin(thetatemp);
        vec3 newvertice1(tx,ty,tz);
        vec3 newvertice2(tx,-ty,tz);
        circle1.push_back(newvertice1);
        circle2.push_back(newvertice2);

    }
    addCylenderNum();
}

void LunaCylender::move(const LMatrix4& rotmat,const LMatrix4& transmat){
    //circle1
    for(int i=0;i<circle1.size();i++){
        vec4 tempvert=vec4(circle1[i][0],circle1[i][1],circle1[i][2],1.0);
        tempvert=rotmat*tempvert;
        tempvert=transmat*tempvert;
        //tempvert=rotmat*tempvert;

        circle1[i]=vec3(tempvert[0],tempvert[1],tempvert[2]);
    }

    //circle2
    for(int i=0;i<circle2.size();i++){
        vec4 tempvert=vec4(circle2[i][0],circle2[i][1],circle2[i][2],1.0);
        //tempvert=transmat*tempvert;
        tempvert=rotmat*tempvert;
        tempvert=transmat*tempvert;
        circle2[i]=vec3(tempvert[0],tempvert[1],tempvert[2]);
    }
}
 void LunaCylender::addCylenderNum(){
    cylendersnum++;
}
