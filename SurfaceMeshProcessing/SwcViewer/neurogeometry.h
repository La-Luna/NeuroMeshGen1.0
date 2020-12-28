#ifndef NEUROGEOMETRY_H
#define NEUROGEOMETRY_H

#include"swcfiledata.h"
#include"neuropoint.h"
#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include<map>
#include"vec3.h"
#include<memory>
using namespace std;

class SwcFileData;
class NeuroGeometry{
public:
    NeuroGeometry(){}
    NeuroGeometry(SwcFileData* originaldata);




public:
    void calPointsInSomaSphere();
    void changePointsRadius();
    void changePointsType();
    void findLargeRadiusPoints();
    void printpointlist(string  vectorname,vector<std::shared_ptr<NeuroPoint>>temp);
    void createNewSwcFIle(string filename,string  vectorname,vector<shared_ptr<NeuroPoint>>temp);
    string neurostruct2string(shared_ptr<NeuroPoint> temp);
    vec3 getSomaPosition(){return soma->pos;}

    void calCylinder();
    void calCylinderTempdata(vector<float>&vertices,vector<unsigned int>& indeces,const vec3& p1,const vec3& p2,double r1,double r2, int level=6);
    void calSphereTempdata(vector<float>&vertices,vector<unsigned int>& indeces,const vec3& p1,double r,int level=12);
    void printverticesdata(float* vertices,int size);
    vector<unsigned int> getFaceIndexArray(int& _indexarraysize){_indexarraysize=this->indexarraysize;return neuritefaceindexarray;}
    vector<unsigned int> getNeuroNodeFaceIndexArray(int& indexarraysize);
    vector<float> getCylinderVerticeArray(int& _vertexsize){this->calCylinder();_vertexsize=this->vertexsize;return neuritecylindervertices;}
    vector<float> getNeuroNodeVerticeArray(int& verticesize);
    vector<float> getlinedata(int & verticesnum);
    virtual void getPointDataonly() = 0;
    virtual void getLinePointsDataonly() = 0;


 public:
    vector<shared_ptr<NeuroPoint>> pointinsomasphere;
    vector<shared_ptr<NeuroPoint>> neuropointslist;
    vector<shared_ptr<NeuroPoint>>newneuropointslist;//包含soma
    shared_ptr<NeuroPoint> soma;
    vector<unsigned int> neuritefaceindexarray;
    vector<unsigned int> neurodefaceindexarray;
    vector<float> neuritecylindervertices;
    vector<float> neurodevertices;



    int indexarraysize;
    int vertexsize;
    int  count=0;


};


#endif
