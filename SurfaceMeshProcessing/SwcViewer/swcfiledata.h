#ifndef SWCFILEDATA_H
#define SWCFILEDATA_H
#include<iostream>
#include"neuropoint.h"
#include"global.h"
#include<vector>
#include<string>
#include<memory>
#include<fstream>
#include<QString>
using namespace std;
class SwcFileData{
public:
    SwcFileData(){cout<<"default SwcFIleData constructor!"<<endl;}
    SwcFileData(const string& filepath,const string& filename);


public:
    vec3 getSomaPosition();
    shared_ptr<NeuroPoint> swcString2Struct(string str);
    void printNeuropointsList();
    void createnewfile(vector<shared_ptr<NeuroPoint>>& temp);
public:
    shared_ptr<NeuroPoint> soma;
    vector<shared_ptr<NeuroPoint>> neuropointslist;
    int countreadpoints=0;
};


#endif
