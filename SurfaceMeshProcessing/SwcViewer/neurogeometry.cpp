#include"neurogeometry.h"
#include"global.h"
#include"lunacylender.h"
NeuroGeometry::NeuroGeometry(SwcFileData* originaldata){
        neuropointslist=originaldata->neuropointslist;
        soma=originaldata->soma;
        //calPointsInSomaSphere();
}
void NeuroGeometry::calPointsInSomaSphere(){
    pointinsomasphere.clear();
    for(int i=0;i<neuropointslist.size();i++){
    //计算每个点到soma pos的距离，小于soma radiusneurostruct2string();NeuroPoint>(*(neuropointslist[i]));
        double distance=caldistance(neuropointslist[i]->pos,soma->pos);
        if(distance<soma->radius){
            shared_ptr<NeuroPoint> piss=make_shared<NeuroPoint>(*(neuropointslist[i]));
            pointinsomasphere.push_back(piss);
        }
    }

    printpointlist("pointsinsomasphere",pointinsomasphere);

}


void NeuroGeometry::changePointsRadius(){
    double _radius=2.0;
    if(!newneuropointslist.empty()){
            for(int i=0;i<newneuropointslist.size();i++){
                newneuropointslist[i]->radius=_radius;
            }
    }
    else{
    //把半径放小
    int _n;
    int _type;
    vec3 _pos;

    int _parentn;
    NeuroPoint temp;
     //soma
    temp=*soma;
    _n=temp.n;_type=temp.type;_pos=temp.pos;_parentn=temp.parentn;
    shared_ptr<NeuroPoint> newpoint=make_shared<NeuroPoint>(_n,_type,_pos,_radius,_parentn);
    newneuropointslist.push_back(newpoint);
    //neuropointslist
    for(int i=0;i<neuropointslist.size();i++){
        temp=*(neuropointslist[i]);
        _n=temp.n;_type=temp.type;_pos=temp.pos;_parentn=temp.parentn;
        //新points放到vector中去
        shared_ptr<NeuroPoint> newpoint=make_shared<NeuroPoint>(_n,_type,_pos,_radius,_parentn);
        newneuropointslist.push_back(newpoint);
    }
    }
}

void  NeuroGeometry::changePointsType(){
    if(newneuropointslist.empty()){
    //把type改为3
    int _n;
    int _type=3;
    vec3 _pos;
    double _radius;
    int _parentn;
    NeuroPoint temp;
     //soma
    temp=*soma;
    _n=temp.n;_radius=temp.radius;_pos=temp.pos;_parentn=temp.parentn;
    shared_ptr<NeuroPoint> newpoint=make_shared<NeuroPoint>(_n,_type,_pos,_radius,_parentn);
    newneuropointslist.push_back(newpoint);
    //neuropointslist
    for(int i=0;i<neuropointslist.size();i++){
        temp=*(neuropointslist[i]);
        _n=temp.n;_radius=temp.radius;_pos=temp.pos;_parentn=temp.parentn;
        //新points放到vector中去
        shared_ptr<NeuroPoint> newpoint=make_shared<NeuroPoint>(_n,_type,_pos,_radius,_parentn);
        newneuropointslist.push_back(newpoint);
    }


    }
    else{
            for(int i=0;i<newneuropointslist.size();i++){
                newneuropointslist[i]->type=3;
            }
    }
}


void NeuroGeometry::findLargeRadiusPoints(){
    cout<<"---------------------the points which radius is large---------------------"<<endl;
    for(int i=0;i<neuropointslist.size();i++){
        if(neuropointslist[i]->radius>5.0){
            cout<<*(neuropointslist[i])<<endl;
        }
    }
}
void NeuroGeometry::printpointlist(string vectorname,vector<shared_ptr<NeuroPoint>>temp){
    cout<<"---------------------print "<<vectorname<<" ---------------------"<<endl;
    cout<<"the list size is"<<temp.size()<<endl;
    for(int i=0;i<temp.size();i++){
        cout<<*(temp[i])<<endl;
    }
    }
void NeuroGeometry::createNewSwcFIle(string filename,string  vectorname,vector<shared_ptr<NeuroPoint>>temp){
    cout<<"---------------------"<<vectorname<<" will be generated into new file called"<<filename<<"---------------------"<<endl;
    ofstream ofs;
    ofs.open(filename.c_str(),ios::out);


    //neuropoints写入，list已有soma
    for(int i=0;i<temp.size();i++){
        string tempstring=neurostruct2string(temp[i]);
        ofs<<tempstring<<endl;;
    }


    ofs.close();
}

string  NeuroGeometry::neurostruct2string(shared_ptr<NeuroPoint> temp){
    string fullstr;
    string _n=int2string(temp->n);
    string _type=int2string(temp->type);
    string _pos=vec3tostring(temp->pos);
    string _radius=double2string(temp->radius);
    string _parentn=int2string(temp->parentn);
    fullstr=fullstr+_n+" "+_type+" "+_pos+" "+_radius+" "+_parentn;
    return fullstr;
}
vector<float> NeuroGeometry::getlinedata(int & verticesnum){
    //modify position according to the position of soma
    vec3 oldsomapos=this->getSomaPosition();
    int linenum=this->neuropointslist.size();
    std::cout<<"the lines' composed points num : "<<linenum<<std::endl;
    int verticesdatanum=linenum*2*6;
    verticesnum=verticesdatanum;
    //float* vertices=new float[verticesdatanum];
    vector<float> vertices(verticesdatanum,0.0);
    int verticeid=0;
    for(int i=0;i<linenum;i++){
        //calcuate pos
        vec3 curpointpos=this->neuropointslist[i]->pos;
        int lastposid=this->neuropointslist[i]->parentn;
        vec3 lastpointpos;
        if(lastposid==1) lastpointpos=oldsomapos;
        else{
             lastpointpos=neuropointslist[lastposid-2]->pos;
        }
        //for(int j=0;j<3;j++){
            curpointpos=vec3(
                        curpointpos[0]-oldsomapos[0],
                        curpointpos[1]-oldsomapos[1],
                        curpointpos[2]-oldsomapos[2]
                        );
            lastpointpos=vec3(
                        lastpointpos[0]-oldsomapos[0],
                        lastpointpos[1]-oldsomapos[1],
                        lastpointpos[2]-oldsomapos[2]
                        );
        //}
        //add pos and color
        vertices[verticeid++]=lastpointpos[0];vertices[verticeid++]=lastpointpos[1];vertices[verticeid++]=lastpointpos[2];
        vertices[verticeid++]=1.0;vertices[verticeid++]=0.0;vertices[verticeid++]=0.0;
        vertices[verticeid++]=curpointpos[0];vertices[verticeid++]=curpointpos[1];vertices[verticeid++]=curpointpos[2];
        vertices[verticeid++]=1.0;vertices[verticeid++]=0.0;vertices[verticeid++]=0.0;


    }
    //print
    //printverticesdata(vertices,verticesdatanum);
    return vertices;

}


    void NeuroGeometry::printverticesdata(float* verticesdata,int size){
        for(int i=0;i<size;i++){
             if(i%6==0)std::cout<<std::endl;
            std::cout<<verticesdata[i]<<" ";

        }
    }

    void NeuroGeometry::calCylinder(){

        //modify position according to the position of soma
        vec3 oldsomapos=this->getSomaPosition();
        int linenum=this->neuropointslist.size();
        int verticeid=0;
        for(int i=0;i<linenum;i++){
            //if(i==0){
            //calcuate pos
            vec3 curpointpos=this->neuropointslist[i]->pos;
            int lastposid=this->neuropointslist[i]->parentn;
            vec3 lastpointpos;
            if(lastposid==1) lastpointpos=oldsomapos;
            else{
                 lastpointpos=neuropointslist[lastposid-2]->pos;
            }
                curpointpos=vec3(
                            curpointpos[0]-oldsomapos[0],
                            curpointpos[1]-oldsomapos[1],
                            curpointpos[2]-oldsomapos[2]
                            );
                lastpointpos=vec3(
                            lastpointpos[0]-oldsomapos[0],
                            lastpointpos[1]-oldsomapos[1],
                            lastpointpos[2]-oldsomapos[2]
                            );

            //for 对每条segement进行画cylinder处理
            //结果合并到this->neuritecylindervertices中去
            calCylinderTempdata(this->neuritecylindervertices,this->neuritefaceindexarray,lastpointpos,curpointpos,0.5,0.5,6);//半径没读！！
            // }//if(i==0)
            //else break;

        }



        //给size_value赋值
        std::cout<<"the cylenders num : "<<LunaCylender::cylendersnum<<std::endl;
        this->vertexsize=this->neuritecylindervertices.size();
        this->indexarraysize=this->neuritefaceindexarray.size();

    }
    vector<float> NeuroGeometry::getNeuroNodeVerticeArray(int& verticesize){
            //得到节点坐标数据
        vec3 oldsomapos=this->getSomaPosition();
        int neurodenums=this->neuropointslist.size();//neuropointslist
        for(int i=0;i<neurodenums;i++){
            //每一个神经元节点
            vec3 curneurode=this->neuropointslist[i]->pos;
            curneurode=curneurode-oldsomapos;
            //都生成一个圆
            this->calSphereTempdata(this->neurodevertices,this->neurodefaceindexarray,curneurode,0.5);


        }
        std::cout<<"the size of unsigned int is : "<<sizeof(unsigned int)<<std::endl;
            //得到vertice数组，index数组
        verticesize=this->neurodevertices.size();
        std::cout<<"the neurodeverticesarrau size is : "<<verticesize<<std::endl;
        return this->neurodevertices;
    }




    void NeuroGeometry::calCylinderTempdata(vector<float>&vertices,vector<unsigned int>& indeces,const vec3& p1,const vec3& p2,double r1,double r2, int level){//默认r1=r2
        //方案1：重新撸个neurite结构

        LunaCylender littlecylender(p1,p2,r1,r2,level);

        //生成index序列

        unsigned int base=vertices.size()/6;
        unsigned int ur,ul,br,bl;
        for(unsigned int i=0;i<level;i++){
            ur=base+i;
            ul=ur+1;
            br=ur+level+1;
            bl=ul+level+1;
            indeces.push_back(ur);indeces.push_back(ul);indeces.push_back(br);
            indeces.push_back(ul);indeces.push_back(br);indeces.push_back(bl);
        }




        //求旋转矩阵，位移矩阵


        vec3 middlepoint(
                    (p1[0]+p2[0])/2.0,
                    (p1[1]+p2[1])/2.0,
                    (p1[2]+p2[2])/2.0
                );

        LMatrix4 mytransform(1.0, 0.0, 0.0, 0.0,//col 1
                   0.0, 1.0, 0.0, 0.0,//col2
                    0.0, 0.0, 1.0, 0.0,//col3
                    middlepoint[0],  middlepoint[1],  middlepoint[2],1.0//col4
        );
        vec3 p1_p2=p1-p2;
        vec3 up=vec3(0.0,1.0,0.0);
        vec3 rotateaxis=up.cross(p1_p2);

        float costheta=p1_p2.dot(up)/(p1_p2.length());

        //float _theta=std::acos(costheta)*(-1);
        //float _costheta=std::cos(_theta);
        //float _sintheta=std::sin(_theta);
        float sintheta=std::sqrt(1-costheta*costheta);
        LMatrix4 myrotate;
        if(rotateaxis==vec3(0.0,0.0,0.0)){}
        else {myrotate=calRotateMatrix(rotateaxis.normal(),costheta,sintheta);}
        littlecylender.move(myrotate,mytransform);

        //cylinder pos1++color1+pos2+color2放进vector中，
        //合到vector中去
        vector<vec3> c1temp=littlecylender.circle1;
        for(int i=0;i<littlecylender.circle1.size();i++){
            vector<float> temp={(float)c1temp[i][0],(float)c1temp[i][1],(float)c1temp[i][2],0.0,1.0,0.0};
            vertices.insert(vertices.end(),temp.begin(),temp.end());
        }
        c1temp=littlecylender.circle2;
        for(int i=0;i<littlecylender.circle1.size();i++){
            vector<float> temp={ (float)c1temp[i][0],(float)c1temp[i][1],(float)c1temp[i][2],0.0,1.0,0.0};
            vertices.insert(vertices.end(),temp.begin(),temp.end());
        }


        //方案2：每个节点处画个球
    }

    void NeuroGeometry::calSphereTempdata(vector<float>&vertices,vector<unsigned int>& indeces,const vec3& p1,double radius,int level){
        vector<float> curspherevertices;
        vector<unsigned int> curindexarray;
        radius=(float)radius;
        int array_size = (level + 1)*(level + 1) * 6;
        //curspherevertices.resize(array_size,0);
        //double dtheta, dphi;
        //dtheta = (2 * PI) / level;
        //dphi = PI / level;
        float dtheta = 360.0 / level;
        float dphi = 180.0 / level;

        float cx, cy, cz;

        for (int i = 0; i <= level; i++){//phi
            for (int j = 0; j <= level; j++){//theta
                float phi_temp =angle2radian(i*dphi) ;
                float theta_temp =angle2radian(j*dtheta) ;
                cx = radius*std::sin(phi_temp)*std::cos(theta_temp);
                cy = radius*std::cos(phi_temp);
                cz = radius*std::sin(phi_temp)*std::sin(theta_temp);
                int curpos = (i*(level+1) + j) * 6;
                curspherevertices.push_back((float)p1[0]+cx);
                curspherevertices.push_back((float)p1[1]+cy);
                curspherevertices.push_back((float)p1[2]+cz);
                curspherevertices.push_back(0.0);
                curspherevertices.push_back(1.0);
                curspherevertices.push_back(0.2);


//                curspherevertices[curpos + 0] = (float)p1[0]+cx;
//                curspherevertices[curpos + 1] = (float)p1[1]+cy;
//                curspherevertices[curpos + 2] = (float)p1[2]+cz;
//                curspherevertices[curpos + 3] = 0.0;
//                curspherevertices[curpos + 4] = 0.0;
//                curspherevertices[curpos + 5] = 1.0;

            }
        }

        //index
        //data
        int indexarray_size = level*level * 2 * 3;
        //curindexarray.resize(indexarray_size);
        unsigned int upl, upr, botl, botr;

        //for ,upl,upr,botl,botr
        unsigned int base=this->neurodevertices.size()/6;
        for (unsigned int i = 0; i < level; i++){
            for (unsigned int j = 0; j < level; j++){
                upr =base+i*(level+1)+j; upl = upr + 1;
                botr =upr + level+1; botl = botr + 1;
                int curpos = (i*level + j) * 6;
                curindexarray.push_back(upr);
                curindexarray.push_back(upl);
                curindexarray.push_back(botr);

                curindexarray.push_back(botr);
                curindexarray.push_back(botl);
                curindexarray.push_back(upl);

//                curindexarray[curpos + 0] = upr;
//                curindexarray[curpos + 1] = upl;
//                curindexarray[curpos + 2] = botr;

//                curindexarray[curpos + 3] = botr;
//                curindexarray[curpos + 4] = botl;
//                curindexarray[curpos + 5] = upl;
            }
        }
        vertices.insert(vertices.end(),curspherevertices.begin(),curspherevertices.end());
        indeces.insert(indeces.end(),curindexarray.begin(),curindexarray.end());
        //count++;
        //std::cout<<count<<" :temp inseces size: "<<indeces.size()<<std::endl;
    }

    vector<unsigned int> NeuroGeometry::getNeuroNodeFaceIndexArray(int& indexarraysize){
        indexarraysize=this->neurodefaceindexarray.size();
       // std::cout<<"the neurodeface index array size is : "<<indexarraysize<<std::endl;
        return this->neurodefaceindexarray;
    }
