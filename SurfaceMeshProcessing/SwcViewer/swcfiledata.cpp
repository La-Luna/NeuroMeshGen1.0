#include"swcfiledata.h"
//仅处理只含有一个soma点的神经元,并且每个points的格式为
//n type x y z radius parents
SwcFileData::SwcFileData(const string& filepath,const string& filename){
    //cout<<"run in SwcFileData(  )"<<endl;
    //auto slash = filepath.find_last_not_of('/');
    auto backslash = filepath.find_last_of('\\');
    string component="/";
    if (backslash != std::string::npos)component = "\\";

    string filefullpath=filepath+component+filename;
    //string filefullpath=filename;
    cout<<"filefullpath:"<<filefullpath<<endl;

    ifstream ifs;
    ifs.open(filefullpath.c_str(),ios::in);
    if(!ifs.is_open()){
        cout<<"read file:"<<filename<<" error!"<<endl;
        return;
    }


    //读出soma到 soma
    //去掉一些##
    //soma确认
    //neuropoints确认
    string buf;

    while(getline(ifs,buf)){

        if(buf[0]!='#'){
                countreadpoints++;
                shared_ptr<NeuroPoint> points_temp;
                points_temp=swcString2Struct(buf);
                //cout<<points_temp<<endl;

                if(points_temp->parentn==-1){
                    soma=points_temp;
                }
                else{
                    neuropointslist.push_back(points_temp);
                }


        };
    }


}


    vec3 SwcFileData:: getSomaPosition(){
        if(soma==NULL){
            return vec3(-1.0,-1.0,-1.0);
        }
        return soma->pos;
    }
     shared_ptr<NeuroPoint> SwcFileData::swcString2Struct(string str){
        // cout<<"run in swcString2Struct(  )"<<endl;
        shared_ptr<NeuroPoint> ans;
         //把swcstring拆到vector<string>中去
        vector<string> dividedswcstring;
        //去掉每行最前面的空格
        int i=0;
        for(;i<str.size();i++){
            if(str[i]!=' ')break;
        }
        //按6个空格去分割字符串
        string temp;
        for(i;i<str.size();i++){
            if(str[i]==' '){
                dividedswcstring.push_back(temp);
                temp.clear();
            }
            else{
                temp+=str[i];
            }
        }
        dividedswcstring.push_back(temp);
        if(dividedswcstring.size()>7){
            cout<<"the swc file is not standard,because the datas of neuropoints is larger than 7!"<<endl;
            return ans;
        }

        //test dividedswcstring
        // for(int i=0;i<dividedswcstring.size();i++){
        //     cout<<dividedswcstring[i]<<" ";
        // }
        // cout<<endl;


         int _n;
         int _type;
         vec3 _pos;
         double _radius;
         int _parentn;
         int  count=7;
         //将vector<string>拆到对应的变量中
         _n=string2int(dividedswcstring[0]);
         _type=string2int(dividedswcstring[1]);
         _pos=string2vec3(dividedswcstring[2],dividedswcstring[3],dividedswcstring[4]);
         _radius=string2double(dividedswcstring[5]);
         _parentn=string2int(dividedswcstring[6]);

         ans=make_shared<NeuroPoint>(_n,_type,_pos,_radius,_parentn);
        return ans;

    }
    void  SwcFileData::printNeuropointsList(){
        cout<<"---------------------the content of soma---------------------"<<endl;
        cout<<*soma<<endl;

        cout<<"The size of neuropointslist is "<<neuropointslist.size()<<endl;
        cout<<"---------------------the content of list---------------------"<<endl;
        for(int i=0;i<neuropointslist.size();i++){
            cout<<*(neuropointslist[i])<<endl;
        }
    }
    void SwcFileData::createnewfile(vector<shared_ptr<NeuroPoint>>& temp){
        return;
    }
