#include"global.h"
#include<iostream>
#include<string>
#include<vector>
#include"vec3.h"
#include<math.h>
#include<algorithm>
#include"neuropoint.h"
#include"lunamatrix.h"

using namespace std;

double angle2radian(double angle){
    return (angle/180)*PI;
}

LMatrix4 calRotateMatrix(const vec3& axis1,float cosA,float sinA){
    float rx = axis1[0];//0
    float ry = axis1[1];//0
    float rz = axis1[2];//1
    //cosA=0.86;sinA=0.5
    float _cosA = 1 - cosA;
    LMatrix4 ans(
        cosA + rx*rx*_cosA, rx*ry*_cosA + rz*sinA, rz*rx*_cosA - ry*sinA, 0,
        rx*ry*_cosA - rz*sinA, cosA + ry*ry*_cosA, rz*ry*_cosA + rx*sinA, 0,
        rx*rz*_cosA + ry*sinA, ry*rz*_cosA - rx*sinA, cosA + rz*rz*_cosA, 0,
        0, 0, 0, 1
        );
    return ans;


}

//-1  或者 1
int string2int(const string& str ){
    int ans=0;
    bool binus=false;
    //把符号提取出来
    int i=0;
    if(str[0]=='-'){binus=true;i++;}
    //int len=str.size();

    for(;i<str.size()&&str[i]!='\r';i++){
        ans=ans*10+str[i]-'0';
    }
    if(binus)return ans*(-1);
    return ans;

}
//-8.3 或者8.3 或者 83 或者0.83
double string2double(const string & str){
//提取负号
bool binus=false;
int i=0;
if(str[i]=='-'){binus=true;i++;};

double fisrtpart=0.0,secondpart=0.0,wholepart=0.0;
//以小数点为分界
bool meetdot=false;
double digits=0.1;
for(;i<str.size();i++){
    if(str[i]=='.'){meetdot=true;continue;}
    if(meetdot){
        secondpart=secondpart+(str[i]-'0')*digits;
        digits*=0.1;
    }
    else{
        fisrtpart=fisrtpart*10.0+str[i]-'0';
    }
}
wholepart=fisrtpart+secondpart;
if(binus)return (-1.0)*wholepart;
return wholepart;
}



vec3 string2vec3(const string& s1,const string& s2,const string&s3){
    vec3 ans(
        string2double(s1),
        string2double(s2),
        string2double(s3)
        );
    return ans;
}

double caldistance(const vec3&a,const vec3& b){
    double temp=(a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]);
    double ans=sqrt(temp);
    return ans;

}

//-1  或者 1  405
string int2string(int temp){
if(temp==0) return "0";
string ans;
if(temp<0){
    ans+='-';
    temp*=(-1);
}
int t=temp;
vector<int> vec;
while(t){
    int  a=t%10;
    vec.push_back(a);
    t/=10;
}
reverse(vec.begin(),vec.end());
for(int i=0;i<vec.size();i++){
    ans+=(vec[i]+'0');
}
return ans;
}

//-8.3 或者8.3 或者 83 或者0.83
//把小数限制为3位！
string double2string(double temp){
    string part1;
    string part2;
    string ans;
    if(temp<0){
        ans+='-';
        temp=temp*(-1);
    }
    int t1=(int)temp;
    double  t2=temp-t1;//0.001  0.01
    part1=int2string(t1);
    ans+=part1;
    ans+='.';
    int mult=10;
    int digit;
    for(int i=0;i<3;i++){
        digit=(int)(t2*mult);
        part2+=digit+'0';
        t2=(t2*mult)-digit;
    }
    ans+=part2;
    return ans;
}

string vec3tostring(const vec3& temp){
    string wholestr;
    string part1=double2string(temp[0]);
    string part2=double2string(temp[1]);
    string part3=double2string(temp[2]);
    wholestr=wholestr+part1+" "+part2+" "+part3;
    return wholestr;
}
std::ostream& operator<<(std::ostream& out,vec3 temp){
    return out<<temp[0]<<" "<<temp[1]<<" "<<temp[2];
}
std::ostream& operator<<( ostream& out,const NeuroPoint& temp){

    return out<<temp.n<<" "<<temp.type<<" "<<temp.pos<<" "<<temp.radius<<" "<<temp.parentn;

}
