#include "neurogeometry2mesh.h"

NeuroGeometry2Mesh::NeuroGeometry2Mesh(SwcFileData* originaldata) :NeuroGeometry(originaldata) {
	this->getPointDataonly();
	this->getLinePointsDataonly();
	this->getSpherePointDataonly();
	this->getCylinderPointDataonly();
}
void NeuroGeometry2Mesh::getPointDataonly() {

	this->posMin = vec3(DBL_MAX, DBL_MAX, DBL_MAX);
	this->posMax = vec3(DBL_MIN, DBL_MIN, DBL_MIN);
	vec3 _somapos = this->getSomaPosition();
	for (auto curpoint : this->neuropointslist) {
		vec3 _point = curpoint->pos - _somapos;
		this->pointsonly.push_back(_point);
		posMin[0] = min(posMin[0], _point[0]);
		posMin[1] = min(posMin[1], _point[1]);
		posMin[2] = min(posMin[2], _point[2]);

		posMax[0] = max(posMax[0], _point[0]);
		posMax[1] = max(posMax[1], _point[1]);
		posMax[2] = max(posMax[2], _point[2]);
	}

}
void NeuroGeometry2Mesh::getLinePointsDataonly() {
	vec3 _somapos = this->getSomaPosition();
	int linenum = this->neuropointslist.size();

	for (int i = 0; i < linenum; i++) {
		vec3 curp = this->neuropointslist[i]->pos;
		vec3 lastp;
		int lastpid = this->neuropointslist[i]->parentn;
		if (lastpid == 1)lastp = vec3(0.0, 0.0, 0.0);
		else
		{
			lastp = this->neuropointslist[lastpid - 2]->pos;
			lastp = lastp - _somapos;
		}
		curp = curp - _somapos;
		this->linepointsonly.push_back(lastp);
		this->linepointsonly.push_back(curp);

	}


 }
void NeuroGeometry2Mesh::getSpherePointDataonly() {
	//soma sphere
	vec3 _somapos(0.0,0.0,0.0);
	
	this->calSphereTempdata(_somapos, soma->radius, 6);


	//every points 2 a sphere
	for (int i = 0; i < this->neuropointslist.size(); i++) {
		this->calSphereTempdata(this->pointsonly[i], this->neuropointslist[i]->radius, 6);
	}
}
void NeuroGeometry2Mesh::getCylinderPointDataonly() {
	//every two points 2 a cylinder
	vec3 _somapos = this->getSomaPosition();
	int linenum = this->neuropointslist.size();

	for (int i = 0; i < linenum; i++) {
		vec3 curp = this->neuropointslist[i]->pos;
		vec3 lastp;
		int lastpid = this->neuropointslist[i]->parentn;
		if (lastpid == 1) {
			lastp = vec3(0.0, 0.0, 0.0);
			curp = curp - _somapos;
			calCylinderTempdata(lastp, curp,soma->radius, neuropointslist[i]->radius,  6);
		}
		else
		{
			lastp = this->neuropointslist[lastpid - 2]->pos;
			lastp = lastp - _somapos;
			curp = curp - _somapos;
			calCylinderTempdata(lastp, curp, neuropointslist[lastpid - 2]->radius, neuropointslist[i]->radius, 6);
		}



	}


}
void NeuroGeometry2Mesh::calSphereTempdata(const vec3& p1, double radius, int level) {
	vector<vec3> curspherepointsvec;
	vector<unsigned int> curindexarray;
	radius = (float)radius;
	int array_size = (level + 1) * (level + 1) * 6;
	float dtheta = 360.0 / level;
	float dphi = 180.0 / level;

	float cx, cy, cz;

	for (int i = 0; i <= level; i++) {//phi
		for (int j = 0; j <= level; j++) {//theta
			float phi_temp = angle2radian(i * dphi);
			float theta_temp = angle2radian(j * dtheta);
			cx = radius * std::sin(phi_temp) * std::cos(theta_temp);
			cy = radius * std::cos(phi_temp);
			cz = radius * std::sin(phi_temp) * std::sin(theta_temp);

			vec3 curpos = vec3((float)p1[0] + cx, (float)p1[1] + cy, (float)p1[2] + cz);
			curspherepointsvec.push_back(curpos);


		}
	}


	int indexarray_size = level * level * 2 * 3;

	unsigned int upl, upr, botl, botr;

	//for ,upl,upr,botl,botr
	unsigned int base = this->neurodevertices.size() / 6;
	for (unsigned int i = 0; i < level; i++) {
		for (unsigned int j = 0; j < level; j++) {
			upr = i * (level + 1) + j; upl = upr + 1;
			botr = upr + level + 1; botl = botr + 1;
			int curpos = (i * level + j) * 6;
			this->spherepointsonly.push_back(curspherepointsvec[upr]);
			this->spherepointsonly.push_back(curspherepointsvec[upl]);
			this->spherepointsonly.push_back(curspherepointsvec[botr]);
 
			this->spherepointsonly.push_back(curspherepointsvec[botr]);
			this->spherepointsonly.push_back(curspherepointsvec[botl]);
			this->spherepointsonly.push_back(curspherepointsvec[upl]);


		}
	}

}
void NeuroGeometry2Mesh::calCylinderTempdata(const vec3& p1, const vec3& p2, double r1, double r2, int level) {

	LunaCylender littlecylender(p1, p2, r1, r2, level);

	//ÇóÐý×ª¾ØÕó£¬Î»ÒÆ¾ØÕó


	vec3 middlepoint(
		(p1[0] + p2[0]) / 2.0,
		(p1[1] + p2[1]) / 2.0,
		(p1[2] + p2[2]) / 2.0
	);

	LMatrix4 mytransform(1.0, 0.0, 0.0, 0.0,//col 1
		0.0, 1.0, 0.0, 0.0,//col2
		0.0, 0.0, 1.0, 0.0,//col3
		middlepoint[0], middlepoint[1], middlepoint[2], 1.0//col4
	);
	vec3 p1_p2 = p1 - p2;
	vec3 up = vec3(0.0, 1.0, 0.0);
	vec3 rotateaxis = up.cross(p1_p2);

	float costheta = p1_p2.dot(up) / (p1_p2.length());

	//float _theta=std::acos(costheta)*(-1);
	//float _costheta=std::cos(_theta);
	//float _sintheta=std::sin(_theta);
	float sintheta = std::sqrt(1 - costheta * costheta);
	LMatrix4 myrotate;
	if (rotateaxis == vec3(0.0, 0.0, 0.0)) {}
	else { myrotate = calRotateMatrix(rotateaxis.normal(), costheta, sintheta); }
	littlecylender.move(myrotate, mytransform);




	//insert into cylinderpointsonly
	unsigned int r, l;
	for (unsigned int i = 0; i < level; i++) {
		r = i;
		l = r + 1;

		this->cylinderpointsonly.push_back(littlecylender.circle1[r]);
		this->cylinderpointsonly.push_back(littlecylender.circle1[l]);
		this->cylinderpointsonly.push_back(littlecylender.circle2[r]);

		this->cylinderpointsonly.push_back(littlecylender.circle2[r]);
		this->cylinderpointsonly.push_back(littlecylender.circle2[l]);
		this->cylinderpointsonly.push_back(littlecylender.circle1[l]);

	}
}