#pragma once
#include "neurogeometry.h"
#include "vec3.h"
#include "swcfiledata.h"
#include "lunacylender.h"

class NeuroGeometry2Mesh : public NeuroGeometry {


public:
	NeuroGeometry2Mesh(SwcFileData* originaldata);
	virtual void getPointDataonly();
	virtual void getLinePointsDataonly();
	void getSpherePointDataonly();
	void getCylinderPointDataonly();
	void calSphereTempdata(const vec3& pos, double radius, int level);
	void calCylinderTempdata(const vec3& p1, const vec3& p2, double r1, double r2, int level);

public:
	vec3 posMin;
	vec3 posMax;
	vector<vec3> pointsonly;
	vector<vec3> linepointsonly;
	vector<vec3> spherepointsonly;
	vector<vec3> cylinderpointsonly;


};