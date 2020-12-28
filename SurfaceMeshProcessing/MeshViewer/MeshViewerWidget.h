#pragma once
#include <QString>
#include "QGLViewerWidget.h"
#include "MeshDefinition.h"
#include "swcfiledata.h"
#include "neurogeometry2mesh.h"

#include <vector>


class MeshViewerWidget  : public QGLViewerWidget
{
	Q_OBJECT
public:
	MeshViewerWidget(QWidget* parent = 0);
	virtual ~MeshViewerWidget(void);
	bool LoadMesh(const std::string & filename);
	bool LoadSwcFile(const std::string& filename);

	//int anamorph_sim(int _argc,char* _argv[]);

	void Clear(void);
	void UpdateMesh(void);
	void UpdateSWC(void);
	bool SaveMesh(const std::string & filename);
	bool ScreenShot(void);
	void SetDrawBoundingBox(bool b);
	void SetDrawBoundary(bool b);
	void EnableLighting(bool b);
	void EnableDoubleSide(bool b);
	void ResetView(void);
	void ViewCenter(void);
	void CopyRotation(void);
	void LoadRotation(void);
	double GetMixedArea(MyMesh::VertexIter v_it);
	void GaussCurvature();
	void MeanCurvature();
	void colorMap(double gaussCur, double& r, double& g, double& b);

public:
	//OpenMesh::VectorT<double,10> gaussCur;
	//std::vector<double> gaussCur;
	//std::vector<double>
		//normalize();

signals:
	void LoadMeshOKSignal(bool, QString);
public slots:
	void PrintMeshInfo(void);
protected:
	virtual void DrawScene(void) override;
	void DrawSceneMesh(void);

private:
	void DrawSWCPoints();
	void DrawSWCLines();
	void DrawSWCCylender();
	void DrawSpheres();

	void DrawPoints(void) const;
	void DrawWireframe(void) const;
	void DrawHiddenLines(void) const;
	void DrawFlatLines(void) const;
	void DrawFlat(void) const;
	void DrawSmooth(void) const;
	void DrawGaussCurve(void) const;
	void DrawBoundingBox(void) const;
	void DrawBoundary(void) const;
protected:

	MyMesh mesh;
	SwcFileData* oldfiledata;
	NeuroGeometry2Mesh* neurogeo;
	//AnaMorph_cellgen am_cellgen;

	QString strMeshFileName;
	QString strMeshBaseName;
	QString strMeshPath;
	MyMesh::Point ptMin;
	MyMesh::Point ptMax;
	bool isEnableLighting;
	bool isTwoSideLighting;
	bool isDrawBoundingBox;
	bool isDrawBoundary;
};
