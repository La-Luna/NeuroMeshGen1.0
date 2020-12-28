#include <QtCore>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>

#include "MeshViewerWidget.h"
//#include "AnaMorph_cellgen.hh"

#include<vector>
#include <utility>
#include <fstream>
#include <map>

//OpenMesh::VProp<std::vector<double>> tempcolor;
std::vector<std::vector<double>> tempcolor;
void normalize(std::vector<double>& t) {
	std::vector<double> ans(t.size());
	//length
	double _lensquar = 0.0;
	for (int i = 0; i < t.size(); i++) {
		_lensquar += t[i] * t[i];
	}
	double _len = sqrt(_lensquar);
	for (int i = 0; i < t.size(); i++) {
		t[i] = t[i] / _len;
	}

}
MeshViewerWidget::MeshViewerWidget(QWidget* parent)
	: QGLViewerWidget(parent),
	ptMin(0.0),
	ptMax(0.0),
	isEnableLighting(true),
	isTwoSideLighting(false),
	isDrawBoundingBox(false),
	isDrawBoundary(false)
{
}

MeshViewerWidget::~MeshViewerWidget(void)
{
}

bool MeshViewerWidget::LoadMesh(const std::string & filename)
{
	Clear();
	bool read_OK;
	//判断一下是什么文件，然后决定调用什么func
	auto _pointpos = filename.find_last_of('.');
	std::string ext = filename.substr(_pointpos + 1);
	//if(ext=="obj") 
		read_OK = MeshTools::ReadMesh(mesh, filename);
	//else {
		//neurogeo = MeshTools::ReadSWC(filename); 
		//ViewAll();
		//this->SetDrawMode(SWCPOINTS);	
	//}

	

	//std::string swcfilepath = "C:/Users/57610/Desktop";
	//std::string swcfilename = "test4.swc";
	//this->oldfiledata = new SwcFileData(swcfilepath, swcfilename);
	//this->neurogeo = new NeuroGeometry2Mesh(oldfiledata);


	std::cout << "Load mesh from file " << filename << std::endl;
	if (read_OK)
	{
		strMeshFileName = QString::fromStdString(filename);
		QFileInfo fi(strMeshFileName);
		strMeshPath = fi.path();
		strMeshBaseName = fi.baseName();
		UpdateMesh();
		update();
		return true;
	}
	return false;
}
//int  MeshViewerWidget::anamorph_sim(int _argc, char* _argv[]) {
//	try {
//		std::srand(0);	// set seed for random number generator to get reproducible outputs
//		AnaMorph_cellgen am_cellgen(_argc, _argv);
//		return am_cellgen.run();
//	}
//	catch (CLAEx& ex) {
//		printf("%s", ex.error_msg.c_str());
//		return EXIT_FAILURE;
//	}
//	catch (...) {
//		printf("ERROR: main(): unhandled exception at top level.\n");
//		return EXIT_FAILURE;
//	}
//}
bool  MeshViewerWidget::LoadSwcFile(const std::string& fname) {

	neurogeo = MeshTools::ReadSWC(fname);
	
	UpdateSWC();
	this->SetDrawMode(SWCPOINTS);
	//update();
	//am_cellgen -i someCell.swc -force-meshing
	//AnaMorph_cellgen(int argc, char *argv[]);
	
	return true;
}
void MeshViewerWidget::Clear(void)
{
	mesh.clear();
}
void MeshViewerWidget::UpdateSWC(void) {
	ptMin[0] = neurogeo->posMin[0];
	ptMin[1] = neurogeo->posMin[1];
	ptMin[2] = neurogeo->posMin[2];

	ptMax[0] = neurogeo->posMax[0];
	ptMax[1] = neurogeo->posMax[1];
	ptMax[2] = neurogeo->posMax[2];


	SetScenePosition((ptMin + ptMax) * 0.5, (ptMin - ptMax).norm() * 0.5);

	std::cout << "Information of the input SWC file:" << std::endl;
	//std::cout << "  [V, E, F] = [" << mesh.n_vertices() << ", " << mesh.n_edges() << ", " << mesh.n_faces() << "]\n";
	std::cout << "the number of points is :  "<<neurogeo->pointsonly.size() << std::endl;
	std::cout << "  BoundingBox:\n";
	std::cout << "  X: [" << ptMin[0] << ", " << ptMax[0] << "]\n";
	std::cout << "  Y: [" << ptMin[1] << ", " << ptMax[1] << "]\n";
	std::cout << "  Z: [" << ptMin[2] << ", " << ptMax[2] << "]\n";
	std::cout << "  Diag length of BBox: " << (ptMax - ptMin).norm() << std::endl;
	//std::cout << "  Edge Length: [" << minlen << ", " << maxlen << "]; AVG: " << avelen / mesh.n_edges() << std::endl;
}
void MeshViewerWidget::UpdateMesh(void)//
{
	mesh.update_normals();
	if (mesh.vertices_empty())
	{
		std::cerr << "ERROR: UpdateMesh() No vertices!" << std::endl;
		return;
	}
	ptMin[0] = ptMin[1] = ptMin[2] = DBL_MAX;
	ptMax[0] = ptMax[1] = ptMax[2] = -DBL_MAX;
	for (const auto& vh : mesh.vertices())
	{
		ptMin.minimize(mesh.point(vh));
		ptMax.maximize(mesh.point(vh));
	}

	double avelen = 0.0;
	double maxlen = 0.0;
	double minlen = DBL_MAX;
	for (const auto& eh : mesh.edges())
	{
		double len = mesh.calc_edge_length(eh);
		maxlen = len > maxlen ? len : maxlen;
		minlen = len < minlen ? len : minlen;
		avelen += len;
	}

	SetScenePosition((ptMin + ptMax)*0.5, (ptMin - ptMax).norm()*0.5);
	std::cout << "Information of the input mesh:" << std::endl;
	std::cout << "  [V, E, F] = [" << mesh.n_vertices() << ", " << mesh.n_edges() << ", " << mesh.n_faces() << "]\n";
	std::cout << "  BoundingBox:\n";
	std::cout << "  X: [" << ptMin[0] << ", " << ptMax[0] << "]\n";
	std::cout << "  Y: [" << ptMin[1] << ", " << ptMax[1] << "]\n";
	std::cout << "  Z: [" << ptMin[2] << ", " << ptMax[2] << "]\n";
	std::cout << "  Diag length of BBox: " << (ptMax - ptMin).norm() << std::endl;
	std::cout << "  Edge Length: [" << minlen << ", " << maxlen << "]; AVG: " << avelen / mesh.n_edges() << std::endl;
}

bool MeshViewerWidget::SaveMesh(const std::string & filename)
{
	return MeshTools::WriteMesh(mesh, filename, DBL_DECIMAL_DIG);
}

bool MeshViewerWidget::ScreenShot()
{
	update();
	QString filename = strMeshPath + "/" + QDateTime::currentDateTime().toString("yyyyMMddHHmmsszzz") + QString(".png");
	QImage image = grabFramebuffer();
	image.save(filename);
	std::cout << "Save screen shot to " << filename.toStdString() << std::endl;
	return true;
}

void MeshViewerWidget::SetDrawBoundingBox(bool b)
{
	isDrawBoundingBox = b;
	update();
}
void MeshViewerWidget::SetDrawBoundary(bool b)
{
	isDrawBoundary = b;
	update();
}
void MeshViewerWidget::EnableLighting(bool b)
{
	isEnableLighting = b;
	update();
}
void MeshViewerWidget::EnableDoubleSide(bool b)
{
	isTwoSideLighting = b;
	update();
}

void MeshViewerWidget::ResetView(void)
{
	ResetModelviewMatrix();
	ViewCenter();
	update();
}

void MeshViewerWidget::ViewCenter(void)
{
	if (!mesh.vertices_empty())
	{
		UpdateMesh();
	}
	update();
}

void MeshViewerWidget::CopyRotation(void)
{
	CopyModelViewMatrix();
}

void MeshViewerWidget::LoadRotation(void)
{
	LoadCopyModelViewMatrix();
	update();
}
double MeshViewerWidget::GetMixedArea(MyMesh::VertexIter v_it) {
	double area = 0.0;
	for (auto vf_it = mesh.vf_iter(*v_it); vf_it.is_valid(); vf_it++) {//寻找vertex的face
		OpenMesh::Vec3d pointTemp[3];
		int k = 0;
		for (auto fv_it = mesh.fv_iter(*vf_it); fv_it.is_valid(); fv_it++) {//寻找face的halfedge
			pointTemp[k] = mesh.point(*fv_it);
			if (k != 0 && (*v_it == *fv_it)) {
				OpenMesh::Vec3d temp;
				temp = pointTemp[0];
				pointTemp[0] = pointTemp[k];
				pointTemp[k] = temp;
			}
			k++;


		}
		double costheta=((pointTemp[1]-pointTemp[0])|(pointTemp[2]-pointTemp[0]))/
						((sqrt((pointTemp[1] - pointTemp[0]) | (pointTemp[1] - pointTemp[0]) ))*
						(sqrt((pointTemp[2] - pointTemp[0]) | (pointTemp[2] - pointTemp[0]))));
		double cosalpha1 = ((pointTemp[0] - pointTemp[1]) | (pointTemp[2] - pointTemp[1]) )/
							((sqrt((pointTemp[0] - pointTemp[1]) | (pointTemp[0] - pointTemp[1]))) *
							(sqrt((pointTemp[2] - pointTemp[1]) | (pointTemp[2] - pointTemp[1]))));
		double cosalpha2=((pointTemp[1] - pointTemp[2]) | (pointTemp[0] - pointTemp[2])) /
							((sqrt((pointTemp[1] - pointTemp[2]) | (pointTemp[1] - pointTemp[2]))) *
							(sqrt((pointTemp[0] - pointTemp[2]) | (pointTemp[0] - pointTemp[2]))));

		if (costheta < 0 || cosalpha1 < 0 || cosalpha2 < 0) {
			double areatemp = 0.5 * sqrt((pointTemp[2]-pointTemp[0])|(pointTemp[2] - pointTemp[0]))
				* sqrt((pointTemp[1] - pointTemp[0]) | (pointTemp[1] - pointTemp[0]))
				* sqrt(1 - costheta * costheta);

				if (costheta < 0) {
					area += 0.5 * areatemp;
				}
				else {
					area += 0.25 * areatemp;
				}
		}
		else {
			double cotalpha1, cotalpha2;
			cotalpha1 = cosalpha1 / sqrt(1-cosalpha1*cosalpha1);
			cotalpha2 = cosalpha2 / sqrt(1 - cosalpha2 * cosalpha2);
			double areatemp = 0.125 * (((pointTemp[2]-pointTemp[0])|(pointTemp[2] - pointTemp[0])) * cotalpha1 + ((pointTemp[1] - pointTemp[0])|(pointTemp[1] - pointTemp[0])) * cotalpha2);
			area += areatemp;
		}


		//std::cout << "area:"<<area<< std::endl;
	}

	return area;
}
void MeshViewerWidget::GaussCurvature() {
	//OpenMesh::VectorT<double, 10> gaussCur;
	std::vector<double> gaussCur;
	gaussCur.clear();
	std::ofstream ofs;
	ofs.open("gaussCur.txt", std::ios::out);
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {//遍历每个点
		double addangle = 0.0;
		for (auto vf_it = mesh.vf_iter(*v_it); vf_it.is_valid(); vf_it++) {//遍历每个点相邻的三角形，将角度保存在addangle
			OpenMesh::Vec3d t[2];
			int _k = 0;
			for (auto fv_it = mesh.fv_iter(*vf_it); fv_it.is_valid();fv_it++) {//每个三角形的点对应向量都取出来
				if (*fv_it != *v_it) {
					t[_k++] = mesh.point(*fv_it) - mesh.point(*v_it);
				}
			}
			double costheta = (t[0] | t[1]) / (sqrt(t[0] | t[0]) * sqrt(t[1] | t[1]));
			double arccostheta = acos(costheta);
			addangle += arccostheta;

		}
		double mixedArea = GetMixedArea(v_it);
		double tempangle = (2 * M_PI - addangle);
		double gaussCurT = tempangle / mixedArea;
		gaussCur.push_back(gaussCurT);
		//ofs << mesh.point(*v_it) << ", ";
		////ofs << gaussCurT << std::endl;
		//ofs << mixedArea<< std::endl;
		

	}

	normalize(gaussCur);
	int _k_ = 0;
	//tempcolor = OpenMesh::VProp<std::vector<double>>(mesh, "tempcolor");
	if (!mesh.has_vertex_colors()) {
		std::cerr << "ERROR:standard vertex property 'color' not available!\n";
		mesh.request_vertex_colors();
	}
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {//循环给每个vertex赋上color
		double r, g, b;
		colorMap(gaussCur[_k_++], r, g, b);
		std::vector<double> t = { r,g,b };
		tempcolor.push_back(t);

		ofs << mesh.point(*v_it) << ", ";
		//ofs << gaussCurT << std::endl;
		ofs << r<<","<<g<<","<<b<< std::endl;
		//tempcolor[*v_it] = t;
		//mesh.set_color(*v_it, OpenMesh::Vec3uc(int(r * 255), int(g * 255), int(b * 255)));
	}
	ofs.close();

}
void MeshViewerWidget::MeanCurvature() {
	std::vector<double> meanCur;
	meanCur.clear();
	std::ofstream ofs;
	ofs.open("meanCur.txt", std::ios::out);
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {//遍历所有顶点
		OpenMesh::Vec3d addij = { 0.0,0.0,0.0 };
		for (auto voh_it = mesh.voh_iter(*v_it); voh_it.is_valid(); voh_it++) {//遍历每一个vertex的out halfedge
			MyMesh::HalfedgeHandle et = *voh_it;
			double cot[2];
			for (int i = 0; i < 2; i++) {
				OpenMesh::Vec3d t[3];
				t[0] = mesh.point(mesh.from_vertex_handle(et));
				t[1] = mesh.point(mesh.to_vertex_handle(et));
				MyMesh::HalfedgeHandle nexte = mesh.next_halfedge_handle(et);
				t[2] = mesh.point(mesh.to_vertex_handle(nexte));
				double cost2 = ((t[0]-t[2])|(t[1]-t[2])) 
								/ (sqrt((t[0] - t[2])|(t[0] - t[2]))*sqrt((t[1] - t[2])|(t[1] - t[2])));
				cot[i] = cost2 / sqrt(1 - cost2 * cost2);
				et = mesh.opposite_halfedge_handle(et);
			}
			OpenMesh::Vec3d vj = mesh.point(mesh.to_vertex_handle(*voh_it));
			OpenMesh::Vec3d vi = mesh.point(*v_it);
			addij += (cot[0]+cot[1]) * (vi - vj);

		}
		double mixArea = GetMixedArea(v_it);
		OpenMesh::Vec3d Ki = 0.5 * addij / mixArea;
		double H = 0.5 * sqrt(Ki | Ki);
		meanCur.push_back(H);
		ofs << mesh.point(*v_it) << ", ";
		ofs <<H<< std::endl;

	}

	ofs.close();
	normalize(meanCur);
	int k = 0;
	if (!mesh.has_vertex_colors()) {
		std::cerr << "ERROR:standard vertex property 'color' not available!\n";
		mesh.request_vertex_colors();
	}
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {
		double r, g, b;
		colorMap(meanCur[k++], r, g, b);
		std::vector<double> temp;
		tempcolor.push_back(temp);
		//mesh.set_color(*v_it, OpenMesh::Vec3uc(int(255 * r), int(g * 255), int(b * 255)));
	}

}
void MeshViewerWidget::colorMap(double gaussCur, double& r, double& g, double& b) {
	const double rone = 0.8;
	const double gone = 1.0;
	const double bone = 1.0;
	double x = (gaussCur < 0.0 ? 0.0 : (gaussCur > 1.0 ? 1.0 : gaussCur));
	if (x < 1.0 / 8.0) {
		r = 0; g = 0;
		b = bone * (0.5 + x / (1.0 / 8.0) * 0.5);
	}
	else if (x < 3.0 / 8.0) {
		r = 0;
		g = gone * (x-1.0/8.0) / (3.0/8.0-1.0/8.0);
		b = bone;
	}
	else if (x < 5.0 / 8.0) {
		r = rone * (x-3.0/8.0) / (5.0/8.0-3.0/8.0);
		g = gone;
		b = bone - (x - 3.0 / 8.0) / (5.0 / 8.0 - 3.0 / 8.0);
	}
	else if (x < 7.0 / 8.0) {
		r = rone;
		g = gone - (x - 5.0 / 8.0) / (7.0 / 8.0 - 5.0 / 8.0);
		b = 0.0;
	}
	else {
		r = rone - (x - 7.0 / 8.0) / (1-7.0 / 8.0)*0.5;
		b = 0;
		g = 0;
	}

}
void MeshViewerWidget::PrintMeshInfo(void)
{
	std::cout << "Mesh Info:\n";
	std::cout << "  [V, E, F] = [" << mesh.n_vertices() << ", " << mesh.n_edges() << ", " << mesh.n_faces() << "]\n";
	std::cout << "  BoundingBox:\n";
	std::cout << "  X: [" << ptMin[0] << ", " << ptMax[0] << "]\n";
	std::cout << "  Y: [" << ptMin[1] << ", " << ptMax[1] << "]\n";
	std::cout << "  Z: [" << ptMin[2] << ", " << ptMax[2] << "]\n";
	std::cout << "  Diag length of BBox: " << (ptMax - ptMin).norm() << std::endl;
}

void MeshViewerWidget::DrawScene(void)
{
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixd(&projectionmatrix[0]);
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixd(&modelviewmatrix[0]);
	//DrawAxis();
	if (isDrawBoundingBox) DrawBoundingBox();
	if (isDrawBoundary) DrawBoundary();
	if (isEnableLighting) glEnable(GL_LIGHTING);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, isTwoSideLighting);
	DrawSceneMesh();
	if (isEnableLighting) glDisable(GL_LIGHTING);
}

void MeshViewerWidget::DrawSceneMesh(void)
{
	if (mesh.n_vertices() == 0&&drawmode>2) { return; }
	if(drawmode>2)SetMaterial();
	switch (drawmode)
	{
	case SWCPOINTS:
		glDisable(GL_LIGHTING);
		DrawSWCPoints();
		break;
	case SWCLINES:
		glDisable(GL_LIGHTING);
		DrawSWCLines();
		break;
	case SWCCYLEN:
		glDisable(GL_LIGHTING);
		DrawSpheres();
		DrawSWCCylender();
		break;
	case POINTS:
		DrawPoints();
		//DrawSWCPoints();
		break;
	case WIREFRAME:
		DrawWireframe();
		break;
	case HIDDENLINES:
		DrawHiddenLines();
		break;
	case FLATLINES:
		DrawFlatLines();
		break;
	case FLAT:
		glColor3d(0.8, 0.8, 0.8);
		DrawFlat();
		break;
	case SMOOTH:
		DrawSmooth();
		break;
	case GAUSSCURV:
		this->GaussCurvature();
		this->DrawGaussCurve();
		//DrawSmooth();

		break;
	default:
		break;
	}
}
void  MeshViewerWidget::DrawSWCPoints() {
	glColor3d(1.0, 0.0, 0.0);
	glPointSize(5);
	glBegin(GL_POINTS);
	for (auto p:this->neurogeo->pointsonly) {
		glVertex3d(p[0], p[1], p[2]);
	}
	glEnd();

}
void MeshViewerWidget::DrawSWCLines() {
	glColor3d(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	int pointlen = this->neurogeo->linepointsonly.size();
	int i;
	for (i = 0; i < pointlen;i+=2 ) {
		vec3 lastp = this->neurogeo->linepointsonly[i];
		vec3 curp = this->neurogeo->linepointsonly[i + 1];
		glVertex3d(lastp[0], lastp[1], lastp[2]);
		glVertex3d(curp[0], curp[1], curp[2]);
	}
	glEnd();

}
void MeshViewerWidget::DrawSWCCylender() {
	glColor3d(0.0, 0.0, 1.0);
	//glBegin(GL_POLYGON_OFFSET_FILL);
	glBegin(GL_LINES);
	for (int i = 0; i < this->neurogeo->cylinderpointsonly.size(); i += 3) {
		vec3 p1 = this->neurogeo->cylinderpointsonly[i];
		vec3 p2 = this->neurogeo->cylinderpointsonly[i + 1];
		vec3 p3 = this->neurogeo->cylinderpointsonly[i + 2];
		glVertex3d(p1[0], p1[1], p1[2]);
		glVertex3d(p2[0], p2[1], p2[2]);

		glVertex3d(p2[0], p2[1], p2[2]);
		glVertex3d(p3[0], p3[1], p3[2]);

		glVertex3d(p3[0], p3[1], p3[2]);
		glVertex3d(p1[0], p1[1], p1[2]);
	}
	glEnd();
}
void MeshViewerWidget::DrawSpheres() {
	glColor3d(1.0, 0.0, 0.0);
	//glBegin(GL_POLYGON_OFFSET_FILL);
	glBegin(GL_LINES);
	for (int i = 0; i < this->neurogeo->spherepointsonly.size(); i += 3) {
		vec3 p1 = this->neurogeo->spherepointsonly[i];
		vec3 p2 = this->neurogeo->spherepointsonly[i+1];
		vec3 p3 = this->neurogeo->spherepointsonly[i+2];
		glVertex3d(p1[0], p1[1], p1[2]);
		glVertex3d(p2[0], p2[1], p2[2]);

		glVertex3d(p2[0], p2[1], p2[2]);
		glVertex3d(p3[0], p3[1], p3[2]);

		glVertex3d(p3[0], p3[1], p3[2]);
		glVertex3d(p1[0], p1[1], p1[2]);
	}
	glEnd();
}
void MeshViewerWidget::DrawPoints(void) const
{ 
	glColor3d(1.0, 0.5, 0.5);
	glPointSize(5);
	glBegin(GL_POINTS);
	for (const auto& vh : mesh.vertices())
	{
		glNormal3dv(mesh.normal(vh).data());
		glVertex3dv(mesh.point(vh).data());
	}
	glEnd();
}

void MeshViewerWidget::DrawWireframe(void) const
{
	glColor3d(0.2, 0.2, 0.2);
	glBegin(GL_LINES);
	for (const auto& eh : mesh.edges())
	{
		auto heh = mesh.halfedge_handle(eh, 0);
		auto vh0 = mesh.from_vertex_handle(heh);
		auto vh1 = mesh.to_vertex_handle(heh);
		glNormal3dv(mesh.normal(vh0).data());
		glVertex3dv(mesh.point(vh0).data());
		glNormal3dv(mesh.normal(vh1).data());
		glVertex3dv(mesh.point(vh1).data());
	}
	glEnd();
}

void MeshViewerWidget::DrawHiddenLines() const
{
	glLineWidth(1.0);
	float backcolor[4];
	glGetFloatv(GL_COLOR_CLEAR_VALUE, backcolor);
	glColor4fv(backcolor);
	glDepthRange(0.01, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	if (glIsEnabled(GL_LIGHTING))
	{
		glDisable(GL_LIGHTING);
		DrawFlat();
		glEnable(GL_LIGHTING);
	}
	else
	{
		DrawFlat();
	}
	glDepthRange(0.0, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glColor3d(.3, .3, .3);
	DrawFlat();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void MeshViewerWidget::DrawFlatLines(void) const
{
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.5f, 2.0f);
	glShadeModel(GL_FLAT);
	//glColor3d(0.8, 0.8, 0.8);
	glColor3d(1.0, 1.0, 1.0);
	DrawFlat();
	glDisable(GL_POLYGON_OFFSET_FILL);
	if (glIsEnabled(GL_LIGHTING))
	{
		glDisable(GL_LIGHTING);
		DrawWireframe();
		glEnable(GL_LIGHTING);
	}
	else
	{
		DrawWireframe();
	}
}

void MeshViewerWidget::DrawFlat(void) const
{
	glBegin(GL_TRIANGLES);
	for (const auto& fh : mesh.faces())
	{
		glNormal3dv(mesh.normal(fh).data());
		for (const auto& fvh : mesh.fv_range(fh))
		{
			glVertex3dv(mesh.point(fvh).data());
		}
	}
	glEnd();
}

void MeshViewerWidget::DrawSmooth(void) const
{
	glColor3d(0.8, 0.8, 0.8);
	glShadeModel(GL_SMOOTH);
	glLoadName(static_cast<GLuint>(mesh.n_vertices()));
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_DOUBLE, 0, mesh.points());
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_DOUBLE, 0, mesh.vertex_normals());
	for (const auto& fh : mesh.faces())
	{
		glBegin(GL_POLYGON);
		for (const auto& fvh : mesh.fv_range(fh))
		{
			glArrayElement(fvh.idx());
		}
		glEnd();
	}
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
}
void MeshViewerWidget::DrawGaussCurve(void) const {


	glPointSize(5);
	glBegin(GL_POINTS);
	MyMesh::VertexIter v_it, v_end(mesh.vertices_end());
	std::vector< std::vector<double>>::iterator color_it;
	for(v_it=mesh.vertices_begin(),
		color_it=tempcolor.begin();v_it!=mesh.vertices_end();v_it++,color_it++)
	
	//for ( auto fh=mesh.faces_begin();fh!=mesh.faces_end();fh++)
	{
		glColor3d((*color_it)[0], (*color_it)[1], (*color_it)[2]);

		glNormal3dv(mesh.normal(*v_it).data());
		glVertex3dv(mesh.point(*v_it).data());


	}
	glEnd();
}
void MeshViewerWidget::DrawBoundingBox(void) const
{
	float linewidth;
	glGetFloatv(GL_LINE_WIDTH, &linewidth);
	glLineWidth(2.0f);
	glColor3d(.3, .7, .3);
	glBegin(GL_LINES);
	for (const auto& i : { 0, 1 })
	{
		for (const auto& j : { 0, 1 })
		{
			for (const auto& k : { 0, 1 })
			{
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(~i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], ~j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], ~k ? ptMin[2] : ptMax[2]);
			}
		}
	}
	glEnd();
	glLineWidth(linewidth);
}

void MeshViewerWidget::DrawBoundary(void) const
{
	float linewidth;
	glGetFloatv(GL_LINE_WIDTH, &linewidth);
	glLineWidth(2.0f);
	glColor3d(0.1, 0.1, 0.1);
	glBegin(GL_LINES);
	for (const auto& eh : mesh.edges())
	{
		if (mesh.is_boundary(eh))
		{
			auto heh = mesh.halfedge_handle(eh, 0);
			auto vh0 = mesh.from_vertex_handle(heh);
			auto vh1 = mesh.to_vertex_handle(heh);
			glNormal3dv(mesh.normal(vh0).data());
			glVertex3dv(mesh.point(vh0).data());
			glNormal3dv(mesh.normal(vh1).data());
			glVertex3dv(mesh.point(vh1).data());
		}
	}
	glEnd();
	glLineWidth(linewidth);
}
