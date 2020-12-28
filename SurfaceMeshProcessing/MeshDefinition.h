#pragma once
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "neurogeometry2mesh.h"
#ifdef _DEBUG
#pragma comment(lib, "OpenMeshCored.lib")
#pragma comment(lib, "OpenMeshToolsd.lib")
#else
#pragma comment(lib, "OpenMeshCore.lib")
#pragma comment(lib, "OpenMeshTools.lib")
#endif
struct MeshTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::Vec3d Point;
	typedef OpenMesh::Vec3d Normal;
	//typedef OpenMesh::Vec3i Color;
	VertexAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
	FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
	EdgeAttributes(OpenMesh::Attributes::Status);
	HalfedgeAttributes(OpenMesh::Attributes::Status);
};
typedef OpenMesh::TriMesh_ArrayKernelT<MeshTraits> MyMesh;

class MeshTools
{
public:
	static NeuroGeometry2Mesh* ReadSWC(const std::string& filename);
	static bool ReadMesh(MyMesh & mesh, const std::string & filename);
	static bool ReadOBJ(MyMesh & mesh, const std::string & filename);
	//static bool ReadOFF(Mesh & mesh, const std::string & filename);
	static bool WriteMesh(const MyMesh & mesh, const std::string & filename, const std::streamsize & precision = 6);
	static bool WriteOBJ(const MyMesh & mesh, const std::string & filename, const std::streamsize & precision = 6);
	static double Area(const MyMesh & mesh);
	static double AverageEdgeLength(const MyMesh & mesh);
	static bool HasBoundary(const MyMesh & mesh);
	static bool HasOneComponent(const MyMesh & mesh);
	static int Genus(const MyMesh & mesh);
	static void BoundingBox(const MyMesh & mesh, MyMesh::Point & bmax, MyMesh::Point & bmin);
	static void Reassign(const MyMesh & mesh1, MyMesh & mesh2);
};
