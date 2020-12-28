#include <OpenMesh/Core/IO/MeshIO.hh>
#include "MeshDefinition.h"
#include "swcfiledata.h"
#include <queue>
#include <iostream>
#include <fstream>
#include <cctype>
NeuroGeometry2Mesh* MeshTools::ReadSWC(const std::string& fname) {
	//std::string fname = "C:/Users/57610/Desktop/test4.swc";
	auto slash = fname.find_last_of('/');
	auto backslash = fname.find_last_of('\\');
	std::string filepath;
	std::string filename;
	std::string ext("swc");
	if (slash != std::string::npos && backslash != std::string::npos) {
		slash = slash > backslash ? slash : backslash;
		filepath = fname.substr(0, slash);
	}
	else if (slash != std::string::npos) {
		filepath = fname.substr(0, slash);
	}
	else if (backslash != std::string::npos) {
		slash = backslash;
		filepath = fname.substr(0, slash);
	}

	filename = fname.substr(slash + 1);
	auto pointpos = filename.find_last_of('.');
	if (pointpos != std::string::npos) {
		ext = filename.substr(pointpos + 1);
		filename = filename.substr(0, pointpos);
	}


	SwcFileData* oldfiledata = new SwcFileData(filepath, filename + '.' + ext);
	NeuroGeometry2Mesh* neuroge = new NeuroGeometry2Mesh(oldfiledata);


	//this->SetDrawMode(SWCPOINTS);


	return neuroge;
}
bool MeshTools::ReadMesh(MyMesh & mesh, const std::string & filename)
{
	//filename	"C:/Users/57610/Desktop/Bigmax_White_OBJ.obj"
	std::string path(".");
	auto slash = filename.find_last_of('/');
	auto backslash = filename.find_last_of('\\');
	if (slash != std::string::npos && backslash != std::string::npos)
	{
		slash = slash > backslash ? slash : backslash;
		path = filename.substr(0, slash);
	}
	else if (backslash != std::string::npos)
	{
		slash = backslash;
		path = filename.substr(0, slash);
	}
	else if (slash != std::string::npos)
	{
		path = filename.substr(0, slash);
	}

	std::string basename = filename.substr(slash + 1);
	auto point = basename.find_last_of('.');
	std::string ext("obj");
	if (point != std::string::npos)
	{
		ext = basename.substr(point + 1);
		basename = basename.substr(0, point);
	}
	std::transform(ext.begin(), ext.end(), ext.begin(), std::tolower);
	if (ext == "obj")
	{
		//return true;
		return ReadOBJ(mesh, path + "/" + basename + "." + ext);
	}
	else
	{
		//std::cout << "Error: the file extension " << ext << " is not supported." << std::endl;
		return OpenMesh::IO::read_mesh(mesh, filename);
	}
}

bool MeshTools::ReadOBJ(MyMesh & mesh, const std::string & filename)
{
	std::ifstream ifs(filename);
	if (!ifs.is_open())
	{
		std::cerr << "Error: cannot open file " << filename << std::endl;
		return false;
	}
	std::string line;
	std::string keyWrd;
	MyMesh::Point::value_type x, y, z;
	std::vector<MyMesh::VertexHandle> vertexHandles;
	std::stringstream stream, lineData, tmp;

	// pass 1: read vertices
	while (ifs && !ifs.eof())
	{
		std::getline(ifs, line);
		if (ifs.bad())
		{
			std::cerr << "  Warning! Could not read file properly!\n";
			return false;
		}

		// Trim Both leading and trailing spaces
		auto start = line.find_first_not_of(" \t\r\n");
		auto end = line.find_last_not_of(" \t\r\n");
		if ((std::string::npos == start) || (std::string::npos == end))
			line = "";
		else
			line = line.substr(start, end - start + 1);

		if (line.size() == 0 || line[0] == '#' || isspace(line[0]))
		{
			continue;
		}
		stream.str(line);
		stream.clear();
		stream >> keyWrd;
		if (keyWrd == "v")
		{
			stream >> x >> y >> z;
			if (!stream.fail())
			{
				vertexHandles.push_back(mesh.add_vertex(MyMesh::Point(x, y, z)));
			}
		}
	}

	// reset stream for second pass
	ifs.clear();
	ifs.seekg(0, std::ios::beg);

	int nCurrentPositions = 0;

	// pass 2: read faces
	while (ifs && !ifs.eof())
	{
		std::getline(ifs, line);
		if (ifs.bad())
		{
			std::cerr << "  Warning! Could not read file properly!\n";
			return false;
		}

		// Trim Both leading and trailing spaces
		auto start = line.find_first_not_of(" \t\r\n");
		auto end = line.find_last_not_of(" \t\r\n");
		line = ((std::string::npos == start) || (std::string::npos == end)) ? "" : line.substr(start, end - start + 1);

		// comment
		if (line.size() == 0 || line[0] == '#' || isspace(line[0]))
		{
			continue;
		}

		stream.str(line);
		stream.clear();
		stream >> keyWrd;

		// track current number of parsed vertex attributes,
		// to allow for OBJs negative indices
		if (keyWrd == "v")
		{
			++nCurrentPositions;
		}

		// faces
		else if (keyWrd == "f")
		{
			int value;

			// read full line after detecting a face
			std::string faceLine;
			std::getline(stream, faceLine);
			lineData.str(faceLine);
			lineData.clear();

			MyMesh::FaceHandle fh;
			std::vector<MyMesh::VertexHandle> faceVertices;

			// work on the line until nothing left to read
			while (!lineData.eof())
			{
				// read one block from the line ( vertex/texCoord/normal )
				std::string vertex;
				lineData >> vertex;

				//get the component (vertex/texCoord/normal)
				size_t found = vertex.find("/");

				// parts are seperated by '/' So if no '/' found its the last component
				if (found != std::string::npos)
				{
					// read the index value
					tmp.str(vertex.substr(0, found));
					tmp.clear();

					// If we get an empty string this property is undefined in the file
					if (vertex.substr(0, found).empty())
					{
						continue;
					}

					// Read current value
					tmp >> value;
				}
				else
				{
					// last component of the vertex, read it.
					tmp.str(vertex);
					tmp.clear();
					tmp >> value;

					// Nothing to read here ( garbage at end of line )
					if (tmp.fail())
					{
						continue;
					}
				}

				if (value < 0) {
					// Calculation of index :
					// -1 is the last vertex in the list
					// As obj counts from 1 and not zero add +1
					value = nCurrentPositions + value + 1;
				}
				// Obj counts from 1 and not zero .. array counts from zero therefore -1
				faceVertices.push_back(MyMesh::VertexHandle(value - 1));
			}

			// note that add_face can possibly triangulate the faces, which is why we have to
			// store the current number of faces first
			size_t n_faces = mesh.n_faces();
			auto endIter = faceVertices.end();
			for (auto iter = faceVertices.begin(); iter != endIter; ++iter)
			{
				endIter = std::remove(iter + 1, endIter, *(iter));
			}
			faceVertices.erase(endIter, faceVertices.end());

			//A minimum of three vertices are required.
			if (faceVertices.size() > 2)
			{
				fh = mesh.add_face(faceVertices);
			}
			std::vector<MyMesh::FaceHandle> newfaces;
			for (size_t i = 0; i < mesh.n_faces() - n_faces; ++i)
			{
				newfaces.push_back(MyMesh::FaceHandle(int(n_faces + i)));
			}
		}
	}
	ifs.close();
	return true;
}

bool MeshTools::WriteMesh(const MyMesh & mesh, const std::string & filename, const std::streamsize & precision)
{
	std::string path(".");
	auto slash = filename.find_last_of('/');
	auto backslash = filename.find_last_of('\\');
	if (slash != std::string::npos && backslash != std::string::npos)
	{
		slash = slash > backslash ? slash : backslash;
		path = filename.substr(0, slash);
	}
	else if (backslash != std::string::npos)
	{
		slash = backslash;
		path = filename.substr(0, slash);
	}
	else if (slash != std::string::npos)
	{
		path = filename.substr(0, slash);
	}

	std::string basename = filename.substr(slash + 1);
	auto point = basename.find_last_of('.');
	std::string ext("obj");
	if (point != std::string::npos)
	{
		ext = basename.substr(point + 1);
		basename = basename.substr(0, point);
	}
	std::transform(ext.begin(), ext.end(), ext.begin(), std::tolower);
	if (ext == "obj")
	{
		return WriteOBJ(mesh, path + "/" + basename + "." + ext, precision);
	}
	else
	{
		//std::cout << "Error: the file extension " << ext << " is not supported." << std::endl;
		return OpenMesh::IO::write_mesh(mesh, filename, OpenMesh::IO::Options::Default, precision);
	}
}

bool MeshTools::WriteOBJ(const MyMesh & mesh, const std::string & filename, const std::streamsize & precision)
{
	std::ofstream ofs(filename);
	if (!ofs.is_open())
	{
		std::cerr << "Error: cannot open file " << filename << std::endl;
		return false;
	}
	ofs.precision(precision);
	for (const auto& vh : mesh.vertices())
	{
		ofs << "v " << mesh.point(vh) << std::endl;
	}
	for (const auto& fh : mesh.faces())
	{
		ofs << "f";
		for (const auto& fvh : mesh.fv_range(fh))
		{
			ofs << " " << fvh.idx() + 1;
		}
		ofs << std::endl;
	}
	ofs.close();
	return true;
}

double MeshTools::Area(const MyMesh & mesh)
{
	double area = 0.0;
	for (const auto & f : mesh.faces())
	{
		auto heh = mesh.halfedge_handle(f);
		const auto & p0 = mesh.point(mesh.from_vertex_handle(heh));
		const auto & p1 = mesh.point(mesh.to_vertex_handle(heh));
		const auto & p2 = mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(heh)));
		area += 0.5 * ((p1 - p0) % (p2 - p0)).norm();
	}
	return area;
}

double MeshTools::AverageEdgeLength(const MyMesh & mesh)
{
	if (mesh.edges_empty()) return 0.0;
	double l = 0.0;
	for (const auto & eh : mesh.edges())
	{
		l += mesh.calc_edge_length(eh);
	}
	l /= mesh.n_edges();
	return l;
}

bool MeshTools::HasBoundary(const MyMesh & mesh)
{
	if (mesh.halfedges_empty()) return false;
	for (const auto & heh : mesh.halfedges())
	{
		if (mesh.is_boundary(heh))
		{
			return true;
		}
	}
	return false;
}

bool MeshTools::HasOneComponent(const MyMesh & mesh)
{
	if (mesh.faces_empty()) return false;
	std::vector<int> visit(mesh.n_faces(), 0);
	visit[0] = 1;
	std::queue<int> q;
	q.push(0);
	while (!q.empty())
	{
		int fid = q.front();
		q.pop();
		for (const auto & ffh : mesh.ff_range(mesh.face_handle(fid)))
		{
			if (!visit[ffh.idx()])
			{
				visit[ffh.idx()] = 1;
				q.push(ffh.idx());
			}
		}
	}
	for (const auto & v : visit)
	{
		if (!v)
		{
			return false;
		}
	}
	return true;
}

int MeshTools::Genus(const MyMesh & mesh)
{
	if (HasBoundary(mesh) || !HasOneComponent(mesh)) return -1;
	return 1 - ((int)mesh.n_vertices() + (int)mesh.n_faces() - (int)mesh.n_edges()) / 2;
}

void MeshTools::BoundingBox(const MyMesh & mesh, MyMesh::Point & bmax, MyMesh::Point & bmin)
{
	if (mesh.vertices_empty()) return;
	bmax = bmin = mesh.point(*mesh.vertices_begin());
	for (const auto vh : mesh.vertices())
	{
		bmax.maximize(mesh.point(vh));
		bmin.minimize(mesh.point(vh));
	}
}

// This function should be used after collapse or split, since the indices of
//   edges may be changed after a local modification. By using this function,
//   the indices of edges are the same as a mesh loaded from obj files.
void MeshTools::Reassign(const MyMesh & mesh1, MyMesh & mesh2)
{
	mesh2.clear();
	for (const auto & vh : mesh1.vertices())
	{
		mesh2.add_vertex(mesh1.point(vh));
	}
	for (const auto & fh : mesh1.faces())
	{
		std::vector<MyMesh::VertexHandle> vhs;
		for (const auto & fvh : mesh1.fv_range(fh))
		{
			vhs.push_back(mesh2.vertex_handle(fvh.idx()));
		}
		mesh2.add_face(vhs);
	}
}
