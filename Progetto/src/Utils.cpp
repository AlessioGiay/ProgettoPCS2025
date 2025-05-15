#include "Utils.hpp"
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

double EdgeLength(const double& x1, const double& y1, const double& z1, const double& x2, const double& y2, const double& z2)
{
	return sqrt(pow(x2 - x1 , 2) + pow(y2 - y1 , 2) + pow(z2 - z1 , 2));	
}
 
namespace PolyhedralLibrary
{
bool ImportVector(const string& path, PolyhedralMesh& mesh)
{
	string answ;
	char p;
	char q;
	char b;
	char c;
	char Id1;
	char Id2;
	char sep;
	bool BestPath = false;
	
	cout << "Inserire ogni valore nella forma p,q,b,c,Id1,Id2 e premere invio. Per gli Id se non si non si vuole il percorso inviare 'n': " << endl;
	cin >> answ;
	
	stringstream ss(answ);
	ss >> p >> sep >> q >> sep >> b >> sep >> c >> sep >> Id1 >> sep >> Id2;
	
	if(Id1 == 'n' && Id2 == 'n')
	{
		PolyhedralChoice(path, mesh, p, q, b, c, BestPath);
	}
	else if(Id1 != 'n' && Id2 != 'n')
	{
		BestPath = true;
		PolyhedralChoice(path, mesh, p, q, b, c, BestPath);
	}
	else
	{
		cerr << "I dati inseriti non sono validi" << endl;
		return false;
	}
	return true;
}
// ***************************************************************************
bool PolyhedralChoice(const string& path, PolyhedralMesh& mesh, const char& p, const char& q, const char& b, const char& c, bool& BestPath)
{
	string addpath;
	string filepath;

	if(p == '3')
	{
		switch(q) 
		{	
			case '3':
				addpath = "/Tetraedro";
				break;
			case '4':
				addpath = "/Ottaedro";
				break;
			case '5':
				addpath = "/Icosaedro";
				break;
			default:
				return false;
		}
		filepath = path + addpath;
		if(!ImportMesh(filepath, mesh))
		{
			return false;
		}
		/*
		if(!(Geodetico(mesh, b, c, BestPath))
		{
			return false;
		}
		if(q == '3')
		{
			if(!(Goldberg(mesh, b, c, BestPath))
			{
				return false;
			}
		}
		*/
		return true;
	}
	return false;
}
// ***************************************************************************
bool ImportMesh(const string& path, PolyhedralMesh& mesh)
{
	cout << "Si IM\n";
	
	if(!ImportCell0Ds(path, mesh))
    {
        cerr << "File Cell0Ds.csv not found" << endl;
        return false;
    }

    if(!ImportCell1Ds(path, mesh))
    {
        cerr << "File Cell1Ds.csv not found" << endl;
        return false;
    }

    if(!ImportCell2Ds(path, mesh))
    {
        cerr << "File Cell2Ds.csv not found" << endl;
        return false;
    }
    return true;
}
// ***************************************************************************
bool ImportCell0Ds(const string& path, PolyhedralMesh& mesh)
{
	cout << "Si IC0\n";
	
	string filePath = path + "/Cell0Ds.csv";
	ifstream file0(filePath);

    if (!file0)
    {
        return false;
    }
    
	list<string> lines;	
    string line;
    while(getline(file0,line))
    {
        lines.push_back(line);
    }

    lines.pop_front();
	mesh.Cell0DsNum = lines.size();
    mesh.Cell0DsCoordinates.reserve(mesh.Cell0DsNum);
    mesh.Cell0DsID.reserve(mesh.Cell0DsNum);
    mesh.Cell0DsMarker.reserve(mesh.Cell0DsNum);
    
    char sep;
    unsigned int Id = 0;
    Vector3d Coordinates;
    
    for (const auto& j : lines)
    {
		stringstream ss(j);
		ss >> Id >> sep >> Coordinates(0) >> sep >> Coordinates(1) >> sep >> Coordinates(2);
		mesh.Cell0DsID.push_back(Id);
		mesh.Cell0DsCoordinates.push_back(Coordinates);
	}
	return true;
}
// ***************************************************************************
bool ImportCell1Ds(const string& path, PolyhedralMesh& mesh)
{
	cout << "Si IC1\n";
	
	string filePath = path + "/Cell1Ds.csv";
	ifstream file1(filePath);

    if (!file1)
    {
        return false;
    }
    
    list<string> lines;
    string line;
    while(getline(file1,line))
    {
        lines.push_back(line);
    }

    lines.pop_front();
    mesh.Cell1DsNum = lines.size();
    mesh.Cell1DsVertices.reserve(mesh.Cell1DsNum);
    mesh.Cell1DsID.reserve(mesh.Cell1DsNum);
    
    Vector2i Vertices;
    char sep;
    unsigned int Id;
    	
	for(const auto& l : lines)
	{
		stringstream ss(l);
		ss >> Id >> sep >> Vertices(0) >> sep >> Vertices(1);
		
		mesh.Cell1DsID.push_back(Id);
		mesh.Cell1DsVertices.push_back(Vertices);
	}
    return true;
}
// ***************************************************************************
bool ImportCell2Ds(const string& path, PolyhedralMesh& mesh)
{
	cout << "Si IC2\n";
	
	string filePath = path + "/Cell2Ds.csv";
	ifstream file2(filePath);

    if (!file2)
    {
        return false;
    }
    
    list<string> lines;
    string line;
    while(getline(file2,line))
    {
        lines.push_back(line);
    }

    lines.pop_front();
    mesh.NumCell2Ds = lines.size();
    mesh.Cell2DsVertices.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsEdges.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsID.reserve(mesh.NumCell2Ds);
    
    char sep;
    unsigned int Id;
    
	for(const auto& l : lines)
	{
		stringstream ss(l);
		ss >> Id >> sep;

		vector<unsigned int> Vertices(3);
		for(unsigned int i = 0; i < 3; i++)
		{
			ss >> Vertices[i] >> sep;
		}
		
		vector<unsigned int> Edges(3);
		for(unsigned int i = 0; i < 3; i++)
		{
			ss >> Edges[i] >> sep;
		}
		
		mesh.Cell2DsVertices.push_back(Vertices);
		mesh.Cell2DsEdges.push_back(Edges);
		mesh.Cell2DsID.push_back(Id);
	}
	return true;
}
/*
// ***************************************************************************
bool ExpPoints(PolyhedralMesh& mesh, const string& FilePath)
{
	cout << "Si ExP\n";
	
	mesh.Points.resize(3, mesh.NumCell0Ds);
	for(size_t i = 0; i < mesh.NumCell0Ds; i++)
	{
		mesh.Points(0,i) = mesh.Cell0DsCoordinates[i][0];
		mesh.Points(1,i) = mesh.Cell0DsCoordinates[i][1];
		mesh.Points(2,i) = mesh.Cell0DsCoordinates[i][2];
	}

	Gedim::UCDUtilities utilities;
	utilities.ExportPoints(FilePath, mesh.Points, {}, {});

	return true;
}
// ***************************************************************************
bool ExpSegments(PolyhedralMesh& mesh, const string& FilePath)
{
	cout << "Si ExS\n";

	mesh.Segments.resize(2, mesh.NumCell1Ds);
	for(size_t i = 0; i < mesh.NumCell1Ds; i++)
	{
		mesh.Segments(0,i) = mesh.Cell1DsVertices[i][0];
		mesh.Segments(1,i) = mesh.Cell1DsVertices[i][1];
	}

	Gedim::UCDUtilities utilities;
	utilities.ExportSegments(FilePath, mesh.Points, mesh.Segments, {}, {}, {});
		
	return true;
}
*/
}