#include "Utils.hpp"
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>

using namespace std;
/*
double EdgeLength(const double& x1, const double& y1, const double& z1, const double& x2, const double& y2, const double& z2)
{
	return sqrt(pow(x2 - x1 , 2) + pow(y2 - y1 , 2) + pow(z2 - z1 , 2));	
}
*/ 
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
		if(!(PolyhedralChoice(path, mesh, p, q, b, c, BestPath)))
		{
			return false;
		}
	}
	else if(Id1 != 'n' && Id2 != 'n')
	{
		BestPath = true;
		if(!(PolyhedralChoice(path, mesh, p, q, b, c, BestPath)))
		{
			return false;
		}
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
				addpath = "/PlatonicSolids/Tetraedro";
				break;
			case '4':
				addpath = "/PlatonicSolids/Ottaedro";
				break;
			case '5':
				addpath = "/PlatonicSolids/Icosaedro";
				break;
			default:
				cerr << "I dati inseriti non sono validi" << endl;
				return false;
		}
		filepath = path + addpath;
		if(!ImportMesh(filepath, mesh))
		{
			return false;
		}
		if(!Output(mesh, path))
		{
			return false;
		}
		/*
		if(!(Geodetico(mesh, b, c, BestPath))
		{
			return false;
		}
		*/
		if(q == '3')
		{
			Goldberg(mesh);
		}
		return true;
	}
	return false;
}
// ***************************************************************************
bool ImportMesh(const string& path, PolyhedralMesh& mesh )
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
	for(const auto& i:mesh.Cell0DsCoordinates)
	{
		for(const auto& j:i)
		{
			cout << "[" << j << "] ";
		}
		cout << endl;
	}
	cout << "*************************\n";
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
	for(const auto& i:mesh.Cell1DsVertices)
	{
		for(const auto& j:i)
		{
			cout << "[" << j << "] ";
		}
		cout << endl;
	}
	cout << "*************************\n";
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
	mesh.Cell2DsNum = lines.size();
	mesh.Cell2DsVertices.reserve(mesh.Cell2DsNum);
	mesh.Cell2DsEdges.reserve(mesh.Cell2DsNum);
	mesh.Cell2DsID.reserve(mesh.Cell2DsNum);
	mesh.Cell2DsCentre.reserve(mesh.Cell2DsNum);
	mesh.Cell2DsArches.reserve(mesh.Cell2DsNum);

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
	for(const auto& i:mesh.Cell2DsVertices)
	{
		for(const auto& j:i)
		{
			cout << "[" << j << "] ";
		}
		cout << endl;
	}
	cout << "*************************\n";
	for(const auto& i:mesh.Cell2DsEdges)
	{
		for(const auto& j:i)
		{
			cout << "[" << j << "] ";
		}
		cout << endl;
	}
	return true;
}
// ***************************************************************************
void Goldberg(PolyhedralMesh& mesh)
{
	for(const auto& i : mesh.Cell2DsVertices)
	{
		Vector3d Centre = Vector3d::Zero();
		
		for(const auto& j : i)
		{
			Centre[0] += mesh.Cell0DsCoordinates[j][0];
			Centre[1] += mesh.Cell0DsCoordinates[j][1];
			Centre[2] += mesh.Cell0DsCoordinates[j][2];
		}
		mesh.Cell2DsCentre.push_back(Centre);
	}
	
	for(unsigned int i = 0; i < mesh.Cell2DsNum; i++)
	{
		for(unsigned int j = i+1; j < mesh.Cell2DsNum; j++)
		{
			unsigned int shared = 0;
			for(int vi = 0; vi < 3; vi++) 
			{
				for(int vj = 0; vj < 3; vj++) 
				{
					if (mesh.Cell2DsVertices[i][vi] == mesh.Cell2DsVertices[j][vj]) 
					{
						shared++;
						break;
					}
				}
			}
			if(shared == 2) 
			{
				Vector2i Indices;
				Indices(0) = i;
				Indices(1) = j;
				mesh.Cell2DsArches.push_back(Indices);
			}
		}
	}

	cout << "******* GOLDBERG ********" << endl;
	for(const auto& i:mesh.Cell2DsCentre)
	{
		for(const auto& j:i)
		{
			cout << "["<< j << "] ";
		}
		cout << endl;
	}
	for(const auto& i:mesh.Cell2DsArches)
	{
		for(const auto& j:i)
		{
			cout << "["<< j << "] ";
		}
		cout << endl;
	}
	
	MatrixXd Points(3, mesh.Cell2DsCentre.size());
	MatrixXi Segments(2, mesh.Cell2DsArches.size());
	
	for(size_t i = 0; i < mesh.Cell2DsCentre.size(); i++)
	{
		Points(0,i) = mesh.Cell2DsCentre[i][0];
		Points(1,i) = mesh.Cell2DsCentre[i][1];
		Points(2,i) = mesh.Cell2DsCentre[i][2];
	}
	
	for(size_t i = 0; i < mesh.Cell2DsArches.size(); i++)
	{
		Segments(0,i) = mesh.Cell2DsArches[i][0];
		Segments(1,i) = mesh.Cell2DsArches[i][1];
	}
	
	string File_0D_Path = "/home/appuser/Data/ProgettoPCS2025/Progetto/Output/Cell0DsGoldberg.inp";
	string File_1D_Path = "/home/appuser/Data/ProgettoPCS2025/Progetto/Output/Cell1DsGoldberg.inp";
	
	Gedim::UCDUtilities utilities;
	utilities.ExportPoints(File_0D_Path, Points, {}, {});
	utilities.ExportSegments(File_1D_Path, Points, Segments, {}, {}, {});

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

bool Output(PolyhedralMesh& mesh, const string& path)
{
	cout << "Si Output\n";
	
	if(!OutputCell0Ds(mesh, path))
    {
        cerr << "File Cell0Ds.txt not found" << endl;
        return false;
    }

    if(!OutputCell1Ds(mesh, path))
    {
        cerr << "File Cell1Ds.txt not found" << endl;
        return false;
    }

    if(!OutputCell2Ds(mesh, path))
    {
        cerr << "File Cell2Ds.txt not found" << endl;
        return false;
    }
/*
	if(!OutputCell3Ds(mesh, path))
    {
        cerr << "File Cell3Ds.txt not found" << endl;
        return false;
    }*/
    return true;
}

bool OutputCell0Ds(PolyhedralMesh& mesh, const string& path)
{
	string filePath = path + "/Output/Cell0Ds.txt";
	ofstream file0(filePath);
	
	if(!(file0))
	{
		return false;
	}
	
	file0 << "Numero di vertici: " << mesh.Cell0DsNum << endl;
	file0 << "Id\tCoordinates" << endl;
	for(size_t i = 0; i < mesh.Cell0DsNum; i++)
	{
		file0 << mesh.Cell0DsID[i] << "\t" << fixed << setprecision(16);
		for(size_t j = 0; j < 3; j++)
		{
			if(mesh.Cell0DsCoordinates[i][j] < 0.0)
			{
				file0 << "[";
			}
			else
			{
				file0 << "[+";
			}
			file0 << mesh.Cell0DsCoordinates[i][j] << "] ";
		}
		file0 << endl;
	}
	
	return true;
}

bool OutputCell1Ds(PolyhedralMesh& mesh, const string& path)
{
	string filePath = path + "/Output/Cell1Ds.txt";
	ofstream file1(filePath);
	
	if(!(file1))
	{
		return false;
	}
	
	cout.unsetf(ios::fixed);
	file1 << "Numero di spigoli: " << mesh.Cell1DsNum << endl;
	file1 << "Id\tVertici" << endl;
	for(size_t i = 0; i < mesh.Cell1DsNum; i++)
	{
		file1 << mesh.Cell1DsID[i] << "\t";
		for(size_t j = 0; j < 2; j++)
		{
			file1 << "[" << mesh.Cell1DsVertices[i][j] << "] ";
		}
		file1 << endl;
	}
	
	return true;
}

bool OutputCell2Ds(PolyhedralMesh& mesh, const string& path)
{
	string filePath = path + "/Output/Cell2Ds.txt";
	ofstream file2(filePath);
	
	if(!(file2))
	{
		return false;
	}

	file2 << "Id\tN. Vertici\tId Vertici\tN. Spigoli\tId Spigoli" << endl;

	for(size_t i = 0; i < mesh.Cell2DsNum; i++)
	{
		file2 << mesh.Cell1DsID[i] << "\t" << "3\t\t";
		for(size_t j = 0; j < 3; j++)
		{
			file2 << "[" << mesh.Cell2DsVertices[i][j] << "]";
		}
		file2 << "\t3\t\t";
		for(size_t j = 0; j < 3; j++)
		{
			file2 << "[" << mesh.Cell2DsEdges[i][j] << "]";
		}
		file2 << endl;
	}

	return true;
}
/*
bool OutputCell3Ds(PolyhedralMesh& mesh, const string& path)
{
	string filePath = path + "/Output/Cell3Ds.txt";
	ofstream file3(filePath);
	
	if(!(file3))
	{
		return false;
	}
}
*/
}