#include "Utils.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <filesystem>
#include <Eigen/Dense>

using namespace std;
namespace fs = filesystem;
 
namespace PolyhedralLibrary
{

bool ImportVector(const string& path, PolyhedralMesh& mesh, PolyhedralData& data)
{
	cout << "\nInserire i valori nella forma:" << endl;
	cout << "  p,q,b,c           (senza percorso)" << endl;
	cout << "  p,q,b,c,Id1,Id2   (con percorso)" << endl;
	cout << "Premere Invio dopo averli scritti tutti su una riga:\n" << endl;

	string answ;
	cin >> answ;
	cout << endl;

	stringstream ss(answ);
	string elmt;
	vector<string> input;

	while (getline(ss, elmt, ',')) 
	{
		input.push_back(elmt);
	}

	if (input.size() != 4 && input.size() != 6) 
	{
		cerr << "Errore: inserire 4 o 6 valori separati da virgola" << endl;
		return false;
	}

	try 
	{
		data.p = stoi(input[0]);
		data.q = stoi(input[1]);
		data.b = stoi(input[2]);
		data.c = stoi(input[3]);

		if (data.p != 3) 
		{
			cerr << "Errore: il valore di p deve essere 3" << endl;
			return false;
		}
		if (data.q != 3 && data.q != 4 && data.q != 5) 
		{
			cerr << "Errore: il valore di q deve essere 3, 4 o 5" << endl;
			return false;
		}
		if (data.b < 0 || data.c < 0)
		{
			cerr << "Errore: b e c devono essere >= 0" << endl;
			return false;
		}

		if (input.size() == 6) 
		{
			const string& id1_str = input[4];
			const string& id2_str = input[5];

			if (all_of(id1_str.begin(), id1_str.end(), ::isdigit) && all_of(id2_str.begin(), id2_str.end(), ::isdigit)) 
			{
				data.Id1 = static_cast<unsigned int>(stoul(id1_str));
				data.Id2 = static_cast<unsigned int>(stoul(id2_str));
				data.BestPath = true;
			}
			else 
			{
				cerr << "Errore: Id1 e Id2 devono essere interi >= 0" << endl;
				return false;
			}
		}
		else
		{
			data.BestPath = false;
		}

		if (!PolyhedralChoice(path, mesh, data))
		{
			return false;
		}
	} 
	catch (const exception& e) 
	{
		cerr << "Errore durante la conversione dei dati: " << e.what() << endl;
		return false;
	}

	return true;
}
// ***************************************************************************
bool PolyhedralChoice(const string& path, PolyhedralMesh& mesh, PolyhedralData& data)
{
	string addpath;
	string filepath;

	switch(data.q) 
	{
		case 3:
			addpath = "/PlatonicSolids/Tetraedro";
			break;
		case 4:
			addpath = "/PlatonicSolids/Ottaedro";
			break;
		case 5:
			addpath = "/PlatonicSolids/Icosaedro";
	}

	filepath = path + addpath;
	if(!ImportMesh(filepath, mesh))
	{
		return false;
	}

	return true;
}
// ***************************************************************************
bool ImportMesh(const string& path, PolyhedralMesh& mesh )
{
	if(!ImportCell0Ds(path, mesh))
	{
		cerr << "File Cell0Ds.csv non trovato" << endl;
		return false;
	}

	if(!ImportCell1Ds(path, mesh))
	{
		cerr << "File Cell1Ds.csv non trovato" << endl;
		return false;
	}

	if(!ImportCell2Ds(path, mesh))
	{
		cerr << "File Cell2Ds.csv non trovato" << endl;
		return false;
	}

	return true;
}
// ***************************************************************************
bool ImportCell0Ds(const string& path, PolyhedralMesh& mesh)
{
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
	
	if (!CheckEdges(mesh)) 
	{
		cerr << "Verifica orientamento fallita: controllare le facce" << endl;
		return false;
	}

	return true;
}
// ***************************************************************************
bool CheckEdges(const PolyhedralMesh& mesh) 
{	
	bool allGood = true;

	for (size_t Id = 0; Id < mesh.Cell2DsNum; Id++) 
	{
		const vector<unsigned int>& faceVerts = mesh.Cell2DsVertices[Id];
		const vector<unsigned int>& faceEdges = mesh.Cell2DsEdges[Id];
		size_t E = faceVerts.size();

		for(size_t i = 0; i < E; ++i) 
		{
			unsigned int v0 = faceVerts[i];
			unsigned int v1 = faceVerts[(i + 1) % E];
			unsigned int EdgeId = faceEdges[i];

			if(EdgeId >= mesh.Cell1DsVertices.size()) 
			{
				cerr << "Errore: edge ID " << EdgeId << " non valido" << endl;
				allGood = false;
				continue;
			}

			unsigned int origin = mesh.Cell1DsVertices[EdgeId][0];
			unsigned int end = mesh.Cell1DsVertices[EdgeId][1];

			bool match = (origin == v0 && end == v1) || (origin == v1 && end == v0);

			if (!match) 
			{
				std::cerr << "Errore nella faccia " << Id << " tra i vertici " << v0 << " -> " << v1 << ", ma il lato " << EdgeId << " connette " << origin << " -> " << end << endl;
				allGood = false;
			}
		}
	}

	return allGood;
}
// ***************************************************************************
bool Output(PolyhedralMesh& mesh, const string& path, PolyhedralData& data)
{
	string folderPath = path + "/Output";
	
	if(!data.Goldby)
	{
		folderPath = folderPath + "/Geodetico";
		try 
		{
			fs::create_directories(folderPath);

			for (const auto& entry : fs::directory_iterator(folderPath)) 
			{
				fs::remove_all(entry.path());
			}
		} 
		catch (const std::exception& e) 
		{
			cerr << "Errore durante la cancellazione: " << e.what() << endl;
			return false;
		}
	}
	else
	{
		folderPath = folderPath + "/Goldberg";
	
		try 
		{
			fs::create_directories(folderPath);

			for (const auto& entry : fs::directory_iterator(folderPath)) 
			{
				fs::remove_all(entry.path());
			}
		} 
		catch (const std::exception& e) 
		{
			cerr << "Errore durante la cancellazione: " << e.what() << endl;
			return false;
		}
	}

	if(!OutputCell0Ds(mesh, folderPath))
	{
		cerr << "Errore nella creazione del file di output Cell0Ds.txt" << endl;
		return false;
	}

	if(!OutputCell1Ds(mesh, folderPath))
	{
		cerr << "Errore nella creazione del file di output Cell1Ds.txt" << endl;
		return false;
	}

	if(!OutputCell2Ds(mesh, folderPath))
	{
		cerr << "Errore nella creazione del file di output Cell2Ds.txt" << endl;
		return false;
	}

	if(!OutputCell3Ds(mesh, folderPath))
	{
		cerr << "Errore nella creazione del file di output Cell3Ds.txt" << endl;
		return false;
	}

	return true;
}
// ***************************************************************************
bool OutputCell0Ds(PolyhedralMesh& mesh, const string& path)
{
	string filePath = path + "/Cell0Ds.txt";
	ofstream file0(filePath);

	if(!(file0))
	{
		return false;
	}

	file0 << "Id\tCoordinates\n" << endl;
	for(size_t i = 0; i < mesh.Cell0DsNum; i++)
	{
		file0 << "[" << mesh.Cell0DsID[i] << "]\t" << fixed << setprecision(16);
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
// ***************************************************************************
bool OutputCell1Ds(PolyhedralMesh& mesh, const string& path)
{
	string filePath = path + "/Cell1Ds.txt";
	ofstream file1(filePath);

	if(!(file1))
	{
		return false;
	}

	cout.unsetf(ios::fixed);

	file1 << "Id\tId Vertici\n" << endl;
	for(size_t i = 0; i < mesh.Cell1DsNum; i++)
	{
		file1 << "[" << mesh.Cell1DsID[i] << "]\t";
		for(size_t j = 0; j < 2; j++)
		{
			file1 << "[" << mesh.Cell1DsVertices[i][j] << "] ";
		}
		file1 << endl;
	}
	return true;
}
// ***************************************************************************
bool OutputCell2Ds(PolyhedralMesh& mesh, const string& path)
{
	string filePath = path + "/Cell2Ds.txt";
	ofstream file2(filePath);

	if(!(file2))
	{
		return false;
	}

	file2 << "Id\tId Vertici\tId Spigoli\n" << endl;

	for(size_t i = 0; i < mesh.Cell2DsNum; i++)
	{
		file2 << "[" << mesh.Cell2DsID[i] << "]\t";
		for(size_t j = 0; j < mesh.Cell2DsVertices[0].size(); j++)
		{
			file2 << "[" << mesh.Cell2DsVertices[i][j] << "]";
		}
		file2 << "\t";
		for(size_t j = 0; j < mesh.Cell2DsEdges[0].size(); j++)
		{
			file2 << "[" << mesh.Cell2DsEdges[i][j] << "]";
		}
		file2 << endl;
	}
	return true;
}
// ***************************************************************************
bool OutputCell3Ds(PolyhedralMesh& mesh, const string& path)
{
	string filePath = path + "/Cell3Ds.txt";
	ofstream file3(filePath);

	if(!(file3))
	{
		return false;
	}

	string type;

	switch(mesh.Cell2DsNum)
	{
		case 4:
			type = "Tetraedro";
			break;
		case 8:
			type = "Ottaedro";
			break;
		case 20:
			type = "Icosaedro";
			break;
		default:
			cerr << "Errore nel numero di facce" << endl;
			return false;
	}

	file3 << "Il solido scelto é: " << type << endl;
	file3 << "\nN. Vertici: " << mesh.Cell0DsNum << "\t\tN. Spigoli: " << mesh.Cell1DsNum << "\t\tN. Facce: " << mesh.Cell2DsNum << endl;
	file3 << "\nId Vertici\tId Spigoli\tId Facce\n" << endl;

	for(size_t i = 0; i < mesh.Cell2DsNum; i++)
	{
		for(size_t j = 0; j < mesh.Cell2DsVertices[0].size(); j++)
		{
			file3 << "[" << mesh.Cell2DsVertices[i][j] << "]";
		}
		file3 << "\t";
		for(size_t j = 0; j < mesh.Cell2DsEdges[0].size(); j++)
		{
			file3 << "[" << mesh.Cell2DsEdges[i][j] << "]";
		}
		file3 << "\t[" << mesh.Cell2DsID[i] << "]" << endl;
	}

	return true;
}
// ***************************************************************************
void Geodetico(PolyhedralMesh& mesh, const string& path)
{	
	MatrixXd Points(3, mesh.Cell0DsNum);
	MatrixXi Segments(2, mesh.Cell1DsNum);

	for(size_t i = 0; i < mesh.Cell0DsNum; i++)
	{
		Points(0,i) = mesh.Cell0DsCoordinates[i][0];
		Points(1,i) = mesh.Cell0DsCoordinates[i][1];
		Points(2,i) = mesh.Cell0DsCoordinates[i][2];
	}

	for(size_t i = 0; i < mesh.Cell1DsNum; i++)
	{
		Segments(0,i) = mesh.Cell1DsVertices[i][0];
		Segments(1,i) = mesh.Cell1DsVertices[i][1];
	}

	string Paraview0Ds = path + "/Output/Geodetico/OriginaleCell0Ds.inp";
	string Paraview1Ds = path + "/Output/Geodetico/OriginaleCell1Ds.inp";

	Gedim::UCDUtilities utilities;
	utilities.ExportPoints(Paraview0Ds, Points, {}, {});
	utilities.ExportSegments(Paraview1Ds, Points, Segments, {}, {}, {});
}
// ***************************************************************************
bool CreateGoldberg(PolyhedralMesh& mesh, PolyhedralMesh& gold)
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
		gold.Cell0DsCoordinates.push_back(Centre);
		gold.Cell0DsID.push_back(gold.Cell0DsNum);
		gold.Cell0DsNum ++;	
	}

	for(unsigned int i = 0; i < mesh.Cell2DsNum; i++)
	{
		for(unsigned int j = i+1; j < mesh.Cell2DsNum; j++)
		{
			unsigned int shared = 0;
			for(unsigned int vi = 0; vi < 3; vi++) 
			{
				for(unsigned int vj = 0; vj < 3; vj++) 
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
				gold.Cell1DsVertices.push_back(Vector2i(i, j));
				gold.Cell1DsID.push_back(gold.Cell1DsNum);
				gold.Cell1DsNum ++;	
			}
		}
	}
	
	if(!(CreateFaces(mesh, gold)))
	{
		return false;
	}

	return true;
}
// ***************************************************************************
void Goldberg(PolyhedralMesh& gold, const string& path)
{
	MatrixXd Points(3, gold.Cell0DsNum);
	MatrixXi Segments(2, gold.Cell1DsNum);

	for(size_t i = 0; i < gold.Cell0DsNum; i++)
	{
		Points(0,i) = gold.Cell0DsCoordinates[i][0];
		Points(1,i) = gold.Cell0DsCoordinates[i][1];
		Points(2,i) = gold.Cell0DsCoordinates[i][2];
	}

	for(size_t i = 0; i < gold.Cell1DsNum; i++)
	{
		Segments(0,i) = gold.Cell1DsVertices[i][0];
		Segments(1,i) = gold.Cell1DsVertices[i][1];
	}

	string Paraview0Ds = path + "/Output/Goldberg/DualeCell0Ds.inp";
	string Paraview1Ds = path + "/Output/Goldberg/DualeCell1Ds.inp";
	
	Gedim::UCDUtilities utilities;
	utilities.ExportPoints(Paraview0Ds, Points, {}, {});
	utilities.ExportSegments(Paraview1Ds, Points, Segments, {}, {}, {});
}
// ***************************************************************************
bool CreateFaces(PolyhedralMesh& mesh, PolyhedralMesh& gold)
{	
	double distMin = 0;
	double dist1 = (mesh.Cell0DsCoordinates[0] - gold.Cell0DsCoordinates[0]).norm();
	double dist2 = (mesh.Cell0DsCoordinates[0] - gold.Cell0DsCoordinates[1]).norm();
	if(dist1 < dist2)
	{
		distMin = dist1;
	}
	else
	{
		distMin = dist2;
	}

	for(unsigned int i = 0; i < mesh.Cell2DsNum; i++)
	{
		gold.Cell2DsID.push_back(i);
		gold.Cell2DsNum++;
		
		vector<unsigned int> FacePoints = {};
			
		for(unsigned int j = 0; j < gold.Cell0DsNum; j++)
		{
			if((mesh.Cell0DsCoordinates[i] - gold.Cell0DsCoordinates[j]).norm() <= distMin * 1.1)
			{
				FacePoints.push_back(j);
			}
		}
		if(FacePoints.size() == 3)
		{
			gold.Cell2DsVertices.push_back(FacePoints);
		}
		else
		{
			cerr << "Errore: la faccia " << i << " ha solo " << FacePoints.size() <<" vertici" << endl;
			return false;
		}
	}
	
	for (unsigned int i = 0; i < gold.Cell2DsVertices.size(); i++)
	{
		vector<unsigned int> FaceEdges;

		for (unsigned int j = 0; j < 3; j++)
		{
			for (unsigned int k = j + 1; k < 3; k++)
			{
				unsigned int v1 = gold.Cell2DsVertices[i][j];
				unsigned int v2 = gold.Cell2DsVertices[i][k];

				for (unsigned int line = 0; line < gold.Cell1DsVertices.size(); line++)
				{
					unsigned int e1 = static_cast<unsigned int>(gold.Cell1DsVertices[line][0]);
					unsigned int e2 = static_cast<unsigned int>(gold.Cell1DsVertices[line][1]);

					if ((e1 == v1 && e2 == v2) || (e1 == v2 && e2 == v1))
					{
						FaceEdges.push_back(line);
						break;
					}
				}
			}
		}

		if (FaceEdges.size() != 3)
		{
			cerr << "Errore: la faccia " << i << " ha solo " << FaceEdges.size() << " lati" << endl;
			return false;
		}

		gold.Cell2DsEdges.push_back(FaceEdges);
	}

	return true;
}
// ***************************************************************************
void CreateTriangulation(const PolyhedralMesh& mesh, PolyhedralData& data, PolyhedralMesh& trg)
{
	Vector3d A;
	Vector3d B;
	Vector3d C;
	
	vector<vector<Vector3d>> Alpha(mesh.Cell2DsNum);

	data.section = data.b + data.c;

	for(size_t count = 0; count < mesh.Cell2DsNum; count++)
	{
		A = mesh.Cell0DsCoordinates[mesh.Cell2DsVertices[count][0]];
		B = mesh.Cell0DsCoordinates[mesh.Cell2DsVertices[count][1]];
		C = mesh.Cell0DsCoordinates[mesh.Cell2DsVertices[count][2]];
		
		for(unsigned int i = 0; i <= data.section; i++)
		{
			for(unsigned int j = 0; j <= data.section - i; j++)
			{
				Vector3d NewPoint;
				unsigned int k = data.section - i - j;

				double I = static_cast<double>(i) / data.section;
				double J = static_cast<double>(j) / data.section;
				double K = static_cast<double>(k) / data.section;

				NewPoint = I*A + J*B + K*C;	

				trg.Cell0DsCoordinates.push_back(NewPoint);
				
				Alpha[count].push_back(NewPoint);
			}
		}
	}

	trg.Cell0DsCoordinates = TrgCleaning(trg);
	trg.Cell1DsVertices = CreateArches(mesh, trg, data, Alpha);	
}
// ***************************************************************************
void Triangulation(PolyhedralMesh& trg, const string& path, PolyhedralData& data)
{
	MatrixXd Points(3, trg.Cell0DsNum);

	for(size_t i = 0; i < trg.Cell0DsNum; i++)
	{
		Points(0,i) = trg.Cell0DsCoordinates[i][0]/(trg.Cell0DsCoordinates[i]).norm();
		Points(1,i) = trg.Cell0DsCoordinates[i][1]/(trg.Cell0DsCoordinates[i]).norm();
		Points(2,i) = trg.Cell0DsCoordinates[i][2]/(trg.Cell0DsCoordinates[i]).norm();
	}

	MatrixXi Segments(2,trg.Cell1DsNum);

	for(size_t i = 0; i < trg.Cell1DsNum; i++)
	{
		Segments(0,i) = trg.Cell1DsVertices[i][0];
		Segments(1,i) = trg.Cell1DsVertices[i][1];
	}
	
	if(data.BestPath)
	{
		ShortPath(trg, data);
	}
	
	VectorXi Materials0Ds(trg.Cell0DsMarker.size());
	for (size_t i = 0; i < trg.Cell0DsMarker.size(); ++i)
	{
		Materials0Ds[i] = static_cast<int>(trg.Cell0DsMarker[i]);
	}
	VectorXi Materials1Ds(trg.Cell1DsMarker.size());
	for (size_t i = 0; i < trg.Cell1DsMarker.size(); ++i)
	{
		Materials1Ds[i] = static_cast<int>(trg.Cell1DsMarker[i]);
	}
	
	string Paraview0Ds;
	string Paraview1Ds;

	if(!data.Goldby)
	{
		Paraview0Ds = path + "/Output/Geodetico/GeodeticoCell0Ds.inp";
		Paraview1Ds = path + "/Output/Geodetico/GeodeticoCell1Ds.inp";
	}
	else
	{
		Paraview0Ds = path + "/Output/Goldberg/GoldbergCell0Ds.inp";
		Paraview1Ds = path + "/Output/Goldberg/GoldbergCell1Ds.inp";
	}

	Gedim::UCDUtilities utilities;
	utilities.ExportPoints(Paraview0Ds, Points, {}, Materials0Ds);
	utilities.ExportSegments(Paraview1Ds, Points, Segments, {}, {}, Materials1Ds);
}
// ***************************************************************************
vector<Vector3d> TrgCleaning(PolyhedralMesh& trg)
{	
	vector<Vector3d> Clean;

	for(unsigned int i = 0; i < trg.Cell0DsCoordinates.size(); i++)
	{
		bool IsDuplicate = false;

		for (const auto& j : Clean)
		{
			if (trg.Cell0DsCoordinates[i] == j)
			{
				IsDuplicate = true;
				break;
			}
		}

		if (!IsDuplicate)
		{
			Clean.push_back(trg.Cell0DsCoordinates[i]);
			trg.Cell0DsID.push_back(trg.Cell0DsNum);
			trg.Cell0DsNum ++;
		}
	}

	return Clean;
}
// ***************************************************************************
vector<Vector2i> CreateArches(const PolyhedralMesh& mesh, PolyhedralMesh& trg, const PolyhedralData& data, vector<vector<Vector3d>>& Alpha)
{
	vector<Vector2i> Arches;
	double ComparisonDistance = (mesh.Cell0DsCoordinates[2] - mesh.Cell0DsCoordinates[1]).norm();
	double expectedLength = ComparisonDistance / data.section;
	double tolerance = expectedLength * 0.1;

	for(unsigned int i = 0; i < trg.Cell0DsCoordinates.size(); ++i)
	{
		for(unsigned int j = i + 1; j < trg.Cell0DsCoordinates.size(); ++j)
		{
			double distance = (trg.Cell0DsCoordinates[i] - trg.Cell0DsCoordinates[j]).norm();

			if (abs(distance - expectedLength) > tolerance)
				continue;

			for(unsigned int riga = 0; riga < Alpha.size(); riga++)
			{
				bool Point1 = false;
				bool Point2 = false;

				for(const auto& punto : Alpha[riga])
				{
					if((trg.Cell0DsCoordinates[i] - punto).norm() <= tolerance)
						Point1 = true;
					if((trg.Cell0DsCoordinates[j] - punto).norm() <= tolerance)
						Point2 = true;
				}

				if(Point1 && Point2)
				{
					Arches.push_back(Vector2i(i, j));
					trg.Cell1DsID.push_back(trg.Cell1DsNum);
					trg.Cell1DsNum++;
					break;
				}
			}
		}
	}

	return Arches;
}
// ***************************************************************************
bool ShortPath(PolyhedralMesh& trg, PolyhedralData& data)
{
	MatrixXi adj = CreateAdjacencyMatrix(trg);

	if (!(data.Id1 < trg.Cell0DsNum))
	{
		cerr << "Errore: l'Id del vertice di inizio percorso non corrisponde a nessun punto" << endl;
		return false;
	}

	if (!(data.Id2 < trg.Cell0DsNum))
	{
		cerr << "Errore: l'Id del vertice di fine percorso non corrisponde a nessun punto" << endl;
		return false;
	}

	unsigned int start = data.Id1;
	unsigned int end = data.Id2;
	unsigned int n = adj.rows();
	
	unsigned int sentinel = trg.Cell0DsNum + 5;

	vector<unsigned int> prev(n, sentinel);
	vector<bool> visited(n, false);
	vector<unsigned int> queue;
	unsigned int front = 0;

	visited[start] = true;
	queue.push_back(start);

	while (front < queue.size())
	{
		unsigned int u = queue[front++];
		if (u == end) break;

		for (unsigned int v = 0; v < n; ++v)
		{
			if (adj(u, v) && !visited[v])
			{
				visited[v] = true;
				prev[v] = u;
				queue.push_back(v);
			}
		}
	}

	vector<unsigned int> path;
	if (prev[end] == sentinel)
	{
		cerr << "Nessun cammino trovato tra i due vertici." << endl;
		return false;
	}

	for (unsigned int at = end; at != sentinel; at = prev[at])
		path.push_back(at);

	reverse(path.begin(), path.end());

	cout <<  "******************************" << endl;
	
	if(!(data.Goldby))
	{
		cout << "\tGeodetico" << endl;
	}
	else
	{
		cout << "\tGoldberg" << endl;
	}
	
	cout << "Distanza : " << (path.size() - 1) * (trg.Cell0DsCoordinates[0]-trg.Cell0DsCoordinates[1]).norm() << endl;

	trg.Cell0DsMarker = vector<unsigned int>(trg.Cell0DsNum, 0);
	trg.Cell1DsMarker = vector<unsigned int>(trg.Cell1DsNum, 0);

	for (unsigned int vid : path)
		trg.Cell0DsMarker[vid] = 1;

	for (size_t i = 0; i < path.size() - 1; ++i)
	{
		unsigned int u = path[i];
		unsigned int v = path[i + 1];
		for (size_t id = 0; id < trg.Cell1DsVertices.size(); id++)
		{
			const auto& edge = trg.Cell1DsVertices[id];
			if ((edge[0] == static_cast<int>(u) && edge[1] == static_cast<int>(v)) ||
				(edge[0] == static_cast<int>(v) && edge[1] == static_cast<int>(u)))
			{
				trg.Cell1DsMarker[id] = 1;
				break;
			}
		}
	}
	
	cout << "Percorso minimo trovato: ";
	for (const auto& i : path)
		cout << i << " ";
	cout << endl << "******************************" << endl;

	return true;
}
// ***************************************************************************
MatrixXi CreateAdjacencyMatrix(PolyhedralMesh& trg)
{	
	MatrixXi Adjacency = MatrixXi::Zero(trg.Cell0DsNum, trg.Cell0DsNum);
	
	for(const auto& i:trg.Cell1DsVertices)
	{
		unsigned int u = i[0];
		unsigned int v = i[1];
		
		Adjacency(u,v) = 1;
		Adjacency(v,u) = 1;
	}
	
	return Adjacency;
}

}