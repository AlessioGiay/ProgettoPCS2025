#include "gtest/gtest.h"
#include "PolyhedralMesh.hpp"
#include "Utils.hpp"
#include <fstream>

using namespace PolyhedralLibrary;
using namespace Eigen;

TEST(UtilsTest, TestImportCell0Ds)
{
	PolyhedralMesh mesh;
	string path = "/home/appuser/Data/ProgettoPCS2025/Progetto/Prove";
	
	ImportCell0Ds(path, mesh);
	
	EXPECT_EQ(mesh.Cell0DsNum, 3);
	
	vector<unsigned int> Id = {0,1,2};
	EXPECT_EQ(mesh.Cell0DsID, Id);
	
	vector<Vector3d> Coo = {
								Vector3d(0.0,0.0,0.0),
								Vector3d(1.0,0.0,0.0),
								Vector3d(0.0,1.0,0.0)
							};
	for (size_t i = 0; i < Coo.size(); i++)
	{
		EXPECT_TRUE(mesh.Cell0DsCoordinates[i].isApprox(Coo[i], 1e-9));
	}
}

TEST(UtilsTest, TestImportCell1Ds)
{
	PolyhedralMesh mesh;
	string path = "/home/appuser/Data/ProgettoPCS2025/Progetto/Prove";
	
	ImportCell1Ds(path, mesh);
	
	EXPECT_EQ(mesh.Cell1DsNum, 3);
	
	vector<unsigned int> Id = {0,1,2};
	EXPECT_EQ(mesh.Cell1DsID, Id);
	
	vector<Vector2i> Ver = {
								Vector2i(0,1),
								Vector2i(1,2),
								Vector2i(2,0)
							};
	EXPECT_EQ(mesh.Cell1DsVertices, Ver);
}

TEST(UtilsTest, TestImportCell2Ds)
{
	PolyhedralMesh mesh;
	string path = "/home/appuser/Data/ProgettoPCS2025/Progetto/Prove";
	
	mesh.Cell1DsVertices = { 
								Vector2i(0,1),
								Vector2i(1,2),
								Vector2i(2,0)
							};
	
	ImportCell2Ds(path, mesh);
	
	EXPECT_EQ(mesh.Cell2DsNum, 1);
	
	vector<unsigned int> Id = {0};
	EXPECT_EQ(mesh.Cell2DsID, Id);
	
	vector<vector<unsigned int>> Ver = {{0,1,2}};
	EXPECT_EQ(mesh.Cell2DsVertices, Ver);
	
	vector<vector<unsigned int>> Edg = {{0,1,2}};
	EXPECT_EQ(mesh.Cell2DsEdges, Edg);
}

TEST(UtilsTest, TestOutputCell0Ds)
{	
	string path = "/home/appuser/Data/ProgettoPCS2025/Progetto/Prove";
	PolyhedralMesh mesh;

	mesh.Cell0DsNum = 3;
	mesh.Cell0DsID = {0, 1, 2};
	mesh.Cell0DsCoordinates = 
	{
		Vector3d(0.0, 0.0, 0.0),
		Vector3d(1.0, 0.0, 0.0),
		Vector3d(0.0, 1.0, 0.0)
	};
	mesh.Cell0DsMarker = {0, 0, 0};

	mesh.Cell1DsNum = 3;
	mesh.Cell1DsID = {0, 1, 2};
	mesh.Cell1DsVertices = 
	{
		Vector2i(0, 1),
		Vector2i(1, 2),
		Vector2i(2, 0)
	};
	mesh.Cell1DsMarker = {0, 0, 0};

	mesh.Cell2DsNum = 1;
	mesh.Cell2DsID = {0};
	mesh.Cell2DsVertices = {{0, 1, 2}};
	mesh.Cell2DsEdges = {{0, 1, 2}};
	
	OutputCell0Ds(mesh, path);
	
	string path1 = path + "/Output/Cell0Ds.txt";
	string path2 = path + "/Cell0Ds.txt";
	
	ifstream RealOutput(path1);
	ifstream ExpetedOutput(path2);
	
	string RealLine;
	string ExpetedLine;
	
	while (getline(ExpetedOutput, ExpetedLine) && getline(RealOutput, RealLine))
	{
		EXPECT_EQ(ExpetedLine, RealLine);
	}
}

TEST(UtilsTest, TestOutputCell1Ds)
{	
	string path = "/home/appuser/Data/ProgettoPCS2025/Progetto/Prove";
	PolyhedralMesh mesh;

	mesh.Cell0DsNum = 3;
	mesh.Cell0DsID = {0, 1, 2};
	mesh.Cell0DsCoordinates = 
	{
		Vector3d(0.0, 0.0, 0.0),
		Vector3d(1.0, 0.0, 0.0),
		Vector3d(0.0, 1.0, 0.0)
	};
	mesh.Cell0DsMarker = {0, 0, 0};

	mesh.Cell1DsNum = 3;
	mesh.Cell1DsID = {0, 1, 2};
	mesh.Cell1DsVertices = 
	{
		Vector2i(0, 1),
		Vector2i(1, 2),
		Vector2i(2, 0)
	};
	mesh.Cell1DsMarker = {0, 0, 0};

	mesh.Cell2DsNum = 1;
	mesh.Cell2DsID = {0};
	mesh.Cell2DsVertices = {{0, 1, 2}};
	mesh.Cell2DsEdges = {{0, 1, 2}};
	
	OutputCell1Ds(mesh, path);
	
	string path1 = path + "/Output/Cell1Ds.txt";
	string path2 = path + "/Cell1Ds.txt";
	
	ifstream RealOutput(path1);
	ifstream ExpetedOutput(path2);
	
	string RealLine;
	string ExpetedLine;
	
	while (getline(ExpetedOutput, ExpetedLine) && getline(RealOutput, RealLine))
	{
		EXPECT_EQ(ExpetedLine, RealLine);
	}
}

TEST(UtilsTest, TestOutputCell2Ds)
{	
	string path = "/home/appuser/Data/ProgettoPCS2025/Progetto/Prove";
	PolyhedralMesh mesh;

	mesh.Cell0DsNum = 3;
	mesh.Cell0DsID = {0, 1, 2};
	mesh.Cell0DsCoordinates = 
	{
		Vector3d(0.0, 0.0, 0.0),
		Vector3d(1.0, 0.0, 0.0),
		Vector3d(0.0, 1.0, 0.0)
	};
	mesh.Cell0DsMarker = {0, 0, 0};

	mesh.Cell1DsNum = 3;
	mesh.Cell1DsID = {0, 1, 2};
	mesh.Cell1DsVertices = 
	{
		Vector2i(0, 1),
		Vector2i(1, 2),
		Vector2i(2, 0)
	};
	mesh.Cell1DsMarker = {0, 0, 0};

	mesh.Cell2DsNum = 1;
	mesh.Cell2DsID = {0};
	mesh.Cell2DsVertices = {{0, 1, 2}};
	mesh.Cell2DsEdges = {{0, 1, 2}};
	
	OutputCell2Ds(mesh, path);
	
	string path1 = path + "/Output/Cell2Ds.txt";
	string path2 = path + "/Cell2Ds.txt";
	
	ifstream RealOutput(path1);
	ifstream ExpetedOutput(path2);
	
	string RealLine;
	string ExpetedLine;
	
	while (getline(ExpetedOutput, ExpetedLine) && getline(RealOutput, RealLine))
	{
		EXPECT_EQ(ExpetedLine, RealLine);
	}
}

TEST(UtilsTest, TestCreateGoldberg)
{
	string path = "/home/appuser/Data/ProgettoPCS2025/Progetto/Prove";
	PolyhedralMesh mesh;
	PolyhedralMesh gold;
	PolyhedralData data;

	mesh.Cell0DsNum = 4;
	mesh.Cell0DsID = {0, 1, 2, 3};
	mesh.Cell0DsCoordinates = 
	{
		Vector3d(1.0000000000000000, 1.0000000000000000, 1.0000000000000000),
		Vector3d(-1.0000000000000000, -1.0000000000000000, 1.0000000000000000),
		Vector3d(-1.0000000000000000, 1.0000000000000000, -1.0000000000000000),
		Vector3d(1.0000000000000000, -1.0000000000000000, -1.0000000000000000)
	};
	mesh.Cell0DsMarker = {0, 0, 0, 0};

	mesh.Cell1DsNum = 6;
	mesh.Cell1DsID = {0, 1, 2, 3, 4, 5};
	mesh.Cell1DsVertices = 
	{
		Vector2i(0, 1),
		Vector2i(0, 2),
		Vector2i(0, 3),
		Vector2i(1, 2),
		Vector2i(1, 3),
		Vector2i(2, 3)
	};
	mesh.Cell1DsMarker = {0, 0, 0, 0, 0, 0};

	mesh.Cell2DsNum = 4;
	mesh.Cell2DsID = {0, 1, 2, 3};
	mesh.Cell2DsVertices = {{0, 1, 2},{0, 2, 3},{0, 3, 1},{1, 2, 3}};
	mesh.Cell2DsEdges = {{0, 3, 1},{1, 5, 2},{2, 4, 0},{3, 5, 4}};
	
	CreateGoldberg(mesh, gold);
	
	PolyhedralMesh exp;

	exp.Cell0DsNum = 4;
	exp.Cell0DsID = {0, 1, 2, 3};
	exp.Cell0DsCoordinates = 
	{
		Vector3d(-1.0000000000000000, 1.0000000000000000, 1.0000000000000000),
		Vector3d(1.0000000000000000, 1.0000000000000000, -1.0000000000000000),
		Vector3d(1.0000000000000000, -1.0000000000000000, 1.0000000000000000),
		Vector3d(-1.0000000000000000, -1.0000000000000000, -1.0000000000000000)
	};
	exp.Cell0DsMarker = {0, 0, 0, 0};

	exp.Cell1DsNum = 6;
	exp.Cell1DsID = {0, 1, 2, 3, 4, 5};
	exp.Cell1DsVertices = 
	{
		Vector2i(0, 1),
		Vector2i(0, 2),
		Vector2i(0, 3),
		Vector2i(1, 2),
		Vector2i(1, 3),
		Vector2i(2, 3)
	};
	exp.Cell1DsMarker = {0, 0, 0, 0, 0, 0};

	exp.Cell2DsNum = 4;
	exp.Cell2DsID = {0, 1, 2, 3};
	exp.Cell2DsVertices = {{0, 1, 2},{0, 2, 3},{0, 1, 3},{1, 2, 3}};
	exp.Cell2DsEdges = {{0, 1, 3},{1, 2, 5},{0, 2, 4},{3, 4, 5}};
	
	for (size_t i = 0; i < exp.Cell0DsNum; i++)
	{
		EXPECT_TRUE(gold.Cell0DsCoordinates[i].isApprox(exp.Cell0DsCoordinates[i], 1e-9));
	}
	
	EXPECT_EQ(gold.Cell0DsNum, exp.Cell0DsNum);
	EXPECT_EQ(gold.Cell0DsID, exp.Cell0DsID);
	EXPECT_EQ(gold.Cell1DsNum, exp.Cell1DsNum);
	EXPECT_EQ(gold.Cell1DsID, exp.Cell1DsID);
	EXPECT_EQ(gold.Cell1DsVertices, exp.Cell1DsVertices);
	EXPECT_EQ(gold.Cell2DsNum, exp.Cell2DsNum);
	EXPECT_EQ(gold.Cell2DsID, exp.Cell2DsID);
	EXPECT_EQ(gold.Cell2DsVertices, exp.Cell2DsVertices);
	EXPECT_EQ(gold.Cell2DsEdges, exp.Cell2DsEdges);
}

TEST(UtilsTest, TestCreateTriangulation)
{
	string path = "/home/appuser/Data/ProgettoPCS2025/Progetto/Prove";
	PolyhedralMesh mesh;

	mesh.Cell0DsNum = 3;
	mesh.Cell0DsID = {0, 1, 2};
	mesh.Cell0DsCoordinates = 
	{
		Vector3d(1.0, 1.0, 1.0),
		Vector3d(-1.0, -1.0, 1.0),
		Vector3d(-1.0, 1.0, -1.0)
	};
	mesh.Cell0DsMarker = {0, 0, 0};

	mesh.Cell1DsNum = 3;
	mesh.Cell1DsID = {0, 1, 2};
	mesh.Cell1DsVertices = 
	{
		Vector2i(0, 1),
		Vector2i(1, 2),
		Vector2i(0, 2)   
	};
	mesh.Cell1DsMarker = {0, 0, 0};

	mesh.Cell2DsNum = 1;
	mesh.Cell2DsID = {0};
	mesh.Cell2DsVertices = {{0, 1, 2}};
	mesh.Cell2DsEdges = {{0, 1, 2}};
	
	PolyhedralData data;
	
	data.BestPath = false;
	data.Goldby = false;
	data.b = 2;
	data.c = 0;
	data.section = 0;
	
	PolyhedralMesh trg;
	
	CreateTriangulation(mesh, data, trg);
	
	PolyhedralMesh exp;

	exp.Cell0DsNum = 6;
	exp.Cell0DsID = {0, 1, 2, 3, 4, 5};
	exp.Cell0DsCoordinates = 
	{
		Vector3d(-1.0, 1.0, -1.0),
		Vector3d(-1.0, 0.0, 0.0),
		Vector3d(-1.0, -1.0, 1.0),
		Vector3d(0.0, 1.0, 0.0),
		Vector3d(0.0, 0.0, 1.0),
		Vector3d(1.0, 1.0, 1.0)
	};
	exp.Cell0DsMarker = {0, 0, 0, 0, 0, 0};

	exp.Cell1DsNum = 9;
	exp.Cell1DsID = {0, 1, 2, 3, 4, 5, 6, 7, 8};
	exp.Cell1DsVertices = 
	{
		Vector2i(0, 1),
		Vector2i(0, 3),
		Vector2i(1, 2),
		Vector2i(1, 3),
		Vector2i(1, 4),
		Vector2i(2, 4),
		Vector2i(3, 4),
		Vector2i(3, 5),
		Vector2i(4, 5)
	};
	exp.Cell1DsMarker = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	
	for (size_t i = 0; i < exp.Cell0DsNum; i++)
	{
		EXPECT_TRUE(trg.Cell0DsCoordinates[i].isApprox(exp.Cell0DsCoordinates[i], 1e-9));
	}
	
	EXPECT_EQ(trg.Cell0DsNum, exp.Cell0DsNum);
	EXPECT_EQ(trg.Cell0DsID, exp.Cell0DsID);
	EXPECT_EQ(trg.Cell1DsNum, exp.Cell1DsNum);
	EXPECT_EQ(trg.Cell1DsID, exp.Cell1DsID);
	EXPECT_EQ(trg.Cell1DsVertices, exp.Cell1DsVertices);
}

TEST(UtilsTest, TestCreateAdjacencyMatrix) 
{
	PolyhedralMesh mesh;

	mesh.Cell0DsNum = 3;
	mesh.Cell0DsID = {1, 2, 3};
	mesh.Cell0DsCoordinates = 
	{
		Vector3d(0.0, 0.0, 0.0),
		Vector3d(1.0, 0.0, 0.0),
		Vector3d(0.0, 1.0, 0.0)
	};
	mesh.Cell0DsMarker = {0, 0, 0};

	mesh.Cell1DsNum = 3;
	mesh.Cell1DsID = {1, 2, 3};
	mesh.Cell1DsVertices = 
	{
		Vector2i(0, 1),
		Vector2i(1, 2),
		Vector2i(2, 0)
	};
	mesh.Cell1DsMarker = {0, 0, 0};

	mesh.Cell2DsNum = 1;
	mesh.Cell2DsID = {1};
	mesh.Cell2DsVertices = {{0, 1, 2}};
	mesh.Cell2DsEdges = {{0, 1, 2}};

	MatrixXi adj = CreateAdjacencyMatrix(mesh);

	ASSERT_EQ(adj.rows(), 3);
	ASSERT_EQ(adj.cols(), 3);
	EXPECT_EQ(adj(0, 1), 1);
	EXPECT_EQ(adj(1, 2), 1);
	EXPECT_EQ(adj(2, 0), 1);
	EXPECT_EQ(adj(0, 2), 1);
	EXPECT_EQ(adj(1, 0), 1);
	EXPECT_EQ(adj(2, 1), 1);
}