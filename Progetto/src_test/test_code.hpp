#include "gtest/gtest.h"
#include "PolyhedralMesh.hpp"
#include "Utils.hpp"

using namespace PolyhedralLibrary;
using namespace Eigen;

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
	
	ifstream RealOutput(path + "/Output/Cell0Ds.txt");
	ifstream ExpetedOutput(path + "/Cell0Ds.txt");
	
	string RealLine;
	string ExpetedLine;
	
	while (getline(ExpetedOutput, ExpetedLine) && getline(RealOutput, RealLine))
	{
		EXPECT_EQ(ExpetedLine, RealLine);
	}
}