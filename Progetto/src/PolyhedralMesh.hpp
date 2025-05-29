#pragma once

#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary 
{

struct PolyhedralMesh
{
	unsigned int Cell0DsNum = 0;
	vector<unsigned int> Cell0DsID = {};
	vector<Vector3d> Cell0DsCoordinates = {};
	vector<unsigned int> Cell0DsMarker = {};

	unsigned int Cell1DsNum = 0;
	vector<unsigned int> Cell1DsID {};
	vector<Vector2i> Cell1DsVertices = {};
	vector<unsigned int> Cell1DsMarker = {};

	unsigned int Cell2DsNum = 0;
	vector<unsigned int> Cell2DsID = {};
	vector<vector<unsigned int>> Cell2DsVertices = {};
	vector<vector<unsigned int>> Cell2DsEdges = {};
};

struct PolyhedralData
{
	int p = 0;
	int q = 0;
	int b = 0;
	int c = 0;
	unsigned int Id1 = 0;
	unsigned int Id2 = 0;
	unsigned int section = 0;
	bool BestPath = false;
	bool Goldby = false;
};

}


/*
	unsigned int Cell3DsNum = 0;
	vector<unsigned int> Cell3DsID = {};
	vector<vector<unsigned int>> Cell3DsVertices = {};
	vector<vector<unsigned int>> Cell3DsEdges = {};
	vector<vector<unsigned int>> Cell3DsFaces = {};
*/