#pragma once

#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;


namespace PolyhedralLibrary {

struct PolyhedralMesh
{
	unsigned int Cell0DsNum = 0;
	vector<unsigned int> Cell0DsID = {};
	vector<unsigned int> Cell0DsMarker = {};
	vector<Vector3d> Cell0DsCoordinates = {};

	unsigned int Cell1DsNum = 0;
	vector<unsigned int> Cell1DsID {};
	vector<unsigned int> Cell1DsMarker = {};
	vector<Vector2i> Cell1DsVertices = {};

	unsigned int Cell2DsNum = 0;
	vector<unsigned int> Cell2DsID = {};
	vector<vector<unsigned int>> Cell2DsVertices = {};
	vector<vector<unsigned int>> Cell2DsEdges = {};
	vector<Vector3d> Cell2DsCentre = {};
	vector<Vector2i> Cell2DsArches = {};
/*
	unsigned int Cell3DsNum = 0;
	vector<unsigned int> Cell3DsID = {};
	vector<vector<unsigned int>> Cell3DsVertices = {};
	vector<vector<unsigned int>> Cell3DsEdges = {};
	vector<vector<unsigned int>> Cell3DsFaces = {};
*/
	const double epsilon = numeric_limits<double>::epsilon();
};

}