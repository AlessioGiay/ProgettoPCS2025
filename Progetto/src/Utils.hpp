#pragma once

#include "PolyhedralMesh.hpp"
#include "UCDUtilities.hpp"

using namespace std;

namespace PolyhedralLibrary
{

bool ImportVector(const string& path, PolyhedralMesh& mesh, PolyhedralData& data);

bool PolyhedralChoice(const string& path, PolyhedralMesh& mesh, PolyhedralData& data);

bool ImportMesh(const string& path, PolyhedralMesh& mesh);

bool ImportCell0Ds(const string& path, PolyhedralMesh& mesh);

bool ImportCell1Ds(const string& path, PolyhedralMesh& mesh);

bool ImportCell2Ds(const string& path, PolyhedralMesh& mesh);

bool Output(PolyhedralMesh& mesh, const string& path);

bool OutputCell0Ds(PolyhedralMesh& mesh, const string& path);

bool OutputCell1Ds(PolyhedralMesh& mesh, const string& path);

bool OutputCell2Ds(PolyhedralMesh& mesh, const string& path);

bool OutputCell3Ds(PolyhedralMesh& mesh, const string& path);

bool Goldberg(PolyhedralMesh& mesh, PolyhedralMesh& gold,  const string& path);

bool CreateFaces(PolyhedralMesh& mesh, PolyhedralMesh& gold);

bool CheckEdges(const PolyhedralMesh& mesh);

bool ShortPath(PolyhedralMesh& trg, PolyhedralData& data);

void Triangulation(const PolyhedralMesh& mesh, PolyhedralData& data, PolyhedralMesh& trg, const string& path,  const bool Goldby);

void Geodetico(PolyhedralMesh& mesh, const string& path);

MatrixXi CreateAdjacencyMatrix(PolyhedralMesh& trg);

vector<Vector3d> TrgCleaning(PolyhedralMesh& trg);

vector<Vector2i> CreateArches(const PolyhedralMesh& mesh, PolyhedralMesh& trg, const PolyhedralData& data, vector<vector<Vector3d>>& Alpha);

}