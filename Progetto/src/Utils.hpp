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

void Triangulation(const PolyhedralMesh& mesh, PolyhedralData& data, PolyhedralMesh& trg, const string& path);

vector<Vector3d> TrgCleaning(PolyhedralMesh& trg);

vector<Vector2i> CreateArches(const PolyhedralMesh& mesh, PolyhedralMesh& trg, const PolyhedralData& data);
}

/*
void Goldberg(PolyhedralMesh& mesh, const string& path, PolyhedralMesh& gold);
*/