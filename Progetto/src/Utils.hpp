#pragma once

#include "PolyhedralMesh.hpp"
#include "UCDUtilities.hpp"
#include <limits>

using namespace std;

namespace PolyhedralLibrary
{
bool ImportVector(const string& path, PolyhedralMesh& mesh);

bool PolyhedralChoice(const string& path, PolyhedralMesh& mesh, const char& p, const char& q,  const char& b, const char& c, bool& BestPath);

bool ImportMesh(const string& path, PolyhedralMesh& mesh);

bool ImportCell0Ds(const string& path, PolyhedralMesh& mesh);

bool ImportCell1Ds(const string& path, PolyhedralMesh& mesh);

bool ImportCell2Ds(const string& path, PolyhedralMesh& mesh);
/*
bool ExpPoints(PolyhedralMesh& mesh, const string& FilePath);

bool ExpSegments(PolyhedralMesh& mesh, const string& FilePath);
*/
}