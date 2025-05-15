#pragma once

#include "PolyhedralMesh.hpp"
#include "UCDUtilities.hpp"
#include <limits>

using namespace std;

namespace PolyhedralLibrary
{
bool ImportCell(const string& path, PolyhedralMesh& mesh);

bool CheckLength(PolyhedralMesh& mesh);

bool CheckAreas(PolyhedralMesh& mesh);

bool ExpPoints(PolyhedralMesh& mesh, const string& FilePath);

bool ExpSegments(PolyhedralMesh& mesh, const string& FilePath);
}