#pragma once

#include "PolyhedrallMesh.hpp"
#include "UCDUtilities.hpp"
#include <limits>

using namespace std;

namespace PolyhedralLibrary
{
bool ImportCell(const string& path, PolygonalMesh& mesh);

bool CheckLength(PolygonalMesh& mesh);

bool CheckAreas(PolygonalMesh& mesh);

bool ExpPoints(PolygonalMesh& mesh, const string& FilePath);

bool ExpSegments(PolygonalMesh& mesh, const string& FilePath);
}