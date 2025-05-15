#include "Utils.hpp"
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

double EdgeLength(const double& x1, const double& y1, const double& z1, const double& x2, const double& y2, const double& z2)
{
	return sqrt(pow(x2 - x1 , 2) + pow(y2 - y1 , 2) + pow(z2 - z1 , 2));	
}

namespace PolyhedralLibrary
{
bool ImportCell(const string& path, PolyhedralMesh& mesh)
{
	return 0;
}
}