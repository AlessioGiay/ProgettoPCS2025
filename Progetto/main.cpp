#include "Utils.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;

int main()
{
	PolyhedralMesh mesh;
	string path = "/home/appuser/Data/ProgettoPCS2025/Progetto";
	string File_0D_Path = "/home/appuser/Data/ProgettoPCS2025/Progetto/Cell0Ds.inp";
	string File_1D_Path = "/home/appuser/Data/ProgettoPCS2025/Progetto/Cell1Ds.inp";

	if(!ImportVector(path, mesh))
	{
		return 1;
	}
/*
	if(!ExpPoints(mesh, File_0D_Path))
	{
		return 1;
	}
	
	if(!ExpSegments(mesh, File_1D_Path))
	{
		return 1;
	}
*/
	cout << "Fatto tutto\n";

	return 0;
}