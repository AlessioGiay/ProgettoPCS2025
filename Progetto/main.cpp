#include "Utils.hpp"
#include "PolyhedralMesh.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;

int main()
{
	PolyhedralMesh mesh;
	PolyhedralMesh gold;
	PolyhedralMesh trg;
	PolyhedralData data;

	string path = "/home/appuser/Data/ProgettoPCS2025/Progetto";
	string File_0D_Path = "/home/appuser/Data/ProgettoPCS2025/Progetto/Cell0Ds.inp";
	string File_1D_Path = "/home/appuser/Data/ProgettoPCS2025/Progetto/Cell1Ds.inp";

	if(!ImportVector(path, mesh, data))
	{
		return 1;
	}

	if(!Output(mesh, path))
	{
		return 1;
	}

//	Triangulation(mesh, data, trg, path);
	if(data.q == 3)
	{
		Goldberg(mesh, path, gold);
	}

	cout << "Fatto tutto\n";

	return 0;
}

/* if(!Geodetico(...))
	{
		return 1;
	}
	
	if(data.q == 3)
	{
		if(!Goldberg(mesh, path, gold))
		{
			return 1;
		}
	}
*/