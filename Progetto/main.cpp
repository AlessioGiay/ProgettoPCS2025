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
	bool Goldby = false;

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

	if(data.b != data.c && (data.b == 0 || data.c == 0) )
	{
		Geodetico(mesh, path);

		Triangulation(mesh, data, trg, path, Goldby);

		if(data.q == 3)
		{
			Goldby = true;
			trg = PolyhedralMesh();
		
			if(!Goldberg(mesh, gold, path))
			{
				return 1;
			}
			else
			{
				Triangulation(gold, data, trg, path, Goldby);
			}
		}
	}

	cout << "Fatto tutto\n";

	return 0;
}

// Problema sull'orientamento dei lati