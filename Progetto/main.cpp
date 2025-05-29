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

	if(!ImportVector(path, mesh, data))
	{
		return 1;
	}

	if(data.b != data.c && (data.b == 0 || data.c == 0) )
	{
		if(!Output(mesh, path, data))
		{
			return 1;
		}
		
		Geodetico(mesh, path);

		CreateTriangulation(mesh, data, trg);
		
		Triangulation(trg, path, data);

		if(data.q == 3)
		{
			data.Goldby = true;
			trg = PolyhedralMesh();
		
			if(!CreateGoldberg(mesh, gold))
			{
				return 1;
			}
			else
			{
				if(!Output(gold, path, data))
				{
					return 1;
				}
				
				Goldberg(gold, path);
				
				CreateTriangulation(gold, data, trg);
				
				Triangulation(trg, path, data);
			}
		}
	}

	return 0;
}