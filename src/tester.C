#include <iostream>
#include <fstream>
#include <subtri.C>

using namespace std;
using namespace lwtv;

main()
{
	int maxIter = 10;
	float triradius = 1;
	float sphRad = 1;
	float resmult = .01;
	subdividedMesh<float> sphere;
	sphere.initTet(sphRad);
	sphere.exampleSubdivide(maxIter, triradius, sphRad, resmult);
	ofstream objFile;
	objFile.open("output.obj");
	objFile<<"# file exported by graded subdivision triangle tester"<<endl;
	objFile<<"# "<<sphere.pointList.size()<<" points\n";
	objFile<<"# "<<(sphere.subtriList.size()*3)<<" vertices\n";
	objFile<<"# "<<sphere.subtriList.size()<<" primitives\n";
	objFile<<"# Bounds: [-1, -1, -1] to [1, 1, 1]"<<endl;
	objFile<<"g\n";
	vector<smPoint<float>>::iterator sphPt;
	for (sphPt = sphere.pointList.begin(); sphPt != sphere.pointList.end(); ++sphPt){
		objFile<<"v "<<sphPt->loc[0]<<" "<<sphPt->loc[1]<<" "<<sphPt->loc[2]<<"\n";
	}	
	objFile<<"g\n";
	vector<subtri<float>>::iterator tri;
	for (tri = sphere.subtriList.begin(); tri != sphere.subtriList.end(); ++tri){
	  if(tri->splitState != FINISHED){
	  	objFile<<"f "<<tri->verts[0].get()+1<<" "<<tri->verts[1].get()+1<<" "<<tri->verts[2].get()+1<<"\n";
	  }
	}
	objFile<<endl;
	objFile.close();
}