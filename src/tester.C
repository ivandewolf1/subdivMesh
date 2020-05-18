#include <iostream>
#include <fstream>
#include <subtri.C>

using namespace std;
using namespace lwtv;

template <typename T>
struct argsStruct
{
	T triradius;
	T sphRad;
	T resmult;
	T maxRatio; 
};

template <typename T>
class testMesh : public subdividedMesh<T>
{
using subdividedMesh<T>::pointList;
using subdividedMesh<T>::subtriList;
using subdividedMesh<T>::newPointsGroup;
bool testSplit(subtri<T> *tri, void *args)  // perhaps pass everything by const &  ?
{
	subtri<T> spit(*tri);
	argsStruct<T> *myArgs = (argsStruct<T>*)args;
	tvec3<T> cenLoc = tri->centroidLoc(pointList);  // cenLoc could be const but operator > does not handle it
 	float mult = 1;
	if(cenLoc[0]>0)mult = myArgs->resmult;
	const T triRadius = tri->maxRad(pointList);
	if(triRadius > (myArgs->triradius*mult)){
  	return true;
	}
	return false;
}

void movePoints(void *args)
{
	argsStruct<T> *myArgs = (argsStruct<T>*)args;
  vector<unsigned int>::iterator newPtIdx;
  for (newPtIdx = newPointsGroup.begin(); newPtIdx != newPointsGroup.end(); ++newPtIdx){
    pointList[*newPtIdx].loc.normalize();
    pointList[*newPtIdx].loc *= myArgs->sphRad;
  }
}
};

main()
{
	int maxIter = 10;
	argsStruct<float> myArgs;
	myArgs.triradius = 4;
	myArgs.sphRad = 1.2;
	myArgs.resmult = .01;
	myArgs.maxRatio = 10;
	testMesh<float> sphere;
	tvec3<float> scale(1.0);
	sphere.initializeTet(myArgs.sphRad, scale);
	sphere.subdivide(maxIter, 10.0, &myArgs);
	/*return 1;
	float triradius = 1.1;
	float sphRad = 1.2;
	float resmult = .01;
	float maxRatio = 10;
	sphere.initTet(sphRad);
	sphere.exampleSubdivide(maxIter, triradius, sphRad, resmult, maxRatio);*/
	ofstream objFile;
	objFile.open("output.obj");
	objFile<<"# file exported by graded subdivision triangle tester"<<endl;
	objFile<<"# "<<sphere.getPointList()->size()<<" points\n";
	objFile<<"# "<<(sphere.getSubtriList()->size()*3)<<" vertices\n";
	objFile<<"# "<<sphere.getSubtriList()->size()<<" primitives\n";
	objFile<<"# Bounds: [-1, -1, -1] to [1, 1, 1]"<<endl;
	objFile<<"g\n";
	vector<smPoint<float>>::iterator sphPt;
	for (sphPt = sphere.getPointList()->begin(); sphPt != sphere.getPointList()->end(); ++sphPt){
		objFile<<"v "<<sphPt->loc[0]<<" "<<sphPt->loc[1]<<" "<<sphPt->loc[2]<<"\n";
	}	
	objFile<<"g\n";
	vector<subtri<float>>::iterator tri;
	for (tri = sphere.getSubtriList()->begin(); tri != sphere.getSubtriList()->end(); ++tri){
	  if(tri->splitState != FINISHED){
	  	objFile<<"f "<<tri->verts[0].get()+1<<" "<<tri->verts[1].get()+1<<" "<<tri->verts[2].get()+1<<"\n";
	  }
	}
	objFile<<endl;
	objFile.close();
}
