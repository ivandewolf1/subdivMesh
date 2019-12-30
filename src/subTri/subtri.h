
#include "tvec3.h"
#include <unordered_set>
#include <algorithm>
#include <vector>

#ifndef __subtri_h__
#define __subtri_h__


namespace subTri {
using namespace std;
using namespace lwtv;
 
template <typename T>
struct triangleCenters
{
  tvec3<T> circumCenter;
  T circumRadius;
  tvec3<T> inCenter;
  T inRadius;
};

template <typename T>
class smPoint
{
 public:
  tvec3<T> loc;
  tvec3<T> N;
  //float u;
  //float v;
  vector <unsigned int> trisIndxs;
  void popIndex(unsigned int index){trisIndxs.erase(std::remove(trisIndxs.begin(), trisIndxs.end(), index),trisIndxs .end());}
};

class listIndex
{
 public:
   listIndex(){exists = false;}
   unsigned int get(){return idx;}
   listIndex& operator = (const unsigned int &val){idx=val; return *this;}
   bool is(){return exists;}
   void setIsnt(){exists = false;}
   void set(unsigned int idxIn){idx = idxIn; exists = true;}
   void set(listIndex liIn){idx = liIn.get(); exists = liIn.is();}
 private:
   unsigned int idx;
   bool exists;
};

 enum splitStates{
   NEW_TRI,
   SPLIT,
   NOT_SPLIT,
   FINISHED
 };
 
 enum flipStates{
   UNKNOWN,
   IS_FLIPPED,
   NOT_FLIPPED
 };
 
template <typename T>
class subtri
{
public:
  subtri();
  ~subtri(){};
  void insertVertNeighbors(vector<subtri> &subtriList, vector<smPoint<T>> &pointList, unordered_set<unsigned int> &splitSet);
  virtual bool testSplit(vector<smPoint<T>> &pointList, double radius, double resmult);
  tvec3<T> centerLoc(vector<smPoint<T>> &pointList);
  double maxRad(vector<smPoint<T>> &pointList);
  void setup(listIndex vert0,listIndex vert1,listIndex vert2, listIndex N0, listIndex N1, listIndex N2, unsigned int myIndex);
  void buildChildren(vector<subtri> &subtriList, vector<smPoint<T>> &pointList);
  void testFlip(vector<subtri> &subtriList, vector<smPoint<T>> &pointList);
  void connectSubtri(unsigned int  parentIndex, int edgeNum, vector<subtri> &subtriList, vector<smPoint<T>> &pointList);

  // utility functions
  triangleCenters<T> computeCenters(unsigned int indexA, unsigned int indexB, unsigned int indexC, vector<smPoint<T>> &pointList);
  
  // accessor functions
  int findIndex(unsigned int meshIndex);
  void setChildren(unsigned int addrIn, unsigned int A, unsigned int B, unsigned int C);
  void setIndex(unsigned int indexIn){index = indexIn;}
  unsigned int vertAddr(int i){return verts[i].get();}

  //data 
  listIndex neighbors[3];
  unsigned int index;
  splitStates splitState;
  listIndex children[3];
  listIndex cenAddr;
  tvec3<T> N;
  tvec3<T> cen;
  
  //private:
  listIndex verts[3];
  flipStates flipped[3];
};

template <typename T>
class subdividedMesh
{
 public:
  subdividedMesh(){}
  virtual ~subdividedMesh(){}
  virtual void exampleSubdivide(int max, double targetRadius, double sphRad, double resmult);
  virtual void exampleTriTest(double targetRadius, double resmult);
  virtual void exampleMovePoints(double radius);
  
  void makePoints();
  void buildChildGroup();
  void prepChildren();
  tvec3<T> incenter(unsigned int indexA, unsigned int indexB, unsigned int indexC, T &radius, tvec3<T> &N);
  tvec3<T> circumcenter(unsigned int indexA, unsigned int indexB, unsigned int indexC, T &radius, tvec3<T> &N);
  virtual void computeNormals();
  void initMesh(double radius);
  void initTet(double radius);

  vector<smPoint<T>> pointList;
  vector<subtri<T>> subtriList;
 private:
  unordered_set<unsigned int> splitSet;
  vector<subtri<T>*> childGroup;

  vector<unsigned int> newPointsGroup;
 };
 
} // end idwTools namespace


#endif
