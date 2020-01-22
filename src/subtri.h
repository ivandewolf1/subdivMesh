
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
  tvec3<T> N;
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
   listIndex& operator = (const unsigned int &val){idx=val; return *this;}
   unsigned int get(){return idx;}
   bool is(){return exists;}
   void set(unsigned int idxIn){idx = idxIn; exists = true;}
   void set(listIndex li){idx = li.get(); exists = li.is();}
 private:
   unsigned int idx;
   bool exists;
};

 enum splitStates{
   NEW_TRI,
   SPLIT,
   FIX_SPLIT,
   NOT_SPLIT,
   FINISHED
 };
 
 enum flipStates{
   UNKNOWN,
   IS_FLIPPED,
   NOT_FLIPPED,
   EDGE
 };
 
template <typename T>
class subtri
{
public:
  subtri();
  ~subtri(){};
  tvec3<T> centroidLoc(vector<smPoint<T>> &pointList);
  double maxRad(vector<smPoint<T>> &pointList);
  void setup(listIndex vert0,listIndex vert1,listIndex vert2, listIndex N0, listIndex N1, listIndex N2, unsigned int myIndex);
  void buildChildren(vector<subtri> &subtriList, vector<smPoint<T>> &pointList);
  void handleBadRatio(vector<subtri> &subtriList, vector<smPoint<T>> &pointList);
  int flipTester(vector<subtri> &subtriList, vector<smPoint<T>> &pointList, triangleCenters<T> centers[4], int i);
  void connectSubtri(unsigned int  parentIndex, int edgeNum, vector<subtri> &subtriList, vector<smPoint<T>> &pointList);
  void badRatioFlip(vector<subtri> &subtriList, vector<smPoint<T>> &pointList, vector<unsigned int> &splitList);
  void unsplitSplitNeighbors(vector<subtri> &subtriList, vector<smPoint<T>> &pointList, vector<unsigned int> &splitList);
  void splitFlip(vector<subtri> &subtriList, vector<smPoint<T>> &pointList);

  // utility functions
  triangleCenters<T> computeCenters(tvec3<T> A, tvec3<T> B, tvec3<T> C);
  // accessor functions
  int findIndex(unsigned int meshIndex);
  void setChildren(unsigned int addrIn, unsigned int A, unsigned int B, unsigned int C);
  void setIndex(unsigned int indexIn){index = indexIn;}
  unsigned int vertAddr(int i){return verts[i].get();}

  //data 
  bool tapped;
  listIndex neighbors[3];
  unsigned int index;
  unsigned int parentIndex;
  splitStates splitState;
  listIndex children[3];
  listIndex cenAddr;
  triangleCenters<T> cen;
  triangleCenters<T> childCen[3];
 
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
  //---------------------------------------------------------------------
  void subdivide(int maxIter, T maxRatio, void *args);
  void splitSort(void *args);
  void flipAndFix(T maxRatio);
  void createGeo();
  void sew();
  virtual bool testSplit(subtri<T> *tri, void *args) = 0;
  virtual void movePoints(void *args) = 0;
  vector<subtri<T>*> newTriGroup;
  //---------------------------------------------------------------------
  
  tvec3<T> incenter(unsigned int indexA, unsigned int indexB, unsigned int indexC, T &radius, tvec3<T> &N);
  tvec3<T> circumcenter(unsigned int indexA, unsigned int indexB, unsigned int indexC, T &radius, tvec3<T> &N);
  void computeNormals();
  void initializeTet(T radius,  tvec3<T> scale);
  vector<smPoint<T>>* getPointList(){return &pointList;}
  vector<subtri<T>>* getSubtriList(){return &subtriList;}

  vector<smPoint<T>> pointList;
  vector<subtri<T>> subtriList;
 //private:
  vector<unsigned int> splitList;
  vector<unsigned int> unSplitList;
  vector<subtri<T>*> childGroup;

  vector<unsigned int> newPointsGroup;
 };
 
} // end subTri namespace


#endif
