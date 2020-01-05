#include "subtri.h"
using namespace subTri;

// to do:
// switch from mean centers to incenters for all computation
// run timing test in subdivide()
// investigate threading
// refactor connectSubtri()
// *** handle edges and holes ***

template <typename T>
void
subdividedMesh<T>::exampleSubdivide(int max, double targetRadius, double sphRad, double resmult)
{
  srand(3.1415926535);//why not
  for(int i = 0; i<max; ++i){
    exampleTriTest(targetRadius, resmult);
    if(splitSet.size() == 0)return;
    prepChildren();// all centers needed to be set by here
    makePoints();
    exampleMovePoints(sphRad);
    buildChildGroup();
    splitSet.clear();// all tri-N needs to be built by this point
    //computeNormals()
  }
}

template <typename T>
void
subdividedMesh<T>::exampleTriTest(double targetRadius, double resmult)
{
  typename vector<subtri<T>*>::iterator tri;
  for (tri = childGroup.begin(); tri != childGroup.end(); ++tri){
    bool split = (*tri)->testSplit(pointList, targetRadius, resmult);
    if(split){
      (*tri)->insertVertNeighbors(subtriList, pointList, splitSet);
    }else{
      if((*tri)->splitState == NEW_TRI) (*tri)->splitState = NOT_SPLIT;
    }
  }
}

// this puts points onto a sphere 
template <typename T>
void
subdividedMesh<T>::exampleMovePoints(double radius)// loops over newPoints
{
  vector<unsigned int>::iterator newPtIdx;
  for (newPtIdx = newPointsGroup.begin(); newPtIdx != newPointsGroup.end(); ++newPtIdx){
  pointList[*newPtIdx].loc.normalize();
  pointList[*newPtIdx].loc *= radius;
   }
}

template <typename T>
void
subdividedMesh<T>::triTest(void *testControls)
{
  typename vector<subtri<T>*>::iterator tri;
  for (tri = childGroup.begin(); tri != childGroup.end(); ++tri){
    bool split = (*tri)->testSplit(pointList, testControls);
    if(split){
      (*tri)->insertVertNeighbors(subtriList, pointList, splitSet);
    }else{
      if((*tri)->splitState == NEW_TRI) (*tri)->splitState = NOT_SPLIT;
    }
  }
}

// this *requires* that all tris have set their N before this is called
// this should be trivial to thread
template <typename T>
void
subdividedMesh<T>::computeNormals()
{
  typename vector<smPoint<T>>::iterator crrPt;
  vector<unsigned int>::iterator crrTri;
  for (crrPt = pointList.begin(); crrPt != pointList.end(); ++crrPt){
    crrPt->N = 0;
    for (crrTri = crrPt->trisIndxs.begin(); crrTri != crrPt->trisIndxs.end(); ++crrTri){
      crrPt->N += subtriList[*crrTri].N;
    }
    crrPt->N.normalize();
  }
}

template <typename T>
void
subdividedMesh<T>::makePoints()
{
  unsigned int crrChildCount = subtriList.size() - (splitSet.size() * 3);
  unsigned int Pcrr = pointList.size();
  pointList.resize(pointList.size() +  splitSet.size());
  newPointsGroup.resize(0);
  unordered_set<unsigned int>::iterator crrSplit;
  for (crrSplit = splitSet.begin(); crrSplit != splitSet.end(); ++crrSplit){
    pointList[Pcrr].loc = subtriList[*crrSplit].centerLoc(pointList);
    newPointsGroup.push_back(Pcrr);
    subtriList[*crrSplit].setChildren(Pcrr, crrChildCount, crrChildCount + 1, crrChildCount + 2);
    crrChildCount += 3;
    Pcrr += 1;
  }
}

// simply loop over children, calling buildChildren
template <typename T>
void
subdividedMesh<T>::buildChildGroup()// loops over splitSet
{
  unordered_set<unsigned int>::iterator splitTri;
  for (splitTri = splitSet.begin(); splitTri != splitSet.end(); ++splitTri){
    subtriList[*splitTri].buildChildren(subtriList, pointList);
  }
}

// create childGroup, loop over childGroup, fill with initial index value
template <typename T>
void
subdividedMesh<T>::prepChildren()// loops over childGroup
{
  unsigned int crrChildCount = subtriList.size(); //assuming no children have been made yet
  childGroup.clear();
  childGroup.resize(splitSet.size() * 3);
  subtriList.resize(subtriList.size() +  (splitSet.size()*3));
  typename vector<subtri<T>*>::iterator crrChild;
  for (crrChild = childGroup.begin(); crrChild != childGroup.end(); ++crrChild){
    *crrChild = &(subtriList[crrChildCount]);
    (*crrChild)->setIndex(crrChildCount);
    crrChildCount++;
  }
  
}

template <typename T>
void
subdividedMesh<T>::initTet(double radius)
{
  tvec3<T> vert0(0.0,1.0,0.0);
  tvec3<T> vert1(-.94281,-1.0/3.0,0.0);
  tvec3<T> vert2(.471405,-1.0/3.0,.816498*.1);
  tvec3<T> vert3(.471405,-1.0/3.0,-.816498*.1);
  vert0.normalize();
  vert1.normalize();
  vert2.normalize();
  vert3.normalize();
  pointList.resize(4);
  pointList[0].loc = (vert0 * radius);
  pointList[1].loc = (vert1 * radius);
  pointList[2].loc = (vert2 * radius);
  pointList[3].loc = (vert3 * radius);
  subtriList.resize(4);
  listIndex indexer0, indexer1, indexer2, indexer3;
  indexer0.set(0);
  indexer1.set(1);
  indexer2.set(2);
  indexer3.set(3);
  subtri<T> cenSource;
  subtriList[0].setup(indexer0,indexer3,indexer2,indexer1,indexer3,indexer2, 0);
  subtriList[0].cen = cenSource.computeCenters(pointList[0].loc,pointList[3].loc,pointList[2].loc);
  subtriList[1].setup(indexer0,indexer1,indexer3,indexer2,indexer3,indexer0, 1);
  subtriList[1].cen = cenSource.computeCenters(pointList[0].loc,pointList[1].loc,pointList[3].loc);
  subtriList[2].setup(indexer0,indexer2,indexer1,indexer0,indexer3,indexer1, 2);
  subtriList[2].cen = cenSource.computeCenters(pointList[0].loc,pointList[2].loc,pointList[1].loc);
  subtriList[3].setup(indexer1,indexer2,indexer3,indexer2,indexer0,indexer1, 3);
  subtriList[3].cen = cenSource.computeCenters(pointList[1].loc,pointList[2].loc,pointList[3].loc);
  pointList[0].trisIndxs.resize(3);
  pointList[0].trisIndxs[0] = 0;
  pointList[0].trisIndxs[1] = 1;
  pointList[0].trisIndxs[2] = 2;
  pointList[1].trisIndxs.resize(3);
  pointList[1].trisIndxs[0] = 1;
  pointList[1].trisIndxs[1] = 2;
  pointList[1].trisIndxs[2] = 3;
  pointList[2].trisIndxs.resize(3);
  pointList[2].trisIndxs[0] = 0;
  pointList[2].trisIndxs[1] = 2;
  pointList[2].trisIndxs[2] = 3;
  pointList[3].trisIndxs.resize(3);
  pointList[3].trisIndxs[0] = 0;
  pointList[3].trisIndxs[1] = 1;
  pointList[3].trisIndxs[2] = 3;
  // fill childGroup
  childGroup.push_back(&(subtriList[0]));
  childGroup.push_back(&(subtriList[1]));
  childGroup.push_back(&(subtriList[2]));
  childGroup.push_back(&(subtriList[3]));
}

//-----------------------------------------------------------------------------------------------------------------------

template <typename T>
subtri<T>::subtri()
{
  neighbors[0] = neighbors[1] = neighbors[2] = 0;
  index = 0;
  splitState = NEW_TRI;
  children[0] = children[1] = children[2] = 0;
  cenAddr = 0;
  verts[0] = verts[1] = verts[2] = 0;
  flipped[0] = flipped[1] = flipped[2] = UNKNOWN;
}

template <typename T>
void
subtri<T>::insertVertNeighbors(vector<subtri> &subtriList, vector<smPoint<T>> &pointList, unordered_set<unsigned int> &splitSet)
{
  for(int i = 0; i < 3; ++i){
    smPoint<T> *crrPt = &pointList[verts[i].get()];
    vector<unsigned int>::iterator indexIter;
    for (indexIter = crrPt->trisIndxs.begin(); indexIter != crrPt->trisIndxs.end(); ++indexIter){
      subtriList[*indexIter].splitState = SPLIT;
      splitSet.insert(*indexIter);
    }
  }
}


// this placeholder tests the maxradius is greater than some value.
template <typename T>
bool
subtri<T>::testSplit(vector<smPoint<T>> &pointList, double targetRadius, double resmult)
{
  
  tvec3<T> cenLoc = centerLoc(pointList);
  float mult = 1;
  if(cenLoc[0]>0)mult = resmult;
  double triRadius = maxRad(pointList);
  if(triRadius > (targetRadius*mult)){
    splitState = SPLIT;
    return true;
  }
  return false;
}

// vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv UNUSED!!!!!!!!!!!!!!!!
// this placeholder tests the maxradius is greater than some value.
template <typename T>
bool
subtri<T>::testSplit(vector<smPoint<T>> &pointList, void *testControls)
{
  float *controls = (float*)testControls;
  cout <<"control zero"<<controls[0]<<endl;
  /*tvec3<T> cenLoc = centerLoc(pointList);
  float mult = 1;
  if(cenLoc[0]>0)mult = *((float*)testData);
  double triRadius = maxRad(pointList);
  if(triRadius > mult){
    splitState = SPLIT;
    return true;
  }*/
  return false;
}

template <typename T>
tvec3<T>
subtri<T>::centerLoc(vector<smPoint<T>> &pointList)
{
  //return cen.inCenter;
  tvec3<T>  avg(0.0, 0.0, 0.0);
  for(int i=0; i<3; ++i){
    avg += pointList[verts[i].get()].loc;
  }
  avg /= 3;
  return avg;
  //cout <<"centers:\n"<<cen.inCenter<<"\n"<<(avg/3)<<" "<<index<<endl;
  //T randVal = (T)rand() / (T)RAND_MAX;
  //cout <<"rv: "<<randVal<<endl;
  //return avg.lerp(cen.inCenter, .5);
}


template <typename T>
double
subtri<T>::maxRad(vector<smPoint<T>> &pointList)
{
  tvec3<T> cenLoc = centerLoc(pointList);
  double maxRad = 0;
  for(int i=0; i<3; ++i){
    double crrRad = (cenLoc - pointList[verts[i].get()].loc).length();
    if(crrRad > maxRad)maxRad = crrRad;
  }
  return maxRad;
}


// simply sets indices as passed in
template <typename T>
void
subtri<T>::setup(listIndex vert0,listIndex vert1,listIndex vert2, listIndex N0, listIndex N1, listIndex N2, unsigned int myIndex)
{
  index = myIndex;
  verts[0].set(vert0);
  verts[1].set(vert1);
  verts[2].set(vert2);
  neighbors[0].set(N0);
  neighbors[1].set(N1);
  neighbors[2].set(N2);
}


// tests for flips, then loops over the 3 children setting each up.
// the flip test is fairly complex,
// but it's the setup part that is most complex.
template <typename T>
void
subtri<T>::buildChildren(vector<subtri> &subtriList, vector<smPoint<T>> &pointList)
{
  testFlip(subtriList, pointList);
  for(int i = 0; i < 3; ++i){
    subtriList[children[i].get()].connectSubtri(index, i, subtriList, pointList);
  }
  for(int i = 0; i < 3; ++i){
    pointList[verts[i].get()].popIndex(index);
  }
  splitState = FINISHED;
}

// tests if the children triangles require a flip or not
// it flips to whichever pair does not contain the worst triangle.
// triangle quality is based on maximizing the ratio of incircle radius to circumcircle radius.
template <typename T>
void
subtri<T>::testFlip(vector<subtri> &subtriList, vector<smPoint<T>> &pointList)
{
  for(int i = 0; i < 3; ++i){
   if(flipped[i] == UNKNOWN){
      if(neighbors[i].is() == false){// this is an edge
        // appending tris was kicked down the road to connectSubtri()...
        tvec3<T> locations[4];
        locations[0] = pointList[verts[i].get()].loc;
        locations[1] = pointList[cenAddr.get()].loc;
        locations[2] = pointList[verts[(i+1)%3].get()].loc;
        locations[3] = (locations[0] + locations[2]) * .5;// midpoint on edge
        triangleCenters<T> centers[4];
        centers[0] = computeCenters(locations[0],locations[1], locations[2]);//unsplit
        centers[1] = computeCenters(locations[0],locations[1], locations[3]);//unsplit
        centers[2] = computeCenters(locations[3],locations[1], locations[2]);//unsplit
        T maxRatio = 0.0;
        int worstIndex = 0;
        for(int centersIndex=0; centersIndex<3; ++centersIndex){
          T ratio = centers[centersIndex].circumRadius / centers[centersIndex].inRadius;
          if(ratio>maxRatio){
            worstIndex = centersIndex;
            maxRatio = ratio;
          }
        }
        if(worstIndex>0){
            flipped[i] = IS_FLIPPED;
          }else{
            flipped[i] = NOT_FLIPPED;
          }
      }else{
        if(subtriList[neighbors[i].get()].splitState != NOT_SPLIT){// if we set our neighbors, this will not be needed
          tvec3<T> locations[4];
          locations[0] = pointList[verts[i].get()].loc;
          locations[1] = pointList[cenAddr.get()].loc;
          locations[2] = pointList[verts[(i+1)%3].get()].loc;
          locations[3] = pointList[subtriList[neighbors[i].get()].cenAddr.get()].loc;
          triangleCenters<T> centers[4];
          centers[0] = computeCenters(locations[0],locations[1], locations[2]);//unflipped
          centers[1] = computeCenters(locations[2],locations[3], locations[0]);//unflipped
          centers[2] = computeCenters(locations[0],locations[1], locations[3]);//flipped
          centers[3] = computeCenters(locations[1],locations[2], locations[3]);//flipped
          T maxRatio = 0.0;
          int worstIndex = 0;
          for(int centersIndex=0; centersIndex<4; ++centersIndex){
            T ratio = centers[centersIndex].circumRadius / centers[centersIndex].inRadius;
            if(ratio>maxRatio){
              worstIndex = centersIndex;
              maxRatio = ratio;
            }
          }
          int neighborindex = subtriList[neighbors[i].get()].findIndex(index);
          if(worstIndex>1){
            flipped[i] = NOT_FLIPPED;
            subtriList[neighbors[i].get()].flipped[neighborindex] = NOT_FLIPPED;
            childCen[i] = centers[0];
            subtriList[neighbors[i].get()].childCen[neighborindex] = centers[1];
          }else{
            flipped[i] = IS_FLIPPED;
            subtriList[neighbors[i].get()].flipped[neighborindex] = IS_FLIPPED;
            childCen[i] = centers[2];
            subtriList[neighbors[i].get()].childCen[neighborindex] = centers[3];
          }
        } // end neighbor is split
      } // end else, i.e. neighbors[i].is()
    } // end if(flipped[i] is unknown)
  } // end for(0-3)
} // end testFlip()

// this is the doozy.
// it is ugly, and probably should be refactored.
template <typename T>
void
subtri<T>::connectSubtri(unsigned int  parentIndex,
						int edgeNum,
						vector<subtri> &subtriList,
					    vector<smPoint<T>> &pointList
						 )
{
  listIndex vert0, vert1, vert2;
  listIndex neighbor0, neighbor1, neighbor2;
  subtri parent = subtriList[parentIndex];
  subtri edgeNeighbor =  subtriList[parent.neighbors[edgeNum].get()];
  subtri backNeighbor =  subtriList[parent.neighbors[(edgeNum+2)%3].get()];
  listIndex edgeIndex, backIndex;
  edgeIndex.set(parent.neighbors[edgeNum]);
  backIndex.set(parent.neighbors[(edgeNum+2)%3]);
  if(parent.flipped[edgeNum] == IS_FLIPPED){
    vert0.set(parent.verts[edgeNum]);
    vert1.set(edgeNeighbor.cenAddr);
    vert2.set(parent.cenAddr);
    int myIndex = edgeNeighbor.findIndex(parentIndex);
    neighbor0.set(edgeNeighbor.children[(myIndex+1)%3]);
    neighbor1.set(edgeNeighbor.children[myIndex]);
  }else{
    vert0.set(parent.verts[edgeNum]);
    vert1.set(parent.verts[(edgeNum+1)%3]);
    vert2.set(parent.cenAddr);
    neighbor0.set(edgeNeighbor.index);
    if(edgeIndex.is()){
      int nghParentIndx = edgeNeighbor.findIndex(parentIndex);
      if(edgeNeighbor.splitState != NOT_SPLIT ){
        neighbor0.set(edgeNeighbor.children[nghParentIndx]);
      }else{
        edgeNeighbor.neighbors[nghParentIndx].set(index);
      }
    }
    neighbor1.set(parent.children[(edgeNum+1)%3]);// should this be [(edgeNum+2)%3] ?
  }
  if(parent.flipped[(edgeNum+2)%3] == IS_FLIPPED){
    int myIndex = backNeighbor.findIndex(parentIndex);
    neighbor2.set(backNeighbor.children[myIndex]);
  }else{
    neighbor2.set(parent.children[(edgeNum+2)%3]);
  }
  unsigned int  childIndex = index;
  subtriList[childIndex].setup(vert0, vert1, vert2, neighbor0, neighbor1, neighbor2, childIndex);
  subtriList[childIndex].cen = parent.childCen[edgeNum];
  pointList[vert0.get()].trisIndxs.push_back(childIndex);
  pointList[vert1.get()].trisIndxs.push_back(childIndex);
  pointList[vert2.get()].trisIndxs.push_back(childIndex);
}

// utilities
template <typename T>
triangleCenters<T>
subtri<T>::computeCenters(tvec3<T> A, tvec3<T> B, tvec3<T> C) 
{
  triangleCenters<T> output;
  tvec3<T> ABmid = (A+B)*.5;
  tvec3<T> ACmid = (A+C)*.5;
  tvec3<T> AB = B - A;
  tvec3<T> AC = C - A;
  tvec3<T> CB = B - C;
  tvec3<T> N = AB.cross(AC).normalized();
  tvec3<T> d1 = AB.cross(N).normalized();
  tvec3<T> d2 = AC.cross(N).normalized();
  tvec3<T> N2 = d2.cross(d1.cross(d2));
  T U, V = (ACmid-ABmid).dot(N2)/d1.dot(N2);
  output.circumCenter = ABmid + (V*d1);
  output.circumRadius = (output.circumCenter-A).length();
  T lengthAB = AB.length();
  T lengthAC = AC.length();
  T lengthCB = CB.length();
  T semiPerimeter = (lengthAB + lengthAC + lengthCB) * .5;
  T len2 = semiPerimeter * (semiPerimeter - lengthAB) * (semiPerimeter - lengthAC) * (semiPerimeter - lengthCB);
  output.inRadius = sqrt(len2)/semiPerimeter;
  U = lengthAB / (semiPerimeter * 2);
  V = lengthAC / (semiPerimeter * 2);
  output.inCenter =  A + (AC * U) + (AB * V);
  output.N = N;
  return output;
}

// accessors
template <typename T>
int
subtri<T>::findIndex(unsigned int meshIndex)
{
  for(int i=0;i<3;++i){
    if(neighbors[i].get() == meshIndex)return i;
  }
  return -1;
}

template <typename T>
void
subtri<T>::setChildren(unsigned int addrIn, unsigned int A, unsigned int B, unsigned int C)
{
  cenAddr.set(addrIn);
  children[0].set(A);
  children[1].set(B);
  children[2].set(C);
}

