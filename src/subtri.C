#include "subtri.h"
using namespace subTri;

template <typename T>
void
printout(int iter, vector<subtri<T>> &subtriList)
{
  return;
  int primnum = 0;
  typename vector<subtri<T>>::iterator tri;
  for (tri = subtriList.begin(); tri != subtriList.end(); ++tri){
    T ratio = tri->cen.circumRadius / tri->cen.inRadius;
    cout<<"iter "<<iter<<" primnum: "<<primnum<<" CR: "<< tri->cen.circumRadius << " IR: "<< tri->cen.inRadius<<" ratio: "<<ratio<<endl;
    //if(primnum == 4)cout<<"iter "<<iter<<" primnum: "<<primnum<<" CR: "<< tri->cen.circumCenter << " IR: "<< tri->cen.inCenter<<" ratio: "<<ratio<<endl;
    //if(primnum == 0)cout<<"iter "<<iter<<" primnum: "<<primnum<<" CR: "<< tri->childCen[0].circumCenter << " IR: "<< tri->childCen[0].inCenter<<" ratio: "<<ratio<<endl;
    ++primnum;
  }
}
/*
    splitSort(args);
    if(splitList.size() == 0)return;
    flipAndFix(maxRatio);
    createGeo();
    movePoints(args);
    sew();
    splitList.clear();
    computeNormals();
*/
template <typename T>
void
subdividedMesh<T>::subdivide(int maxIter, T maxRatio, void *args)
{
  cout<<"maxratio: "<<maxRatio<<endl;
  for(int i = 0; i<maxIter; ++i){
    //cout<<"subdivide 1"<<endl;
    splitSort(args, i);
    //cout<<"subdivide 2"<<endl;
    if(splitList.size() == 0)return;
    //cout<<"subdivide 3"<<endl;
    flipTest(maxRatio);// this can append neighbors to the list if a tri is badratio
    //cout<<"subdivide 4"<<endl;
    createGeo();//resizeGeoArrays();
    //cout<<"subdivide 5"<<endl;
    movePoints(args);
    //cout<<"subdivide 6"<<endl;
    //splitBadNeighbors()// happens *after* move points.... 
    //cout<<"subdivide 7"<<endl;
    sew(maxRatio);//attachTriangles();
    //cout<<"subdivide 8"<<endl;
    //flipBadTris();
    //cout<<"subdivide 9"<<endl;
    splitList.clear();
    //cout<<"subdivide 10"<<endl;
    //computeNormals();
  }
}

/*template <typename T>
void
subdividedMesh<T>::subdivide(int maxIter, T maxRatio, void *args)
{
  int print = 0;
  if(print)cout<<"subdivide 0"<<endl;
  for(int i = 0; i<maxIter; ++i){
     if(print)cout<<"subdivide 1 INDEX: "<<i<<endl;
     printout(i+51, subtriList);
    splitSort(args, i);
     printout(i+52, subtriList);
     if(print)cout<<"subdivide 2"<<endl;
    if(splitList.size() == 0)return;
     if(print)cout<<"subdivide 3"<<endl;
    //flipAndFix(maxRatio);
     flipTest(maxRatio);
     printout(i+53, subtriList);
     if(print)cout<<"subdivide 4"<<endl;
    createGeo();
     printout(i+54, subtriList);
     if(print)cout<<"subdivide 5"<<endl;
    movePoints(args);
     printout(i+55, subtriList);
     if(print)cout<<"subdivide 6"<<endl;
     printout(i, subtriList);
    sew();
     printout(i+100, subtriList);
     if(print)cout<<"subdivide 7"<<endl;
    splitList.clear();
     if(print)cout<<"subdivide 8"<<endl;
    computeNormals();
    //printout(i);
  }
  //cout<<"subdivide 9"<<endl;
}*/

template <typename T>
void
subdividedMesh<T>::splitSort(void *args, int iter)
{
  //typename vector<subtri<T>*>::iterator tri;
  typename vector<unsigned int>::iterator triIdx;
  //bool dump = false;
  //for (tri = newTriGroup.begin(); tri != newTriGroup.end(); ++tri){
  for (triIdx = newTriList.begin(); triIdx != newTriList.end(); ++triIdx){
    if((subtriList[*triIdx].cen.circumRadius / subtriList[*triIdx].cen.inRadius) > 10){
      for(int i=0; i<3;++i){
        if(subtriList[subtriList[*triIdx].neighbors[i].get()].splitState == NOT_SPLIT){
          subtriList[subtriList[*triIdx].neighbors[i].get()].splitState = SPLIT;
          splitList.push_back(subtriList[*triIdx].neighbors[i].get());
        }
      }
    }
  }


  for (triIdx = newTriList.begin(); triIdx != newTriList.end(); ++triIdx){
    //if(*triIdx == 166)cout<<"~~~~~~~~~~~~~   I GOT 166!!!!!!!!!!!!!!"<<endl;
    if(subtriList[*triIdx].splitState != FINISHED){
      //if(*triIdx == 166)cout<<" testing 166 "<<endl;
      if(subtriList[*triIdx].splitState != NEW_TRI)cout<<" got a not new, "<<(*triIdx)<<", state: "<<subtriList[*triIdx].splitState<<endl;
      bool split = testSplit(&(subtriList[*triIdx]), args);

      if(split){
        //if(*triIdx == 166)cout<<" 166 is split"<<endl;
        //if(subtriList[*triIdx].index == 163)cout<<"+++++++++++++ setting 163 to SPLIT"<<endl;
        subtriList[*triIdx].splitState = SPLIT;
        splitList.push_back(*triIdx);
      }else{
        //if(*triIdx == 166)cout<<"  166  not split"<<endl;
        //if(*triIdx == 163)cout<<"++++++++++++++ setting 163 to NOT_SPLIT "<<NOT_SPLIT<<endl;
        subtriList[*triIdx].setSplitState(NOT_SPLIT);
        //if(*triIdx == 163)cout<<"163 splitstate: "<< subtriList[163].splitState<<endl;
        //if(*tri == 163)(*tri)->splitState = NOT_SPLIT;
        //if(*tri == 163)cout<<"163 splitstate after: "<< subtriList[163].splitState<<endl;
        //if(*tri == 163)dump=true;
      }
    }
  }
  //if(dump){
  //  for(unsigned int i=0;i<subtriList.size();++i){
  //    if(subtriList[i].index == 163)cout<<"+++++++ I got a 163 at "<<i<<endl;
  //  }
  //}
}


template <typename T>
void
subdividedMesh<T>::flipTest(T maxRatio)
{
  unordered_set<unsigned int> additionalSplit;
  vector<unsigned int>::iterator crrSplit;
  for (crrSplit = splitList.begin(); crrSplit != splitList.end(); ++crrSplit){
    subtriList[*crrSplit].splitFlip(subtriList, pointList, additionalSplit, maxRatio, true);
  }

  unordered_set<unsigned int>::iterator addedSplit;
  unordered_set<unsigned int> loopSplit;
  //for(int i=0;i<2;++i){
    for (addedSplit = additionalSplit.begin(); addedSplit != additionalSplit.end(); ++addedSplit){
      subtriList[*addedSplit].splitFlip(subtriList, pointList, loopSplit, maxRatio, false);
      splitList.push_back(*addedSplit);
    }
    additionalSplit = loopSplit;
  //}
}

template <typename T>
void
subdividedMesh<T>::createGeo()
{
  // save the current end index for the newtrigroup creation
  unsigned int newTriCounter = subtriList.size(); 
  // create new triangles VVVVVVVVVVVVV
  subtriList.resize(subtriList.size() +  (splitList.size()*3));
  //create and fill the newtrigroup 
  //newTriGroup.clear();
  //newTriGroup.resize(splitList.size() * 3);
  newTriList.clear();
  newTriList.resize(splitList.size() * 3);
  //typename vector<subtri<T>*>::iterator crrNewTri;
  typename vector<unsigned int>::iterator crrNewTri;
  for (crrNewTri = newTriList.begin(); crrNewTri != newTriList.end(); ++crrNewTri){
    *crrNewTri = newTriCounter;
    subtriList[*crrNewTri].setIndex(newTriCounter);
    newTriCounter++;
  }

  unsigned int crrChildCount = subtriList.size() - (splitList.size() * 3);
  unsigned int Pcrr = pointList.size();
  // create new points ................
  pointList.resize(pointList.size() +  splitList.size());
  // create and fill newPointsGroup
  newPointsGroup.resize(0);
  // iterate over the splitList, setting:
  // 1. location of the next new point
  // 2. current split triangle's cenAddr
  // 3. current split triangle's  children indices
  vector<unsigned int>::iterator crrSplit;
  for (crrSplit = splitList.begin(); crrSplit != splitList.end(); ++crrSplit){
    //if(Pcrr==0)cout<<"!!!!!!!!!!!!! zero point"<<endl;
    pointList[Pcrr].loc = subtriList[*crrSplit].centroidLoc(pointList);
    newPointsGroup.push_back(Pcrr);
    subtriList[*crrSplit].setChildren(Pcrr, crrChildCount, crrChildCount + 1, crrChildCount + 2);
    crrChildCount += 3;
    Pcrr += 1;
  }
}


template <typename T>
void
subdividedMesh<T>::sew(T maxRatio)
{
  //cout<<"sew buildChildren"<<endl;
  vector<unsigned int> badRatioList;
  vector<unsigned int>::iterator splitTri;
  for (splitTri = splitList.begin(); splitTri != splitList.end(); ++splitTri){
    if(subtriList[*splitTri].splitState != FINISHED)
      subtriList[*splitTri].buildChildren(subtriList, pointList, newTriList, badRatioList, maxRatio);
  }

  //typename vector<subtri<T>>::iterator loopTri;
  //for (loopTri = subtriList.begin(); loopTri != subtriList.end(); ++loopTri){
    //if(loopTri->splitState == NEW_TRI)cout<<"before, newtri index: "<<loopTri->index<<endl;
  //}
  //return;
  /*cout<<"sew resize"<<endl;
  unsigned int newTriCounter = subtriList.size();
  subtriList.resize(subtriList.size() +  (badRatioList.size()*2));
  for(unsigned int newIndex = newTriCounter; newIndex < subtriList.size(); ++newIndex){
    subtriList[newIndex].setIndex(newIndex);
    //subtriList[newIndex].splitState = FINISHED;
    //newTriGroup.push_back(&subtriList[newIndex]);
  }*/
  /*typename vector<subtri<T>*>::iterator crrNewTri = subtriList.end();
  subtriList.resize(subtriList.size() +  (badRatioList.size()*2));
  for (; crrNewTri != subtriList.end(); ++crrNewTri){
    (*crrNewTri)->setIndex(newTriCounter);
    newTriCounter++;
  }*/
  //for (loopTri = subtriList.begin(); loopTri != subtriList.end(); ++loopTri){
    //if(loopTri->splitState == NEW_TRI)cout<<"after, newtri index: "<<loopTri->index<<endl;
  //}
  return;
  //cout<<"sew fixBadRatio"<<endl;
  unsigned int newTriCounter = subtriList.size();
  vector<unsigned int>::iterator ratioTriIdx;
  for(ratioTriIdx = badRatioList.begin(); ratioTriIdx != badRatioList.end(); ++ratioTriIdx){
    subtriList[*ratioTriIdx].fixBadRatio(subtriList, pointList, newTriCounter, newTriList);
    newTriCounter += 2;
  }
  //cout<<"sew finished"<<endl;
  //badRatioList.clear();
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
      crrPt->N += subtriList[*crrTri].cen.N;
    }
    //crrPt->N /= crrPt->trisIndxs.size();
  }
}

template <typename T>
void
subdividedMesh<T>::initializeTet(T radius, tvec3<T> scale)
{
  tvec3<T> vert0(0.0,1.0,0.0);
  tvec3<T> vert1(-.94281,-1.0/3.0,0.0);
  tvec3<T> vert2(.471405,-1.0/3.0,.816498);
  tvec3<T> vert3(.471405,-1.0/3.0,-.816498);
  vert0 *= scale;
  vert1 *= scale;
  vert2 *= scale;
  vert3 *= scale;
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
  subtriList[0].setup(indexer0,indexer3,indexer2,indexer1,indexer3,indexer2, 0, -1, subtriList);
  subtriList[0].cen = cenSource.computeCenters(pointList[0].loc,pointList[3].loc,pointList[2].loc);
  subtriList[1].setup(indexer0,indexer1,indexer3,indexer2,indexer3,indexer0, 1, -1, subtriList);
  subtriList[1].cen = cenSource.computeCenters(pointList[0].loc,pointList[1].loc,pointList[3].loc);
  subtriList[2].setup(indexer0,indexer2,indexer1,indexer0,indexer3,indexer1, 2, -1, subtriList);
  subtriList[2].cen = cenSource.computeCenters(pointList[0].loc,pointList[2].loc,pointList[1].loc);
  subtriList[3].setup(indexer1,indexer2,indexer3,indexer2,indexer0,indexer1, 3, -1, subtriList);
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
  newTriList.push_back(0);
  newTriList.push_back(1);
  newTriList.push_back(2);
  newTriList.push_back(3);
  computeNormals();
  printout(50, subtriList);
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
  tapped = false;
}

template <typename T>
void
subtri<T>::splitFlip(vector<subtri> &subtriList, vector<smPoint<T>> &pointList, unordered_set<unsigned int> &additionalSplit, T maxRatio, bool addAdditional)
{
  bool builtBadRatio = false;
  //if(index==112)cout<<"++++ splitFlip"<<endl;
  for(int i = 0; i < 3; ++i){
    if(flipped[i] == UNKNOWN){
      if(neighbors[i].is() == false){// this is an edge
        // edge case
      }else{// NOT an edge
        //if(index==112)cout<<"++++ splitFlip 1 "<< neighbors[i].get()<<" "<<subtriList[neighbors[i].get()].splitState<<endl;
        //if(subtriList[neighbors[i].get()].splitState == NEW_TRI)cout<<"    SPLIT FLIP NEW_TRI !!! "<<neighbors[i].get()<<endl;
        if((subtriList[neighbors[i].get()].splitState != NOT_SPLIT)){// if we set our neighbors, this will not be needed
        //if((subtriList[neighbors[i].get()].splitState == SPLIT) || (subtriList[neighbors[i].get()].splitState == FIX_SPLIT)){// if we set our neighbors, this will not be needed
          array<triangleCenters<T>,4> centers;
          int worstIndex = flipTester(subtriList, pointList, centers, i);
          int neighborIndex = subtriList[neighbors[i].get()].findIndex(index);
          
          if(neighbors[i].get() == 309)cout<<" neighbor of 309 index: "<<index<<" worstindex: "<<worstIndex<<endl;
          if(worstIndex>1){
            //if(index == 0)printout(-13, subtriList);
            flipped[i] = NOT_FLIPPED;
            subtriList[neighbors[i].get()].flipped[neighborIndex] = NOT_FLIPPED;
            childCen[i] = centers[0];
            subtriList[neighbors[i].get()].childCen[neighborIndex] = centers[1];
            if((centers[0].circumRadius / centers[0].inRadius) > maxRatio) builtBadRatio = true;
            if((centers[1].circumRadius / centers[1].inRadius) > maxRatio) builtBadRatio = true;
          }else{
            //if(index==112)cout<<"++++ splitFlip 2 "<<endl;
            //if(index == 0)printout(-14, subtriList);
            flipped[i] = IS_FLIPPED;
            subtriList[neighbors[i].get()].flipped[neighborIndex] = IS_FLIPPED;
            childCen[i] = centers[2];
            subtriList[neighbors[i].get()].childCen[neighborIndex] = centers[3];
            if((centers[2].circumRadius / centers[2].inRadius) > maxRatio) builtBadRatio = true;
            if((centers[3].circumRadius / centers[3].inRadius) > maxRatio) builtBadRatio = true;
          }
        } // end neighbor is split
        else
        {
          //if(neighbors[i].get() == 309)cout<<" neighbor of 309 index: "<<index<<" not split "<<addAdditional<<endl;
          childCen[i] = computeCenters(pointList[verts[i].get()].loc, centroidLoc(pointList), pointList[verts[(i+1)%3].get()].loc);
          if((childCen[i].circumRadius / childCen[i].inRadius) > maxRatio){
            if(addAdditional){
              additionalSplit.insert(neighbors[i].get());
              //cout<<" additionalSplit insert "<<neighbors[i].get()<<endl;
              subtriList[neighbors[i].get()].splitState=SPLIT;
            }
          }
        }
      } // end else, i.e. neighbors[i].is()
    } // end if(flipped[i] is unknown)
  } // end for(0-2)
  //if(!addAdditional)return;
  /*if(builtBadRatio){
    for(int i = 0; i < 3; ++i){
      if(neighbors[i].is())
        if(subtriList[neighbors[i].get()].splitState == NOT_SPLIT){
          if(0){
            additionalSplit.insert(neighbors[i].get());
          }else{
            cout<<"BADRATIO index: "<< index<<" NEIGHBOR: "<<neighbors[i].get()<<endl;
          }
        }
    }
  }*/
}

template <typename T>
tvec3<T>
subtri<T>::centroidLoc(vector<smPoint<T>> &pointList)
{
  tvec3<T>  centroid(0.0, 0.0, 0.0);
  for(int i=0; i<3; ++i){
    centroid += pointList[verts[i].get()].loc;
  }
  centroid /= 3;
  return centroid;
}


template <typename T>
double
subtri<T>::maxRad(vector<smPoint<T>> &pointList)
{
  tvec3<T> cenLoc = centroidLoc(pointList);
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
subtri<T>::setup(listIndex vert0,listIndex vert1,listIndex vert2, listIndex N0, listIndex N1, listIndex N2, unsigned int myIndex, unsigned int parentIndex,  vector<subtri> &subtriList)
{
  index = myIndex;
  verts[0].set(vert0);
  verts[1].set(vert1);
  verts[2].set(vert2);
  neighbors[0].set(N0);
  neighbors[1].set(N1);
  neighbors[2].set(N2);
  int myNeighborNum = subtriList[N0.get()].findIndex(parentIndex);
  if(myNeighborNum != -1)subtriList[N0.get()].neighbors[myNeighborNum].set(myIndex);
  myNeighborNum = subtriList[N1.get()].findIndex(parentIndex);
  if(myNeighborNum != -1)subtriList[N1.get()].neighbors[myNeighborNum].set(myIndex);
  myNeighborNum = subtriList[N2.get()].findIndex(parentIndex);
  if(myNeighborNum != -1)subtriList[N2.get()].neighbors[myNeighborNum].set(myIndex);
}


// tests for flips, then loops over the 3 children setting each up.
// the flip test is fairly complex,
// but it's the setup part that is most complex.
template <typename T>
void
subtri<T>::buildChildren(vector<subtri> &subtriList, vector<smPoint<T>> &pointList, vector<unsigned int> &newTriList, vector<unsigned int> &badRatioList, T maxRatio)
{
  for(int i = 0; i < 3; ++i){
     //if(children[i].get() == 4)printout(-3, subtriList);
    subtriList[children[i].get()].connectSubtri(index, i, subtriList, pointList, newTriList, badRatioList, maxRatio);
     //if(children[i].get() == 4)printout(-4, subtriList);
  }
  for(int i = 0; i < 3; ++i){
    pointList[verts[i].get()].popIndex(index);
  }
  splitState = FINISHED;
}


template <typename T>
void
subtri<T>::get4Locations(vector<subtri> &subtriList, vector<smPoint<T>> &pointList, int i, array<tvec3<T>,4> locations)
{
  locations[0] = pointList[verts[i].get()].loc;
  locations[1] = centroidLoc(pointList);
  locations[2] = pointList[verts[(i+1)%3].get()].loc;
  locations[3] = subtriList[neighbors[i].get()].centroidLoc(pointList);
}

template <typename T>
int
subtri<T>::test4Locations(array<triangleCenters<T>,4> &centers, array<tvec3<T>,4> locations)
{
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
  return worstIndex;
}


template <typename T>
int
subtri<T>::flipTester(vector<subtri> &subtriList, vector<smPoint<T>> &pointList, array<triangleCenters<T>,4> centers, int i)
{
  array<tvec3<T>,4> locations;
  get4Locations(subtriList, pointList, i, locations);
  return test4Locations(centers, locations);
}


// this is the doozy.
// it is ugly, and probably should be refactored.
template <typename T>
void
subtri<T>::connectSubtri(unsigned int  parentIndex,
						int edgeNum,
						vector<subtri> &subtriList,
					    vector<smPoint<T>> &pointList,
              vector<unsigned int> &newTriList,
              vector<unsigned int> &badRatioList,
              T maxRatio
						 )
{

  listIndex vert0, vert1, vert2;
  listIndex neighbor0, neighbor1, neighbor2;
  subtri parent = subtriList[parentIndex];
  subtri edgeNeighbor =  subtriList[parent.neighbors[edgeNum].get()];
  //if(parentIndex==112)cout<<"112: en: "<<parent.neighbors[edgeNum].get()<<" "<<index<<endl;
  subtri backNeighbor =  subtriList[parent.neighbors[(edgeNum+2)%3].get()];
  listIndex edgeIndex, backIndex;
  edgeIndex.set(parent.neighbors[edgeNum]);
  backIndex.set(parent.neighbors[(edgeNum+2)%3]);
  //if(parentIndex==112)cout<<"112: pi: "<<
  //if(index==213)cout<<"  213 some info "<<parent.index<<" "<<edgeNum<<" "<<edgeNeighbor.index<<endl;
  if(parent.flipped[edgeNum] == IS_FLIPPED){
    vert0.set(parent.verts[edgeNum]);
    //if(index==213)cout<<" set 213 vert1 A to "<<edgeNeighbor.cenAddr.get()<<" enidx: "<<edgeNeighbor.index<<endl;
    vert1.set(edgeNeighbor.cenAddr);
    vert2.set(parent.cenAddr);
    int myIndex = edgeNeighbor.findIndex(parentIndex);
    neighbor0.set(edgeNeighbor.children[(myIndex+1)%3]);
    neighbor1.set(edgeNeighbor.children[myIndex]);
    //if(parentIndex==112)cout<<"112: A n0: "<<neighbor0.get()<<" mi: "<<myIndex<<endl;
    //if(parentIndex==112)cout<<"112: A n1: "<<neighbor1.get()<<endl;
    //if(parentIndex==112)cout<<"112: A en d: "<<edgeNeighbor.index<<" "<<edgeNeighbor.children[2].get()<<endl;
  }else{
    vert0.set(parent.verts[edgeNum]);
    //if(index==213)cout<<" set 213 vert1 B to "<<parent.verts[(edgeNum+1)%3].get()<<" pidx: "<<parent.index<<endl;
    vert1.set(parent.verts[(edgeNum+1)%3]);
    vert2.set(parent.cenAddr);
    neighbor0.set(edgeNeighbor.index);
    //if(parentIndex==112)cout<<"112: B n0: "<<neighbor0.get()<<endl;
    if(edgeIndex.is()){
      int nghParentIndx = edgeNeighbor.findIndex(parentIndex);
      if(edgeNeighbor.splitState != NOT_SPLIT ){
        neighbor0.set(edgeNeighbor.children[nghParentIndx]);
        //if(parentIndex==112)cout<<"112: C n0: "<<neighbor0.get()<<endl;
      }else{
        subtriList[edgeNeighbor.index].neighbors[nghParentIndx].set(index);
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

  subtriList[index].cen = computeCenters(pointList[vert0.get()].loc, pointList[vert1.get()].loc, pointList[vert2.get()].loc);
  T ratio = subtriList[index].cen.circumRadius / subtriList[index].cen.inRadius;
  if(ratio>maxRatio) badRatioList.push_back(index);
  //if(parentIndex==112)cout<<"112: v: "<<vert0.get()<<" "<<vert1.get()<<" "<<vert2.get()<<" n: "<<neighbor0.get()<<" "<<neighbor1.get()<<" "<<neighbor2.get()<<" i: "<<index<<" pi: "<<parentIndex<<endl;
  //if(index==213)cout<<">>>>>>>>> 213: v: "<<vert0.get()<<" "<<vert1.get()<<" "<<vert2.get()<<" n: "<<neighbor0.get()<<" "<<neighbor1.get()<<" "<<neighbor2.get()<<" i: "<<index<<" pi: "<<parentIndex<<endl;
  //if(index==2657)cout<<">>>>>>>>> 2657: v: "<<vert0.get()<<" "<<vert1.get()<<" "<<vert2.get()<<" n: "<<neighbor0.get()<<" "<<neighbor1.get()<<" "<<neighbor2.get()<<" i: "<<index<<" pi: "<<parentIndex<<endl;
  if(index==2533)cout<<">>>>>>>>> 2533: v: "<<vert0.get()<<" "<<vert1.get()<<" "<<vert2.get()<<" n: "<<neighbor0.get()<<" "<<neighbor1.get()<<" "<<neighbor2.get()<<" i: "<<index<<" pi: "<<parentIndex<<endl;
  subtriList[index].setup(vert0, vert1, vert2, neighbor0, neighbor1, neighbor2, index, parentIndex, subtriList);
  subtriList[index].parentIndex = parentIndex;
  pointList[vert0.get()].trisIndxs.push_back(index);
  pointList[vert1.get()].trisIndxs.push_back(index);
  pointList[vert2.get()].trisIndxs.push_back(index);
}


// uuuugly uuugly 
template <typename T>
void
subtri<T>::fixBadRatio(vector<subtri> &subtriList, vector<smPoint<T>> &pointList, unsigned int newTriIdx, vector<unsigned int> &newTriList) 
{
  cout<<" fixBadRatio 0 index: --------- "<<index<<endl;
  if(splitState == FINISHED)return;
  //if(newTriIdx==166)cout<<"~~~~~~~~~~~~~~ huh  "<<subtriList[newTriIdx].splitState<<endl;

  //cout<<"sew resize"<<endl;
  unsigned int newTriCounter = subtriList.size();
  subtriList.resize(subtriList.size() + 2);
  for(unsigned int newIndex = newTriCounter; newIndex < subtriList.size(); ++newIndex){
    subtriList[newIndex].setIndex(newIndex);
  }


  array<tvec3<T>,4> locations;
  T bestRatio = 100000;
  int bestNeighbor, nIndx;
  unsigned int farvert;
  /*for(int i = 0; i < 3; ++i){
    nIndx = subtriList[neighbors[i].get()].findIndex(index);
    farvert = subtriList[neighbors[i].get()].verts[(nIndx+2)%3].get();
    locations[0] = pointList[verts[(i+0)%3].get()].loc;
    locations[1] = pointList[farvert].loc;
    locations[2] = pointList[verts[(i+1)%3].get()].loc;
    locations[3] = pointList[verts[(i+2)%3].get()].loc;
    triangleCenters<T> centers0 = computeCenters(locations[1],locations[2], locations[3]);
    triangleCenters<T> centers1 = computeCenters(locations[3],locations[0], locations[1]);
    T ratio0 = centers0.circumRadius / centers0.inRadius;
    T ratio1 = centers1.circumRadius / centers1.inRadius;
    T worstOfTwo = (ratio0>ratio1 ? ratio0 : ratio1);
    if(worstOfTwo < bestRatio){
      bestRatio = worstOfTwo;
      bestNeighbor = i;
    }
  }*/
  T largestVal = 0;
  for(int i = 0; i < 3; ++i){
    if(subtriList[neighbors[i].get()].cen.inRadius > largestVal){
      largestVal = subtriList[neighbors[i].get()].cen.inRadius;
      bestNeighbor = i;
    }
  }
  /*if(pointList.size() > 455){
    cout<<
    return;
  }*/
  //if(newTriIdx==2657)cout<<" ~~~~~~~~~~~~~~~~ handling 2657, making: "<<newTriIdx<<" and: "<<(newTriIdx+1)<<endl;
  cout<<" fixBadRatio index: "<<index<<" parentIndex: "<<parentIndex<<" verts: "<<verts[0].get()<<" "<<verts[1].get()<<" "<<verts[2].get()<<endl;
  nIndx = subtriList[neighbors[bestNeighbor].get()].findIndex(index);
  //farvert = subtriList[neighbors[bestNeighbor].get()].verts[(nIndx+2)%3].get();
  listIndex vIndex0  = verts[(bestNeighbor+0)%3];
  listIndex vIndex1 = subtriList[neighbors[bestNeighbor].get()].verts[(nIndx+2)%3];
  listIndex vIndex2 = verts[(bestNeighbor+1)%3];
  listIndex vIndex3 = verts[(bestNeighbor+2)%3];
  listIndex nIndex0 = neighbors[(bestNeighbor+2)%3];
  listIndex nIndex1 = subtriList[neighbors[bestNeighbor].get()].neighbors[(nIndx+1)%3];
  listIndex nIndex2 = subtriList[neighbors[bestNeighbor].get()].neighbors[(nIndx+2)%3];
  listIndex nIndex3 = neighbors[(bestNeighbor+1)%3];

  cout<<" I got " <<bestNeighbor<<" "<<neighbors[bestNeighbor].get()<<endl;
  cout<<" triA: "<<newTriIdx<<" "<<vIndex3.get()<<" "<<vIndex0.get()<<" "<<vIndex1.get()<<" n: "<<nIndex0.get()<<" "<<nIndex1.get()<<" newguy "<<endl;
  cout<<" triB: "<<(newTriIdx+1)<<" "<<vIndex3.get()<<" "<<vIndex1.get()<<" "<<vIndex2.get()<<" n: newguy "<<nIndex2.get()<<" "<<nIndex3.get()<<endl;

  listIndex new0, new1;
  new0.set(newTriIdx);
  new1.set(newTriIdx + 1);

  //cout<<" fixBadRatio 2"<<endl;
  subtriList[newTriIdx].setup(vIndex3, vIndex0, vIndex1, nIndex0, nIndex1, new1, newTriIdx, -1, subtriList);
  subtriList[newTriIdx].cen = computeCenters(pointList[vIndex3.get()].loc,pointList[vIndex0.get()].loc,pointList[vIndex1.get()].loc);
  pointList[vIndex3.get()].popIndex(index);
  pointList[vIndex3.get()].trisIndxs.push_back(newTriIdx);
  pointList[vIndex0.get()].popIndex(index);
  pointList[vIndex0.get()].trisIndxs.push_back(newTriIdx);
  pointList[vIndex1.get()].popIndex(index);
  pointList[vIndex1.get()].trisIndxs.push_back(newTriIdx);
  newTriList.push_back(newTriIdx);
  //if(newTriIdx==2657)cout<<"~~~~~~~~~~~~~~ pushing 2657 onto the list 0 "<<subtriList[newTriIdx].splitState<<endl;
  //if(newTriIdx==3974)cout<<"~~~~~~~~~~~~~~ pushing 3974 onto the list 0 "<<subtriList[newTriIdx].splitState<<" index: "<<index<<endl;

  //cout<<" fixBadRatio 3"<<endl;
  subtriList[newTriIdx+1].setup(vIndex3, vIndex1, vIndex2, new0, nIndex2, nIndex3, newTriIdx+1, -1, subtriList);
  subtriList[newTriIdx+1].cen = computeCenters(pointList[vIndex3.get()].loc,pointList[vIndex1.get()].loc,pointList[vIndex2.get()].loc);
  pointList[vIndex3.get()].popIndex(neighbors[bestNeighbor].get());
  pointList[vIndex3.get()].trisIndxs.push_back(newTriIdx+1);
  pointList[vIndex1.get()].popIndex(neighbors[bestNeighbor].get());
  pointList[vIndex1.get()].trisIndxs.push_back(newTriIdx+1);
  pointList[vIndex2.get()].popIndex(neighbors[bestNeighbor].get());
  pointList[vIndex2.get()].trisIndxs.push_back(newTriIdx+1);
  newTriList.push_back(newTriIdx+1);
  //if((newTriIdx+1)==2657)cout<<"~~~~~~~~~~~~~~ pushing 2657 onto the list 1 "<<subtriList[newTriIdx].splitState<<endl;
  //if((newTriIdx+1)==3974)cout<<"~~~~~~~~~~~~~~ pushing 3974 onto the list 1 "<<subtriList[newTriIdx].splitState<<endl;

  //cout<<" fixBadRatio 4"<<endl;
  subtriList[index].splitState = FINISHED;
  cout<<"  setting finished: "<<neighbors[bestNeighbor].get()<<endl;
  subtriList[neighbors[bestNeighbor].get()].splitState = FINISHED;

  //cout<<" fixBadRatio 5"<<endl;
  unsigned int neighbor0, neighbor1, neighbor2, neighbor3;
  neighbor0 = neighbors[(bestNeighbor+2)%3].get();
  neighbor1 = subtriList[neighbors[bestNeighbor].get()].neighbors[(nIndx+1)%3].get();
  neighbor2 = subtriList[neighbors[bestNeighbor].get()].neighbors[(nIndx+2)%3].get();
  neighbor3 = neighbors[(bestNeighbor+1)%3].get();

  //cout<<" fixBadRatio 6 index: "<<index<<" "<<subtriList.size()<< endl;
  unsigned int n0Indx, n1Indx, n2Indx, n3Indx;
  n0Indx = subtriList[neighbor0].findIndex(index);
  cout<<" n0Indx: "<<neighbor0<<" "<<index<<" "<<parentIndex<<" "<<n0Indx<<" "<<subtriList[neighbor0].neighbors[0].get()<<" "<<subtriList[neighbor0].neighbors[1].get()<<" "<<subtriList[neighbor0].neighbors[2].get()<<endl;

  n1Indx = subtriList[neighbor1].findIndex(neighbors[bestNeighbor].get());
  n2Indx = subtriList[neighbor2].findIndex(neighbors[bestNeighbor].get());
  n3Indx = subtriList[neighbor3].findIndex(index);

  cout<<" fixBadRatio 7 neigh: "<<neighbor0<<" "<<neighbor1<<" "<<neighbor2<<" "<<neighbor3<<" Nindx: "<<n0Indx<<" "<<n1Indx<<" "<<n2Indx<<" "<<n3Indx<<endl;
  /*if(pointList.size() > 455){
    cout<<"   oh shit... "<<index<<" "<<parentIndex<<endl;
    return;
  }*/
  subtriList[neighbor0].neighbors[n0Indx] = newTriIdx;
  subtriList[neighbor1].neighbors[n1Indx] = newTriIdx;
  subtriList[neighbor2].neighbors[n2Indx] = newTriIdx+1;
  subtriList[neighbor3].neighbors[n3Indx] = newTriIdx+1;

  cout<<" fixBadRatio 8 (finished)"<<endl;
  //unsigned int anotherIdx = subtriList[neighbors[(bestNeighbor+2)%3].get()].findIndex(index);
  //subtriList[neighbors[(bestNeighbor+2)%3].get()].neighbors[anotherIdx] = newTriIdx;
  //anotherIdx = subtriList[neighbors[bestNeighbor].get()].neighbors[(nIndx+1)%3].findIndex(neighbors[bestNeighbor].get());
  //subtriList[neighbors[bestNeighbor].get()].neighbors[anotherIdx] = newTriIdx;
  //anotherIdx = subtriList[neighbors[bestNeighbor].get()].neighbors[(nIndx+2)%3].findIndex(neighbors[bestNeighbor].get());
  //anotherIdx = subtriList[neighbors[(bestNeighbor+2)%3].get()].findIndex(index);
  //subtriList[neighbors[bestNeighbor].get()].neighbors[(nIndx+1)%3] = anotherIdx;

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



