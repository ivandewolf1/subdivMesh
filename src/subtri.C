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
    //cout<<"iter "<<iter<<" primnum: "<<primnum<<" CR: "<< tri->cen.circumRadius << " IR: "<< tri->cen.inRadius<<" ratio: "<<ratio<<endl;
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
  int print = 0;
  if(print)cout<<"subdivide 0"<<endl;
  for(int i = 0; i<maxIter; ++i){
     if(print)cout<<"subdivide 1 INDEX: "<<i<<endl;
     printout(i+51, subtriList);
    splitSort(args);
     printout(i+52, subtriList);
     if(print)cout<<"subdivide 2"<<endl;
    if(splitList.size() == 0)return;
     if(print)cout<<"subdivide 3"<<endl;
    flipAndFix(maxRatio);
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
}

template <typename T>
void
subdividedMesh<T>::splitSort(void *args)
{
  typename vector<subtri<T>*>::iterator tri;
  for (tri = newTriGroup.begin(); tri != newTriGroup.end(); ++tri){
    bool split = testSplit(*tri, args);
    if(split){
      (*tri)->splitState = SPLIT;
      splitList.push_back((*tri)->index);
    }else{
      (*tri)->splitState = NOT_SPLIT;
      unSplitList.push_back((*tri)->index);
    }
  }
}

template <typename T>
void
subdividedMesh<T>::flipAndFix(T maxRatio)
{
  int print = 0;

  vector<unsigned int>::iterator crrUnSplit;
  for (crrUnSplit = unSplitList.begin(); crrUnSplit != unSplitList.end(); ++crrUnSplit){
    /*if(*crrUnSplit == 892){
      print = 0;
    }else{
      print = 0;
    }*/
    //cout<<"flipAndFix 1 cr: "<< subtriList[*crrUnSplit].cen.circumRadius << " ir: "<< subtriList[*crrUnSplit].cen.inRadius<<endl;
    T ratio = subtriList[*crrUnSplit].cen.circumRadius / subtriList[*crrUnSplit].cen.inRadius;
    if(ratio > maxRatio && false){
      //cout<<"hit badRatioFlip "<<endl;
      subtriList[*crrUnSplit].badRatioFlip(subtriList, pointList, splitList);
    }else{
      subtriList[*crrUnSplit].unsplitSplitNeighbors(subtriList, pointList, splitList);
    }
    printout(62, subtriList);
  }  
  vector<unsigned int>::iterator crrSplit;
  for (crrSplit = splitList.begin(); crrSplit != splitList.end(); ++crrSplit){
    subtriList[*crrSplit].splitFlip(subtriList, pointList);
  }
  unSplitList.clear();
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
  newTriGroup.clear();
  newTriGroup.resize(splitList.size() * 3);
  typename vector<subtri<T>*>::iterator crrNewTri;
  for (crrNewTri = newTriGroup.begin(); crrNewTri != newTriGroup.end(); ++crrNewTri){
    *crrNewTri = &(subtriList[newTriCounter]);
    (*crrNewTri)->setIndex(newTriCounter);
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
    if(Pcrr==0)cout<<"!!!!!!!!!!!!! zero point"<<endl;
    pointList[Pcrr].loc = subtriList[*crrSplit].centroidLoc(pointList);
    newPointsGroup.push_back(Pcrr);
    subtriList[*crrSplit].setChildren(Pcrr, crrChildCount, crrChildCount + 1, crrChildCount + 2);
    crrChildCount += 3;
    Pcrr += 1;
  }
}


template <typename T>
void
subdividedMesh<T>::sew()
{
  vector<unsigned int>::iterator splitTri;
  for (splitTri = splitList.begin(); splitTri != splitList.end(); ++splitTri){
     if((*splitTri) == 4)printout(-1, subtriList);
    subtriList[*splitTri].buildChildren(subtriList, pointList, newTriGroup);
     if((*splitTri) == 4)printout(-2, subtriList);
  }


  // postpass, handle ratioNeighborGroup here 


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
  newTriGroup.push_back(&(subtriList[0]));
  newTriGroup.push_back(&(subtriList[1]));
  newTriGroup.push_back(&(subtriList[2]));
  newTriGroup.push_back(&(subtriList[3]));
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

  twoBadRatioNeighbors = true;
}

template <typename T>
void
subtri<T>::badRatioFlip(vector<subtri> &subtriList, vector<smPoint<T>> &pointList, vector<unsigned int> &splitList)
{
  cout<<"\nbadRatioFlip 0"<<endl;
  triangleCenters<T> centers[4];
  int nFlipped = 0;
  for(int i=0; i<3; ++i){
    if(1){//flipped[i] == UNKNOWN){
      int worstIndex = flipTester(subtriList, pointList, centers, i);
      int neighborIndex = subtriList[neighbors[i].get()].findIndex(index);
      cout<<index<<" "<<i<<" "<<subtriList[neighbors[i].get()].splitState<<endl;
      /*if(worstIndex>1){
        flipped[i] = NOT_FLIPPED;
        subtriList[neighbors[i].get()].flipped[neighborIndex] = NOT_FLIPPED;
        childCen[i] = centers[0];
        subtriList[neighbors[i].get()].childCen[neighborIndex] = centers[1];
      }else{// if flipped
        nFlipped += 1;
        flipped[i] = IS_FLIPPED;
        if((subtriList[neighbors[i].get()].splitState != SPLIT) && (subtriList[neighbors[i].get()].splitState != FIX_SPLIT)){
          subtriList[neighbors[i].get()].splitState = FIX_SPLIT;
        }
        subtriList[neighbors[i].get()].flipped[neighborIndex] = IS_FLIPPED;
        childCen[i] = centers[2];
        subtriList[neighbors[i].get()].childCen[neighborIndex] = centers[3];
        splitList.push_back(neighbors[i].get());
      }// end else, flipped */
    }
    if(children[i].get() == 4)printout(-9, subtriList);
  }// end for i in 0-2
  /*if(nFlipped > 0){// if nothing improves it, don't split it.
    splitState = FIX_SPLIT;
    splitList.push_back(index);
  }*/
  cout<<"badRatioFlip 10"<<endl;
}

template <typename T>
void
subtri<T>::unsplitSplitNeighbors(vector<subtri> &subtriList, vector<smPoint<T>> &pointList, vector<unsigned int> &splitList)
{
  if(splitState == SPLIT)return;
  triangleCenters<T> centers[4];
  int nFlipped = 0;
  for(int i=0; i<3;++i){
    if(subtriList[neighbors[i].get()].splitState == SPLIT){
      int worstIndex = flipTester(subtriList, pointList, centers, i);// <<<<<<<<<<<<<<<<< AWRY
      int neighborIndex = subtriList[neighbors[i].get()].findIndex(index);
      if(worstIndex>1){
        printout(-11, subtriList);
        flipped[i] = NOT_FLIPPED;
        subtriList[neighbors[i].get()].flipped[neighborIndex] = NOT_FLIPPED;
        childCen[i] = centers[0];
        subtriList[neighbors[i].get()].childCen[neighborIndex] = centers[1];
      }else{// if flipped
        printout(-12, subtriList);
        nFlipped += 1;
        flipped[i] = IS_FLIPPED;
        subtriList[neighbors[i].get()].flipped[neighborIndex] = IS_FLIPPED;
        childCen[i] = centers[2];
        subtriList[neighbors[i].get()].childCen[neighborIndex] = centers[3];
        splitState = FIX_SPLIT;
        printout(-13, subtriList);
      }// end else, flipped 
      printout(-14, subtriList);
    }// end if neighbor split
    printout(-15, subtriList);
  }// end of for i in 0-2
  if(nFlipped > 0){// if nothing improves it, don't split it.
    splitList.push_back(index);
  }
}

template <typename T>
void
subtri<T>::splitFlip(vector<subtri> &subtriList, vector<smPoint<T>> &pointList)
{
  for(int i = 0; i < 3; ++i){
    if(index==326)cout<<"checkit"<<endl;
    if(flipped[i] == UNKNOWN){
      if(index==326)cout<<"unknown"<<endl;
      if(neighbors[i].is() == false){// this is an edge
        // edge case
      }else{// NOT an edge
        if(subtriList[neighbors[i].get()].splitState != NOT_SPLIT){// if we set our neighbors, this will not be needed
          triangleCenters<T> centers[4];
          int worstIndex = flipTester(subtriList, pointList, centers, i);
          int neighborIndex = subtriList[neighbors[i].get()].findIndex(index);
          
          if(worstIndex>1){
            if(index == 0)printout(-13, subtriList);
            flipped[i] = NOT_FLIPPED;
            subtriList[neighbors[i].get()].flipped[neighborIndex] = NOT_FLIPPED;
            childCen[i] = centers[0];
            subtriList[neighbors[i].get()].childCen[neighborIndex] = centers[1];
          }else{
            if(index == 0)printout(-14, subtriList);
            flipped[i] = IS_FLIPPED;
            subtriList[neighbors[i].get()].flipped[neighborIndex] = IS_FLIPPED;
            childCen[i] = centers[2];
            subtriList[neighbors[i].get()].childCen[neighborIndex] = centers[3];
          }
        } // end neighbor is split
      } // end else, i.e. neighbors[i].is()
    } // end if(flipped[i] is unknown)
  } // end for(0-2)
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
subtri<T>::buildChildren(vector<subtri> &subtriList, vector<smPoint<T>> &pointList, vector<subtri<T>*> &newTriGroup)
{
  for(int i = 0; i < 3; ++i){
     if(children[i].get() == 4)printout(-3, subtriList);
    subtriList[children[i].get()].connectSubtri(index, i, subtriList, pointList, newTriGroup);
     if(children[i].get() == 4)printout(-4, subtriList);
  }
  for(int i = 0; i < 3; ++i){
    pointList[verts[i].get()].popIndex(index);
  }
  splitState = FINISHED;
}


template <typename T>
int
subtri<T>::flipTester(vector<subtri> &subtriList, vector<smPoint<T>> &pointList, triangleCenters<T> centers[4], int i)
{
  tvec3<T> locations[4];
  locations[0] = pointList[verts[i].get()].loc;
  locations[1] = centroidLoc(pointList);
  locations[2] = pointList[verts[(i+1)%3].get()].loc;
  locations[3] = subtriList[neighbors[i].get()].centroidLoc(pointList);
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


// this is the doozy.
// it is ugly, and probably should be refactored.
template <typename T>
void
subtri<T>::connectSubtri(unsigned int  parentIndex,
						int edgeNum,
						vector<subtri> &subtriList,
					    vector<smPoint<T>> &pointList,
              vector<subtri<T>*> &newTriGroup
						 )
{
  int special = 283;
  bool print = false;
  listIndex vert0, vert1, vert2;
  listIndex neighbor0, neighbor1, neighbor2;
  subtri parent = subtriList[parentIndex];
  subtri edgeNeighbor =  subtriList[parent.neighbors[edgeNum].get()];
  subtri backNeighbor =  subtriList[parent.neighbors[(edgeNum+2)%3].get()];
  listIndex edgeIndex, backIndex;
  edgeIndex.set(parent.neighbors[edgeNum]);
  backIndex.set(parent.neighbors[(edgeNum+2)%3]);
  if(index==special)cout<<"special 0 "<<edgeNeighbor.index<<endl;
  if(parent.flipped[edgeNum] == IS_FLIPPED){
    vert0.set(parent.verts[edgeNum]);
    vert1.set(edgeNeighbor.cenAddr);
    if(index==special)cout<<"special 1 "<<vert1.get()<<endl;
    //if(edgeNeighbor.cenAddr.get() == 0)cout<<"cenaddr0 idx: "<<edgeNeighbor.index<<" index: "<<index<<endl;
    vert2.set(parent.cenAddr);
    int myIndex = edgeNeighbor.findIndex(parentIndex);
    neighbor0.set(edgeNeighbor.children[(myIndex+1)%3]);
    neighbor1.set(edgeNeighbor.children[myIndex]);
  }else{
    vert0.set(parent.verts[edgeNum]);
    vert1.set(parent.verts[(edgeNum+1)%3]);
    if(index==special)cout<<"special 2 "<<vert1.get()<<endl;
    vert2.set(parent.cenAddr);
    neighbor0.set(edgeNeighbor.index);
    if(edgeIndex.is()){
      int nghParentIndx = edgeNeighbor.findIndex(parentIndex);
      if(edgeNeighbor.splitState != NOT_SPLIT ){
        neighbor0.set(edgeNeighbor.children[nghParentIndx]);
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

  //unsigned int  childIndex = index;
  //subtriList[index].setup(vert0, vert1, vert2, neighbor0, neighbor1, neighbor2, index);
  //subtriList[childIndex].cen = parent.childCen[edgeNum];
  subtriList[index].cen = computeCenters(pointList[vert0.get()].loc, pointList[vert1.get()].loc, pointList[vert2.get()].loc);
  T ratio = subtriList[index].cen.circumRadius / subtriList[index].cen.inRadius;
  if(ratio>10){
    if(edgeNeighbor.splitState==NOT_SPLIT){


      // this is where the ratioNeighborGroup should get filled



        cout<<"connectSubtri bad ratio "<<index<<endl;
        //cout<<vert0.get()<<" "<<vert1.get()<<" "<<vert2.get()<<endl;
      //cout<<edgeNeighbor.index<<" "<<neighbor0.get()<<endl;
      int myIndex = edgeNeighbor.findIndex(parentIndex);
      int farvert = (myIndex+2)%3;
      //cout<<myIndex<<" "<<edgeNeighbor.verts[farvert].get()<<endl;
      listIndex nVert0, nVert1, nVert2, nNgh0, nNgh1, nNgh2;
      nVert0.set(vert0);
      nVert1.set(edgeNeighbor.verts[farvert]);
      nVert2.set(vert2);
      nNgh0.set((myIndex+1)%3);
      nNgh1.set(index);
      nNgh2.set(neighbor2);
        //cout<< myIndex<<" "<<nVert0.get()<<" "<<nVert1.get()<<" "<<nVert2.get()<<" "<<endl;
        //cout<<edgeNeighbor.neighbors[myIndex].get()<<endl;

      subtri<T> newTri;
      unsigned int newIndex = subtriList.size();
      subtriList.push_back(newTri);
      subtriList[newIndex].cen = computeCenters(pointList[nVert0.get()].loc, pointList[nVert1.get()].loc, pointList[nVert2.get()].loc);
      if(index==special)cout<<"special badratio new cen "<<endl;
      subtriList[newIndex].setup(nVert0, nVert1, nVert2, nNgh0, nNgh1, nNgh2, newIndex, edgeNeighbor.index, subtriList);
        cout<<"setupnew: "<<nVert0.get()<<" "<<nVert1.get()<<" "<<nVert2.get()<<" "<<nNgh0.get()<<" "<<nNgh1.get()<<" "<<nNgh2.get()<<"          -> "<<newIndex<<endl;
      subtriList[newIndex].parentIndex = edgeNeighbor.index;
      pointList[nVert0.get()].trisIndxs.push_back(newIndex);
      pointList[nVert1.get()].trisIndxs.push_back(newIndex);
      pointList[nVert2.get()].trisIndxs.push_back(newIndex);
      subtriList[newIndex].splitState = NOT_SPLIT;
      newTriGroup.push_back(&(subtriList[newIndex]));

      vert0.set(edgeNeighbor.verts[farvert]);
      neighbor0.set(edgeNeighbor.neighbors[(myIndex+2)%3]);
      neighbor2.set(newIndex);
      subtriList[index].cen = computeCenters(pointList[vert0.get()].loc, pointList[vert1.get()].loc, pointList[vert2.get()].loc);
        cout<<"setup: "<<vert0.get()<<" "<<vert1.get()<<" "<<vert2.get()<<" "<<neighbor0.get()<<" "<<neighbor1.get()<<" "<<neighbor2.get()<<" "<<index<<" "<<parentIndex<<endl;
    }
  }
  if(edgeNeighbor.cenAddr.get() == 0){
    //cout<<"cenaddr0 idx: "<<edgeNeighbor.index<<" index: "<<index<<endl;
    cout<<vert0.get()<<" "<<vert1.get()<<" "<< vert2.get()<<" "<< neighbor0.get()<<" "<< neighbor1.get()<<" "<< neighbor2.get()<<" "<< index<<" "<< parentIndex<<" "<<edgeNeighbor.index<<endl;
  }
  subtriList[index].setup(vert0, vert1, vert2, neighbor0, neighbor1, neighbor2, index, parentIndex, subtriList);
  subtriList[index].parentIndex = parentIndex;
  pointList[vert0.get()].trisIndxs.push_back(index);
  pointList[vert1.get()].trisIndxs.push_back(index);
  pointList[vert2.get()].trisIndxs.push_back(index);


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

