#include "math.h"
#include "stdlib.h"
#include "MashmBool.h"

/**
 * @brief Take a two-dimensional rectagular domain of elements and map the elements to processes
 *
 * Lexicographic ordering
 *
 * N N+1 
 * 0   1 2 ... N-1
 */
int decomp2dRectNumElements(int m, int n, int rank, int totalNumRanks) {
  int totalNumElements;
  int floorElementsPerRank;
  int ceilElementsPerRank;
  int numberOfElements;
  int numCeil;

  totalNumElements = m*n;
  floorElementsPerRank = (int) floor(((double)totalNumElements)/totalNumRanks);
  ceilElementsPerRank = (int) ceil(((double)totalNumElements)/totalNumRanks);

  numCeil = totalNumElements % totalNumRanks;
  if (rank < numCeil) {
    numberOfElements = ceilElementsPerRank;
  }
  else {
    numberOfElements = floorElementsPerRank;
  }
  return numberOfElements;
}


void decomp2dRectGetElements(int m, int n, int rank, int totalNumRanks, int* elemIndices) {
  int totalNumElements = m*n;
  int floorElementsPerRank;
  int ceilElementsPerRank;
  int numberOfElements;
  int numCeil;
  int offset;
  int i;

  floorElementsPerRank = (int) floor(((double)totalNumElements)/totalNumRanks);
  ceilElementsPerRank = (int) ceil(((double)totalNumElements)/totalNumRanks);

  numCeil = totalNumElements % totalNumRanks;

  if (rank < numCeil) {
    numberOfElements = ceilElementsPerRank;
    offset = rank*ceilElementsPerRank;
  }
  else {
    numberOfElements = floorElementsPerRank;
    offset = numCeil*ceilElementsPerRank + (rank - numCeil)*floorElementsPerRank;
  }

  for (i = 0; i < numberOfElements; i++) {
    elemIndices[i] = offset + i;
  }
}

int decomp2dGetOwnerRank(int mIndex, int nIndex, int totalNumRanks, int m, int n) {
  int totalNumElements;
  int floorElementsPerRank;
  int ceilElementsPerRank;
  int numberOfElements;
  int numCeil;
  int i;
  int elemIndex, rankCeil;

  totalNumElements = m*n;
  floorElementsPerRank = (int) floor(((double)totalNumElements)/totalNumRanks);
  ceilElementsPerRank = (int) ceil(((double)totalNumElements)/totalNumRanks);

  numCeil = totalNumElements % totalNumRanks;

  elemIndex = nIndex + n*mIndex;
  rankCeil = 0;
  for (i = 0; i < numCeil; i++) {
    rankCeil = rankCeil + ceilElementsPerRank;
    if (elemIndex < rankCeil) {
      printf("map1: nIndex, mIndex %d, %d owned by %d\n", mIndex, nIndex, i);
      return i;
    }
  }

  for (i = numCeil; i < totalNumRanks; i++) {
    rankCeil = rankCeil + floorElementsPerRank;
    if (elemIndex < rankCeil) {
      printf("map2: nIndex, mIndex %d, %d owned by %d\n", mIndex, nIndex, i);
      return i;
    }
    
  }
  return -1;
}

void decomp2dCreateGraph(int m, int n, int rank, int totalNumRanks, int* numElements, int** elements, int** neighbors, int** msgSizes, int* numNeighbors) {
  //int numElements;
  //int* elements;
  int neighborsMaxSize;
  int* neighborsMax;
  int* msgSizesMax;
  //int* neighbors;
  //int* msgSizes;
  //int numNeighbors;
  int i, j, k, iNeigh;
  int mIndex, nIndex, msgSize, neighborRank;
  MashmBool found;
  int mIndexOther, nIndexOther;
  int neighborCounter;
  int numberOfElements;

  /* Determine the number of elements */
  numberOfElements = decomp2dRectNumElements(m, n, rank, totalNumRanks);
  *numElements = numberOfElements;
  printf("Rank %d numElements %d\n", rank, *numElements);

  *elements = (int*) malloc(sizeof(int)*(numberOfElements));

  printf("Rank %d after allocating numElements %d\n", rank, *numElements);
  printf("Rank %d after allocating numElements %d, ptrVal = %p\n", rank, *numElements, *elements);
  if ( (*elements) == NULL ) {
    printf("Error 123: Rank %d, malloc failed\n", rank);
  }
  for (i = 0; i < numberOfElements; i++) {
    (*elements)[i] = 0;
  }
  for (i = 0; i < numberOfElements; i++) {
    printf("Test 456: Rank %d, Elements %i = %d\n", rank, i, (*elements)[i]);
  }
  decomp2dRectGetElements(m, n, rank, totalNumRanks, *elements);
  return;
  /* Allocate the maximum possible number of elements */
  neighborsMaxSize = 8*(*numElements);
  neighborsMax = (int*) malloc(sizeof(int)*neighborsMaxSize);
  msgSizesMax = (int*) malloc(sizeof(int)*neighborsMaxSize);

  /* Iterate through neighbors and add the communication */
  neighborCounter = 0;
  for (i = 0; i < *numElements; i++) {
    mIndex = (*elements)[i] /n;
    nIndex = (*elements)[i] % n;

    for (j = -1; j < 2; j++) {
      mIndexOther = mIndex + j;
      if (mIndexOther < 0 || mIndexOther >= m) {
        continue;
      }
      for (k = -1; k < 2; k++) {
        nIndexOther = nIndex + k;
        if (nIndexOther < 0 || nIndexOther >= n) {
          continue;
        }

        neighborRank = decomp2dGetOwnerRank(mIndexOther, nIndexOther, totalNumRanks, m, n);

        if (neighborRank == rank) {
          continue;
        }

        if (k == 0 || j == 0) {
          msgSize = 4;
        }
        else {
          msgSize = 1;
        }

        /* Iterate through the list of neighbors to see if it has already been added
         * update the msgSize as you go
         */
        found = false;
        for (iNeigh = 0; iNeigh < neighborCounter; iNeigh++) {
          if (neighborsMax[iNeigh] == neighborRank) {
            msgSizesMax[iNeigh] += msgSize;
            found = true;
          }
        }
        if (! found) {
          neighborsMax[neighborCounter] = neighborRank;
          msgSizesMax[neighborCounter] = msgSize;
          neighborCounter += 1;
        }
      }
    }
  }
   
  *numNeighbors = neighborCounter;
  /* Count the number of unique entries */
  *neighbors = (int*) malloc(sizeof(int)*(*numNeighbors));
  *msgSizes = (int*) malloc(sizeof(int)*(*numNeighbors));
  for (i = 0; i < (*numNeighbors); i++) {
    (*neighbors)[i] = neighborsMax[i];
    (*msgSizes)[i] = msgSizesMax[i];
  }
}

void  decomp2dDestroyGraph(int** neighbors, int** msgSizes) {
  free(*neighbors);
  free(*msgSizes);
}

void decomp2dGetCommGraph(int m, int n, int rank, int totalNumRanks, int* commProcIds, int* msgSizes) {

}
