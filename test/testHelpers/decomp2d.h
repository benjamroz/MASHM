#include "math.h"
#include "stdlib.h"
#include "MashmBool.h"

#include "cmake_fortran_c_interface.h"

#define decomp2dCreateGraph FCI_GLOBAL(decomp2dcreategraph,DECOMP2DCREATEGRAPH)
#define decomp2dDestroyGraph FCI_GLOBAL(decomp2ddestroygraph,DECOMP2DDESTROYGRAPH)
/**
 * @brief Take a two-dimensional rectagular domain of elements and map the elements to processes
 *
 * Lexicographic ordering
 *
 * N N+1 
 * 0   1 2 ... N-1
 */
int decomp2dRectNumElements(int m, int n, int rank, int totalNumRanks); 

void decomp2dRectGetElements(int m, int n, int rank, int totalNumRanks, int* elemIndices); 

int decomp2dGetOwnerRank(int mIndex, int nIndex, int totalNumRanks, int m, int n); 

void decomp2dCreateGraph(int m, int n, int rank, int totalNumRanks, int* numElements, int** elements, int** neighbors, int** msgSizes, int* numNeighbors); 

void  decomp2dDestroyGraph(int** neighbors, int** msgSizes); 

void decomp2dGetCommGraph(int m, int n, int rank, int totalNumRanks, int* commProcIds, int* msgSizes);
