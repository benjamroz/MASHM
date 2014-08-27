#ifndef MASHM_H
#define MASHM_H

#include "mpi.h"

#include "MashmBool.h"
#include "MashmDefs.h"
#include "MashmFortranBindings.h"
#include "cmake_fortran_c_interface.h"

/**
 */
/** @struct Mashm
 *  @brief An opaque pointer to be used by the user to use the Mashm library
 *  @var Mashm::p
 *  Member 'p' contains a pointer to a private (internal) struct
 */
typedef struct {
  struct MashmPrivate* p;
} Mashm;

/**
 *  \memberof Mashm
 *  @brief Initialize the Mashm type
 *  @param in_mashm Pointer to Mashm type 
 *  @param in_comm MPI communicator
 */
void MashmInit(Mashm* in_mashm, MPI_Comm in_comm);

/**
 *  \memberof Mashm
 *  @brief Binding for initializing from Fortran
 *  @param in_mashm Mashm object
 *  @param f_comm Fortran MPI communicator
 *
 *  Special binding for initializing using a Fortran MPI Communicator.
 *  This gets mapped to mashmInit in Fortran.
 */
void MashmInitF2C(Mashm* in_mashm, MPI_Fint f_comm); 

/**
 *  \memberof Mashm
 *  @brief Returns the MPI communicator handed to the Mashm object
 *  @param Mashm object
 *  @return MPI communicator
 */
MPI_Comm MashmGetComm(const Mashm in_mashm);

/**
 *  \memberof Mashm
 *  @brief Return the number of MPI processes in the communicator
 *  @param Mashm object
 *  @return Number of MPI processes
 */

int MashmGetSize(const Mashm in_mashm);

/**
 *  \memberof Mashm
 *  @brief Get the rank of the calling process
 *  @param in_mashm Mashm object
 *  @return Rank of the calling process
 */
int MashmGetRank(const Mashm in_mashm);

/**
 *  \memberof Mashm
 *  @brief Set details about an MPI communication (symmetric)
 *  @param in_mashm Mashm object
 *  @param commIndex Index of the message to set
 *  @param pairRank Rank of the communication partner
 *  @param msgSize Size of the message
 */
void MashmSetComm(Mashm in_mashm, int commIndex, int pairRank, int msgSize);

/**
 *  \memberof Mashm
 *  @brief Allocate location for message data and set up pointers
 *  @param in_mashm Mashm object
 */
void MashmCommFinish(Mashm in_mashm);

/**
 *  \memberof Mashm
 *  @brief Print information about the Mashm object
 *  @param in_mashm Mashm object
 */
void MashmPrintInfo(const Mashm in_mashm);

/**
 *  \memberof Mashm
 *  @brief Determine if the message is on node or off node
 *  @param in_mashm Mashm object
 *  @param msgIndex
 *  @return True or False
 *  True if rank on shared memory node
 *  False if the rank is on shared memory node
 */
MashmBool MashmIsMsgOnNode(Mashm in_mashm, int msgIndex);

/**
 *  \memberof Mashm
 *  @brief Set the Communication method to use
 *  @param in_mashm Mashm object
 *  @param commType Communication type
 */
void MashmSetCommMethod(Mashm in_mashm, MashmCommType commType);

/**
 *  \memberof Mashm
 *  @brief Set the number of messages for the calling process
 *  @param in_mashm Mashm object
 *  @param numNeighbors Number of messages
 */
void MashmSetNumComms(Mashm in_mashm, int numMessages);

/**
 *  \memberof Mashm
 *  @brief Return the communication method
 *  @param in_mashm Mashm object
 *  @return The communication method
 */
MashmCommType MashmGetCommMethod(Mashm in_mashm);

/**
 *  \memberof Mashm
 *  @brief Determine if a process rank is on a shared memory node
 *  @param in_mashm Mashm object
 *  @param pairRank Message pair rank
 *  @return Whether the specified rank is on the same shared memory node as the calling process
 */
MashmBool MashmIsIntraNodeRank(Mashm in_mashm, int pairRank);

/**
 *  \memberof Mashm
 *  @brief Return the buffer pointer for a message
 *  @param in_mashm Mashm object
 *  @param msgIndex Index of the message
 *  @param sendReceive Whether to return the send or receive buffer pointer
 *  @return A pointer to the location of a send or receive buffer message
 */
double* MashmGetBufferPointer(Mashm in_mashm, int msgIndex, MashmSendReceive sendReceive);

/**
 *  \memberof Mashm
 *  @brief Retire (set to null) the buffer pointer for a message
 *  @param in_mashm Mashm object
 *  @param bufPtr The pointer to a message buffer
 */
void MashmRetireBufferPointer(Mashm in_mashm, double** bufPtr);

/**
 *  \memberof Mashm
 *  @brief Return the buffer pointer for a certain process
 *  @param in_mashm Mashm object
 *  @param destRank The rank of the message pair rank
 *  @param sendReceive Whether to return the send or receive buffer pointer
 *  @return A pointer to the location of a send or receive buffer for the message to the specified rank
 */
double* MashmGetBufferPointerForDest(Mashm in_mashm, int destRank, MashmSendReceive sendReceive);

/**
 *  \memberof Mashm
 *  @brief Begin the internode communication
 *  @param in_mashm Mashm object
 */
void MashmInterNodeCommBegin(Mashm in_mashm);

/**
 *  \memberof Mashm
 *  @brief Finish the internode communication
 *  @param in_mashm Mashm object
 */
void MashmInterNodeCommEnd(Mashm in_mashm);

/**
 *  \memberof Mashm
 *  @brief Begin the intranode communication
 *  @param in_mashm Mashm object
 */
void MashmIntraNodeCommBegin(Mashm in_mashm);

/**
 *  \memberof Mashm
 *  @brief Finish the intranode communication
 *  @param in_mashm Mashm object
 */
void MashmIntraNodeCommEnd(Mashm in_mashm);

/**
 *  \memberof Mashm
 *  @brief Destroy the Mashm object
 *  @param in_mashm Mashm object
 */
void MashmDestroy(Mashm* in_mashm);

#endif 
