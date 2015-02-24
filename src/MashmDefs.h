#include "MashmConfig.h"

#ifndef MASHM_DEFS_H
#define MASHM_DEFS_H

typedef enum { MASHM_COMM_STANDARD, MASHM_COMM_INTRA_MSG, MASHM_COMM_INTRA_SHARED, MASHM_COMM_MIN_AGG } MashmCommType;

typedef enum { MASHM_SEND, MASHM_RECEIVE } MashmSendReceive;

typedef enum { MASHM_MIN_AGG_ROUND_ROBIN, MASHM_MIN_AGG_ROOT_PROC } MashmMinAggType;

#endif
