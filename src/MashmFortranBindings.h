#ifndef MASHM_FORT_BIND_H
#define MASHM_FORT_BIND_H

#define MashmInitF2C FCI_GLOBAL(mashminit,MASHMINIT)
//#define MashmIsMsgOnNode FCI_GLOBAL(mashmismsgonnodec,MASHMISMSGONNODEC)
//#define MashmSetComm FCI_GLOBAL(mashmsetcommc,MASHMSETCOMMc)
#define MashmGetBufferPointer FCI_GLOBAL(mashmgetbufferpointerc,MASHMGETBUFFERPOINTERC)
//!#define MashmGetBufferPointer2 FCI_GLOBAL(mashmgetbufferpointer2c,MASHMGETBUFFERPOINTER2C)
#define MashmRetireBufferPointer FCI_GLOBAL(mashmretirebufferpointerc,MASHMRETIREBUFFERPOINTERC)

#if 0
#define MashmGetComm FCI_GLOBAL(mashmgetcomm,MASHMGETCOMM)

#define MashmGetSize FCI_GLOBAL(mashmgetsize,MASHMGETSIZE)

#define MashmGetRank FCI_GLOBAL(mashmgetrank,MashmGetRank)

#define MashmGetCommRank FCI_GLOBAL(mashmgetcommrank,MASHMGETCOMMRANK)

#define MashmGetCommSize FCI_GLOBAL(mashmgetcommsize,MASHMGETCOMMSIZE)

#define MashmCommFinish FCI_GLOBAL(mashmcommfinish,MASHMCOMMFINISH)

#define MashmPrintInfo FCI_GLOBAL(mashmprintinfo,MASHMPRINTINFO)


#define MashmSetCommMethod FCI_GLOBAL(mashmsetcommmethod,MASHMSETCOMMMETHOD)

#define MashmSetNumComms FCI_GLOBAL(mashmsetnumcomms,MASHMSETNUMCOMMS)

#define MashmGetCommMethod FCI_GLOBAL(mashmgetcommmethod,MASHMGETCOMMMETHOD)

#define MashmNumMpiMsgs FCI_GLOBAL(mashmnummpimsgs,MASHMNUMMPIMSGS)

#define MashmNumIntraNodeMsgs FCI_GLOBAL(mashmnumintranodemsgs,MASHMNUMINTRANODEMSGS)

#define MashmIsIntraNodeRank FCI_GLOBAL(mashmisintranoderank,MASHMISINTRANODERANK)

#define MashmCalcMsgBufferSize FCI_GLOBAL(mashmcalcmsgbuffersize,MASHMCALCMSGBUFFERSIZE)

#define MashmCalcNumConnectedNodes FCI_GLOBAL(mashmcalcnumconnectednodes,MASHMCALCNUMCONNECTEDNODES)

#define MashmSetupStandardComm FCI_GLOBAL(mashmsetupstandardcomm,MASHMSETUPSTANDARDCOMM)

#define MashmGetBufferPointerForDest FCI_GLOBAL(mashmgetbufferpointerfordest,MASHMGETBUFFERPOINTERFORDEST)

#define MashmInterNodeCommBegin FCI_GLOBAL(mashminternodecommbegin,MASHMINTERNODECOMMBEGIN)

#define MashmInterNodeCommEnd FCI_GLOBAL(mashminternodecommend,MASHMINTERNODECOMMEND)

#define MashmIntraNodeCommBegin FCI_GLOBAL(mashmintranodecommbegin,MASHMINTRANODECOMMBEGIN)

#define MashmIntraNodeCommEnd FCI_GLOBAL(mashmintranodecommend,MASHMINTRANODECOMMEND)

#define MashmDestroy FCI_GLOBAL(mashmdestroy,MASHMDESTROY)

#define MashmPrintCommCollection FCI_GLOBAL(mashmprintcommcollection,MASHMPRINTCOMMCOLLECTION)
#endif
#endif
