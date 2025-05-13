#ifndef VOLE_SENDER_TENSOR_H
#define VOLE_SENDER_TENSOR_H


#include "Tools/VopeSender.h"
#include <cryptoTools/Network/Session.h>
#include <cryptoTools/Network/IOService.h>
#include <cryptoTools/Common/BitVector.h>
#include <libOTe_Tests/Common.h>
#include <coproto/Socket/BufferingSocket.h>
#include <coproto/Socket/Socket.h>
#include <libOTe/config.h>
#include "Tools/Coeff128.h"
#include <libOTe/Tools/CoeffCtx.h>
#include <libOTe/TwoChooseOne/SoftSpokenOT/SoftSpokenMalOtExt.h>

using namespace osuCrypto;




template<typename G, typename Ctx, u64 row>
void VoleSender_TensorVector(std::vector<std::vector<G>> &ret, std::vector<std::vector<G>> A, MultType type, bool debug, bool mal, Socket & ch)
{
   Ctx mctx;
   std::array<G,row> delta;
   SilentVopeSender<std::array<G,row>, G, CoeffCtxArrayPrime<G, row>> send;
   using VecF = typename Ctx::template Vec<std::array<G,row>>;
   u64 column=A[0].size();
   u64 rightMatrixcolumn=ret[0].size();
   VecF temp(rightMatrixcolumn);
   mctx.zero(temp.begin(),temp.end());
  PRNG prng(sysRandomSeed());
  for(u64 i=0;i<column;i++) 
  {
      for(u64 k=0;k<row;k++)
      {
        delta[k]=A[k][i];  
       // std::cout<<delta[k]<<std::endl; 
      }
       // std::cout<<"i="<<i<<"j="<<j<<std::endl; 
        VecF b(rightMatrixcolumn);
        cp::sync_wait(
        macoro::when_all_ready(send.silentSend(delta, b, prng, ch)));
        
      for(u64 h=0; h<rightMatrixcolumn; h++) 
      {
        for(u64 k=0;k<row;k++)
        {
          mctx.plus(temp[h][k], temp[h][k],b[h][k]);  
        // std::cout<<delta[k]<<std::endl; 
        }
      }
      //std::cout<<"i="<<i<<std::endl;
  }
    for(u64 i=0;i<rightMatrixcolumn;i++) 
    {
      for(u64 k=0;k<row;k++)
      {
        ret[k][i]=temp[i][k];  
       // std::cout<<delta[k]<<std::endl; 
      }
    }
  }

  #endif