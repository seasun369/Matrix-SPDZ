#include "Tools/VopeReceiver.h"
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


template<typename G, typename Ctx, u64 leftMatrixrow>
void VoleReceiver_TensorVector(std::vector<std::vector<G>> &ret, std::vector<std::vector<G>> &C, block seed, MultType type, bool debug, bool mal, Socket & ch)
{
   Ctx mctx;
   //block seed(100,0);
   PRNG prng(seed);
   PRNG prng1(sysRandomSeed());
   u64 column=C[0].size();
   u64 row=C.size();
   SilentVopeReceiver<std::array<G,leftMatrixrow>, G, CoeffCtxArrayPrime<G,leftMatrixrow>> recv;
   using VecF = typename Ctx::template Vec<std::array<G,leftMatrixrow>>;
   using VecG = typename Ctx::template Vec<G>;
   VecF temp(column);
   mctx.zero(temp.begin(),temp.end());
   std::vector<block> seeds(row);
  // PRNG prng(sysRandomSeed());
  for(u64 i=0;i<row;i++) 
  {
    seeds[i]= prng.get();
    PRNG setupSeed(seeds[i]);
       // std::cout<<"i="<<i<<"j="<<j<<std::endl; 
    VecF a(column);
    VecG c(column);
    cp::sync_wait(
    macoro::when_all_ready(recv.silentReceive(c, a, prng1, setupSeed, ch)));
    // std::cout<<"i="<<i<<std::endl; 
    mctx.copy(c.begin(),c.end(),C[i].begin());
      for(u64 h=0;h<leftMatrixrow;h++) 
      {
        for(u64 k=0;k<column;k++)
        {
          mctx.plus(temp[k][h],temp[k][h],a[k][h]);  
        // std::cout<<delta[k]<<std::endl; 
        }
      }

  }
  for(u64 i=0;i<leftMatrixrow;i++) 
    {
      for(u64 k=0;k<column;k++)
      {
        ret[i][k]=temp[k][i];  
       // std::cout<<delta[k]<<std::endl; 
      }
    }
    //std::cout<<ret[0][0]<<std::endl;

}