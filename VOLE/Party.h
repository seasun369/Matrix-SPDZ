
#include <coproto/Socket/Socket.h>
#include <libOTe/config.h>
#include <cryptoTools/Common/block.h>
using namespace osuCrypto;

#ifndef _PARTY
#define _PARTY

template<typename T, typename Ctx>
class Party
{
   public: 
   u64            Party_Id;
   std::vector<T> MatrixKey;
   T              SpdzKey;
   block          seed;
   Socket         chl;
   Socket         chll;    
   

   void configure(u64 party, u64 Maxlength, Socket &ch1, Socket &ch2, block commonseed)
   {
      Party_Id=party;
      seed=commonseed;
      PRNG prng1(sysRandomSeed());
      Ctx mctx;
      mctx.resize(MatrixKey, Maxlength);
      for(u64 i=0;i<MatrixKey.size();i++)
      {
        mctx.fromBlock(MatrixKey[i], prng1.get<block>());
       // std::cout<<"matrix="<<MatrixKey[i]<<"i="<<i<<std::endl;
      }
        mctx.fromBlock(SpdzKey, prng1.get<block>());
        chl=ch1;   
        chll=ch2;   
   }
   block random_coin()
   {
     PRNG prng(seed);
     seed=prng.get();
     return seed;
   }

};
#endif