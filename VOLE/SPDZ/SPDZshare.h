#include <libOTe/config.h>
#include <cryptoTools/Common/block.h>



#ifndef _SPDZ_SHARE
#define _SPDZ_SHARE

using namespace osuCrypto;

template<typename T, typename Ctx>
class SpdzShare
{
   public: 
   std::vector<T>              share;
   std::vector<T>              shareMac;
   block                       shareSeed;

   void configure(u64 batchSize)
   {
      Ctx ctx;
      share.resize(batchSize);
      shareMac.resize(batchSize);
      ctx.zero(share.begin(),share.end());
      ctx.zero(shareMac.begin(),shareMac.end());
   }
};

#endif

