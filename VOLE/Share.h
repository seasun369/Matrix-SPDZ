#include <libOTe/config.h>
#include <cryptoTools/Common/block.h>


#ifndef _SHARE
#define _SHARE

using namespace osuCrypto;


template<typename T, typename Ctx>
class Share
{
   public: 
   u64                         Row;
   u64                         Column;
   std::vector<std::vector<T>> MatrixShare;   //the share of the secret matrix
   block                       MatrixSeed;    // the seed to generate each row of the matrixshare
   std::vector<T>              leftMatrixShareMac;  // the share of the mac of the secret matrix
   std::vector<T>              rightMatrixShareMac;
   u64                         Sharetype;
   


   void configure(u64 RowLength, u64 ColumnLength)
   {
      Row=RowLength;
      Column=ColumnLength;
      Ctx mctx;
      mctx.resize(MatrixShare, Row);
      mctx.resize(rightMatrixShareMac, Row);
      mctx.resize(leftMatrixShareMac, Column);
      mctx.zero(leftMatrixShareMac.begin(),leftMatrixShareMac.end());
      PRNG prng(sysRandomSeed());
      MatrixSeed=prng.get();
      for(u64 i=0; i<Row;i++)
      {
        mctx.resize(MatrixShare[i], Column);
        mctx.zero(MatrixShare[i].begin(), MatrixShare[i].end());
      }
   }

};
#endif