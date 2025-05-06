#include <libOTe_Tests/Common.h>
#include "libOTe/TwoChooseOne/SoftSpokenOT/SoftSpokenShOtExt.h"
#include <coproto/Socket/BufferingSocket.h>
#include <coproto/Socket/Socket.h>
#include "cryptoTools/Common/BitVector.h"
#include <libOTe/config.h>
#include "Tools/Coeff128.h"
//#include <libOTe/Tools/CoeffCtx.h>
#include "Party.h"
#include "Share.h"
#include "Authentication.h"
#include "MatrixTranspose.h"
#include "Localcomputation.h"
#include "MatrixTripple.h"
//#include "SPDZ/BatchOLE.h"
//#include "SPDZ/SPDZprotocol.h"
//#include "SPDZ/SPDZshare.h"
#define Alice 1
#define Bob 2
using namespace osuCrypto;





int main()
{

  u64 Maxlength=10000;
  CoeffCtxIntegerPrime ctx;
  const u64 column=128;
  const u64 row=128;
  const u64 rightcolumn=128;
  u64 batchSize=5000;
  u64 batch=10;
  auto chl1 = cp::LocalAsyncSocket::makePair();
  auto chl2 = cp::LocalAsyncSocket::makePair();
  Party<u64, CoeffCtxIntegerPrime> party1, party2;
  PRNG prng(sysRandomSeed());
  block commonseed(0,0);
  party1.configure(Alice, Maxlength, chl1[0], chl2[0], commonseed);
  party2.configure(Bob, Maxlength, chl1[1],chl2[1], commonseed);
  std::vector<Share<u64, CoeffCtxIntegerPrime>> LeftMatrix1(batch), LeftMatrix2(batch), LeftMatrix1_t(batch), LeftMatrix2_t(batch);
  std::vector<Share<u64, CoeffCtxIntegerPrime>> RandomMatrix1(batch), RandomMatrix2(batch), RandomMatrix1_t(batch), RandomMatrix2_t(batch);
  std::vector<Share<u64, CoeffCtxIntegerPrime>> RightMatrix1(batch), RightMatrix2(batch);
  std::vector<Share<u64, CoeffCtxIntegerPrime>> MatrixProduct1(batch), MatrixProduct2(batch);
  Share<u64, CoeffCtxIntegerPrime> matrixproduct1, matrixproduct2; 
  Share<u64, CoeffCtxIntegerPrime> sacrifice1, sacrifice2; 
  SoftSpokenMalOtSender baseSend;
  SoftSpokenMalOtReceiver baseRecv;
 // SpdzShare<u64, CoeffCtxIntegerPrime> shareA1, shareA2, shareB1, shareB2, shareC1, shareC2;
  std::vector<u64> b(batchSize),d(batchSize),a(batchSize),c(batchSize);
 // std::vector<double> a(batchSize),c(batchSize);
  u64 sum0,sum1;
  double sum;
  //std::cout<<pow<<std::endl;
 /* std::cout<<a[0]<<std::endl;
  std::cout<<b[0]<<std::endl;
  convert_to_field<u64, CoeffCtxIntegerPrime>(a,b,precision);
    std::cout<<a[0]<<std::endl;
    std::cout<<b[0]<<std::endl;
  convert_to_field<u64, CoeffCtxIntegerPrime>(c,d,precision);
  */

  std::thread t0([&]() 
  {
     //Inner_Product_Sender<u64, CoeffCtxIntegerPrime>(b, chl1[0], sum0);
     //Batch_OLE_Sender<u64, CoeffCtxIntegerPrime>(d, c, baseSend, chl1[0]);
     generate_matrix_triples<u64,CoeffCtxIntegerPrime,row>(party1, LeftMatrix1, RightMatrix1, MatrixProduct1,  column, rightcolumn, DefaultMultType);
    // generate_Matrix_and_MAC_with_seed<u64,CoeffCtxIntegerPrime>(party1, LeftMatrix1[0],  row,  column, DefaultMultType);
     //generate_Matrix_and_MAC_with_seed<u64,CoeffCtxIntegerPrime>(party1, LeftMatrix1[1],  row,  column, DefaultMultType);
     generate_matrix_transpose<u64,CoeffCtxIntegerPrime>(party1, LeftMatrix1, LeftMatrix1_t, RandomMatrix1, RandomMatrix1_t, row, column, rightcolumn,  DefaultMultType);
    std::cout<<"party1 Bytes send="<<(party1.chl.bytesSent()+party1.chll.bytesSent())/(float)(1024*1024)<<"MB"<<std::endl;
  });
  
  std::thread t1([&]()
  {
      //Batch_OLE_Receiver<u64, CoeffCtxIntegerPrime>(b, a, baseRecv, chl1[1]);
     //generate_Beaver_tripple(party2, shareA2, shareB2, shareC2,batchSize,DefaultMultType);
     //generate_spdz_share_with_mac(party2, share2, batchSize, DefaultMultType);
     //authentication_spdz(party2,share2);
     generate_matrix_triples<u64,CoeffCtxIntegerPrime,row>(party2, LeftMatrix2, RightMatrix2, MatrixProduct2,  column, rightcolumn, DefaultMultType);
     //generate_Matrix_and_MAC_with_seed<u64,CoeffCtxIntegerPrime>(party2, LeftMatrix2[0],  row,  column, DefaultMultType);
     //generate_Matrix_and_MAC_with_seed<u64,CoeffCtxIntegerPrime>(party2, LeftMatrix2[1],  row,  column, DefaultMultType);
     generate_matrix_transpose<u64,CoeffCtxIntegerPrime>(party2, LeftMatrix2, LeftMatrix2_t, RandomMatrix2, RandomMatrix2_t, row, column, rightcolumn,  DefaultMultType);
     //std::cout<<"party2 Bytes send="<<party2.chl.bytesSent()/(float)(1024*1024)+party2.chll.bytesSent()/(float)(1024*1024)<<"MB"<<std::endl;
     //Inner_Product_Receiver<u64, CoeffCtxIntegerPrime>(d, chl1[1], sum1);
   // for(u64 i=0;i<c.size();i++)
     // {
     //   ctx.minus(c[i], 0, c[i]);
     // }
     //MatrixTripple<u64,CoeffCtxIntegerPrime,row>(party2, LeftMatrix2, RightMatrix2, MatrixProduct2,  column, rightcolumn, DefaultMultType);
    //std::cout<<"Here2"<<std::endl;
     //generate_Matrix_and_MAC_with_seed(party2, sacrifice2,  row,  column, DefaultMultType);
     //authentication(party2, sacrifice2, LeftMatrix2);
  });
  t0.join();
  t1.join();
  //convert_to_float_element<u64, CoeffCtxIntegerPrime>(sum0, sum, 0);
  //std::cout<<sum<<std::endl;
 /*
  u64 temp1, temp2;
  ctx.plus(temp1, shareA1.share[1], shareA2.share[1]);
  ctx.plus(temp2, shareB1.share[1], shareB2.share[1]);
  ctx.mul(temp1,temp1,temp2);
  ctx.plus(temp2,shareC1.share[1],shareC2.share[1]);
  std::cout<<"temp1="<<temp1<<std::endl;
  std::cout<<"temp2="<<temp2<<std::endl;
  ctx.plus(temp1,party1.SpdzKey,party2.SpdzKey);
  ctx.mul(temp1,temp1,temp2);
  std::cout<<"temp1="<<temp1<<std::endl;
  ctx.plus(temp2,shareC1.shareMac[1],shareC2.shareMac[1]);
  std::cout<<"temp2="<<temp2<<std::endl;
  */
  
  
  /*    
  u64 temp;
  u64 spdzKey;
  ctx.plus(spdzKey, party1.SpdzKey, party2.SpdzKey);
  for(u64 i=0;i<batchSize;i++)
  {
    ctx.plus(temp, share1.share[i], share2.share[i]);
    ctx.mul(temp, temp, spdzKey);
    std::cout<<temp<<std::endl;
    ctx.plus(temp, share1.shareMac[i], share2.shareMac[i]);
    std::cout<<temp<<std::endl;
  }
  */
  //std::vector<Share<u64, CoeffCtxIntegerPrime>> share1(5); 
  //std::vector<Share<u64, CoeffCtxIntegerPrime>> share2(5); 
/*
  std::thread t0([&]() 
  {
    
    generate_Matrix_and_MAC_with_seed(party1, sacrifice1,  column,  rightcolumn, DefaultMultType);
    for(u64 i=0;i<5;i++)
    {
      generate_Matrix_and_MAC_with_seed(party1, share1[i],  column,  rightcolumn, DefaultMultType);
    }
    //  authentication(party1, sacrifice1, share1);
  });
  
  std::thread t1([&]()
  {
    //generate_Matrix_and_MAC_with_seed(party2, sacrifice2,  column,  rightcolumn, DefaultMultType);
    for(u64 i=0;i<5;i++)
    {
      generate_Matrix_and_MAC_with_seed(party2, share2[i],  column,  rightcolumn, DefaultMultType);
    }
      //authentication(party2, sacrifice2, share2);
  });
  
  
 
  t0.join();
  t1.join();  
  
  */
 /*
   std::cout<<"before check"<<std::endl;
   ctx.vector_plus(party1.MatrixKey, party1.MatrixKey, party2.MatrixKey);
   ctx.matrix_plus(MatrixProduct1[0].MatrixShare, MatrixProduct1[0].MatrixShare, MatrixProduct2[0].MatrixShare);
   //ctx.vector_plus(LeftMatrix1[0].leftMatrixShareMac,LeftMatrix1[0].leftMatrixShareMac, LeftMatrix2[0].leftMatrixShareMac);
   std::cout<<"matrix add clear"<<std::endl;
   //ctx.matrix_plus(share3.MatrixShare, share3.MatrixShare, share4.MatrixShare);
   u64 sum=0;
   u64 temp;
  for(u64 i=0; i<row;i++)
  {
    
    ctx.mul(temp, MatrixProduct1[0].MatrixShare[i][0], party1.MatrixKey[i]);
    ctx.plus(sum, sum, temp);
  }
      ctx.plus(temp, MatrixProduct1[0].leftMatrixShareMac[0], MatrixProduct2[0].leftMatrixShareMac[0]);
      std::cout<<sum<<std::endl;
      std::cout<<temp<<std::endl; 
*/
/*
ctx.matrix_plus(MatrixProduct1[0].MatrixShare, MatrixProduct1[0].MatrixShare, MatrixProduct2[0].MatrixShare);
ctx.matrix_plus(LeftMatrix1[0].MatrixShare, LeftMatrix1[0].MatrixShare, LeftMatrix2[0].MatrixShare);
ctx.matrix_plus(RightMatrix1[0].MatrixShare, RightMatrix1[0].MatrixShare, RightMatrix2[0].MatrixShare);
ctx.matrix_product(MatrixProduct2[0].MatrixShare, LeftMatrix1[0].MatrixShare, RightMatrix1[0].MatrixShare);
std::cout<<MatrixProduct1[0].MatrixShare[1][1]<<std::endl;
std::cout<<MatrixProduct2[0].MatrixShare[1][1]<<std::endl;
*/


      
 // std::vector<u64> globalkey(row);
 
 /*
  for(u64 j=0;j<column;j++)
  {
    u64 sum=0;
    for(u64 i=0;i<row;i++)
    {
      u64 temp;
      ctx.plus(globalkey[i], party1.MatrixKey[i],party2.MatrixKey[i]);
      ctx.plus(temp,share2.MatrixShare[i][j], share1.MatrixShare[i][j]);
      ctx.mul(temp,temp, globalkey[i]);
      ctx.plus(sum,sum,temp);
    }
      u64 mac;
      ctx.plus(mac, share1.leftMatrixShareMac[j], share2.leftMatrixShareMac[j]);
      std::cout<<mac<<std::endl;
      std::cout<<sum<<std::endl; 
      
  }
  */
  
      //ctx.matrix_plus(share1.MatrixShare, share1.MatrixShare, share2.MatrixShare);
      //ctx.matrix_scaler(share1.MatrixShare, share1.MatrixShare[0][0], share2.MatrixShare);

  /*
   PRNG prng(sysRandomSeed());
   CoeffCtxIntegerPrime ctx;
   const u64 leftMatrixrow=128;
   const u64 rightMatrixcolumn=1024;
   const u64 middleMatrix=128;
   std::vector<std::vector<u64>> A(leftMatrixrow);
   for(u64 i=0;i<A.size();i++)
   {
     ctx.resize(A[i], middleMatrix);
     for(u64 j=0;j<middleMatrix;j++)
     {
        ctx.fromBlock(A[i][j],prng.get<block>());
     }
   }
   std::vector<std::vector<u64>> C(middleMatrix);
   for(u64 i=0;i<C.size();i++)
   {
     ctx.resize(C[i], rightMatrixcolumn);
   }
   std::vector<std::vector<u64>> ret1(leftMatrixrow), ret2(leftMatrixrow);
   for(u64 i=0;i<ret1.size();i++)
   {
     ctx.resize(ret1[i], rightMatrixcolumn);
     ctx.resize(ret2[i], rightMatrixcolumn);
   }
   auto chls = cp::LocalAsyncSocket::makePair();
  std::thread t0([&]() 
  {
    
   VoleSender_TensorVector<u64,CoeffCtxIntegerPrime, leftMatrixrow>(ret1, A, DefaultMultType, false, false, chls[0]);

  });
  
  std::thread t1([&]()
  {
   VoleReciver_TensorVector<u64,CoeffCtxIntegerPrime,leftMatrixrow>(ret2, C, prng, DefaultMultType, false, false, chls[1]);
  });
  t0.join();
  t1.join();

  for(u64 k=0;k<A.size();k++)
  {
    u64 sum=0;
    u64 temp;
    for(u64 i=0;i<C.size();i++)
    {
      ctx.mul(temp, A[k][i], C[i][0]);
      ctx.plus(sum, sum, temp);
      //std::cout<<sum<<std::endl;
    }
    ctx.minus(temp, ret2[k][0],ret1[k][0]);
    std::cout<<ctx.eq(temp,sum)<<std::endl;
    std::cout<<"k="<<k<<"  sum="<<sum<<std::endl;
    std::cout<<"minus="<<temp<<std::endl;
  }
  */
/*
  u64 Maxlength=100;
  CoeffCtxIntegerPrime ctx;
  auto chls = cp::LocalAsyncSocket::makePair();
  Party<u64, CoeffCtxIntegerPrime> party1, party2;
  party1.configure(Alice, Maxlength, chls[0]);
  party2.configure(Bob, Maxlength, chls[1]);
  */
 //CoeffCtxIntegerPrime ctx;
  /*std::vector<Share<u64, CoeffCtxIntegerPrime>> share;
  share.resize(10);
  for(u64 i=0;i<10;i++)
  {
    share[i].configure(100,120);

  }
 
 std::cout<<share[0].MatrixSeed<<std::endl;
 std::cout<<share[0].leftMatrixShareMac.size()<<std::endl;
 std::cout<<share.end()-share.begin()<<std::endl;
  //std::cout<<party1.MatrixKey[0]<<std::endl;
  //std::cout<<party2.MatrixKey[50]<<std::endl;
*/
}