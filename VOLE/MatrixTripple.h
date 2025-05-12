

//#include "Authentication.h"
#include <libOTe_Tests/Common.h>
#include <coproto/Socket/BufferingSocket.h>
#include "libOTe/Tools/Coproto.h"
#include <coproto/Socket/Socket.h>
#include <libOTe/config.h>
#include "Tools/Coeff128.h"
#include <libOTe/Tools/CoeffCtx.h>
//#include "Party.h"
//#include "Share.h"
//#include "Localcomputation.h"
#define Alice 1
#define Bob 2
using namespace osuCrypto;


template <typename G, typename Ctx, u64 leftMatrixRow>
void generate_matrix_triples(class Party<G,Ctx> party, std::vector<Share<G,Ctx>> &leftMatrix,std::vector<Share<G,Ctx>> &rightMatrix, 
std::vector<Share<G,Ctx>> &productMatrix, u64 leftMatrixColumn, u64 rightMatrixColumn, MultType type)
{

  MACORO_TRY{	
			if(leftMatrix.size()!=rightMatrix.size() || productMatrix.size()!=rightMatrix.size())
				throw std::runtime_error("The matrix tripple size does not match." LOCATION);
			} MACORO_CATCH(eptr) {
				std::rethrow_exception(eptr);
			}
 // std::cout<<"pass size computation"<<std::endl;
  u64 batchsize=leftMatrix.size();
  productMatrix.resize(batchsize);
  Ctx ctx;
  std::vector<Share<G, Ctx>> B(batchsize), C(batchsize),Matrixproduct(batchsize);
  //std::vector<Share<G, Ctx>> MaskMatrix(2*batchsize);
  Share<G, Ctx> sacrifice1, sacrifice2,sacrifice3,sacrifice4,sacrifice5;
  std::vector<std::vector<G>> TempMacA(1), TempMacB(1);
  ctx.resize(TempMacA[0],leftMatrixColumn);
  ctx.resize(TempMacB[0],rightMatrixColumn);
  for(u64 i=0;i<batchsize;i++)
  {
    //generate the matrix sharing <A[i]>, <B1[i]>, <B2[i]> and <mask1[i]>, <mask2[i]>, [C1[i]], [C2[i]] 
    Matrixproduct[i].configure(leftMatrixRow,rightMatrixColumn);
    productMatrix[i].configure(leftMatrixRow,rightMatrixColumn);
    //leftMatrix[i].configure(leftMatrixRow, leftMatrixColumn);
    //rightMatrix[i].configure(leftMatrixColumn, rightMatrixColumn);
    //B[i].configure(leftMatrixColumn, rightMatrixColumn);
    C[i].configure(leftMatrixColumn, rightMatrixColumn);
    //C[i].configure(leftMatrixColumn, rightMatrixColumn);
    generate_Matrix_and_MAC_with_seed(party, leftMatrix[i], leftMatrixRow, leftMatrixColumn, type);
    generate_Matrix_and_MAC_with_seed(party, rightMatrix[i],leftMatrixColumn, rightMatrixColumn, type);
    generate_Matrix_and_MAC_with_seed(party, B[i], leftMatrixColumn, rightMatrixColumn, type);
    //generate the additive sharing of C1[i]=A[i]B1[i] and C2[i]=A[i]B2[i]
    generate_Matrix_product_share<u64, Ctx, leftMatrixRow>(party, leftMatrix[i].MatrixShare, rightMatrixColumn, rightMatrix[i].MatrixSeed,  productMatrix[i].MatrixShare, type);
    generate_Matrix_product_share<u64, Ctx, leftMatrixRow>(party, leftMatrix[i].MatrixShare, rightMatrixColumn, B[i].MatrixSeed, Matrixproduct[i].MatrixShare,  type);
    //generate the mac of [C1[i]] and [C2[i]] by computing (party1_MAC(A[i]+party_2MAC(A[i]))(party1_share(B1[i])+party2_Share(B1[i])) 
    ctx.copy(leftMatrix[i].leftMatrixShareMac.begin(),leftMatrix[i].leftMatrixShareMac.begin()+leftMatrixColumn, TempMacA[0].begin());
    generate_Matrix_product_share<u64, Ctx, 1>(party, TempMacA, rightMatrixColumn, B[i].MatrixSeed, TempMacB,  type);
    ctx.vector_plus(Matrixproduct[i].leftMatrixShareMac, Matrixproduct[i].leftMatrixShareMac, TempMacB[0]);
    generate_Matrix_product_share<u64, Ctx, 1>(party, TempMacA, rightMatrixColumn, rightMatrix[i].MatrixSeed, TempMacB,  type);
    ctx.vector_plus(productMatrix[i].leftMatrixShareMac, productMatrix[i].leftMatrixShareMac, TempMacB[0]);

  }


   
    // using sacrifice1,2,3,4 to get authenticated sharing <A[i]>, <B1[i]>, <B2[i]> and <mask1[i]>, <mask2[i]>
    generate_Matrix_and_MAC_with_seed(party, sacrifice1, leftMatrixRow, leftMatrixColumn, type);
    generate_Matrix_and_MAC_with_seed(party, sacrifice2, leftMatrixColumn, rightMatrixColumn, type);
    generate_Matrix_and_MAC_with_seed(party, sacrifice3, leftMatrixColumn, rightMatrixColumn, type);
    generate_Matrix_and_MAC_with_seed(party, sacrifice4, leftMatrixRow,rightMatrixColumn, type);
    generate_Matrix_and_MAC_with_seed(party, sacrifice5, leftMatrixRow,rightMatrixColumn, type);
   

    authentication(party, sacrifice1, leftMatrix);
    authentication(party, sacrifice2, B);
    authentication(party, sacrifice3, rightMatrix);
    authentication(party, sacrifice4, productMatrix);
    authentication(party, sacrifice5, Matrixproduct);
    
    std::vector<std::vector<G>> recvShare(leftMatrixColumn);
    Share<G,Ctx> TempShare;
    TempShare.configure(leftMatrixRow, rightMatrixColumn);
    
    for(u64 i=0;i<leftMatrixColumn; i++)
    {
      ctx.resize(recvShare[i], rightMatrixColumn);
    }
    
    for(u64 i=0;i<batchsize;i++)
  {
    //using authenticated sharing <mask[i]> to obtain <C1[i]> and <C2[i]>
    G publicCoin;
    //Compute <B1>-x<B2> and open it to D=B1-xB2
    // Then local compute <C1>-x<C2>-<A>D which should be open to 0;
    ctx.fromBlock(publicCoin, party.random_coin());
    local_scaler(B[i], publicCoin, B[i], ctx);
    local_minus(B[i], rightMatrix[i], B[i], ctx);
    open_matrix(party,B[i].MatrixShare,recvShare);
    //ctx.matrix_plus(recvShare,recvShare, B[i].MatrixShare);
    check_mac(party, B[i].leftMatrixShareMac, recvShare);
    local_matrix_product(TempShare, leftMatrix[i], recvShare, ctx);
    local_scaler(Matrixproduct[i], publicCoin, Matrixproduct[i], ctx);
    local_minus(Matrixproduct[i], productMatrix[i], Matrixproduct[i], ctx);
    local_minus(Matrixproduct[i], Matrixproduct[i], TempShare,ctx);
    open_matrix(party, Matrixproduct[i].MatrixShare, TempShare.MatrixShare);
    check_mac(party, Matrixproduct[i].leftMatrixShareMac, TempShare.MatrixShare);
    //ctx.matrix_plus(TempShare.MatrixShare, TempShare.MatrixShare,Matrixproduct[i].MatrixShare);
    for(u64 j=0;j<TempShare.MatrixShare.size();j++)
    {
      for(u64 k=0;k<TempShare.MatrixShare[0].size();k++)
      {
          MACORO_TRY{
        
        if (TempShare.MatrixShare[j][k]!= 0)
          throw std::runtime_error("Shares are corrupted. Abort the protocol" LOCATION);

        } MACORO_CATCH(eptr) {
          std::rethrow_exception(eptr);
        }
      }
    }

  }
   
}

template<typename G, typename Ctx, u64 leftMatrixRow>
void generate_Matrix_product_share(class Party<G,Ctx> party, std::vector<std::vector<G>> leftMatrix, u64 rightMatrixColumn, block MatrixSeed, std::vector<std::vector<G>> &ret, MultType type)
{

  u64 leftMatrixColumn=leftMatrix[0].size();
   //std::cout<<"resize 0 clear"<<std::endl;
  std::vector<std::vector<G>> TempMatrix(leftMatrixColumn);
  std::vector<std::vector<G>> ret1(leftMatrixRow), ret2(leftMatrixRow);
  Ctx ctx;
  ctx.resize(ret,leftMatrixRow);
  //std::cout<<"resize 1 clear"<<std::endl;
  for(u64 i=0;i<leftMatrixRow;i++)
  {
    ctx.resize(ret[i], rightMatrixColumn);
    ctx.resize(ret1[i], rightMatrixColumn);
    ctx.resize(ret2[i], rightMatrixColumn);
  }
  for(u64 j=0;j<leftMatrixColumn;j++)
  {
    ctx.resize(TempMatrix[j], rightMatrixColumn);
  }
  //wait for vole
  if(party.Party_Id==Alice)
   {
     VoleSender_TensorVector<G, Ctx, leftMatrixRow>(ret1, leftMatrix, type, false, false, party.chl);
     VoleReceiver_TensorVector<G, Ctx, leftMatrixRow>(ret2,  TempMatrix, MatrixSeed, type,false, false, party.chll);
   }
   else if(party.Party_Id==Bob)
   {
      VoleReceiver_TensorVector<G, Ctx, leftMatrixRow>(ret2,  TempMatrix, MatrixSeed,type,false, false, party.chl);
      VoleSender_TensorVector<G, Ctx,  leftMatrixRow>(ret1, leftMatrix, type, false, false, party.chll);
   }
    //std::cout<<"rightmatrix="<<TempMatrix[0][0]<<std::endl;
     for(u64 i=0;i<leftMatrixRow;i++)
     {
        ctx.vector_matrix_product(ret[i], leftMatrix[i], TempMatrix);
     }
      //  std::cout<<"vector matrix product clear"<<std::endl;
        ctx.matrix_plus(ret, ret, ret2);
        ctx.matrix_minus(ret, ret, ret1);
}
