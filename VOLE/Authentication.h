#include "VoleReceiver_TensorVector.h"
#include "VoleSender_TensorVector.h"
#include <libOTe_Tests/Common.h>
#include <coproto/Socket/BufferingSocket.h>
#include "libOTe/Tools/Coproto.h"
#include <coproto/Socket/Socket.h>
#include <libOTe/config.h>
#include "Tools/Coeff128.h"
#include <libOTe/Tools/CoeffCtx.h>
#include <libOTe/TwoChooseOne/SoftSpokenOT/SoftSpokenMalOtExt.h>
#include "cryptoTools/Crypto/Commit.h"
//#nclude "Party.h"
//#include "Share.h"
#define Alice 1
#define Bob 2
using namespace osuCrypto;


template<typename G, typename Ctx>
void generate_Matrix_and_MAC_with_seed(class Party<G,Ctx> party, Share<G, Ctx> &share, u64 row, u64 column, MultType type)
{
    Ctx ctx;
    std::vector<std::vector<G>> ret1(1), ret2(1), matrixKey(1);
    PRNG prng(sysRandomSeed());
    share.MatrixSeed=prng.get();
    ctx.resize(ret1[0],column);
    ctx.resize(ret2[0],column);
    ctx.resize(matrixKey[0],row);
    share.configure(row,column);
    ctx.copy(party.MatrixKey.begin(),party.MatrixKey.begin()+row, matrixKey[0].begin());
   if(party.Party_Id==Alice)
   {
       VoleSender_TensorVector<G, Ctx,1>(ret1, matrixKey, type, false, false, party.chl);
       VoleReceiver_TensorVector<G, Ctx, 1>(ret2, share.MatrixShare, share.MatrixSeed,type,false, false, party.chll);
   }
   else if(party.Party_Id==Bob)
   {
      VoleReceiver_TensorVector<G, Ctx, 1>(ret2, share.MatrixShare, share.MatrixSeed,type,false, false, party.chl); 
      VoleSender_TensorVector<G, Ctx,1>(ret1, matrixKey, type, false, false, party.chll);
      
   } 
      
      ctx.vector_matrix_product(share.leftMatrixShareMac, matrixKey[0], share.MatrixShare);
      for(u64 i=0;i<column;i++)
      {
          ctx.plus(share.leftMatrixShareMac[i], share.leftMatrixShareMac[i], ret2[0][i]);
          ctx.minus(share.leftMatrixShareMac[i], share.leftMatrixShareMac[i], ret1[0][i]);
          //std::cout<<share.leftMatrixShareMac.size()<<std::endl;
      }
   
}

template<typename G, typename Ctx>
void generate_MAC_with_authenticated_sharing(class Party<G,Ctx> party, Share<G, Ctx> &share, Share<G, Ctx> authenticated_random_sharing, MultType type)
{
    Ctx ctx;
    u64 row=share.MatrixShare.size();
    u64 column=share.MatrixShare[0].size();
    Share<G, Ctx> MatrixForOpenning, TempMatrix; 
    std::vector<G> matrixKey(row);
    MatrixForOpenning.configure(row, column);
    TempMatrix.configure(row, column);
    share.leftMatrixShareMac.resize(column);
    ctx.copy(party.MatrixKey.begin(),party.MatrixKey.begin()+row, matrixKey.begin());
    //locally compute [M]+[R]=[M+R]
    ctx.matrix_plus(MatrixForOpenning.MatrixShare,share.MatrixShare,authenticated_random_sharing.MatrixShare);
    open_matrix(party, MatrixForOpenning.MatrixShare, TempMatrix.MatrixShare);
      //ctx.matrix_plus(MatrixForOpenning,MatrixForOpenning,TempMatrix);
    //send the share of [M+R] to reconstruct M+R and generate the new mac [v(M+R)]-[vR]
    ctx.vector_matrix_product(share.leftMatrixShareMac,matrixKey,TempMatrix.MatrixShare);
    ctx.vector_minus(share.leftMatrixShareMac, share.leftMatrixShareMac, authenticated_random_sharing.leftMatrixShareMac);
}

template<typename G, typename Ctx>
void authentication(class Party<G,Ctx> party, Share<G, Ctx> &sacrifice, const std::vector<Share<G, Ctx>> share_with_mac)
{
   Ctx ctx;
   u64 numberofShare=share_with_mac.size();
   u64 row=sacrifice.MatrixShare.size();
   u64 column=sacrifice.MatrixShare[0].size();
   //std::cout<<"here"<<std::endl;
   std::vector<G> matrixKey(row);
   Share<G, Ctx> ShareTemp;
   ctx.copy(party.MatrixKey.begin(),party.MatrixKey.begin()+row, matrixKey.begin());
   ShareTemp.configure(row,column);
   for(u64 i=0;i<numberofShare;i++)
   {
      G publicCoin;
      ctx.fromBlock(publicCoin, party.random_coin());
      local_scaler(ShareTemp, publicCoin, share_with_mac[i],ctx);
      local_plus(sacrifice, ShareTemp, sacrifice,ctx);
   }
  
    open_matrix(party,sacrifice.MatrixShare,ShareTemp.MatrixShare);
    check_mac(party, sacrifice.leftMatrixShareMac, ShareTemp.MatrixShare);
}


template<typename G, typename Ctx>
void open_matrix(class Party<G,Ctx> party, std::vector<std::vector<G>> &send_share, std::vector<std::vector<G>> &receive_share)
{
   Ctx ctx;
   u64 row=send_share.size();
   u64 column=send_share[0].size();
   std::vector<G> buff(row*column),buff1(row*column);
   std::vector<u8> sendbuff,recvbuff;
   sendbuff.resize(buff.size() * ctx.template byteSize<G>());
   recvbuff.resize(buff1.size() * ctx.template byteSize<G>());
   u64 k=0;
   for(u64 i=0;i<row;i++)
   {
      for(u64 j=0;j<column;j++)
      {
         buff[k]=send_share[i][j];
         k=k+1;
      }

   }
   if(party.Party_Id==Alice)
    {
      ctx.serialize(buff.begin(), buff.end(), sendbuff.begin());
      cp::sync_wait(macoro::when_all_ready(party.chl.send(sendbuff)));
      cp::sync_wait(macoro::when_all_ready(party.chll.recv(recvbuff)));
      ctx.deserialize(recvbuff.begin(), recvbuff.end(), buff1.begin());
    }
    else if(party.Party_Id==Bob)
    {
      ctx.serialize(buff.begin(), buff.end(), sendbuff.begin());
      cp::sync_wait(macoro::when_all_ready(party.chl.recv(recvbuff)));
      cp::sync_wait(macoro::when_all_ready(party.chll.send(sendbuff)));
      ctx.deserialize(recvbuff.begin(), recvbuff.end(), buff1.begin());
    }
    k=0; 
    for(u64 i=0;i<row;i++)
    {
      for(u64 j=0;j<column;j++)
      {
        ctx.plus(receive_share[i][j], send_share[i][j], buff1[k]);
        k=k+1;
      }
    }
}


template<typename G, typename Ctx>
void open_mac_with_commit(class Party<G,Ctx> party, std::vector<G> &send_mac, std::vector<G> &receive_mac)
{
   Ctx ctx;
   u64 row=send_mac.size();
   std::vector<G> buff(row),buff1(row);
   std::vector<u8> sendbuff,recvbuff;
   sendbuff.resize(buff.size() * ctx.template byteSize<G>());
   recvbuff.resize(buff1.size() * ctx.template byteSize<G>());
   Commit mycomm,theircomm,recvcomm;
   for(u64 k=0;k<row;k++)
   {
         buff[k]=send_mac[k];
   }
   ctx.serialize(buff.begin(), buff.end(), sendbuff.begin());
   mycomm=Commit(sendbuff);
   if(party.Party_Id==Alice)
    {
      cp::sync_wait(macoro::when_all_ready(party.chl.send(std::move(mycomm)), party.chll.recv(recvcomm)));
      cp::sync_wait(macoro::when_all_ready(party.chl.send(sendbuff), party.chll.recv(recvbuff)));
    }
    else if(party.Party_Id==Bob)
    {
      cp::sync_wait(macoro::when_all_ready(party.chll.send(std::move(mycomm)), party.chl.recv(recvcomm)));
      cp::sync_wait(macoro::when_all_ready(party.chl.recv(recvbuff), party.chll.send(sendbuff)));
    }
      ctx.deserialize(recvbuff.begin(), recvbuff.end(), buff1.begin());
    for(u64 k=0;k<row;k++)
    {
        receive_mac[k]=buff1[k];
    }

    theircomm=Commit(recvbuff);
     MACORO_TRY{
			if (theircomm!= recvcomm)
      {

				throw std::runtime_error("Shares are corrupted. Abort the protocol" LOCATION);
      }

			} MACORO_CATCH(eptr) {
				std::rethrow_exception(eptr);
			}
}


template<typename G, typename Ctx>
void check_mac(class Party<G,Ctx> party, std::vector<G> &mac, std::vector<std::vector<G>> &opened_matrix)
{
   Ctx ctx;
   std::vector<G> recv_mac(mac.size());
   std::vector<G> matrixKey(opened_matrix.size());
   ctx.copy(party.MatrixKey.begin(),party.MatrixKey.begin()+opened_matrix.size(), matrixKey.begin());
   ctx.vector_matrix_product(recv_mac, matrixKey, opened_matrix);
   ctx.vector_minus(mac,mac, recv_mac);
   open_mac_with_commit(party,mac,recv_mac);
    MACORO_TRY{
          for(u64 i=0;i<mac.size();i++)
          {
            G temp;
            ctx.plus(temp,recv_mac[i],mac[i]);
          
            //std::cout<<temp<<std::endl;
            if (temp!= 0)
            {
            //std::cout<<"temp"<<temp<<std::endl;

              throw std::runtime_error("Shares are corrupted. Abort the protocol" LOCATION);
            }      
          }
        } MACORO_CATCH(eptr) {
                        std::rethrow_exception(eptr);
                      }
}