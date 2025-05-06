#include <libOTe_Tests/Common.h>
#include <coproto/Socket/BufferingSocket.h>
#include "libOTe/Tools/Coproto.h"
#include <coproto/Socket/Socket.h>
#include <libOTe/config.h>
//#include "Tools/Coeff128.h"
#include <libOTe/Tools/CoeffCtx.h>
#include<libOTe/Vole/Silent/SilentVoleSender.h>
#include<libOTe/Vole/Silent/SilentVoleReceiver.h>
#include "SPDZshare.h"
#include "SPDZ_LocalComputation.h"
#include "cryptoTools/Crypto/Commit.h"
#include "libOTe/TwoChooseOne/SoftSpokenOT/SoftSpokenShOtExt.h"
#define Alice 1
#define Bob 2
using namespace osuCrypto;




template<typename G, typename Ctx>
void generate_spdz_share_with_mac(class Party<G,Ctx> party, SpdzShare<G, Ctx> &share, u64 batchSize, MultType type)
{
  
  using VecG = typename Ctx::template Vec<G>;
  VecG tempShare(batchSize), tempMac1(batchSize),tempMac2(batchSize);
  Ctx ctx;
  SilentVoleReceiver<G, G, Ctx> recv;
  SilentVoleSender<G, G, Ctx> send;
  share.configure(batchSize);
  recv.mMultType = type;
  send.mMultType = type;
  PRNG prng(sysRandomSeed());
  share.shareSeed=prng.get();
  PRNG prng1(share.shareSeed);
  //std::cout<<"generate spdz with mac set up clear"<<std::endl;
  if(party.Party_Id==Alice)
  {
   cp::sync_wait(
        macoro::when_all_ready(send.silentSend(party.SpdzKey, tempMac1, prng, party.chl), 
        recv.silentReceive(tempShare, tempMac2, prng1, party.chll)));
  }
  else if(party.Party_Id==Bob)
  {
    cp::sync_wait(
        macoro::when_all_ready(send.silentSend(party.SpdzKey, tempMac1, prng, party.chll), 
        recv.silentReceive(tempShare, tempMac2, prng1, party.chl)));
  }
 // std::cout<<"communication clear"<<std::endl;
  for(u64 i=0;i<batchSize;i++)
  {
   ctx.mul(share.shareMac[i], party.SpdzKey, tempShare[i]);
   ctx.minus(share.shareMac[i], share.shareMac[i], tempMac1[i]);
   ctx.plus(share.shareMac[i], share.shareMac[i], tempMac2[i]);
   share.share[i]=tempShare[i];
  }
}


template<typename G, typename Ctx>
void authentication_spdz(class Party<G,Ctx> party, SpdzShare<G, Ctx> &share)
{
   u64 batchSize=share.share.size();
   G publicCoin, temp;
   Ctx ctx;
   for(u64 i=0;i<batchSize-1;i++)
   {
     ctx.fromBlock(publicCoin, party.random_coin());
     ctx.mul(temp, publicCoin, share.share[i]);
     ctx.plus(share.share[batchSize-1], temp, share.share[batchSize-1]);
     ctx.mul(temp, publicCoin, share.shareMac[i]);
     ctx.plus(share.shareMac[batchSize-1], temp, share.shareMac[batchSize-1]);
   }
   open_share(party, share.share[batchSize-1], temp);
   ctx.mul(temp,  temp,  party.SpdzKey);
   ctx.minus(temp, share.shareMac[batchSize-1],temp); 
   open_spdz_mac(party, temp);
   share.share.resize(batchSize-1);
   share.shareMac.resize(batchSize-1);
     
}

template<typename G, typename Ctx>
void open_share(class Party<G,Ctx> party, G &sendshare, G &recvshare)
{
   Ctx ctx;
   std::vector<u8> sendbuff,recvbuff;
   std::vector<G> buff(1);
   buff[0]=sendshare;
   sendbuff.resize(ctx.template byteSize<G>());
   recvbuff.resize(ctx.template byteSize<G>());
   ctx.serialize(buff.begin(), buff.end(), sendbuff.begin());
   //std::cout<<"open share set up clear"<<std::endl;
    if(party.Party_Id==Alice)
    {
      
      cp::sync_wait(macoro::when_all_ready(party.chl.send(sendbuff), party.chll.recv(recvbuff)));
    }
    else if(party.Party_Id==Bob)
    {
     cp::sync_wait(macoro::when_all_ready(party.chll.send(sendbuff), party.chl.recv(recvbuff)));
    }
    //std::cout<<"communication clear"<<std::endl;
    ctx.deserialize(recvbuff.begin(), recvbuff.end(), &recvshare);
    ctx.plus(recvshare,recvshare,sendshare);
}

template<typename G, typename Ctx>
void open_batch_share(class Party<G,Ctx> party, std::vector<G> &sendshare, std::vector<G> &recvshare)
{
   Ctx ctx;
   std::vector<u8> sendbuff,recvbuff;
   u64 batchSize=sendshare.size();
   //std::vector<G> buff(1);
   //buff[0]=sendshare;
   sendbuff.resize(batchSize*ctx.template byteSize<G>());
   recvbuff.resize(batchSize*ctx.template byteSize<G>());
   ctx.serialize(sendshare.begin(), sendshare.end(), sendbuff.begin());
   //std::cout<<"open share set up clear"<<std::endl;
    if(party.Party_Id==Alice)
    {
      
      cp::sync_wait(macoro::when_all_ready(party.chl.send(sendbuff), party.chll.recv(recvbuff)));
    }
    else if(party.Party_Id==Bob)
    {
     cp::sync_wait(macoro::when_all_ready(party.chll.send(sendbuff), party.chl.recv(recvbuff)));
    }
    //std::cout<<"communication clear"<<std::endl;
    ctx.deserialize(recvbuff.begin(), recvbuff.end(), recvshare.begin());
    ctx.vector_plus(recvshare,recvshare,sendshare);
}

template<typename G, typename Ctx>
void open_spdz_mac(class Party<G,Ctx> party, G &sendshare)
{
   std::vector<u8> sendbuff,recvbuff;
   std::vector<G> buff(1);
   G recvshare;
   Ctx ctx;
   sendbuff.resize(ctx.template byteSize<G>());
   recvbuff.resize(ctx.template byteSize<G>());
   buff[0]=sendshare;
   ctx.serialize(buff.begin(), buff.end(), sendbuff.begin());
   Commit mycomm,recvcomm,theircomm;
   mycomm=Commit(sendbuff);
   if(party.Party_Id==Alice)
    {
      cp::sync_wait(macoro::when_all_ready(party.chl.send(std::move(mycomm)), party.chll.recv(recvcomm)));
      cp::sync_wait(macoro::when_all_ready(party.chl.send(sendbuff), party.chll.recv(recvbuff)));
      
    }
    else if(party.Party_Id==Bob)
    {
       cp::sync_wait(macoro::when_all_ready(party.chll.send(std::move(mycomm)), party.chl.recv(recvcomm)));
       cp::sync_wait(macoro::when_all_ready(party.chll.send(sendbuff), party.chl.recv(recvbuff))); 
    }
    ctx.deserialize(recvbuff.begin(), recvbuff.end(), &recvshare);
    
    G temp;
    ctx.plus(temp,sendshare, recvshare);

    MACORO_TRY{
       theircomm=Commit(recvbuff);
			if (temp==0 && theircomm==recvcomm)
      {

      }
      else{
				throw std::runtime_error("Shares are corrupted. Abort the protocol" LOCATION);
      }
			} MACORO_CATCH(eptr) {
				std::rethrow_exception(eptr);
			}
    }

template<typename G, typename Ctx>
void generate_Beaver_tripple(class Party<G,Ctx> party, SpdzShare<G, Ctx> &shareA, SpdzShare<G, Ctx> &shareB, 
SpdzShare<G, Ctx> &shareC, u64 batchSize, MultType type)
{
  u64 tau=4;
  Ctx ctx;
  std::vector<SpdzShare<G,Ctx>> A(tau), C(tau);
  G publicCoin;
  SoftSpokenMalOtSender baseSend;
  SoftSpokenMalOtReceiver baseRecv;
  shareB.configure(batchSize+1);
  generate_spdz_share_with_mac(party,shareB,batchSize+1,type);
  for(u64 i=0;i<tau;i++)
  {
    generate_spdz_share_with_mac(party,A[i],batchSize+1,type);
    C[i].configure(batchSize+1);
  }
  if(party.Party_Id==Alice)
  {
    for(u64 i=0;i<tau;i++)
    {
      std::vector<G> tempShare(batchSize+1);
      Batch_OLE_Sender<G,Ctx>(A[i].share, C[i].share, baseSend, party.chl);
      Batch_OLE_Receiver<G,Ctx>(shareB.share, tempShare, baseRecv, party.chll);
      for(u64 j=0;j<batchSize+1;j++)
      {
        G temp;
        ctx.mul(temp, A[i].share[j],shareB.share[j]);
        ctx.minus(C[i].share[j], C[i].share[j],tempShare[j]);
        ctx.plus(C[i].share[j], C[i].share[j], temp);
      }
      Batch_OLE_Sender<G,Ctx>(A[i].shareMac, C[i].shareMac, baseSend, party.chl);
      Batch_OLE_Receiver<G,Ctx>(shareB.share, tempShare, baseRecv, party.chll);
      for(u64 j=0;j<batchSize+1;j++)
      {
        G temp;
        ctx.mul(temp, A[i].shareMac[j],shareB.share[j]);
        ctx.minus(C[i].shareMac[j], C[i].shareMac[j],tempShare[j]);
        ctx.plus(C[i].shareMac[j], C[i].shareMac[j], temp);
      }

    }
  }
  else if(party.Party_Id==Bob)
  {
    for(u64 i=0;i<tau;i++)
    {
      std::vector<G> tempShare(batchSize+1);
      Batch_OLE_Receiver<G,Ctx>(shareB.share, tempShare, baseRecv, party.chl);
      Batch_OLE_Sender<G,Ctx>(A[i].share, C[i].share, baseSend, party.chll);
      
      for(u64 j=0;j<batchSize+1;j++)
      {
        G temp;
        ctx.mul(temp, A[i].share[j],shareB.share[j]);
        ctx.minus(C[i].share[j], C[i].share[j],tempShare[j]);
        ctx.plus(C[i].share[j], C[i].share[j], temp);
      }
      Batch_OLE_Receiver<G,Ctx>(shareB.share, tempShare, baseRecv, party.chl);
      Batch_OLE_Sender<G,Ctx>(A[i].shareMac, C[i].shareMac, baseSend, party.chll);
      for(u64 j=0;j<batchSize+1;j++)
      {
        G temp;
        ctx.mul(temp, A[i].shareMac[j],shareB.share[j]);
        ctx.minus(C[i].shareMac[j], C[i].shareMac[j],tempShare[j]);
        ctx.plus(C[i].shareMac[j], C[i].shareMac[j], temp);
      }

    }
  }
      SpdzShare<G, Ctx> TempShareA, shareAA, shareCC;
      TempShareA.configure(batchSize+1);
      shareA.configure(batchSize+1);
      shareC.configure(batchSize+1);
      shareAA.configure(batchSize+1);
      shareCC.configure(batchSize+1);
      for(u64 i=0;i<tau;i++)
      {
        ctx.fromBlock(publicCoin, party.random_coin());
        SPDZ_local_scaler(TempShareA, publicCoin, A[i],ctx);
        SPDZ_local_plus(shareA, shareA, TempShareA,ctx);
        SPDZ_local_scaler(TempShareA, publicCoin, C[i],ctx);
        SPDZ_local_plus(shareC, shareC, TempShareA,ctx);
      }
      for(u64 i=0;i<tau;i++)
      {
        ctx.fromBlock(publicCoin, party.random_coin());
        SPDZ_local_scaler(TempShareA, publicCoin, A[i],ctx);
        SPDZ_local_plus(shareAA, shareAA, TempShareA,ctx);
        SPDZ_local_scaler(TempShareA, publicCoin, C[i],ctx);
        SPDZ_local_plus(shareCC, shareCC, TempShareA,ctx);
      }
     // std::cout<<"authentication start"<<std::endl;
      authentication_spdz(party, shareA);
      authentication_spdz(party, shareAA);
      authentication_spdz(party, shareB);
      authentication_spdz(party, shareC);
      authentication_spdz(party, shareCC);
      ctx.fromBlock(publicCoin, party.random_coin());
      TempShareA.configure(batchSize);
      SPDZ_local_scaler(TempShareA, publicCoin, shareA,ctx);
      SPDZ_local_minus(shareAA, shareAA, TempShareA,ctx);
      SPDZ_local_scaler(TempShareA, publicCoin, shareC,ctx);
      SPDZ_local_minus(shareCC, shareCC, TempShareA,ctx);
      open_batch_share(party, shareAA.share,TempShareA.share);
      G tempMac=0;
      //std::cout<<"authentication is over"<<std::endl;
      for(u64 i=0;i<batchSize;i++)
      {
        G temp;
        ctx.mul(temp,party.SpdzKey,TempShareA.share[i]);
        ctx.minus(temp, shareAA.shareMac[i],temp);
        ctx.fromBlock(publicCoin, party.random_coin());
        ctx.mul(temp,publicCoin,temp);
        ctx.plus(tempMac,tempMac,temp);
        ctx.mul(TempShareA.shareMac[i], TempShareA.share[i], shareB.shareMac[i]);
        ctx.mul(TempShareA.share[i], TempShareA.share[i], shareB.share[i]);
      }
        //std::cout<<"open spdz mac starts"<<std::endl;
        //std::cout<<"tempMac="<<tempMac<<std::endl;
        //std::cout<<"party_ID"<<party.Party_Id<<std::endl;
        open_spdz_mac(party, tempMac);
        //std::cout<<"open spdz mac over"<<std::endl;
        SPDZ_local_minus(TempShareA, shareCC, TempShareA,ctx);
        std::vector<G> recvShareA(batchSize);
        open_batch_share(party, TempShareA.share, recvShareA);
        tempMac=0;
        //std::cout<<"open batch share is over"<<std::endl;
     for(u64 i=0;i<batchSize;i++)
      {
        G temp;
        ctx.fromBlock(publicCoin, party.random_coin());
        ctx.mul(temp,publicCoin,TempShareA.shareMac[i]);
        ctx.plus(tempMac,tempMac,temp);

        MACORO_TRY{
        //std::cout<<temp<<std::endl;
        if (recvShareA[i]==0)
        {

        }
        else{
        //std::cout<<"temp="<<temp<<"  party_id="<<party.Party_Id<<std::endl;
          throw std::runtime_error("Shares are corrupted. Abort the protocol" LOCATION);
        }
        } MACORO_CATCH(eptr) {
          std::rethrow_exception(eptr);
        }
        
      } 
      open_spdz_mac(party, tempMac);
}