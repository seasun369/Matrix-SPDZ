#include <libOTe_Tests/Common.h>
#include "libOTe/TwoChooseOne/SoftSpokenOT/SoftSpokenShOtExt.h"
#include <coproto/Socket/BufferingSocket.h>
#include <coproto/Socket/Socket.h>
#include "cryptoTools/Common/BitVector.h"
#include <libOTe/config.h>
//#include "Tools/Coeff128.h"


using namespace osuCrypto;


/* 

This protocol realizes the batch implementation of oblivious linear evaluations. 
In particular, the input vector of sender is b and the input vector of receiver is d. 
The output a and c satisfies that a=b*d-c. The idea is to compute the OLE with correlated OT. 
For example, write delta as binary representation (delta[i]) and to obliviously compute delta*x, one 
compute the delta[i]*x with the help of OT. 

*/


template<typename G, typename Ctx>
void Batch_OLE_Receiver(std::vector<G> b, std::vector<G> & a, OtReceiver& ot, Socket & chl)
{
      Ctx ctx;
      BitVector bv;
      PRNG prng(sysRandomSeed());
      std::vector<G> msg;
      auto buffer = std::vector<u8>{};
      u64 fieldSize=ctx.template bitSize<G>();
      // each bit of bv serves as the choice bits.
      for(u64 i=0;i<b.size();i++)
      {
        bv.append(ctx.binaryDecomposition(b[i]));
      }
      auto otMsg = AlignedVector<block>{ };
		  otMsg.resize(bv.size());
      cp::sync_wait(ot.receive(bv, otMsg, prng, chl));
      ctx.zero(a.begin(), a.end());
      u64 bufferSize=a.size() * ctx.template byteSize<G>()*fieldSize;
      buffer.resize(bufferSize);
			cp::sync_wait(chl.recv(buffer));
			ctx.resize(msg, a.size()*fieldSize);
			ctx.deserialize(buffer.begin(), buffer.end(), msg.begin());
      u64 k=0;
      for(u64 i=0;i<a.size();i++)
      {
           G temp;
           for(u64 j=0;j<fieldSize;j++)
           {
              prng.SetSeed(otMsg[k]);
              ctx.fromBlock(temp, prng.get<block>());
              if(bv[k])
              {
                ctx.plus(temp, msg[k], temp);
              }
                ctx.plus(a[i], a[i], temp);
                k=k+1;
           }
      }

}

template<typename G, typename Ctx>
void Batch_OLE_Sender(std::vector<G> d, std::vector<G> & c, OtSender& ot, Socket & chl)
{
      PRNG prng(sysRandomSeed());
      Ctx ctx; 
      auto otMsg = AlignedVector<std::array<block, 2>>{};
      auto buffer = std::vector<u8>{};
      std::vector<G> temp,msg;
      u64 fieldSize=ctx.template bitSize<G>();
	  otMsg.resize(c.size()*fieldSize);
	  cp::sync_wait(ot.send(otMsg, prng, chl));  
    ctx.zero(c.begin(), c.end());
	  ctx.resize(msg, otMsg.size());
	  ctx.resize(temp, 2);
      u64 k=0;
      for (u64 i = 0; i < c.size(); ++i)
	{
        for(u64 j=0;j<fieldSize;j++)
        {
          prng.SetSeed(otMsg[k][0]);
          ctx.powerOfTwo(temp[1], j);
          ctx.fromBlock(msg[k], prng.get<block>());
          ctx.plus(c[i], c[i], msg[k]);
          ctx.mul(temp[0], temp[1], d[i]);
          ctx.minus(msg[k], msg[k], temp[0]);
          prng.SetSeed(otMsg[k][1]);
          ctx.fromBlock(temp[0], prng.get<block>());
          ctx.minus(msg[k], msg[k], temp[0]);
          k=k+1;
        }
      }
      
      buffer.resize(msg.size() * ctx.template byteSize<G>());
			ctx.serialize(msg.begin(), msg.end(), buffer.begin());
			cp::sync_wait(chl.send(std::move(buffer)));

}

template<typename G, typename Ctx>
void Inner_Product_Receiver(std::vector<G> b, Socket & chl, u64 &sum)
{
  SoftSpokenMalOtReceiver baseRecv;
  std::vector<G> a(b.size());
  Ctx ctx;
  Batch_OLE_Receiver<G,Ctx>(b, a, baseRecv, chl);
  u64 temp=0;
  for(u64 i=0;i<a.size();i++)
  {
    ctx.plus(temp,temp,a[i]);
  }
    cp::sync_wait(macoro::when_all_ready(chl.send(std::move(temp))));
    cp::sync_wait(macoro::when_all_ready(chl.recv(sum)));
    ctx.minus(sum,sum,temp);
}

template<typename G, typename Ctx>
void Inner_Product_Sender(std::vector<G> d, Socket & chl, u64 &sum)
{
  SoftSpokenMalOtSender baseSend;
  std::vector<G> c(d.size());
  Ctx ctx;
  Batch_OLE_Sender<G,Ctx>(d, c, baseSend, chl);
  u64 temp=0;
  for(u64 i=0;i<c.size();i++)
  {
    ctx.plus(temp,temp,c[i]);
  }
  cp::sync_wait(macoro::when_all_ready(chl.recv(sum)));
  cp::sync_wait(macoro::when_all_ready(chl.send(std::move(temp))));
  ctx.minus(sum,temp,sum);
}

template<typename G, typename Ctx>
void convert_to_field(std::vector<double> & a, std::vector<G> &output, u64 precision)
{
   Ctx ctx;
   output.resize(a.size());
   u64 pow=1<<precision;
   for(u64 i=0;i<a.size();i++)
   {
     output[i]=std::floor(a[i]*pow);
     ctx.plus(output[i], output[i],0);
   }
}

template<typename G, typename Ctx>
void convert_to_float_vector(std::vector<G> & a, std::vector<double> &output, u64 precision)
{
   Ctx ctx;
   output.resize(a.size());
   std::vector<__int128_t> temp(a.size());
   __int128_t prime=ctx.Prime(); 
   u64 pow=1<<precision;
   for(u64 i=0;i<a.size();i++)
   {
     if(a[i]>prime/2)
     {
       temp[i]=a[i]-prime;
     }
     else
     {
       temp[i]=a[i];
     }
     output[i]=double(temp[i])/pow;
   }
}

template<typename G, typename Ctx>
void convert_to_float_element(G & a, double &output, u64 precision)
{
   Ctx ctx;
   __int128_t temp;
   __int128_t prime=ctx.Prime(); 
   u64 pow=1<<precision;
     if(a>prime/2)
     {
       temp=a-prime;
     }
     else
     {
       temp=a;
     }
   output=double(temp)/pow;
}
