// © 2022 Visa.
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <libOTe/config.h>


#include "cryptoTools/Common/Defines.h"
#include "cryptoTools/Common/Timer.h"
#include "cryptoTools/Crypto/PRNG.h"
#include "libOTe/Tools/Coproto.h"
#include "libOTe/TwoChooseOne/OTExtInterface.h"
#include "Coeff128.h"

using namespace osuCrypto;

	template <
		typename F,
		typename G,
		typename CoeffCtx
	>
	class NoisyVoleReceiver_TensorVector : public TimerAdapter
	{
	public:
		using VecF = typename CoeffCtx::template Vec<F>;
		using VecG = typename CoeffCtx::template Vec<G>;


		// for chosen c, compute a such htat
		//
		//  a = b + c * delta
		//
		template<typename VecG, typename VecF>
		task<> receive(VecG& c, VecF& a, PRNG& prng,
			OtSender& ot, Socket& chl, CoeffCtx ctx)
		{
			MACORO_TRY{

			auto otMsg = AlignedUnVector<std::array<block, 2>>{};

			setTimePoint("NoisyVoleReceiver.ot.begin");
			otMsg.resize(ctx.template bitSize<F>());
			co_await(ot.send(otMsg, prng, chl));

			setTimePoint("NoisyVoleReceiver.ot.end");

			co_await(receive(c, a, prng, otMsg, chl, ctx));

			} MACORO_CATCH(eptr) {
				if (!chl.closed()) co_await chl.close();
				std::rethrow_exception(eptr);
			}
		}

		// for chosen c, compute a such htat
		//
		//  a = b + c * delta
		//
		template<typename VecG, typename VecF>
		task<> receive(VecG& c, VecF& a, PRNG& _,
			span<std::array<block, 2>> otMsg,
			Socket& chl, CoeffCtx ctx)
		{
			MACORO_TRY{
			auto buff = std::vector<u8>{};
			auto msg = VecG{};
			auto temp = VecG{};
			auto prng = PRNG{};

			if (c.size() != a.size())
				throw RTE_LOC;
			if (a.size() == 0)
				throw RTE_LOC;

			setTimePoint("NoisyVoleReceiver.begin");

			ctx.zero(a.begin(), a.end());
			ctx.resize(msg, otMsg.size() * a.size());
			ctx.resize(temp, 2);
			u64 batchsize=sizeof(a[0][0])*8*a.size();
            //std::cout<< "otMSg.size="<<otMsg.size()<<"  c.size="<<c.size()<<std::endl;
            //std::cout<< "msg.size="<<msg.size()<<"  a.size="<<a.size()<<std::endl;
			//std::cout<<"  ctx.size="<<sizeof(a[0][0])*8<<std::endl;
			//u64 batchsize=8*sizeof(a[0][0]);
			u64 fieldSize=8*sizeof(a[0][0]);
			for (size_t i = 0, k = 0; i < otMsg.size(); ++i)
			{
				prng.SetSeed(otMsg[i][0], a.size());
                
				// t1 = 2^i
				ctx.powerOfTwo(temp[1], (i % fieldSize));
                u64 current_index=k/batchsize;
				//std::cout << temp[1] << std::endl;
				for (size_t j = 0; j < c.size(); ++j, ++k)
				{
					// msg[i,j] = otMsg[i,j,0]
					ctx.fromBlock(msg[k], prng.get<block>());
					
				
					//ctx.zero(msg.begin() + k, msg.begin() + k + 1);
					//std::cout << "m" << i << ",0 = " << ctx.str(msg[k]) << std::endl;

					// a[j] += otMsg[i,j,0]
					ctx.plus(a[j][current_index], a[j][current_index], msg[k]);

				
					//std::cout << "z = " << ctx.str(a[j]) << std::endl;

					// temp = 2^i * c[j]
					ctx.mul(temp[0], temp[1], c[j]);

					
					//std::cout << "2^i y = " << ctx.str(temp[0]) << std::endl;
                     
					// msg[i,j] = otMsg[i,j,0] + 2^i * c[j]
					ctx.minus(msg[k], msg[k], temp[0]);

					
					//std::cout << "m" << i << ",0 + 2^i y = " << ctx.str(msg[k]) << std::endl;
				}

				k -= c.size();
				prng.SetSeed(otMsg[i][1], a.size());

				for (size_t j = 0; j < c.size(); ++j, ++k)
				{
					// temp = otMsg[i,j,1]
					ctx.fromBlock(temp[0], prng.get<block>());
					//ctx.zero(temp.begin(), temp.begin() + 1);
					//std::cout << "m" << i << ",1 = " << ctx.str(temp[0]) << std::endl;

					// enc one message under the OT msg.
					// msg[i,j] = (otMsg[i,j,0] + 2^i * c[j]) - otMsg[i,j,1]
					ctx.minus(msg[k], msg[k], temp[0]);
					//std::cout << "m" << i << ",0 + 2^i y - m" << i << ",1 = " << ctx.str(msg[k]) << std::endl << std::endl;
				}
			}

          // std::cout<<"computation is over"<<std::endl;
           //std::cout<<"buffersize="<<msg.size() * ctx.template byteSize<G>()<<std::endl;
			buff.resize(msg.size() * ctx.template byteSize<G>());
			ctx.serialize(msg.begin(), msg.end(), buff.begin());

           //std::cout<<"buffer resize recv"<<std::endl;
			co_await(chl.send(std::move(buff)));
			setTimePoint("NoisyVoleReceiver.done");

			} MACORO_CATCH(eptr) {
				if (!chl.closed()) co_await chl.close();
				std::rethrow_exception(eptr);
			}
		}
       
	};
	  

  // namespace osuCrypto

