using namespace osuCrypto;


/*


*/


template<typename G, typename Ctx>
void SPDZ_local_plus(SpdzShare<G,Ctx> &ret, const SpdzShare<G,Ctx> &shareA, const SpdzShare<G,Ctx> &shareB, Ctx ctx)
{
    ctx.vector_plus(ret.share, shareA.share, shareB.share);
    ctx.vector_plus(ret.shareMac, shareA.shareMac, shareB.shareMac);
}


template<typename G, typename Ctx>
void SPDZ_local_minus(SpdzShare<G,Ctx> &ret, const SpdzShare<G,Ctx> &shareA, const SpdzShare<G,Ctx> &shareB, Ctx ctx)
{
    ctx.vector_minus(ret.share, shareA.share, shareB.share);
    ctx.vector_minus(ret.shareMac, shareA.shareMac, shareB.shareMac);
}

template<typename G, typename Ctx>
void SPDZ_local_scaler(SpdzShare<G,Ctx> &ret, const G x, const SpdzShare<G,Ctx> &shareA, Ctx ctx)
{
    ctx.vector_scaler(ret.share, x, shareA.share);
    ctx.vector_scaler(ret.shareMac, x, shareA.shareMac);
}

