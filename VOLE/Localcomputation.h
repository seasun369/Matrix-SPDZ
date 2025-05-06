

using namespace osuCrypto;

template<typename G, typename Ctx>
void local_plus(Share<G,Ctx> &ret, const Share<G,Ctx> &shareA, const Share<G,Ctx> &shareB, Ctx ctx)
{
    ctx.matrix_plus(ret.MatrixShare, shareA.MatrixShare, shareB.MatrixShare);
    ctx.vector_plus(ret.leftMatrixShareMac, shareA.leftMatrixShareMac, shareB.leftMatrixShareMac);
}

template<typename G, typename Ctx>
void local_minus(Share<G,Ctx> &ret, const Share<G,Ctx> &shareA, const Share<G,Ctx> &shareB, Ctx ctx)
{
   ctx.matrix_minus(ret.MatrixShare, shareA.MatrixShare, shareB.MatrixShare);
   ctx.vector_minus(ret.leftMatrixShareMac, shareA.leftMatrixShareMac, shareB.leftMatrixShareMac);
}

template<typename G, typename Ctx>
void local_scaler(Share<G,Ctx> &ret, const G x, const Share<G,Ctx> &shareA, Ctx ctx)
{
    ctx.matrix_scaler(ret.MatrixShare, x, shareA.MatrixShare);
    ctx.vector_scaler(ret.leftMatrixShareMac, x, shareA.leftMatrixShareMac);
}

template<typename G, typename Ctx>
void local_matrix_product(Share<G,Ctx> &ret,  const Share<G,Ctx> &shareA, const std::vector<std::vector<G>> publicMatrix, Ctx ctx)
{
    ctx.matrix_product(ret.MatrixShare, shareA.MatrixShare, publicMatrix);
    ctx.vector_matrix_product(ret.leftMatrixShareMac, shareA.leftMatrixShareMac, publicMatrix);
}