#define Alice 1
#define Bob 2
using namespace osuCrypto;


template<typename G, typename Ctx>
void generate_matrix_transpose(class Party<G,Ctx> party, std::vector<Share<G,Ctx>> &A, std::vector<Share<G,Ctx>> &A_t, 
std::vector<Share<G,Ctx>> &R, std::vector<Share<G,Ctx>> &R_t,  u64 leftMatrixRow, u64 leftMatrixColumn, u64 rightMatrixColumn, 
MultType type)
{
    Ctx ctx;
    u64 batch_size=A.size()+1;
    A.resize(batch_size);
    std::vector<Share<G, Ctx>> authenticated_sharing1(batch_size),  authenticated_sharing2(batch_size);
    Share<G, Ctx> sacrifice1, sacrifice2, sacrifice3, sacrifice4;
    A_t.resize(batch_size);
    R.resize(batch_size);
    R_t.resize(batch_size);
    generate_Matrix_and_MAC_with_seed(party, A[batch_size-1], leftMatrixRow, leftMatrixColumn, type);
    //generate <A^T>, <R>, <R^T>
    for(u64 i=0; i<batch_size; i++)
    {
        A_t[i].configure(leftMatrixColumn, leftMatrixRow);
        R_t[i].configure(rightMatrixColumn,leftMatrixRow);
        generate_Matrix_and_MAC_with_seed(party, R[i], leftMatrixRow, rightMatrixColumn, type);
        generate_Matrix_and_MAC_with_seed(party, authenticated_sharing1[i], leftMatrixColumn,leftMatrixRow, type);
        generate_Matrix_and_MAC_with_seed(party, authenticated_sharing2[i], rightMatrixColumn, leftMatrixRow, type);
        ctx.transpose(R_t[i].MatrixShare, R[i].MatrixShare);
        ctx.transpose(A_t[i].MatrixShare, A[i].MatrixShare);
    }
        generate_Matrix_and_MAC_with_seed(party, sacrifice1, leftMatrixRow, rightMatrixColumn, type);
        generate_Matrix_and_MAC_with_seed(party, sacrifice2, leftMatrixColumn,leftMatrixRow, type);
        generate_Matrix_and_MAC_with_seed(party, sacrifice3, rightMatrixColumn, leftMatrixRow, type);
        generate_Matrix_and_MAC_with_seed(party, sacrifice4, leftMatrixRow, leftMatrixColumn, type);
        authentication(party, sacrifice1, R);
        authentication(party, sacrifice2, authenticated_sharing1);
        authentication(party, sacrifice3, authenticated_sharing2);
        authentication(party, sacrifice4, A);
        for(u64 i=0; i<batch_size; i++)
        {
            generate_MAC_with_authenticated_sharing(party, A_t[i], authenticated_sharing1[i], type);
            generate_MAC_with_authenticated_sharing(party, R_t[i], authenticated_sharing2[i], type);     
        }   
        generate_Matrix_and_MAC_with_seed(party, sacrifice1, leftMatrixColumn,leftMatrixRow, type);
        generate_Matrix_and_MAC_with_seed(party, sacrifice2, rightMatrixColumn, leftMatrixRow, type);
        authentication(party, sacrifice1, A_t);
        authentication(party, sacrifice2, R_t);
        //check whether <A> and <A^T>, <R> and <R^T> are the transpose of each other
        for(u64 i=0; i<batch_size-1; i++)
        {
            G public_coin;
            ctx.fromBlock(public_coin, party.random_coin());
            sacrifice1.configure(leftMatrixRow, leftMatrixColumn);
            local_scaler(sacrifice1, public_coin, A[i], ctx);
            local_plus(A[batch_size-1],sacrifice1, A[batch_size-1],ctx);
            sacrifice2.configure(leftMatrixColumn,leftMatrixRow);
            local_scaler(sacrifice2, public_coin, A_t[i], ctx);
            local_plus(A_t[batch_size-1],sacrifice2, A_t[batch_size-1],ctx); 
            sacrifice3.configure(leftMatrixRow,rightMatrixColumn);
            local_scaler(sacrifice3, public_coin, R[i], ctx);
            local_plus(R[batch_size-1],sacrifice3, R[batch_size-1],ctx);  
            sacrifice4.configure(rightMatrixColumn,leftMatrixRow);
            local_scaler(sacrifice4, public_coin, R_t[i], ctx);
            local_plus(R_t[batch_size-1],sacrifice4, R_t[batch_size-1],ctx);  
        }   
           //take the linear combination of M1=\sum r_i<A[i]> and M_2=\sum r_i<A[i]^T> to check if 
           //  M1=M2^T. Do the same thing for <R[i]> and <R[i]^T>
            open_matrix(party, A[batch_size-1].MatrixShare, sacrifice1.MatrixShare);
            check_mac(party, A[batch_size-1].leftMatrixShareMac, sacrifice1.MatrixShare);
            open_matrix(party, A_t[batch_size-1].MatrixShare, sacrifice2.MatrixShare);
            check_mac(party, A_t[batch_size-1].leftMatrixShareMac, sacrifice2.MatrixShare);
            open_matrix(party, R[batch_size-1].MatrixShare, sacrifice3.MatrixShare);
            check_mac(party, R[batch_size-1].leftMatrixShareMac, sacrifice3.MatrixShare);
            open_matrix(party, R_t[batch_size-1].MatrixShare, sacrifice4.MatrixShare);
            check_mac(party, R_t[batch_size-1].leftMatrixShareMac, sacrifice4.MatrixShare);
            for(u64 i=0;i<sacrifice1.MatrixShare.size();i++)
            {
                for(u64 j=0;j<sacrifice1.MatrixShare[0].size();j++)
                {
                    MACORO_TRY{
        
                        if (sacrifice1.MatrixShare[i][j]!=sacrifice2.MatrixShare[j][i])
                          throw std::runtime_error("Shares are corrupted. Abort the protocol" LOCATION);
                
                        } MACORO_CATCH(eptr) {
                          std::rethrow_exception(eptr);
                        }
                }
            }

            for(u64 i=0;i<sacrifice3.MatrixShare.size();i++)
            {
                for(u64 j=0;j<sacrifice3.MatrixShare[0].size();j++)
                {
                    MACORO_TRY{
        
                        if (sacrifice3.MatrixShare[i][j]!=sacrifice4.MatrixShare[j][i])
                          throw std::runtime_error("Shares are corrupted. Abort the protocol" LOCATION);
                
                        } MACORO_CATCH(eptr) {
                          std::rethrow_exception(eptr);
                        }
                }
            }
            A.resize(batch_size-1);
            A_t.resize(batch_size-1);
            R.resize(batch_size-1);
            R_t.resize(batch_size-1);
}