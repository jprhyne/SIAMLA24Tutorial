/**
 * This file will create data files to be used for figures in the presentation.
 * We will be taking advantage of the templating features by running it for single and double precision
 */
#include <cstring>
#include <cstdint>
#include "testutils.hpp"
#include <tlapack/lapack/lacpy.hpp>
#include <tlapack/lapack/lange.hpp>
#include <tlapack/lapack/gehrd.hpp>
#include <tlapack/lapack/unghr.hpp>
#include <tlapack/lapack/getri.hpp>
#include <tlapack/lapack/getrf.hpp>
#include <tlapack/blas/gemm.hpp>

#include <tlapack/lapack/eispack_hqr.hpp>
#include <tlapack/lapack/eispack_schurToEigen.hpp>
#include <tlapack/lapack/multishift_qr.hpp>
#include <sys/time.h>

// Helper function to more easily grab the current time
// used to facilitate timing executions
// Stolen from: https://stackoverflow.com/a/17440673 as well as previous C projects.
double getTime()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
}
using namespace tlapack;
int main( int argc, char **argv)
{
    /**********************************************************************************************
     * Single precision section                                                                   *
     **********************************************************************************************/
    using T1 = float;
    using idx_t = std::size_t;
    using complex_t1 = std::complex<T1>;
    Create<legacyMatrix<T1,idx_t,Layout::ColMajor>> new_matrix1;
    idx_t i,n;
    n = -1;
    // Do some input parsing to determine if we give the size of the matrix we want to work with.
    for(i = 1; i < argc; ++i){
        if( strcmp( *(argv + i), "-n") == 0) {
            n  = atoi( *(argv + i + 1) );
            i++;
        }
    }

    if (n == -1)
        n = 30;


    // Create matrices
    std::vector<T1> A1_; auto A1 = new_matrix1( A1_, n, n);
    std::vector<T1> H1_; auto H1 = new_matrix1( H1_, n, n);
    std::vector<T1> Hs1_; auto Hs1 = new_matrix1( Hs1_, n, n);
    std::vector<T1> Q1_; auto Q1 = new_matrix1( Q1_, n, n);
    std::vector<T1> Qs1_; auto Qs1 = new_matrix1( Qs1_, n, n);
    std::vector<complex_t1> s1(n);
    //Populate A and U with random numbers
    rand_generator gen;
    idx_t seed = 17519;
    gen.seed(seed);
    for(idx_t i = 0; i < n; i++) {
        for(idx_t j = 0; j < n; j++) {
            T1 val = rand_helper<T1>(gen);
            A1(i,j) = val; 
            H1(i,j) = val; 
        }
    }
    std::vector<T1> tau1(n);
    T1 zero1 = T1(0);
    T1 one1 = T1(1);
    gehrd(0, n, H1, tau1);
    lacpy(Uplo::General, H1, Q1);
    lacpy(Uplo::General, H1, Hs1);
    unghr(0, n, Q1, tau1);
    lacpy(Uplo::General, Q1, Qs1);
    idx_t ilo = 0;
    idx_t igh = n;
    for (idx_t i = 1; i < n; i++)
        for (idx_t j = 0; j < i - 1; j++)
            H1(i,j) = zero1;
    double timeHqr = -getTime();
    //auto normA = lange(tlapack::frob_norm, A);
    T1 norm1 = zero1;
    auto retCode = tlapack::eispack_hqr(true, true, ilo, igh, H1, s1, Q1, norm1);
    timeHqr += getTime();
    std::cout << "Single precision\n";
    // Print out the time for computing both the eigenvalues, T, and the schur vectors
    std::cout << "TQE:" << timeHqr << "\n";
    // copy over everything back to the original matrices
    lacpy(Uplo::General, Hs1, H1);
    lacpy(Uplo::General, Qs1, Q1);
    timeHqr = -getTime();
    retCode = tlapack::eispack_hqr(true, false, ilo, igh, H1, s1, Q1, norm1);
    timeHqr += getTime();
    std::cout << "TE:" << timeHqr << "\n";
    // copy over everything back to the original matrices
    lacpy(Uplo::General, Hs1, H1);
    lacpy(Uplo::General, Qs1, Q1);
    timeHqr = -getTime();
    retCode = tlapack::eispack_hqr(false, false, ilo, igh, H1, s1, Q1, norm1);
    timeHqr += getTime();
    std::cout << "E:" << timeHqr << "\n";
    /**********************************************************************************************
     * Double precision section                                                                   *
     **********************************************************************************************/
    using T2 = double;
    using idx_t = std::size_t;
    using complex_t2 = std::complex<T2>;
    Create<legacyMatrix<T2,idx_t,Layout::ColMajor>> new_matrix2;
    // Create matrices
    std::vector<T2> A2_; auto A2 = new_matrix2( A2_, n, n);
    std::vector<T2> H2_; auto H2 = new_matrix2( H2_, n, n);
    std::vector<T2> Hs2_; auto Hs2 = new_matrix2( Hs2_, n, n);
    std::vector<T2> Q2_; auto Q2 = new_matrix2( Q2_, n, n);
    std::vector<T2> Qs2_; auto Qs2 = new_matrix2( Qs2_, n, n);
    std::vector<complex_t2> s2(n);
    //Populate A and U with random numbers
    gen.seed(seed);
    for(idx_t i = 0; i < n; i++) {
        for(idx_t j = 0; j < n; j++) {
            T2 val = rand_helper<T2>(gen);
            A2(i,j) = val; 
            H2(i,j) = val; 
        }
    }
    std::vector<T2> tau2(n);
    T2 zero2 = T2(0);
    T2 one2 = T2(1);
    gehrd(0, n, H2, tau2);
    lacpy(Uplo::General, H2, Q2);
    lacpy(Uplo::General, H2, Hs2);
    unghr(0, n, Q2, tau2);
    lacpy(Uplo::General, Q2, Qs2);
    for (idx_t i = 1; i < n; i++)
        for (idx_t j = 0; j < i - 1; j++)
            H2(i,j) = zero2;
    timeHqr = -getTime();
    //auto normA = lange(tlapack::frob_norm, A);
    T2 norm2 = zero2;
    retCode = tlapack::eispack_hqr(true, true, ilo, igh, H2, s2, Q2, norm2);
    timeHqr += getTime();
    std::cout << "Double precision\n";
    // Print out the time for computing both the eigenvalues, T, and the schur vectors
    std::cout << "TQE:" << timeHqr << "\n";
    // copy over everything back to the original matrices
    lacpy(Uplo::General, Hs2, H2);
    lacpy(Uplo::General, Qs2, Q2);
    timeHqr = -getTime();
    retCode = tlapack::eispack_hqr(true, false, ilo, igh, H2, s2, Q2, norm2);
    timeHqr += getTime();
    std::cout << "TE:" << timeHqr << "\n";
    // copy over everything back to the original matrices
    lacpy(Uplo::General, Hs2, H2);
    lacpy(Uplo::General, Qs2, Q2);
    timeHqr = -getTime();
    retCode = tlapack::eispack_hqr(false, false, ilo, igh, H2, s2, Q2, norm2);
    timeHqr += getTime();
    std::cout << "E:" << timeHqr << "\n";
}
