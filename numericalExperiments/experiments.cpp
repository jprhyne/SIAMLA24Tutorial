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
template <typename matrix_t>
void computeTimes(std::size_t n)
{
    using T = type_t<matrix_t>;
    using idx_t = size_type<matrix_t>;
    using real_t = real_type<T>;
    using complex_t = std::complex<T>;
    // Create functors
    Create<matrix_t> new_matrix;
    // Create matrices
    std::vector<T> A_; auto A = new_matrix( A_, n, n);
    std::vector<T> H_; auto H = new_matrix( H_, n, n);
    std::vector<T> Hs_; auto Hs = new_matrix( Hs_, n, n);
    std::vector<T> Q_; auto Q = new_matrix( Q_, n, n);
    std::vector<T> Qs_; auto Qs = new_matrix( Qs_, n, n);
    std::vector<complex_t> s(n);
    //Populate A and U with random numbers
    rand_generator gen;
    idx_t seed = 17519;
    gen.seed(seed);
    for(idx_t i = 0; i < n; i++) {
        for(idx_t j = 0; j < n; j++) {
            T val = rand_helper<T>(gen);
            A(i,j) = val; 
            H(i,j) = val; 
        }
    }
    std::vector<T> tau1(n);
    T zeroT = T(0);
    real_t zero = real_t(0);
    gehrd(0, n, H, tau1);
    lacpy(Uplo::General, H, Q);
    lacpy(Uplo::General, H, Hs);
    unghr(0, n, Q, tau1);
    lacpy(Uplo::General, Q, Qs);
    idx_t ilo = 0;
    idx_t igh = n;
    for (idx_t i = 1; i < n; i++)
        for (idx_t j = 0; j < i - 1; j++)
            H(i,j) = zeroT;
    double timeHqr = -getTime();
    real_t norm1 = zero;
    auto retCode = tlapack::eispack_hqr(true, true, ilo, igh, H, s, Q, norm1);
    timeHqr += getTime();
    // Print out the time for computing the eigenvalues, T, and the schur vectors
    std::cout << "TQE:" << timeHqr << "\n";
    // copy over everything back to the original matrices
    lacpy(Uplo::General, Hs, H);
    lacpy(Uplo::General, Qs, Q);
    timeHqr = -getTime();
    retCode = tlapack::eispack_hqr(true, false, ilo, igh, H, s, Q, norm1);
    timeHqr += getTime();
    // Print out the time for computing the eigenvalues and T
    std::cout << "TE:" << timeHqr << "\n";
    // copy over everything back to the original matrices
    lacpy(Uplo::General, Hs, H);
    lacpy(Uplo::General, Qs, Q);
    timeHqr = -getTime();
    retCode = tlapack::eispack_hqr(false, false, ilo, igh, H, s, Q, norm1);
    timeHqr += getTime();
    // Print out the time for computing just the eigenvalues
    std::cout << "E:" << timeHqr << "\n";
    return;
}

int main( int argc, char **argv)
{
    std::size_t i, n;
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
    /**********************************************************************************************\
    |* Single precision section                                                                   *|
    \**********************************************************************************************/
    std::cout << "Single precision" << std::endl;
    computeTimes<legacyMatrix<float,std::size_t,Layout::ColMajor>>(n);
    /**********************************************************************************************\
    |* Double precision section                                                                   *|
    \**********************************************************************************************/
    std::cout << "Double precision" << std::endl;
    computeTimes<legacyMatrix<double,std::size_t,Layout::ColMajor>>(n);
    // Other precisions can/will follow
}
