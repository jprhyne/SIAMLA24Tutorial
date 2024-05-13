#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

int main(int argc, char *argv[]) {
    // integer variables
    int info, lda, ldq, m, n, k, lwork, nb, i, j;
    int workQuery = -1;
    // double variables
    double *A, *Q, *As, *tau, *work, *T, *wr, *wi=NULL;
    double normA, tmp;
    double perform_refL, elapsed_refL, norm_orth_1, norm_repres_1;
    double zero = 0;
    double one = 1;
    double negOne = -1;
    // struct for help with timing
    struct timeval tp;
    // character variables
    char aChar = 'A';
    char iChar = 'I';
    char sChar = 'S';

    // Dummy value that is used when calling fortran
    // functions that take characters from C
    size_t dummy = 0;

    bool timesOnly = false;

    n = 30;
    for(i = 1; i < argc; ++i){
        if( strcmp( *(argv + i), "-ldq") == 0) {
            ldq  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-lda") == 0) {
            lda  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-n") == 0) {
            n  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-nb") == 0) {
            nb  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-t") == 0) {
            timesOnly = true;
        }
    }

    if( lda < 0 ) lda = n;
    if( ldq < 0 ) ldq = n;
    m = n;
    k = n;

    // allocate memory for the matrices and vectors that we need to use
    A =   (double *) malloc(lda * k * sizeof(double));
    As =  (double *) malloc(lda * k * sizeof(double));
    Q =   (double *) malloc(ldq * k * sizeof(double));
    tau = (double *) malloc(k * sizeof(double));
    wr =  (double *) malloc(k * sizeof(double));
    wi =  (double *) malloc(k * sizeof(double));
    T =   (double *) malloc(n * n * sizeof(double));

    // Generate A as a random matrix
    for (i = 0; i < lda * k; i++)
        A[i] = (double) rand() / (double) (RAND_MAX) - 0.5e+00;
    // Store random data inside Q to ensure that we do not assume anything about Q
    for (i = 0; i < ldq * k; i++)
        Q[i] = (double) rand() / (double) (RAND_MAX) - 0.5e+00;

    // Store a copy of A inside As
    dlacpy_(&aChar, &m, &k, A, &lda, As, &lda, dummy);
    // Find the norm of A for use in later accuracy computations
    normA = 0.0;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            tmp = A[i + j*m];
            normA += tmp*tmp;
        }
    }
    normA = sqrt(normA);
    int ilo = 1;
    int ihi = n;
    // Create the work array to do workspace queries
    work = (double *) malloc(sizeof(double));
    // Determine how much workspace is needed for our operations
    // Check dgeqrf first
    dgehrd_(&n, &ilo, &ihi, A, &lda, tau, work, &workQuery, &info );
    lwork = work[0];

    dhseqr_ref_(&sChar, &iChar, &n, &ilo, &ihi, A, &n, wr, wi, Q, &n, work, &workQuery, &info);
    if (work[0] > lwork)
        lwork = work[0];

    // reallocate work to be of the right size
    work = (double *) realloc(work, lwork * sizeof(double));

    // Call dgeqrf first
    dgehrd_(&n, &ilo, &ihi, A, &lda, tau, work, &lwork, &info );

    // Copy A into Q for use with dorgkr
    dlacpy_(&aChar, &m, &k, A, &lda, Q, &ldq, dummy);

    // Take the current time for use with timing dorg2r
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

    dhseqr_ref_(&sChar, &iChar, &n, &ilo, &ihi, Q, &n, wr, wi, A, &n, work, &lwork, &info);

    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

    printf("Ref: %10.10e\n", elapsed_refL);
    free(Q);
    free(As);
    free(A);
    free(tau);
    free(work);
    free(T);
}
