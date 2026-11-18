#include "utils_math.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/******************************************************************************/

void zero_out(const int n, double* v){

    memset(v, 0, n*sizeof(double));
}

void izero_out(const int n, int* v){

    memset(v, 0, n*sizeof(int));
}

void zzero_out(const int n, double complex* v){

    memset(v, 0, n*sizeof(double complex));
}

/******************************************************************************/

/* Check if two floating point numbers of type double are close within */
/* tolerance.                                                          */
bool isclose(const double a, const double b,
        const double rel_tol, const double abs_tol){

    if (a == b){
        return true;
    }

    if ( (!isfinite(a)) || (!isfinite(b)) ){
        return false;
    }
    else {
        double diff = fabs(b-a);
        bool ret = ( diff <= fabs(rel_tol*b)
                    || diff <= fabs(rel_tol*a)
                    || diff <= abs_tol );
        return ret;
    }
}

/******************************************************************************/

/* Check if two arrays of type double are elementwise close within tolerance. */
bool allclose(const int n, const double* a, const double* b,
        const double rel_tol, const double abs_tol){

    bool ret;
    bool tmp;

    ret = true;
    for (int i=0; i<n; ++i){
        tmp = isclose(a[i], b[i], rel_tol, abs_tol);
        if (!tmp){
            ret = false;
            break;
        }
    }
    return ret;
}

/******************************************************************************/

/* Kronecker delta function, 0 <= i,j <= 2 */
double kron_delta(const int i, const int j){
    if (i == j){
        return 1.0;
    }
    else
    {
        return 0.0;
    }
}

/******************************************************************************/

/* Unnormalized sinc function, f(x) = sin(x)/x */
double sinc(const double x){
    // At |x|=1e-4, the first omitted term from the Taylor expansion is 
    // of the order of 1e-18, comfortably below double-precision rounding error
    // near 1.
    if (fabs(x) < 1e-4) {
        const double x2 = x * x;
        return 1.0 - x2 / 6.0;
    }
    else {
        return sin(x) / x;
    }
}

/******************************************************************************/

/* Calculates the mean of a vector based on Neumaier summation */
double vec_mean(const int n, const double* v){

    double sum = 0.0;
    double c = 0.0;
    double t, x;

    for (int i=0; i < n; ++i){
        x = v[i];
        t = sum + x;
        if (fabs(sum) >= fabs(x)){
            c += (sum - t) + x;
        }
        else {
            c += (x - t) + sum;
        }
        sum = t;
    }
    sum += c;
    return sum/n ;
}

/******************************************************************************/

/* Calculates the 2-norm of a vector */
double norm(const int n, const double* a){

    /* scaling factor */
    double scale = 0.0;
    for (int i=0; i<n; ++i){
        scale = fmax(scale, fabs(a[i]));
    }
    if (scale == 0.0){
        return 0.0; /* norm of a zero vector is zero */
    } 
    /* calculate the square of the scaled norm */
    double snorm = 0.0;
    double tmp = 0.0;
    for (int i=0; i<n; ++i){
        tmp = a[i]/scale;
        snorm += (tmp*tmp);
    }
    /* unscale and return */
    return scale*sqrt(snorm);
}

/******************************************************************************/

/* Normalizes a vector in-place */
int normalize(const int n, double* a){

    /* scaling factor */
    double scale = 0.0;
    for (int i=0; i<n; ++i){
        scale = fmax(scale, fabs(a[i]));
    }
    if (scale == 0.0){
        return EXIT_FAILURE; /* cannot normalize a zero vector */
    } 

    /* scale the vector */
    for (int i=0; i<n; ++i){
        a[i] /= scale;
    }

    double snorm = sqrt(dot(n, a, a)); /* scaled norm */

    /* divide by the scaled norm */
    for (int i=0; i<n; ++i){
        a[i] /= snorm;
    }
    return EXIT_SUCCESS;
}

/******************************************************************************/

double dot(const int n, const double* a, const double* b){

    double ret;

    ret = 0.0;
    for (int i=0; i<n; ++i){
        ret += a[i]*b[i];
    }
    return ret;
}

/******************************************************************************/

void cross(const double a[3], const double b[3], double c[3]){

    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

/******************************************************************************/

void cross_mat(const double a[3], double* ax){

    const int n = 3;

    ax[n*0+0] =   0.0;
    ax[n*0+1] = -a[2];
    ax[n*0+2] =  a[1];
    ax[n*1+0] =  a[2];
    ax[n*1+1] =   0.0;
    ax[n*1+2] = -a[0];
    ax[n*2+0] = -a[1];
    ax[n*2+1] =  a[0];
    ax[n*2+2] =   0.0;
}

/******************************************************************************/

double vec_angle_3d(const double a[3], const double b[3]){

    double a_[3], b_[3], c[3], s;

    s = fmax( fabs(a[0]), fabs(a[1]) );
    s = fmax( s, fabs(a[2]) );
    a_[0] = a[0]/s; a_[1] = a[1]/s; a_[2] = a[2]/s;

    s = fmax( fabs(b[0]), fabs(b[1]) );
    s = fmax( s, fabs(b[2]) );
    b_[0] = b[0]/s; b_[1] = b[1]/s; b_[2] = b[2]/s;

    cross(a_, b_, c);
    double c_mag = norm(3, c);

    double ab = dot(3, a_, b_);
    return atan2(c_mag, ab);
}

/******************************************************************************/

double uvec_angle(const int n, const double* u, const double* v, 
                  double* const wrk){

    for(int i=0; i<n; ++i){
        wrk[i] = u[i] + v[i];
    }
    double s_mag = norm(n, wrk);

    for(int i=0; i<n; ++i){
        wrk[i] = u[i] - v[i];
    }
    double d_mag = norm(n, wrk);

    return 2.0*atan2(d_mag, s_mag);
}

/******************************************************************************/

/* Calculates a.(b x c) */
double scalar_triple_product(const double a[3], const double b[3],
        const double c[3]){

    double d[3];
    double ret;

    cross(b, c, d);
    ret = a[0]*d[0] + a[1]*d[1] + a[2]*d[2]; /* a.d */
    return ret;
}

/******************************************************************************/

/* Calculates a x (b x c) */
void vector_triple_product(const double a[3], const double b[3],
        const double c[3], double d[3]){

    double ac;
    double ab;

    ac = a[0]*c[0] + a[1]*c[1] + a[2]*c[2]; /* a.c */
    ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; /* a.b */

    /* a x (b x c) = b (a.c) - c (a.b) */
    d[0] = b[0]*ac - c[0]*ab;
    d[1] = b[1]*ac - c[1]*ab;
    d[2] = b[2]*ac - c[2]*ab;
}

/******************************************************************************/

void outer(const int na, const int nb, const double* a,
        const double* b, double* c){

    for (int i=0; i<na; ++i){
        for (int j=0; j<nb; ++j){
            c[i*nb+j] = a[i]*b[j];
        }
    }
}

/******************************************************************************/

/* See David Simpson's fortran code at        */
/* http://www.davidgsimpson.com/software.html */
/* for formulas for 5x5 and 6x6 determinants. */
/* For a general NxN matrix do an LU decomp.  */

double det(const int n, const double* A){

    double ret = 0.0;

    switch (n) {
        case 2 :
            ret = A[0]*A[3] - A[1]*A[2];
            break;
        case 3 :
            ret = A[0]*A[4]*A[8] - A[0]*A[5]*A[7]
                - A[1]*A[3]*A[8] + A[1]*A[5]*A[6]
                + A[2]*A[3]*A[7] - A[2]*A[4]*A[6] ;
            break;
        case 4 :
            ret = A[0] * (A[5] * (A[10]*A[15] - A[11]*A[14])
                        + A[6] * (A[11]*A[13] - A[ 9]*A[15])
                        + A[7] * (A[ 9]*A[14] - A[10]*A[13]))
                - A[1] * (A[4] * (A[10]*A[15] - A[11]*A[14])
                        + A[6] * (A[11]*A[12] - A[ 8]*A[15])
                        + A[7] * (A[ 8]*A[14] - A[10]*A[12]))
                + A[2] * (A[4] * (A[ 9]*A[15] - A[11]*A[13])
                        + A[5] * (A[11]*A[12] - A[ 8]*A[15])
                        + A[7] * (A[ 8]*A[13] - A[ 9]*A[12]))
                - A[3] * (A[4] * (A[ 9]*A[14] - A[10]*A[13])
                        + A[5] * (A[10]*A[12] - A[ 8]*A[14])
                        + A[6] * (A[ 8]*A[13] - A[ 9]*A[12])) ;
            break;
    }
    return ret;
}

/******************************************************************************/

/* Returns sum of the elements along the diagonal of an n x n square matrix */
double trace(const int n, const double* A){

    double ret = 0.0;

    for (int i=0; i<n; ++i){
        ret += A[i*n+i];
    }
    return ret;
}

/******************************************************************************/

/* Calculates the mean of a matrix along dimension dim */
void mat_mean(const int nrows, const int ncols, const double* A,
        const int dim, double* out){

    if (dim == 0){
        for (int j=0; j < ncols; ++j){
            out[j] = 0.0;
        }
        for (int i=0; i < nrows; ++i){
            for (int j=0; j < ncols; ++j){
                out[j] += A[ncols*i+j];
            }
        }
        for (int j=0; j < ncols; ++j){
            out[j] /= nrows;
        }
    }
    else if (dim == 1){
        for (int i=0; i < nrows; ++i){
            out[i] = 0.0;
        }
        for (int i=0; i < nrows; ++i){
            for (int j=0; j < ncols; ++j){
                out[i] += A[ncols*i+j];
            }
        }
        for (int i=0; i < nrows; ++i){
            out[i] /= ncols;
        }
    }
}

/******************************************************************************/

/* Generates an identity matrix of size n x n */
void eye(const int n, double* a){

    /* zero out all elements */
    for (int i=0; i<n; ++i){
        for (int j=0; j<n; ++j){
            a[i*n+j] = 0.0;
        }
    }
    /* set diagonal elements to 1.0 */
    for (int i=0; i<n; ++i){
        a[i*n+i] = 1.0;
    }
}

/******************************************************************************/

/* Performs a matrix vector multiplication: c_i = A_ij b_j */
void mv_mult(const int nrows, const int ncols, const double* A,
        const bool transA, const double* b, double* c){

    // If transA is true, the lengths of b must be equal to nrows and c must be
    // equal to at least ncols.
    // After applying transA, A is interpreted as an m x n matrix.
    int m = transA ? ncols : nrows;
    int n = transA ? nrows : ncols;

    // Row and column strides
    const int ar = transA ? 1 : n;
    const int ac = transA ? m : 1;

    zero_out(m, c); // Zero out the output buffer

    for (int i=0; i<m; ++i){
        for (int j=0; j<n; ++j){
            c[i] += A[ar*i + ac*j] * b[j];
        }
    }
}

/******************************************************************************/

/* Performs a matrix matrix multiplication: C_ij = A_ik B_kj */
void mm_mult(const int nrows_A, const int ncols_A, const double* A,
        const bool transA, const int nrows_B, const int ncols_B,
        const double* B, const bool transB, double* C){

    // On entry, A is nrows_A x ncols_A and B is nrows_B x ncols_B
    // After transposing, A is m x n, B is n x p, C is m x p
    int m = transA ? ncols_A : nrows_A;
    int nA = transA ? nrows_A : ncols_A;

    int nB = transB ? ncols_B : nrows_B;
    int p = transB ? nrows_B : ncols_B;

    // Compatibility condition is nA==nB
    if (nA != nB){
        printf("Incompatible matrices."
               " `nA`(= %d) must equal `nB`(=%d).\n", nA, nB);
        abort();
    }
    int n = nA; //equal to nB

    // Row and column strides
    const int ar = transA ? 1 : n;
    const int ac = transA ? m : 1;
    const int br = transB ? 1 : p;
    const int bc = transB ? n : 1; 
    
    zero_out(m*p, C); // Zero out the output buffer

    for (int i=0; i<m; ++i){
        for (int j=0; j<p; ++j){
            for (int k=0; k<n; ++k){
                C[p*i+j] += A[ar*i + ac*k] * B[br*k + bc*j];
            }
        }
    }
}

/******************************************************************************/

/* Transposes a matrix: B = transpose(A) */
void transpose(const int nrows, const int ncols, const double* A, double* B){

    for (int i=0; i<nrows; ++i){
        for (int j=0; j<ncols; ++j){
            B[nrows*j+i] = A[ncols*i+j];
        }
    }
}

/******************************************************************************/

/* Adds a square matrix and its transpose in place: Aij = Aij + Aji */
void add_transpose(const int n, double* A){

    /* Dealing with the upper triangular part (including the diagonal) */
    for (int i=0; i < n; ++i){
        for (int j=i; j < n; ++j){
            A[n*i+j] += A[n*j+i];
        }
    }
    /* Dealing with the strictly lower triangular part */
    for (int i=0; i < n; ++i){
        for (int j=0; j < i; ++j){
            A[n*i+j] = A[n*j+i];
        }
    }
}

/******************************************************************************/

/* Multiplies a matrix with its transpose: B = A * A^T */
void mult_transpose(const int nrows, const int ncols, const double* A,
        double* B){

    int m = nrows;
    int n = ncols;

    /* A is m x n, A^T is n x m, so B is m x m */
    for (int i=0; i<m; ++i){
        for (int j=0; j<m; ++j){
            B[m*i+j] = 0.0;
        }
    }

    for (int i=0; i<m; ++i){
        for (int j=0; j<m; ++j){
            for (int k=0; k<n; ++k){
                B[m*i+j] += A[n*i+k]*A[n*j+k];
            }
        }
    }
}

/******************************************************************************/
