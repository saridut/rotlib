/** @file */
#ifndef UTILS_MATH_H
#define UTILS_MATH_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <complex.h>

/**
 * @defgroup util_math Math utilities
 * @{
 */

/** @{ */

/*
 *  @brief Zeros out a block of memory of a given length.
 *
 *  `zero_out`, `izero_out`, and `zzero_out` operate on blocks of type
 *  `double`, `int`, and `double complex`, respectively.
 *
 * @param[in] n Number of elements in `v`. Must be > 0.
 * @param[out] v Writable storage for at least `n` elements of the
 * corresponding type.
 */
void zero_out(const int n, double* v);
void izero_out(const int n, int* v);
void zzero_out(const int n, double complex* v);

/** @} */

/**
 * @brief Checks whether two floating-point numbers are close within relative
 * or absolute tolerance.
 *
 * The comparison is symmetric. The absolute difference is compared with
 * `abs_tol` and with `rel_tol` times the magnitude of either operand. Meeting
 * any one of these criteria is sufficient.
 *
 * Equal values, including infinities of the same sign, are close. Unequal
 * non-finite values and comparisons involving NaN return `false`.
 *
 * @param[in] a A floating-point number.
 * @param[in] b Another floating-point number.
 * @param[in] rel_tol Relative tolerance. Must be >= 0.
 * @param[in] abs_tol Absolute tolerance. Must be >= 0.
 * @returns `true` if `a` and `b` are close, otherwise `false`.
 */
bool isclose(const double a, const double b,
        const double rel_tol, const double abs_tol);

/**
 * @brief Checks whether corresponding elements of two vectors are close
 * within relative or absolute tolerance.
 *
 * Each pair of elements is compared using isclose(). The input vectors may
 * overlap.
 *
 * @param[in] n Number of elements to compare in `a` and `b`. Must be > 0.
 * @param[in] a A vector of length `n`.
 * @param[in] b Another vector of length `n`.
 * @param[in] rel_tol Relative tolerance for each pair of elements. Must be >= 0.
 * @param[in] abs_tol Absolute tolerance for each pair of elements. Must be >= 0.
 * @returns `true` if every pair of corresponding elements is close,
 * otherwise `false`.
 */
bool allclose(const int n, const double* a, const double* b,
        const double rel_tol, const double abs_tol);

/**
 * @brief Calculates the Kronecker delta for indices `i` and `j`.
 *
 * Index ranges are not checked.
 *
 * @param[in] i Index in [0, 2].
 * @param[in] j Index in [0, 2].
 * @returns `1.0` if `i` equals `j`, otherwise `0.0`.
 */
double kron_delta(const int i, const int j);

/**
 * @brief Calculates the unnormalized sinc function, `sin(x)/x`, with its
 * value at zero defined as `1.0`.
 *
 * @param[in] x Finite argument in radian.
 * @returns Unnormalized sinc of `x`.
 */
double sinc(const double x);

/**
 * @brief Calculates the arithmetic mean of the elements of a vector.
 *
 * Uses Neumaier compensated summation to reduce rounding error. Intermediate
 * sums may overflow even when the mean is representable.
 *
 * @param[in] n Number of elements in `v`. Must be > 0.
 * @param[in] v Vector of length `n` with finite elements.
 * @returns Arithmetic mean of the elements of `v`.
 */
double vec_mean(const int n, const double* v);

/**
 * @brief Calculates the Euclidean (2-)norm of a vector.
 *
 * Uses maximum-component scaling before accumulating squares. The returned
 * norm may overflow if it exceeds the range of `double`.
 *
 * @param[in] n Number of elements in `a`. Must be > 0.
 * @param[in] a Vector of length `n` with finite elements.
 * @returns Euclidean norm of `a`, or `0.0` if `a` is a zero vector.
 */
double norm(const int n, const double* a);

/**
 * @brief Normalizes a vector in-place. In case of a zero vector, `a` is not
 * changed and `EXIT_FAILURE` is returned.
 *
 * @param[in] n Number of elements in `a`. Must be > 0.
 * @param[in,out] a `in:` A vector of length `n` with finite elements.\n
 * `out:` The normalized vector if `a` is non-zero, otherwise unchanged.
 * @returns `EXIT_FAILURE` if `a` is a zero vector, else `EXIT_SUCCESS`.
 */
int normalized(const int n, double* a);

/**
 * @brief Calculates the dot product of two vectors.
 *
 * The input vectors may overlap.
 *
 * @param[in] n Number of elements in `a` and `b`. Must be > 0.
 * @param[in] a A vector of length `n` with finite elements.
 * @param[in] b Another vector of length `n` with finite elements.
 * @returns Dot product of `a` and `b`.
 */
double dot(const int n, const double* a, const double* b);

/**
 * @brief Calculates the cross product `a x b` and puts the result in `c`.
 *
 * `c` must not overlap with `a` or `b`. The input vectors may overlap.
 *
 * @param[in] a A 3-vector with finite elements.
 * @param[in] b Another 3-vector with finite elements.
 * @param[out] c Cross product of `a` and `b`. Must provide writable storage
 * for at least 3 doubles.
 */
void cross(const double a[3], const double b[3], double c[3]);

/**
 * @brief Calculates the skew-symmetric cross product matrix of a vector.
 * Multiplying this matrix by a column vector `b` gives `a x b`.
 *
 * `ax` must not overlap with `a`.
 *
 * @param[in] a A 3-vector with finite elements.
 * @param[out] ax Cross product matrix of size `3 x 3`. Must provide writable
 * storage for at least 9 doubles.
 */
void cross_mat(const double a[3], double* ax);

/**
 * @brief Calculates the angle between two 3D vectors (not necessarily unit
 * vectors).
 *
 * @param[in] a (3,) vector
 * @param[in] b (3,) vector
 * @returns Angle between `a` and `b` in [0, *pi*].
 */
double vec_angle_3d(const double* a, const double* b);

/**
 * @brief Calculates the angle between two unit vectors.
 *
 * @param[in] n Number of elements in `u`. Must be equal to that in `v`.
 * @param[in] u Unit vector
 * @param[in] v Unit vector
 * @param[in] wrk Workspace vector of the same length as `u` (or `v`). It will
 * be overwritten and must not overlap with `u` or `v`.
 * @returns Angle between `u` and `v` in [0, *pi*].
 */
double uvec_angle(const int n, const double* u, const double* v,
        double* const wrk);

/**
 * @brief Calculates the scalar triple product `a . (b x c)`.
 *
 * The input vectors may overlap.
 *
 * @param[in] a A 3-vector with finite elements.
 * @param[in] b Another 3-vector with finite elements.
 * @param[in] c Another 3-vector with finite elements.
 * @returns Scalar triple product of `a`, `b`, and `c`.
 */
double scalar_triple_product(const double a[3], const double b[3],
        const double c[3]);

/**
 * @brief Calculates the vector triple product `a x (b x c)` and puts the
 * result in `d`.
 *
 * The input vectors may overlap. `d` may coincide exactly with any input
 * vector, but must not partially overlap an input vector.
 *
 * @param[in] a A 3-vector with finite elements.
 * @param[in] b Another 3-vector with finite elements.
 * @param[in] c Another 3-vector with finite elements.
 * @param[out] d Vector triple product. Must provide writable storage for
 * at least 3 doubles.
 */
void vector_triple_product(const double a[3], const double b[3],
        const double c[3], double d[3]);

/**
 * @brief Calculates the outer product of two vectors and puts the result in
 * the matrix `c`, with entries `c_ij = a_i * b_j`.
 *
 * `c` must not overlap with `a` or `b`. The input vectors may overlap.
 *
 * @param[in] na Number of elements in `a`. Must be > 0.
 * @param[in] nb Number of elements in `b`. Must be > 0.
 * @param[in] a A vector of length `na` with finite elements.
 * @param[in] b A vector of length `nb` with finite elements.
 * @param[out] c Outer product matrix of size `na x nb`. Must provide
 * writable storage for at least `na*nb` doubles.
 */
void outer(const int na, const int nb, const double* a,
        const double* b, double* c);

/**
 * @brief Calculates the determinant of a square matrix of order 2, 3, or 4.
 *
 * @param[in] n Number of rows and columns of `A`. Must be 2, 3, or 4.
 * @param[in] A Matrix of size `n x n` with finite elements.
 * @returns Determinant of `A`.
 */
double det(const int n, const double* A);

/**
 * @brief Calculates the trace of a square matrix, i.e. the sum of its
 * diagonal elements.
 *
 * @param[in] n Number of rows and columns of `A`. Must be > 0.
 * @param[in] A Matrix of size `n x n` with finite elements.
 * @returns Trace of `A`.
 */
double trace(const int n, const double* A);

/**
 * @brief Calculates the arithmetic mean of a matrix along dimension `dim`.
 *
 * If `dim` is 0, averages over rows to give one mean per column. If `dim` is
 * 1, averages over columns to give one mean per row.
 *
 * `out` must not overlap with `A`.
 *
 * @param[in] nrows Number of rows of `A`. Must be > 0.
 * @param[in] ncols Number of columns of `A`. Must be > 0.
 * @param[in] A Matrix of size `nrows x ncols` with finite elements.
 * @param[in] dim Dimension along which to average. Must be 0 or 1.
 * @param[out] out Vector of means. Must provide writable storage for at least
 * `ncols` doubles if `dim` is 0, or `nrows` doubles if `dim` is 1.
 */
void mat_mean(const int nrows, const int ncols, const double* A,
        const int dim, double* out);

/**
 * @brief Generates an identity matrix of size `n x n`.
 *
 * @param[in] n Number of rows and columns of `a`. Must be > 0.
 * @param[out] a Identity matrix. Must provide writable storage for at least
 * `n*n` doubles.
 */
void eye(const int n, double* a);

/**
 * @brief Multiplies matrix `A` and vector `b` according to one of the following
 * forms: (1) Ab or (2) A^Tb and puts the result in `c`.
 *
 * Dimensions of `A` describe stored matrices, e.g. `nrows` and
 * `ncols` are the dimensions of `A` before transposition. Even if `transA` is
 * `true`, the stored matrix is not changed, rather different strides are used
 * to handle the transposed case.
 *
 *  If `transA` is `false`, `b` must be of size `ncols` and `c` of size 
 * `nrows`. If `transA` is `true`, `b` must be of size `nrows` and `c` of size
 * `ncols`. 
 *
 * `c` must not overlap with `A` or `b` and must provide writeable storage for
 * at least `nrows` (or `ncols`) doubles, depending on `transA`.
 *
 * @param[in] nrows Number of rows of `A`. Must be > 0.
 * @param[in] ncols Number of columns of `A`. Must be > 0.
 * @param[in] A Matrix of size `nrows x ncols`.
 * @param[in] transA `true` if A should be transposed, `false` otherwise.
 * @param[in] b Vector of size `ncols` if `transA` is `false`, otherwise of
 * size `nrows`.
 * @param[out] c Result of the multiplication. Must be of size `nrows` if
 * `transA` is `false`, otherwise of size `ncols`.
 */
void mv_mult(const int nrows, const int ncols, const double* A,
        const bool transA, const double* b, double* c);

/**
 * @brief Multiplies two matrices `A` and `B` according to one of the following
 * forms: (1) AB, (2) A^TB, (3) AB^T, or (4) A^TB^T and puts the result in `C`.
 *
 * Dimensions of `A` and `B` describe stored matrices, e.g. `nrows_A` and
 * `ncols_A` are the dimensions of `A` before transposition. Even if `transA` is
 * `true`, the stored matrix is not changed, rather different strides are used
 * to handle the transposed case.
 *
 * After applying `transA`, `A` is interpreted to be an `m x nA` array.
 * Similarly, after applying `transB`, `B` is viewed as an `nB x p` array.
 * Compatilibity requires that `nA = nB`. Calls with non-compatible matrices
 * will print an error message and abort.
 *
 * `C` must not overlap with `A` or `B` and must provide writeable storage for
 * at least `m*p` doubles, where `m` and `p` depend on the transpose flags.
 *
 * @param[in] nrows_A Number of rows of `A`. Must be > 0.
 * @param[in] ncols_A Number of columns of `A`. Must be > 0.
 * @param[in] A Matrix of size `nrows_A x ncols_A`.
 * @param[in] transA `true` if A should be transposed, `false` otherwise.
 * @param[in] nrows_B Number of rows of `B`. Must be > 0.
 * @param[in] ncols_B Number of columns of `B`. Must be > 0.
 * @param[in] B Matrix of size `nrows_B x ncols_B`.
 * @param[in] transB `true` if B should be transposed, `false` otherwise.
 * @param[out] C Result of the multiplication.
 */
void mm_mult(const int nrows_A, const int ncols_A, const double* A,
        const bool transA, const int nrows_B, const int ncols_B,
        const double* B, const bool transB, double* C);

/**
 * @brief Transposes matrix `A` and puts the result in `B`.
 *
 * `B` must not overlap with `A`.
 *
 * @param[in] nrows Number of rows of `A`. Must be > 0.
 * @param[in] ncols Number of columns of `A`. Must be > 0.
 * @param[in] A Matrix of size `nrows x ncols` with finite elements.
 * @param[out] B Transposed matrix of size `ncols x nrows`. Must provide
 * writable storage for at least `nrows*ncols` doubles.
 */
void transpose(const int nrows, const int ncols, const double* A, double* B);

/**
 * @brief Adds a square matrix and its transpose in-place.
 *
 * The result is symmetric, and the diagonal elements are doubled.
 *
 * @param[in] n Number of rows and columns of `A`. Must be > 0.
 * @param[in,out] A `in:` A matrix of size `n x n` with finite elements.\n
 * `out:` The sum of the original matrix and its transpose.
 */
void add_transpose(const int n, double* A);

/**
 * @brief Multiplies a matrix by its transpose and puts the result `A A^T`
 * in `B`.
 *
 * `B` must not overlap with `A`.
 *
 * @param[in] nrows Number of rows of `A`. Must be > 0.
 * @param[in] ncols Number of columns of `A`. Must be > 0.
 * @param[in] A Matrix of size `nrows x ncols` with finite elements.
 * @param[out] B Product matrix of size `nrows x nrows`. Must provide writable
 * storage for at least `nrows*nrows` doubles.
 */
void mult_transpose(const int nrows, const int ncols, const double* A,
        double* B);

/**
 * @}
 */

#ifdef __cplusplus
}
#endif

#endif // UTILS_MATH_H
