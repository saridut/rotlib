#include "rotlib.h"
#include "utils_math.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Euler angle sequence: XYZ (world). First rotation about X, second rotation   */
/* about Y, and the third rotation about Z axis of the world(i.e. fixed) frame. */


/******************************************************************************/
                       /* INTERNAL ROUTINES */
/******************************************************************************/

static void rotmat_XYZ(const double* euler, double* rotmat){

    const int n = 3;
    double phi, theta, psi;
    double sin_phi, sin_theta, sin_psi;
    double cos_phi, cos_theta, cos_psi;

    phi = euler[0]; theta = euler[1]; psi = euler[2];

    sin_phi = sin(phi); sin_theta = sin(theta); sin_psi = sin(psi);
    cos_phi = cos(phi); cos_theta = cos(theta); cos_psi = cos(psi);

    rotmat[n*0+0] = cos_theta*cos_psi;
    rotmat[n*0+1] = sin_phi*sin_theta*cos_psi - cos_phi*sin_psi;
    rotmat[n*0+2] = cos_phi*sin_theta*cos_psi + sin_phi*sin_psi;
    rotmat[n*1+0] = cos_theta*sin_psi;
    rotmat[n*1+1] = sin_psi*sin_theta*sin_phi + cos_phi*cos_psi;
    rotmat[n*1+2] = cos_phi*sin_theta*sin_psi - sin_phi*cos_psi;
    rotmat[n*2+0] = -sin_theta;
    rotmat[n*2+1] = sin_phi*cos_theta;
    rotmat[n*2+2] = cos_phi*cos_theta;
}

/******************************************************************************/

static void factor_rotmat_XYZ(const double* rotmat, const double tol, 
        double* euler){

    const int n = 3;
    double phi, theta, psi;

    /* Recover the middle angle theta */
    double s = -rotmat[n*2+0];
    double c = hypot(rotmat[n*0+0], rotmat[n*1+0]); /* = abs(cos theta) */
    theta = atan2(s,c); /* Returns in [-pi/2, pi/2] */

    /* Calculate phi (accounting for gimbal lock) */
    if (c > tol){
        phi = atan2(rotmat[n*2+1], rotmat[n*2+2]);
    }
    else {
        phi = 0.0;
    }
    /* Calculate psi using phi */
    double sphi = sin(phi);
    double cphi = cos(phi);
    psi = atan2(sphi*rotmat[n*0+2] - cphi*rotmat[n*0+1], 
                cphi*rotmat[n*1+1] - sphi*rotmat[n*1+2]);

    euler[0] = phi; euler[1] = theta; euler[2] = psi;
}

/******************************************************************************/

static void rotmat_XZY(const double* euler, double* rotmat){

    const int n = 3;
    double phi, theta, psi;
    double sin_phi, sin_theta, sin_psi;
    double cos_phi, cos_theta, cos_psi;

    phi = euler[0]; theta = euler[1]; psi = euler[2];

    sin_phi = sin(phi); sin_theta = sin(theta); sin_psi = sin(psi);
    cos_phi = cos(phi); cos_theta = cos(theta); cos_psi = cos(psi);

    rotmat[n*0+0] = cos_theta*cos_psi;
    rotmat[n*0+1] = sin_phi*sin_theta - cos_phi*cos_theta*sin_psi;
    rotmat[n*0+2] = cos_phi*sin_theta + sin_phi*cos_theta*sin_psi;
    rotmat[n*1+0] = sin_psi;
    rotmat[n*1+1] = cos_phi*cos_psi;
    rotmat[n*1+2] = -sin_phi*cos_psi;
    rotmat[n*2+0] = -sin_theta*cos_psi;
    rotmat[n*2+1] = sin_phi*cos_theta + cos_phi*sin_theta*sin_psi;
    rotmat[n*2+2] = cos_phi*cos_theta - sin_phi*sin_theta*sin_psi;
}

/******************************************************************************/

static void factor_rotmat_XZY(const double* rotmat, const double tol,
        double* euler){

    const int n = 3;
    double phi, theta, psi;

    /* Recover the middle angle psi */
    double s = rotmat[n*1+0];
    double c = hypot(rotmat[n*0+0], rotmat[n*2+0]); /* = abs(cos psi) */
    psi = atan2(s,c); /* Returns in [-pi/2, pi/2] */

    /* Calculate phi (accounting for gimbal lock) */
    if (c > tol){
        phi = atan2(-rotmat[n*1+2], rotmat[n*1+1]);
    }
    else {
        phi = 0.0;
    }
    /* Calculate theta using phi */
    double sphi = sin(phi);
    double cphi = cos(phi);
    theta = atan2(sphi*rotmat[n*0+1] + cphi*rotmat[n*0+2], 
                  sphi*rotmat[n*2+1] + cphi*rotmat[n*2+2]);

    euler[0] = phi; euler[1] = theta; euler[2] = psi;
}

/******************************************************************************/

static void rotmat_YXZ(const double* euler, double* rotmat){

    const int n = 3;
    double phi, theta, psi;
    double sin_phi, sin_theta, sin_psi;
    double cos_phi, cos_theta, cos_psi;

    phi = euler[0]; theta = euler[1]; psi = euler[2];

    sin_phi = sin(phi); sin_theta = sin(theta); sin_psi = sin(psi);
    cos_phi = cos(phi); cos_theta = cos(theta); cos_psi = cos(psi);

    rotmat[n*0+0] = cos_theta*cos_psi - sin_phi*sin_theta*sin_psi;
    rotmat[n*0+1] = -cos_phi*sin_psi;
    rotmat[n*0+2] = sin_theta*cos_psi + sin_phi*cos_theta*sin_psi;
    rotmat[n*1+0] = sin_phi*sin_theta*cos_psi + cos_theta*sin_psi;
    rotmat[n*1+1] = cos_phi*cos_psi;
    rotmat[n*1+2] = sin_theta*sin_psi - sin_phi*cos_theta*cos_psi;
    rotmat[n*2+0] = -cos_phi*sin_theta;
    rotmat[n*2+1] = sin_phi;
    rotmat[n*2+2] = cos_phi*cos_theta;
}

/******************************************************************************/

static void factor_rotmat_YXZ(const double* rotmat, const double tol, 
        double* euler){

    const int n = 3;
    double phi, theta, psi;

    /* Recover the middle angle phi */
    double s = rotmat[n*2+1];
    double c = hypot(rotmat[n*0+1], rotmat[n*1+1]); /* = abs(cos phi) */
    phi = atan2(s,c); /* Returns in [-pi/2, pi/2] */
    
    /* Calculate theta (accounting for gimbal lock) */
    if (c > tol){
        theta = atan2(-rotmat[n*2+0], rotmat[n*2+2]);
    }
    else {
        theta = 0.0;
    }
    /* Calculate psi using theta */
    double stheta = sin(theta);
    double ctheta = cos(theta);
    psi = atan2(ctheta*rotmat[n*1+0] + stheta*rotmat[n*1+2], 
                ctheta*rotmat[n*0+0] + stheta*rotmat[n*0+2]);

    euler[0] = phi; euler[1] = theta; euler[2] = psi;
}

/******************************************************************************/

static void rotmat_YZX(const double* euler, double* rotmat){

    const int n = 3;
    double phi, theta, psi;
    double sin_phi, sin_theta, sin_psi;
    double cos_phi, cos_theta, cos_psi;

    phi = euler[0]; theta = euler[1]; psi = euler[2];

    sin_phi = sin(phi); sin_theta = sin(theta); sin_psi = sin(psi);
    cos_phi = cos(phi); cos_theta = cos(theta); cos_psi = cos(psi);

    rotmat[n*0+0] = cos_theta*cos_psi;
    rotmat[n*0+1] = -sin_psi;
    rotmat[n*0+2] = sin_theta*cos_psi;
    rotmat[n*1+0] = sin_phi*sin_theta + cos_phi*cos_theta*sin_psi;
    rotmat[n*1+1] = cos_phi*cos_psi;
    rotmat[n*1+2] = cos_phi*sin_theta*sin_psi - sin_phi*cos_theta;
    rotmat[n*2+0] = sin_phi*cos_theta*sin_psi - cos_phi*sin_theta;
    rotmat[n*2+1] = sin_phi*cos_psi;
    rotmat[n*2+2] = sin_phi*sin_theta*sin_psi + cos_phi*cos_theta;
}

/******************************************************************************/

static void factor_rotmat_YZX(const double* rotmat, const double tol,
        double* euler){

    const int n = 3;
    double phi, theta, psi;

    /* Recover the middle angle psi */
    double s = -rotmat[n*0+1];
    double c = hypot(rotmat[n*1+1], rotmat[n*2+1]); /* = abs(cos psi) */
    psi = atan2(s,c); /* Returns in [-pi/2, pi/2] */
    
    /* Calculate theta (accounting for gimbal lock) */
    if (c > tol){
        theta = atan2(rotmat[n*0+2], rotmat[n*0+0]);
    }
    else {
        theta = 0.0;
    }
    /* Calculate phi using theta */
    double stheta = sin(theta);
    double ctheta = cos(theta);
    phi = atan2(stheta*rotmat[n*1+0] - ctheta*rotmat[n*1+2], 
                ctheta*rotmat[n*2+2] - stheta*rotmat[n*2+0]);

    euler[0] = phi; euler[1] = theta; euler[2] = psi;
}

/******************************************************************************/

static void rotmat_ZXY(const double* euler, double* rotmat){

    const int n = 3;
    double phi, theta, psi;
    double sin_phi, sin_theta, sin_psi;
    double cos_phi, cos_theta, cos_psi;

    phi = euler[0]; theta = euler[1]; psi = euler[2];

    sin_phi = sin(phi); sin_theta = sin(theta); sin_psi = sin(psi);
    cos_phi = cos(phi); cos_theta = cos(theta); cos_psi = cos(psi);

    rotmat[n*0+0] = cos_theta*cos_psi + sin_phi*sin_theta*sin_psi;
    rotmat[n*0+1] = sin_phi*sin_theta*cos_psi - cos_theta*sin_psi;
    rotmat[n*0+2] = cos_phi*sin_theta;
    rotmat[n*1+0] = cos_phi*sin_psi;
    rotmat[n*1+1] = cos_phi*cos_psi;
    rotmat[n*1+2] = -sin_phi;
    rotmat[n*2+0] = sin_phi*cos_theta*sin_psi - sin_theta*cos_psi;
    rotmat[n*2+1] = sin_phi*cos_theta*cos_psi + sin_theta*sin_psi;
    rotmat[n*2+2] = cos_phi*cos_theta;
}

/******************************************************************************/

static void factor_rotmat_ZXY(const double* rotmat, const double tol,
        double* euler){

    const int n = 3;
    double phi, theta, psi;

    /* Recover the middle angle phi */
    double s = -rotmat[n*1+2];
    double c = hypot(rotmat[n*0+2], rotmat[n*2+2]); /* = abs(cos phi) */
    phi = atan2(s,c); /* Returns in [-pi/2, pi/2] */
    
    /* Calculate psi (accounting for gimbal lock) */
    if (c > tol){
        psi = atan2(rotmat[n*1+0], rotmat[n*1+1]);
    }
    else {
        psi = 0.0;
    }
    /* Calculate theta using psi */
    double spsi = sin(psi);
    double cpsi = cos(psi);
    theta = atan2(spsi*rotmat[n*2+1] - cpsi*rotmat[n*2+0], 
                  cpsi*rotmat[n*0+0] - spsi*rotmat[n*0+1]);

    euler[0] = phi; euler[1] = theta; euler[2] = psi;
}

/******************************************************************************/

static void rotmat_ZYX(const double* euler, double* rotmat){

    const int n = 3;
    double phi, theta, psi;
    double sin_phi, sin_theta, sin_psi;
    double cos_phi, cos_theta, cos_psi;

    phi = euler[0]; theta = euler[1]; psi = euler[2];

    sin_phi = sin(phi); sin_theta = sin(theta); sin_psi = sin(psi);
    cos_phi = cos(phi); cos_theta = cos(theta); cos_psi = cos(psi);

    rotmat[n*0+0] = cos_theta*cos_psi;
    rotmat[n*0+1] = -cos_theta*sin_psi;
    rotmat[n*0+2] = sin_theta;
    rotmat[n*1+0] = sin_phi*sin_theta*cos_psi + cos_phi*sin_psi;
    rotmat[n*1+1] = cos_phi*cos_psi - sin_phi*sin_theta*sin_psi;
    rotmat[n*1+2] = -sin_phi*cos_theta;
    rotmat[n*2+0] = sin_phi*sin_psi - cos_phi*sin_theta*cos_psi;
    rotmat[n*2+1] = sin_phi*cos_psi + cos_phi*sin_theta*sin_psi;
    rotmat[n*2+2] = cos_phi*cos_theta;
}

/******************************************************************************/

static void factor_rotmat_ZYX(const double* rotmat, const double tol,
        double* euler){

    const int n = 3;
    double phi, theta, psi;

    /* Recover the middle angle theta */
    double s = rotmat[n*0+2];
    double c = hypot(rotmat[n*1+2], rotmat[n*2+2]); /* = abs(cos theta) */
    theta = atan2(s,c); /* Returns in [-pi/2, pi/2] */
    
    /* Calculate psi (accounting for gimbal lock) */
    if (c > tol){
        psi = atan2(-rotmat[n*0+1], rotmat[n*0+0]);
    }
    else {
        psi = 0.0;
    }
    /* Calculate phi using psi */
    double spsi = sin(psi);
    double cpsi = cos(psi);
    phi = atan2(spsi*rotmat[n*2+0] + cpsi*rotmat[n*2+1], 
                spsi*rotmat[n*1+0] + cpsi*rotmat[n*1+1]);

    euler[0] = phi; euler[1] = theta; euler[2] = psi;
}

/******************************************************************************/

static void rotmat_euler(const double* euler, enum EulangSeq seq, const bool world,
        double* rotmat){

    double euler_[3];
    double rotmat_[9]; /* 3 x 3 matrix as an 1-D array */

    if (!world){
        for (int i=0; i<3; ++i){
            euler_[i] = -euler[i];
        }
    }
    else {
        for (int i=0; i<3; ++i){
            euler_[i] = euler[i];
        }
    }

    switch (seq){
        case ES_XYZ :
            rotmat_XYZ(euler_, rotmat_);
            break;
        case ES_XZY :
            rotmat_XZY(euler_, rotmat_);
            break;
        case ES_YXZ :
            rotmat_YXZ(euler_, rotmat_);
            break;
        case ES_YZX :
            rotmat_YZX(euler_, rotmat_);
            break;
        case ES_ZXY :
            rotmat_ZXY(euler_, rotmat_);
            break;
        case ES_ZYX :
            rotmat_ZYX(euler_, rotmat_);
            break;
    }

    if (!world) {
        transpose(3, 3, rotmat_, rotmat);
    }
    else {
        memcpy(rotmat, rotmat_, 9*sizeof(double));
    }
}

/******************************************************************************/

static void factor_rotmat(const double* rotmat, enum EulangSeq seq,
        const bool world, const double tol, double* euler){

    double rotmat_[9]; /* 3 x 3 matrix as an 1-D array */

    if (!world) {
        transpose(3, 3, rotmat, rotmat_);
    }
    else {
        memcpy(rotmat_, rotmat, 9*sizeof(double));
    }

    switch (seq){
        case ES_XYZ :
            factor_rotmat_XYZ(rotmat_, tol, euler);
            break;
        case ES_XZY :
            factor_rotmat_XZY(rotmat_, tol, euler);
            break;
        case ES_YXZ :
            factor_rotmat_YXZ(rotmat_, tol, euler);
            break;
        case ES_YZX :
            factor_rotmat_YZX(rotmat_, tol, euler);
            break;
        case ES_ZXY :
            factor_rotmat_ZXY(rotmat_, tol, euler);
            break;
        case ES_ZYX :
            factor_rotmat_ZYX(rotmat_, tol, euler);
            break;
    }

    if (!world){
        for (int i=0; i<3; ++i){
            euler[i] = -euler[i];
        }
    }
}

/******************************************************************************/

static void rotate_vectors (const int n, const double* v, const double* rotmat,
        double* vrot) {

    mm_mult(n, 3, v, false, 3, 3, rotmat, true, vrot);
}

/******************************************************************************/

static void shift_vectors (const int n, const double* v, const double* shiftmat,
        double* vshift){

    mm_mult(n, 3, v, false, 3, 3, shiftmat, true, vshift);
}

/******************************************************************************/

static void shift_tensor2 (const double* A, const double* shiftmat,
        double* Ashift){

    const int n=3;

    zero_out(n*n, Ashift); // Zero out the output buffer

    for (int i=0; i<n; ++i){
        for (int j=0; j<n; ++j){
            for (int p=0; p<n; ++p){
                for (int q=0; q<n; ++q){
                    Ashift[n*i+j] += shiftmat[n*i+p]*shiftmat[n*j+q]*A[n*p+q];
                }
            }
        }
    }

}

/******************************************************************************/

static void shift_tensor3 (const double* A, const double* shiftmat,
        double* Ashift){

    const int n = 3;

    zero_out(n*n*n, Ashift); // Zero out the output buffer

    for (int i=0; i<n; ++i){
        for (int j=0; j<n; ++j){
            for (int k=0; k<n; ++k){
                for (int p=0; p<n; ++p){
                    for (int q=0; q<n; ++q){
                        for (int r=0; r<n; ++r){
                            Ashift[n*n*i+n*j+k] 
                                += shiftmat[n*i+p]*shiftmat[n*j+q]*shiftmat[n*k+r]
                                *A[n*n*p+n*q+r];
                        }
                    }
                }
            }
        }
    }

}

/******************************************************************************/
                            /* UTILITY FUNCTIONS */
/******************************************************************************/

void align(const int m, const double* v, const int n, const double* old, 
        const double* new, double* valigned){

    double sin_angle, cos_angle, angle;
    double axis[3], z_old[3], z_new[3];
    double mat[9], axes_old[9], axes_new[9];
    int estat;

    if (n == 1) {
        //Angle between old and new
        cos_angle = dot(3, old, new);
        cross(old, new, axis);
        sin_angle = norm(3, axis); // always non-negative
        // angle will be in [0, pi] as sin_angle >= 0
        angle = atan2(sin_angle, cos_angle);
        // try normalizing the axis
        estat = normalize(3, axis);
        if (estat == EXIT_SUCCESS){
            // Non-collinear vectors: standard case
            aa_rotate_vectors(m, v, axis, angle, valigned);
        } 
        else {
            // Collinear vectors: parallel or anti-parallel
            if (cos_angle > 0.0){
                // angle = 0; no rotation necessary
                memcpy(valigned, v, 3*m*sizeof(double));
            }
            else {
                //angle = pi: Rotation axis not unique. Pick a vector by
                //crossing old with a cartesian direction least parallel to it.  
                //
                //1. Find cartesian direction least parallel to old.
                int min_idx = fabs(old[0]) < fabs(old[1]) ? 0 : 1;
                min_idx = fabs(old[min_idx]) < fabs(old[2]) ? min_idx : 2;

                //2. Create that basis vector
                double vtmp[3] = {0.0, 0.0, 0.0};
                vtmp[min_idx] = 1.0;

                //3. Cross old with vtmp and normalize to get the rotation axis
                cross(old, vtmp, axis); //axis will always have +ve norm
                aa_fix(axis, &angle, true);

                //Rotate
                aa_rotate_vectors(m, v, axis, angle, valigned);
            }
        }
    }
    else if (n == 2) {
        //The third vector orthogonal to `old`. 
        cross(old, old+3, z_old);
        estat = normalize(3, z_old); // Norm always +ve
        //The third vector orthogonal to `new`
        cross(new, new+3, z_new);
        estat = normalize(3, z_new); // Norm always +ve
        //The `old` and new `frames`
        memcpy(axes_old, old, 6*sizeof(double));
        memcpy(axes_old+6, z_old, 3*sizeof(double));

        memcpy(axes_new, new, 6*sizeof(double));
        memcpy(axes_new+6, z_new, 3*sizeof(double));

        mm_mult(3, 3, axes_new, true, 3, 3, axes_old, false, mat);
        mm_mult(m, 3, v, false, 3, 3, mat, true, valigned);
    }
    else if (n == 3) {
        mm_mult(3, 3, new, true, 3, 3, old, false, mat);
        mm_mult(m, 3, v, false, 3, 3, mat, true, valigned);
    }
    else {
        printf("`n`(= %d) must be 1, 2, or 3.\n", n);
        abort();
    }
}

/******************************************************************************/

bool mat_is_rotmat(const double* mat, const double tol){

    double tmp[9]; /* 3 x 3 matrix as 1-D array */
    double maxval;

    double determinant = det(3, mat);
    bool det_is_one = fabs(determinant-1.0) <= tol;

    // tmp = mat @ mat^T
    mm_mult(3, 3, mat, false, 3, 3, mat, true, tmp);
    //Subtract the identity matrix
    for(int i=0; i<3; ++i){
        tmp[3*i+i] -= 1.0;
    }
    //Find maximum absolute value of tmp
    maxval = 0.0;
    for(int i=0; i<9; ++i){
        maxval = fmax(maxval, fabs(tmp[i]));
    }

    bool is_orthogonal = maxval <= tol;

    return (is_orthogonal && det_is_one);
}

/******************************************************************************/
                               /* QUATERNIONS */
/******************************************************************************/

void quat_canonicalized(double q[4]){

    for(int i=0; i<4; ++i){
        if (q[i] > 0.0){
            break;
        }
        else if (q[i] < 0.0){
            q[0] = -q[0]; q[1] = -q[1]; q[2] = -q[2]; q[3] = -q[3];
            break;
        }
    }
}

/******************************************************************************/

void quat_rand(double q[4]){

    double u0 = (double) rand() / (double) RAND_MAX;
    double u1 = (double) rand() / (double) RAND_MAX;
    double u2 = (double) rand() / (double) RAND_MAX;

    double sqrt_u0 = sqrt(u0);
    double sqrt_1mu0 = sqrt(1.0 - u0);

    double cos_u1 = cos(2.0*M_PI*u1);
    double sin_u1 = sin(2.0*M_PI*u1);
    double cos_u2 = cos(2.0*M_PI*u2);
    double sin_u2 = sin(2.0*M_PI*u2);

    q[0] = sqrt_1mu0*cos_u1;
    q[1] = sqrt_1mu0*sin_u1;
    q[2] = sqrt_u0*sin_u2;
    q[3] = sqrt_u0*cos_u2;

    quat_canonicalized(q);
}

/******************************************************************************/

void quat_identity(double q[4]){

    q[0] = 1.0; q[1] = 0.0; q[2] = 0.0; q[3] = 0.0;
}

/******************************************************************************/

/* Performs in-place conjugation of a quaternion */
void quat_conjugated(double q[4]){

    q[1] = -q[1]; q[2] = -q[2]; q[3] = -q[3];
}

/******************************************************************************/

/* Performs in-place inversion of a quaternion */
int quat_inverted(double q[4]){

    double q_scaled[4];

    // Scale the quaternion
    double scale = 0.0;
    for(int i=0; i<4; ++i){
        scale = fmax(scale, fabs(q[i]));
    }
    if (scale == 0.0){
        printf("Cannot invert a zero quaternion.");
        return EXIT_FAILURE;
    }
    for(int i=0; i<4; ++i){
        q_scaled[i] = q[i]/scale;
    }

    double d = dot(4, q_scaled, q_scaled);
    quat_conjugated(q_scaled);

    for(int i=0; i<4; ++i){
        q[i] = (q_scaled[i]/d)/scale;
    }
    return EXIT_SUCCESS;
}

/******************************************************************************/

/* Performs in-place normalization of a quaternion */
int quat_normalized(double q[4], const bool canonicalize){

    int estat = normalize(4, q);
    if (estat == EXIT_FAILURE){
        return estat;
    }
    if (canonicalize){
        quat_canonicalized(q);
    }
    return EXIT_SUCCESS;
}

/******************************************************************************/

/* Returns true if a quaternion is normalized within tolerance*/
bool quat_is_normalized(const double q[4], const double tol){

    return fabs(norm(4,q) - 1.0) <= tol;
}

/******************************************************************************/

/* Calculates the product of two quaternions */
void quat_prod(const double p[4], const double q[4], double pq[4]){

    pq[0] = p[0]*q[0] - p[1]*q[1] - p[2]*q[2] - p[3]*q[3];
    pq[1] = p[1]*q[0] + p[0]*q[1] - p[3]*q[2] + p[2]*q[3];
    pq[2] = p[2]*q[0] + p[3]*q[1] + p[0]*q[2] - p[1]*q[3];
    pq[3] = p[3]*q[0] - p[2]*q[1] + p[1]*q[2] + p[0]*q[3];
}

/******************************************************************************/

double quat_angle_between(const double p[4], const double q[4]){

    double q_[4] = {q[0], q[1], q[2], q[3]};
    double tmp[4];
    double cos_angle = dot(4, p, q_);
    if (cos_angle < 0.0){
        // Negate q_ to ensure both p and q_ lie on the same hemisphere
        q_[0] = -q_[0]; q_[1] = -q_[1]; q_[2] = -q_[2]; q_[3] = -q_[3];
    }

    for(int i=0; i<4; ++i){
        tmp[i] = p[i] - q_[i];
    }
    double d_mag = norm(4, tmp);

    for(int i=0; i<4; ++i){
        tmp[i] = p[i] + q_[i];
    }
    double s_mag = norm(4, tmp);

    return 4.0*atan2(d_mag, s_mag);
}

/******************************************************************************/

void quat_interpolate(const double q1[4], const double q2[4], const double t,
        double q[4]){

    double q2_[4] = {q2[0], q2[1], q2[2], q2[3]};
    double tmp[4];
    double cos_angle = dot(4, q1, q2_);
    if (cos_angle < 0.0){
        // Negate q2_ to ensure both q1 and q2_ lie on the same hemisphere
        q2_[0] = -q2_[0]; q2_[1] = -q2_[1]; q2_[2] = -q2_[2]; q2_[3] = -q2_[3];
    }

    for(int i=0; i<4; ++i){
        tmp[i] = q1[i] - q2_[i];
    }
    double d_mag = norm(4, tmp);

    for(int i=0; i<4; ++i){
        tmp[i] = q1[i] + q2_[i];
    }
    double s_mag = norm(4, tmp);

    double theta = 2.0*atan2(d_mag, s_mag);

    // SLERP
    double c = sinc(theta);
    double a = (1.0 - t) * sinc ( (1.0-t)*theta ) / c;
    double b = t * sinc(t*theta) / c;

    for (int i=0; i<4; ++i){
        q[i] = a*q1[i] + b*q2_[i];
    }
    quat_normalized(q, false);
}

/******************************************************************************/

void quat_deriv_to_angvel_mat(const double q[4], const bool world, double* mat){

    /* mat is 3 x 4 */
    const int n = 4;

    if (world){
        mat[n*0+0] = -q[1]; mat[n*0+1] =  q[0]; mat[n*0+2] = -q[3]; mat[n*0+3] =  q[2];
        mat[n*1+0] = -q[2]; mat[n*1+1] =  q[3]; mat[n*1+2] =  q[0]; mat[n*1+3] = -q[1];
        mat[n*2+0] = -q[3]; mat[n*2+1] = -q[2]; mat[n*2+2] =  q[1]; mat[n*2+3] =  q[0];
    }
    else
    {
        mat[n*0+0] = -q[1]; mat[n*0+1] =  q[0]; mat[n*0+2] =  q[3]; mat[n*0+3] = -q[2];
        mat[n*1+0] = -q[2]; mat[n*1+1] = -q[3]; mat[n*1+2] =  q[0]; mat[n*1+3] =  q[1];
        mat[n*2+0] = -q[3]; mat[n*2+1] =  q[2]; mat[n*2+2] = -q[1]; mat[n*2+3] =  q[0];
    }

    for (int i=0; i<3; ++i){
        for (int j=0; j<4; ++j){
            mat[4*i+j] *= 2.0;
        }
    }
}

/******************************************************************************/

void quat_deriv_to_angvel(const double q[4], const double qdot[4], 
        const bool world, double angvel[3]){

    double mat[12]; /* 3 x 4 matrix as 1-D array */

    quat_deriv_to_angvel_mat(q, world, mat);
    mv_mult(3, 4, mat, false, qdot, angvel);
}

/******************************************************************************/

void quat_deriv_from_angvel_mat(const double q[4], const bool world, 
        double* mat){

    /* mat is 4 x 3 */
    const int n = 3;

    if (world){
        mat[n*0+0] = -q[1]; mat[n*0+1] = -q[2]; mat[n*0+2] = -q[3];
        mat[n*1+0] =  q[0]; mat[n*1+1] =  q[3]; mat[n*1+2] = -q[2];
        mat[n*2+0] = -q[3]; mat[n*2+1] =  q[0]; mat[n*2+2] =  q[1];
        mat[n*3+0] =  q[2]; mat[n*3+1] = -q[1]; mat[n*3+2] =  q[0];
    }
    else
    {
        mat[n*0+0] = -q[1]; mat[n*0+1] = -q[2]; mat[n*0+2] = -q[3];
        mat[n*1+0] =  q[0]; mat[n*1+1] = -q[3]; mat[n*1+2] =  q[2];
        mat[n*2+0] =  q[3]; mat[n*2+1] =  q[0]; mat[n*2+2] = -q[1];
        mat[n*3+0] = -q[2]; mat[n*3+1] =  q[1]; mat[n*3+2] =  q[0];
    }

    for (int i=0; i<4; ++i){
        for (int j=0; j<3; ++j){
            mat[3*i+j] *= 0.5;
        }
    }
}

/******************************************************************************/

void quat_deriv_from_angvel(const double q[4], const double angvel[3],
        const bool world, double qdot[4]){

    double mat[12]; /* 4 x 3 matrix as 1-D array */

    quat_deriv_from_angvel_mat(q, world, mat);
    mv_mult(4, 3, mat, false, angvel, qdot);
}

/******************************************************************************/

void quat_rotmat(const double q[4], double* rotmat){

    const int n = 3;

    double q0sq = q[0]*q[0];
    double q1sq = q[1]*q[1];
    double q2sq = q[2]*q[2];
    double q3sq = q[3]*q[3];
    double q0q1 = q[0]*q[1];
    double q0q2 = q[0]*q[2];
    double q0q3 = q[0]*q[3];
    double q1q2 = q[1]*q[2];
    double q1q3 = q[1]*q[3];
    double q2q3 = q[2]*q[3];

    rotmat[n*0+0] = 2*(q0sq + q1sq) - 1.0;
    rotmat[n*0+1] = 2*(q1q2 - q0q3);
    rotmat[n*0+2] = 2*(q1q3 + q0q2);
    rotmat[n*1+0] = 2*(q1q2 + q0q3);
    rotmat[n*1+1] = 2*(q0sq + q2sq) - 1.0;
    rotmat[n*1+2] = 2*(q2q3 - q0q1);
    rotmat[n*2+0] = 2*(q1q3 - q0q2);
    rotmat[n*2+1] = 2*(q2q3 + q0q1);
    rotmat[n*2+2] = 2*(q0sq + q3sq) - 1.0;
}

/******************************************************************************/

void quat_shiftmat(const double q[4], const bool forward, double* shiftmat){

    double conj_q[4];

    if (forward){
        memcpy(conj_q, q, 4*sizeof(double));
        quat_conjugated(conj_q);
        quat_rotmat(conj_q, shiftmat);
    }
    else {
        quat_rotmat(q, shiftmat);
    }
}

/******************************************************************************/

void quat_rotate_vectors(const int n, const double* v, const double q[4],
        double* vrot){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    quat_rotmat(q, rotmat);
    rotate_vectors(n, v, rotmat, vrot);
}

/******************************************************************************/

void quat_shift_vectors(const int n, const double* v, const double q[4],
        const bool forward, double* vshift){

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    quat_shiftmat(q, forward, shiftmat);
    shift_vectors(n, v, shiftmat, vshift);
}

/******************************************************************************/

void quat_shift_tensor2(const double* A, const double q[4],
        const bool forward, double* Ashift){

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    quat_shiftmat(q, forward, shiftmat);
    shift_tensor2(A, shiftmat, Ashift);
}

/******************************************************************************/

void quat_shift_tensor3(const double* A, const double q[4],
        const bool forward, double* Ashift){

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    quat_shiftmat(q, forward, shiftmat);
    shift_tensor3(A, shiftmat, Ashift);
}

/******************************************************************************/

void quat_to_aa(const double q[4], double axis[3], double* angle){

    double p[3], q_[4];

    memcpy(q_, q, 4*sizeof(double));
    if (q_[0] < 0.0){
        q_[0] = -q_[0]; q_[1] = -q_[1]; q_[2] = -q_[2]; q_[3] = -q_[3];
    }

    double s = 0.0;
    for (int i=1; i<4; ++i){
        // Only the vector part
        s = fmax(s, fabs(q_[i]));
    }
    if (s == 0.0){
        // Identity rotation (angle = 0)
        axis[0] = 1.0; axis[1] = 0.0; axis[2] = 0.0;
        *angle = 0.0;
    }
    else {
        // Usual case
        p[0] = q_[1]/s; p[1] = q_[2]/s; p[2] = q_[3]/s; 
        double r = norm(3, p);
        axis[0] = p[0]/r; axis[1] = p[1]/r; axis[2] = p[2]/r;
        *angle = 2.0*atan2(s*r, q_[0]);
    }
}

/******************************************************************************/

void quat_to_dcm(const double q[4], double* dcm){

    quat_shiftmat(q, true, dcm);
}

/******************************************************************************/

void quat_to_euler(const double q[4], enum EulangSeq seq, const bool world,
        const double tol, double euler[3]){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    quat_rotmat(q, rotmat);
    factor_rotmat(rotmat, seq, world, tol, euler);
}

/******************************************************************************/
                            /* AXIS  ANGLE */
/******************************************************************************/

/* Returns angle in [0, pi) */
int aa_fix(double axis[3], double* angle, const bool normalize_axis){

    int estat; // Returned exit status

    // Check if axis is a zero vector
    if ((axis[0]==0.0) && (axis[1]==0.0) && (axis[2]==0.0)){
        return EXIT_FAILURE;
    }

    if (normalize_axis) {
        estat = normalize(3, axis);
        if (estat == EXIT_FAILURE) {return estat;}
    }
    // Get principal signed angle in [-pi, pi]
    *angle = atan2(sin(*angle), cos(*angle));

    // Make the angle non-negative
    if (*angle < 0.0){
        *angle = - *angle;
        axis[0] = -axis[0]; axis[1] = -axis[1]; axis[2] = -axis[2]; 
    }
    return EXIT_SUCCESS;
}

/******************************************************************************/

/* Generates a random pair of axis-angle. */
void aa_rand(double axis[3], double* angle){

    double q[4];

    quat_rand(q);
    quat_to_aa(q, axis, angle);
}

/******************************************************************************/

void aa_rotmat(const double axis[3], const double angle, double* rotmat) {

    const int n = 3;
    double sine, sine_2, cosine, icos;

    sine = sin(angle);
    cosine = cos(angle);
    sine_2 = sin(angle/2.0);
    icos = 2.0*sine_2*sine_2; // equals 1.0 - cosine

    rotmat[n*0+0] = axis[0]*axis[0]*icos + cosine;
    rotmat[n*0+1] = axis[0]*axis[1]*icos - axis[2]*sine;
    rotmat[n*0+2] = axis[0]*axis[2]*icos + axis[1]*sine;
    rotmat[n*1+0] = axis[0]*axis[1]*icos + axis[2]*sine;
    rotmat[n*1+1] = axis[1]*axis[1]*icos + cosine;
    rotmat[n*1+2] = axis[1]*axis[2]*icos - axis[0]*sine;
    rotmat[n*2+0] = axis[2]*axis[0]*icos - axis[1]*sine;
    rotmat[n*2+1] = axis[1]*axis[2]*icos + axis[0]*sine;
    rotmat[n*2+2] = axis[2]*axis[2]*icos + cosine;
}

/******************************************************************************/

void aa_from_rotmat(const double* rotmat, double axis[3], double* angle){

    const int n = 3;
    double s;

    // Cos of the angle from the trace
    double tr = trace(3, rotmat);
    double cos_angle = (tr - 1.0)/2.0;
    //Clamp cos_angle to [-1, 1]
    cos_angle = fmax(-1.0, fmin(1.0, cos_angle));

    // Sine of the angle from the anti-symmetric vector w
    double w[] = { 0.5*(rotmat[n*2+1] - rotmat[n*1+2]),
                   0.5*(rotmat[n*0+2] - rotmat[n*2+0]),
                   0.5*(rotmat[n*1+0] - rotmat[n*0+1]) };
    double sin_angle = norm(3, w);

    // Calculate the angle
    *angle = atan2(sin_angle, cos_angle);

    if (cos_angle >= 0.0){
        // Angle lies in [0, pi/2]
        if (sin_angle > 0.0){
            axis[0] = w[0]/sin_angle;
            axis[1] = w[1]/sin_angle;
            axis[2] = w[2]/sin_angle;
        }
        else {
            // angle = 0, axis non-unique
            *angle = 0.0;
            axis[0] = 1.0; axis[1] = 0.0; axis[2] = 0.0;
        }
    }
    else {
        // Find the largest entry on the diagonal of rotmat
        int k = rotmat[n*0+0] > rotmat[n*1+1] ? 0 : 1;
        k = rotmat[n*k+k] > rotmat[n*2+2] ? k : 2;
        switch (k){
            case 0:
                axis[0] = sqrt((rotmat[n*0+0]-cos_angle)/(1.0-cos_angle));
                s = 1.0/(2.0*(1.0-cos_angle)*axis[0]);
                axis[1] = s*(rotmat[n*1+0]+rotmat[n*0+1]);
                axis[2] = s*(rotmat[n*2+0]+rotmat[n*0+2]);
                break;
            case 1:
                axis[1] = sqrt((rotmat[n*1+1]-cos_angle)/(1.0-cos_angle));
                s = 1.0/(2.0*(1.0-cos_angle)*axis[1]);
                axis[0] = s*(rotmat[n*0+1]+rotmat[n*1+0]);
                axis[2] = s*(rotmat[n*2+1]+rotmat[n*1+2]);
                break;
            case 2:
                axis[2] = sqrt((rotmat[n*2+2]-cos_angle)/(1.0-cos_angle));
                s = 1.0/(2*(1.0-cos_angle)*axis[2]);
                axis[0] = s*(rotmat[n*0+2]+rotmat[n*2+0]);
                axis[1] = s*(rotmat[n*1+2]+rotmat[n*2+1]);
                break;
        }
    }
    // Resolve axis aign
    if ( dot(3, axis, w) < 0.0 ){
        axis[0] = -axis[0]; axis[1] = -axis[1]; axis[2] = -axis[2];
    }

    aa_fix(axis, angle, true);
}

/******************************************************************************/

void aa_shiftmat(const double axis[3], const double angle,
        const bool forward, double* shiftmat){

    double axis_[3];
    double rotmat[9]; /* 3 x 3 matrix as 1-D array */

    for (int i=0; i<3; ++i){
        axis_[i] = -axis[i];
    }

    aa_rotmat(axis_, angle, rotmat);

    if (!forward) {
        transpose(3, 3, rotmat, shiftmat);
    }
    else {
        memcpy(shiftmat, rotmat, 9*sizeof(double));
    }
}

/******************************************************************************/

/* Rotates vectors about axis by angle. v is a pointer to an array of
 * vectors, i.e. a n x 3 array. */
void aa_rotate_vectors(const int n, const double* v, const double axis[3],
        const double angle, double* vrot){

    double rotmat[9]; /* 3 x 3 matrix as 1-D array */

    aa_rotmat(axis, angle, rotmat);
    rotate_vectors(n, v, rotmat, vrot);
}

/******************************************************************************/

void aa_shift_vectors(const int n, const double* v, const double axis[3], 
        const double angle, const bool forward, double* vshift){

    double shiftmat[9]; /* 3 x 3 matrix as 1-D array */

    aa_shiftmat(axis, angle, forward, shiftmat);
    shift_vectors (n, v, shiftmat, vshift);
}

/******************************************************************************/

void aa_shift_tensor2(const double* A, const double axis[3],
        const double angle, const bool forward, double* Ashift){

    double shiftmat[9]; /* 3 x 3 matrix as 1-D array */

    aa_shiftmat(axis, angle, forward, shiftmat);
    shift_tensor2(A, shiftmat, Ashift);
}

/******************************************************************************/

void aa_shift_tensor3(const double* A, const double axis[3],
        const double angle, const bool forward, double* Ashift){

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    aa_shiftmat(axis, angle, forward, shiftmat);
    shift_tensor3(A, shiftmat, Ashift);
}

/******************************************************************************/

void aa_to_quat(const double axis[3], const double angle, double q[4]){

    double half_angle = 0.5*angle;
    double sn = sin(half_angle);
    
    q[0] = cos(half_angle);
    q[1] = sn*axis[0];
    q[2] = sn*axis[1];
    q[3] = sn*axis[2];

    quat_normalized(q, true);
}

/******************************************************************************/

void aa_to_dcm(const double axis[3], const double angle,
        double* dcm){

    aa_shiftmat(axis, angle, true, dcm);
}

/******************************************************************************/

void aa_to_euler(const double axis[3], const double angle, enum EulangSeq seq,
        const bool world, const double tol, double euler[3]){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    aa_rotmat(axis, angle, rotmat);
    factor_rotmat(rotmat, seq, world, tol, euler);
}

/******************************************************************************/
                         /* DIRECTION COSINE MATRIX */
/******************************************************************************/

bool mat_is_dcm(const double* mat, const double tol){

    return mat_is_rotmat(mat, tol);
}

/******************************************************************************/

void dcm_from_axes(const double* A, const double* B, double* dcm) {

    mm_mult(3, 3, B, false, 3, 3, A, true, dcm);
}

/******************************************************************************/

void dcm_rotmat(const double* dcm, double* rotmat){

    transpose(3, 3, dcm, rotmat);
}

/******************************************************************************/

void dcm_shiftmat(const double* dcm, const bool forward, double* shiftmat){

    if (!forward) {
        transpose(3, 3, dcm, shiftmat);
    }
    else {
        memcpy(shiftmat, dcm, 9*sizeof(double));
    }
}

/******************************************************************************/

void dcm_rotate_vectors(const int n, const double* v, const double* dcm,
        double* vrot){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    dcm_rotmat(dcm, rotmat);
    rotate_vectors (n, v, rotmat, vrot);
}

/******************************************************************************/

void dcm_shift_vectors(const int n, const double* v, const double* dcm,
        const bool forward, double* vshift){

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    dcm_shiftmat(dcm, forward, shiftmat);
    shift_vectors (n, v, shiftmat, vshift);
}

/******************************************************************************/

void dcm_shift_tensor2(const double* A, const double* dcm, const bool forward,
        double* Ashift){

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    dcm_shiftmat(dcm, forward, shiftmat);
    shift_tensor2(A, shiftmat, Ashift);
}

/******************************************************************************/

void dcm_shift_tensor3(const double* A, const double* dcm, const bool forward,
        double* Ashift){

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    dcm_shiftmat(dcm, forward, shiftmat);
    shift_tensor3(A, shiftmat, Ashift);
}

/******************************************************************************/

void dcm_to_quat(const double* dcm, double q[4]){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */
    double axis[3];
    double angle = 0.0;

    dcm_rotmat(dcm, rotmat);
    aa_from_rotmat(rotmat, axis, &angle);
    aa_to_quat(axis, angle, q);
}

/******************************************************************************/

void dcm_to_aa(const double* dcm, double axis[3], double* angle) {

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    dcm_rotmat(dcm, rotmat);
    aa_from_rotmat(rotmat, axis, angle);
}

/******************************************************************************/

void dcm_to_euler(const double* dcm, enum EulangSeq seq, const bool world,
        const double tol, double euler[3]){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    dcm_rotmat(dcm, rotmat);
    factor_rotmat(rotmat, seq, world, tol, euler);
}

/******************************************************************************/
                               /* EULER ANGLES */
/******************************************************************************/

void euler_rotmat(const double euler[3], enum EulangSeq seq, const bool world,
        double* rotmat){

    rotmat_euler(euler, seq, world, rotmat);
}

/******************************************************************************/

void euler_shiftmat(const double euler[3], enum EulangSeq seq, const bool world,
        const bool forward, double* shiftmat){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    rotmat_euler(euler, seq, world, rotmat);

    if (forward){
        transpose(3, 3, rotmat, shiftmat);
    }
    else {
        memcpy(shiftmat, rotmat, 9*sizeof(double));
    }
}

/******************************************************************************/

void euler_rotate_vectors(const int n, const double* v, const double euler[3],
        enum EulangSeq seq, const bool world, double* vrot) {

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    euler_rotmat(euler, seq, world, rotmat);
    rotate_vectors(n, v, rotmat, vrot);
}

/******************************************************************************/

void euler_shift_vectors(const int n, const double* v, const double euler[3],
        enum EulangSeq seq, const bool world, const bool forward,
        double* vshift) {

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    euler_shiftmat(euler, seq, world, forward, shiftmat);
    shift_vectors (n, v, shiftmat, vshift);
}

/******************************************************************************/

void euler_shift_tensor2(const double* A, const double euler[3],
        enum EulangSeq seq, const bool world, const bool forward,
        double* Ashift) {

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    euler_shiftmat(euler, seq, world, forward, shiftmat);
    shift_tensor2(A, shiftmat, Ashift);
}

/******************************************************************************/

void euler_shift_tensor3(const double* A, const double euler[3],
        enum EulangSeq seq, const bool world, const bool forward,
        double* Ashift) {

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    euler_shiftmat(euler, seq, world, forward, shiftmat);
    shift_tensor3(A, shiftmat, Ashift);
}

/******************************************************************************/

void euler_to_quat(const double euler[3], enum EulangSeq seq, const bool world,
        double q[4]){

    double axis[3];
    double angle;

    euler_to_aa(euler, seq, world, axis, &angle);
    aa_to_quat(axis, angle, q);
}

/******************************************************************************/

void euler_to_aa(const double euler[3], enum EulangSeq seq, 
        const bool world, double axis[3], double* angle){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    euler_rotmat(euler, seq, world, rotmat);
    aa_from_rotmat(rotmat, axis, angle);
}

/******************************************************************************/

void euler_to_dcm(const double euler[3], enum EulangSeq seq, const bool world,
        double* dcm){

    euler_shiftmat(euler, seq, world, true, dcm);
}

/******************************************************************************/

void euler_to_euler(const double euler[3], enum EulangSeq seq, const bool world,
        enum EulangSeq to_seq, const bool to_world, const double tol,
        double to_euler[3]){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    euler_rotmat(euler, seq, world, rotmat);
    factor_rotmat(rotmat, to_seq, to_world, tol, to_euler);
}

/******************************************************************************/
