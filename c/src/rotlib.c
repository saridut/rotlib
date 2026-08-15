#include "rotlib.h"
#include "utils_math.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Euler angle sequence: XYZ (world). First rotation about X, second rotation   */
/* about Y, and the third rotation about Z axis of the world(i.e. fixed) frame. */
/* This is the same as the sequence used in Blender.                            */
/* In contrast, the XYZ sequence is understood in the Aerospace community as:   */
/* First rotation about Z-axis, second rotation about Y-axis, and the third     */
/* rotation about X-axis of the body frame.                                     */


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

static void factor_rotmat_XYZ(const double* rotmat, double* euler){

    const int n = 3;
    double phi, theta, psi;

    if (rotmat[n*2+0] < 1.0){
        if (rotmat[n*2+0] > -1.0){
            theta = asin(-rotmat[n*2+0]);
            psi = atan2(rotmat[n*1+0], rotmat[n*0+0]);
            phi = atan2(rotmat[n*2+1], rotmat[n*2+2]);
        }
        else {
            /* Not unique: phi - psi = atan2(-rotmat[1][2], rotmat[1][1]) */
            theta = M_PI_2;
            psi = -atan2(-rotmat[n*1+2], rotmat[n*1+1]);
            phi = 0.0;
        }
    }
    else {
        /* Not unique: phi + psi = atan2(-rotmat[1][2], rotmat[1][1]) */
        phi = 0.0;
        theta = -M_PI_2;
        psi = atan2(-rotmat[n*1+2], rotmat[n*1+1]);
    }

    euler[0] = phi; euler[1] = theta; euler[2] = psi;
}

/******************************************************************************/

static void factor_rotmat_XZY(const double* rotmat, double* euler){

    const int n = 3;
    double phi, theta, psi;

    if (rotmat[n*1+0] < 1.0){
        if (rotmat[n*1+0] > -1.0){
            phi = atan2(-rotmat[n*1+2], rotmat[n*1+1]);
            theta = atan2(-rotmat[n*2+0], rotmat[n*0+0]);
            psi = asin(rotmat[n*1+0]);
        }
        else {
            /* Not unique: phi - theta = atan2(rotmat[2][1], rotmat[2][2]) */
            phi = 0.0;
            theta = -atan2(rotmat[n*2+1], rotmat[n*2+2]);
            psi = -M_PI_2;
        }
    }
    else {
        /* Not unique: phi + theta = atan2(rotmat[2][1], rotmat[2][2]) */
        phi = 0.0;
        theta = atan2(rotmat[n*2+1], rotmat[n*2+2]);
        psi = M_PI_2;
    }

    euler[0] = phi; euler[1] = theta; euler[2] = psi;
}

/******************************************************************************/

static void factor_rotmat_YXZ(const double* rotmat, double* euler){

    const int n = 3;
    double phi, theta, psi;

    if (rotmat[n*2+1] < 1.0) {
        if (rotmat[n*2+1] > -1.0) {
            phi = asin(rotmat[n*2+1]);
            theta = atan2(-rotmat[n*2+0], rotmat[n*2+2]);
            psi = atan2(-rotmat[n*0+1], rotmat[n*1+1]);
        }
        else {
            /* Not unique: theta - psi = atan2(rotmat[0][2], rotmat[0][0]) */
            phi = -M_PI_2;
            theta = 0.0;
            psi = -atan2(rotmat[n*0+2], rotmat[n*0+0]);
        }
    }
    else {
        /* Not unique: theta + psi = atan2(rotmat[0][2], rotmat[0][0]) */
        phi = M_PI_2;
        theta = 0.0;
        psi = atan2(rotmat[n*0+2], rotmat[n*0+0]);
    }

    euler[0] = phi; euler[1] = theta; euler[2] = psi;
}

/******************************************************************************/

static void factor_rotmat_YZX(const double* rotmat, double* euler){

    const int n = 3;
    double phi, theta, psi;

    if (rotmat[n*0+1] < 1.0) {
        if (rotmat[n*0+1] > -1.0) {
            phi = atan2(rotmat[n*2+1], rotmat[n*1+1]);
            theta = atan2(rotmat[n*0+2], rotmat[n*0+0]);
            psi = asin(-rotmat[n*0+1]);
        }
        else {
            /* Not unique: theta - phi = atan2(-rotmat[2][0], rotmat[2][2]) */
            phi = -atan2(-rotmat[n*2+0], rotmat[n*2+2]);
            theta = 0.0;
            psi = M_PI_2;
        }
    }
    else {
        /* Not unique: theta + phi = atan2(-rotmat[2][0], rotmat[2][2]) */
        phi = atan2(-rotmat[n*2+0], rotmat[n*2+2]);
        theta = 0.0;
        psi = -M_PI_2;
    }

    euler[0] = phi; euler[1] = theta; euler[2] = psi;
}

/******************************************************************************/

static void factor_rotmat_ZXY(const double* rotmat, double* euler){

    const int n = 3;
    double phi, theta, psi;

    if (rotmat[n*1+2] < 1.0) {
        if (rotmat[n*1+2] > -1.0) {
            phi = asin(-rotmat[n*1+2]);
            theta = atan2(rotmat[n*0+2], rotmat[n*2+2]);
            psi = atan2(rotmat[n*1+0], rotmat[n*1+1]);
        }
        else {
            /* Not unique: psi - theta = atan2(-rotmat[0][1], rotmat[0][0]) */
            phi = M_PI_2;
            theta = -atan2(-rotmat[n*0+1], rotmat[n*0+0]);
            psi = 0.0;
        }
    }
    else {
        /* Not unique: psi + theta = atan2(-rotmat[0][1], rotmat[0][0]) */
        phi = -M_PI_2;
        theta = atan2(-rotmat[n*0+1], rotmat[n*0+0]);
        psi = 0.0;
    }

    euler[0] = phi; euler[1] = theta; euler[2] = psi;
}

/******************************************************************************/

static void factor_rotmat_ZYX(const double* rotmat, double* euler){

    const int n = 3;
    double phi, theta, psi;

    if (rotmat[n*0+2] < 1.0) {
        if (rotmat[n*0+2] > -1.0) {
            phi = atan2(-rotmat[n*1+2], rotmat[n*2+2]);
            theta = asin(rotmat[n*0+2]);
            psi = atan2(-rotmat[n*0+1], rotmat[n*0+0]);
        }
        else {
            /* Not unique: psi - phi = atan2(rotmat[1][0], rotmat[1][1]) */
            phi = -atan2(rotmat[n*1+0], rotmat[n*1+1]);
            theta = -M_PI_2;
            psi = 0.0;
        }
    }
    else {
        /* Not unique: psi + phi = atan2(rotmat[1][0], rotmat[1][1]) */
        phi = atan2(rotmat[n*1+0], rotmat[n*1+1]);
        theta = M_PI_2;
        psi = 0.0;
    }

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
        const bool world, double* euler){

    double rotmat_[9]; /* 3 x 3 matrix as an 1-D array */

    if (!world) {
        transpose(3, 3, rotmat, rotmat_);
    }
    else {
        memcpy(rotmat_, rotmat, 9*sizeof(double));
    }

    switch (seq){
        case ES_XYZ :
            factor_rotmat_XYZ(rotmat_, euler);
            break;
        case ES_XZY :
            factor_rotmat_XZY(rotmat_, euler);
            break;
        case ES_YXZ :
            factor_rotmat_YXZ(rotmat_, euler);
            break;
        case ES_YZX :
            factor_rotmat_YZX(rotmat_, euler);
            break;
        case ES_ZXY :
            factor_rotmat_ZXY(rotmat_, euler);
            break;
        case ES_ZYX :
            factor_rotmat_ZYX(rotmat_, euler);
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

    double rotmat_trans[9]; /* 3 x 3 array as an 1-D array*/ 

    transpose(3, 3, rotmat, rotmat_trans);
    mm_mult(n, 3, v, 3, rotmat_trans, vrot);
}

/******************************************************************************/

static void shift_vectors (const int n, const double* v, const double* shiftmat,
        const bool forward, double* vshift){

    double shiftmat_trans[9]; /* 3 x 3 array as an 1-D array*/ 

    transpose(3,3, shiftmat, shiftmat_trans);
    mm_mult(n, 3, v, 3, shiftmat_trans, vshift);
}

/******************************************************************************/

static void shift_tensor2 (const double* A, const double* shiftmat,
        const bool forward, double* Ashift){

    const int n=3;

    for (int i=0; i<n; ++i){
        for (int j=0; j<n; ++j){
            Ashift[n*i+j] = 0.0;
        }
    }

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
        const bool forward, double* Ashift){

    const int n = 3;

    for (int i=0; i<n; ++i){
        for (int j=0; j<n; ++j){
            for (int k=0; k<n; ++k){
                Ashift[n*n*i+n*j+k] = 0.0;
            }
        }
    }

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

void align(const int n, const double* v, const int m, const double* old, 
        const double* new, double* valigned){

    double cos_angle, angle;
    double axis[3], z_old[3], z_new[3];
    double dcm[9], axes_old[9], axes_new[9];

    if (m == 1) {
        //Angle between old and new
        cos_angle = dot(3, old, new);

        if ( (cos_angle > 1.0) || (cos_angle < -1.0) ){
            if (isclose(fabs(cos_angle), 1.0, 1e-9, 1e-15)){
                cos_angle = copysign(1.0, cos_angle);
            } else {
                printf("Cos of angle out-of-range.\n");
                abort();
            }
        }
        if (isclose(cos_angle, 1.0, 1e-9, 1e-15)){
            //Angle = zero: No rotation needed
            memcpy(valigned, v, 3*n*sizeof(double));
        } else if (isclose(cos_angle, -1.0, 1e-9, 1e-15)){
            //Angle = pi: Rotation axis not unique. We flip the
            //vectors `v` to align.
            memcpy(valigned, v, 3*n*sizeof(double));
            for (int i=0; i<3*n; ++i){
                *valigned = -*valigned;
                ++valigned;
            }
        } else {
            angle = acos(cos_angle);
            //Axis of rotation
            cross(old, new, axis);
            unitize(3, axis);
            aa_rotate_vectors(n, v, axis, angle, valigned);
        }
    } else if (m == 2) {
        //The third vector orthogonal to `old`
        cross(old, old+3, z_old);
        unitize(3, z_old);
        //The third vector orthogonal to `new`
        cross(new, new+3, z_new);
        unitize(3, z_new);
        //The `old` and new `frames`
        memcpy(axes_old, old, 6*sizeof(double));
        memcpy(axes_old+6, z_old, 3*sizeof(double));

        memcpy(axes_new, new, 6*sizeof(double));
        memcpy(axes_new+6, z_new, 3*sizeof(double));

        dcm_from_axes(axes_old, axes_new, dcm);
        dcm_rotate_vectors(n, v, dcm, valigned);
    } else if (m == 3) {
        dcm_from_axes(old, new, dcm);
        dcm_rotate_vectors(n, v, dcm, valigned);
    } else {
        printf("`m`(= %d) must be 1, 2, or 3.\n", m);
        abort();
    }
}

/******************************************************************************/

bool mat_is_rotmat(const double* mat){

    double identity[9];  /* 3 x 3 matrix as 1-D array */
    double mt[9];        /* 3 x 3 matrix as 1-D array */
    double mmt[9];       /* 3 x 3 matrix as 1-D array */

    double determinant = det(3, mat);
    bool det_is_one = isclose(determinant, 1.0, 1e-12, 1e-12);

    eye(3, identity);

    transpose(3, 3, mat, mt);
    mm_mult(3, 3, mat, 3, mt, mmt);

    bool is_orthogonal = allclose(9, mmt, identity, 1e-12, 1e-12);

    return (is_orthogonal && det_is_one);
}

/******************************************************************************/
                               /* QUATERNIONS */
/******************************************************************************/

void quat_rand(double q[4]){

    double axis[3];
    double angle;

    aa_rand(axis, &angle);
    aa_to_quat(axis, angle, q);
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
void quat_inverted(double q[4]){

    quat_conjugated(q);
    quat_normalized(q);
}

/******************************************************************************/

/* Performs in-place normalization of a quaternion */
void quat_normalized(double q[4]){

    double nrm = norm(4, q);

    for (int i=0; i<4; ++i){
        q[i] /= nrm;
    }
}

/******************************************************************************/

/* Returns true if a quaternion is normalized */
bool quat_is_normalized(const double q[4]){

    double nrm = norm(4, q);

    if (isclose(nrm, 1.0, 1e-14, 0.0)){
        return true;
    }
    else {
        return false;
    }
}

/******************************************************************************/

/* Calculates the product of two quaternions */
void quat_prod(const double p[4], const double q[4], double pq[4],
        bool normalize){

    pq[0] = p[0]*q[0] - p[1]*q[1] -p [2]*q[2] -p [3]*q[3];
    pq[1] = p[1]*q[0] + p[0]*q[1] -p [3]*q[2] +p [2]*q[3];
    pq[2] = p[2]*q[0] + p[3]*q[1] +p [0]*q[2] -p [1]*q[3];
    pq[3] = p[3]*q[0] - p[2]*q[1] +p [1]*q[2] +p [0]*q[3];

    if (normalize) quat_normalized(pq);
}

/******************************************************************************/

double quat_angle_between(const double p[4], const double q[4]){

    return acos(dot(4, p, q));
}

/******************************************************************************/

void quat_interpolate(const double q1[4], const double q2[4], const double t,
        double q[4]){

    double theta = quat_angle_between(q1, q2);
    double sin_theta = sin(theta);
    double sin_t_theta = sin(t*theta);
    double sin_it_theta = sin((1.0-t)*theta);

    for (int i=0; i<4; ++i){
        q[i] = ( sin_it_theta*q1[i] + sin_t_theta*q2[i] )/sin_theta;
    }
}

/******************************************************************************/

void quat_deriv_to_angvel_mat(const double q[4], double* mat){

    /* mat is 3 x 4 */
    const int n = 4;

    mat[n*0+0] = -q[1]; mat[n*0+1] =  q[0]; mat[n*0+2] = -q[3]; mat[n*0+3] =  q[2];
    mat[n*1+0] = -q[2]; mat[n*1+1] =  q[3]; mat[n*1+2] =  q[0]; mat[n*1+3] = -q[1];
    mat[n*2+0] = -q[3]; mat[n*2+1] = -q[2]; mat[n*2+2] =  q[1]; mat[n*2+3] =  q[0];

    for (int i=0; i<3; ++i){
        for (int j=0; j<4; ++j){
            mat[4*i+j] *= 2.0;
        }
    }
}

/******************************************************************************/

void quat_deriv_to_angvel(const double q[4], const double qdot[4], 
        double angvel[3]){

    double mat[12]; /* 3 x 4 matrix as 1-D array */

    quat_deriv_to_angvel_mat(q, mat);
    mv_mult(3, 4, mat, qdot, angvel);
}

/******************************************************************************/

void quat_deriv_from_angvel_mat(const double q[4], double* mat){

    /* mat is 4 x 3 */
    const int n = 3;

    mat[n*0+0] = -q[1]; mat[n*0+1] = -q[2]; mat[n*0+2] = -q[3];
    mat[n*1+0] =  q[0]; mat[n*1+1] =  q[3]; mat[n*1+2] = -q[2];
    mat[n*2+0] = -q[3]; mat[n*2+1] =  q[0]; mat[n*2+2] =  q[1];
    mat[n*3+0] =  q[2]; mat[n*3+1] = -q[1]; mat[n*3+2] =  q[0];

    for (int i=0; i<4; ++i){
        for (int j=0; j<3; ++j){
            mat[3*i+j] *= 0.5;
        }
    }
}

/******************************************************************************/

void quat_deriv_from_angvel(const double q[4], const double angvel[3],
        double qdot[4]){

    double mat[12]; /* 4 x 3 matrix as 1-D array */

    quat_deriv_from_angvel_mat(q, mat);
    mv_mult(4, 3, mat, angvel, qdot);
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

void quat_shiftmat(const double q[4], const bool forward,
        double* shiftmat){

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
    shift_vectors(n, v, shiftmat, forward, vshift);
}

/******************************************************************************/

void quat_shift_tensor2(const double* A, const double q[4],
        const bool forward, double* Ashift){

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    quat_shiftmat(q, forward, shiftmat);
    shift_tensor2(A, shiftmat, forward, Ashift);
}

/******************************************************************************/

void quat_shift_tensor3(const double* A, const double q[4],
        const bool forward, double* Ashift){

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    quat_shiftmat(q, forward, shiftmat);
    shift_tensor3(A, shiftmat, forward, Ashift);
}

/******************************************************************************/

void quat_to_aa(const double q[4], double axis[3], double* angle){

    double sine = sqrt(1.0-q[0]*q[0]);

    *angle = 2*acos(q[0]);

    if (*angle > 0.0) {
        if (*angle < M_PI){
            axis[0] = q[1]/sine;
            axis[1] = q[2]/sine;
            axis[2] = q[3]/sine;
        }
        else {
            double rotmat[9];  /* 3 x 3 matrix as 1-D array */
            quat_rotmat(q, rotmat);
            aa_from_rotmat(rotmat, axis, angle);
        }
    }
    else {
        axis[0] = 1.0; axis[1] = 0.0; axis[2] = 0.0;
    }
    aa_fix(axis, angle, true);
}

/******************************************************************************/

void quat_to_dcm(const double q[4], double* dcm){

    quat_shiftmat(q, true, dcm);
}

/******************************************************************************/

void quat_to_euler(const double q[4], enum EulangSeq seq, const bool world,
        double euler[3]){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    quat_rotmat(q, rotmat);
    factor_rotmat(rotmat, seq, world, euler);
}

/******************************************************************************/
                            /* AXIS  ANGLE */
/******************************************************************************/

/* Returns angle in [0, pi) */
void aa_fix(double axis[3], double* angle, const bool normalize){

    if (normalize) {
        double nrm = norm(3, axis);
        if (!isclose(nrm, 1.0, 1e-14, 1e-14)) {
            for (int i=0; i<3; ++i) {
                axis[i] /= nrm;
            }
        }
    }

    *angle = fmod(*angle, 2*M_PI);

    if (*angle < 0.0) {
        *angle = -*angle;
        for (int i=0; i<3; ++i) {
            axis[i] = -axis[i];
        }
    }
    else if (*angle > M_PI) {
        *angle = 2*M_PI - *angle;
        for (int i=0; i<3; ++i) {
            axis[i] = -axis[i];
        }
    }
}

/******************************************************************************/

/* Generates a random pair of axis-angle. The axis is a random vector from */
/* the surface of a unit sphere. Algorithm from Allen & Tildesley p. 349.  */
void aa_rand(double axis[3], double* angle){

    double zeta1, zeta2, zetasq, rt;

    /* Generate angle: A uniform random number from [0.0, 2*pi) */
    *angle = 2.0*M_PI*rand()/RAND_MAX;

    while (true) {
        /* Generate two uniform random numbers from [-1, 1) */
        zeta1 = 2.0*rand()/RAND_MAX - 1.0;
        zeta2 = 2.0*rand()/RAND_MAX - 1.0;

        zeta1 = 2.0*zeta1 - 1.0;
        zeta2 = 2.0*zeta2 - 1.0;
        zetasq = zeta1*zeta1 + zeta2*zeta2;
        if (zetasq <= 1.0) break;
    }

    rt = sqrt(1.0-zetasq);

    axis[0] = 2.0*zeta1*rt ;
    axis[1] = 2.0*zeta2*rt ;
    axis[2] = 1.0 - 2.0*zetasq ;

    aa_fix(axis, angle, true);
}

/******************************************************************************/

void aa_rotmat(const double axis[3], const double angle,
        double* rotmat) {

    const int n = 3;
    double sine, cosine, icos;

    sine = sin(angle);
    cosine = cos(angle);
    icos = 1.0 - cosine;

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
    double tr, s;
    double tmp;
    int k;

    tr = trace(3, rotmat);
    *angle = acos(0.5*(tr-1.0));

    if (*angle > 0) {
        if (*angle < M_PI) {
            axis[0] = rotmat[n*2+1] - rotmat[n*1+2];
            axis[1] = rotmat[n*0+2] - rotmat[n*2+0];
            axis[2] = rotmat[n*1+0] - rotmat[n*0+1];
        }
        else {
            /* Find the location of the largest entry in the diagonal of rotmat */
            k = 0; tmp = rotmat[n*0+0];
            if (rotmat[n*1+1] > rotmat[n*0+0]) {
                tmp = rotmat[n*1+1];
                k = 1;
            }
            if (rotmat[n*2+2] > tmp){
                k = 2;
            }

            switch (k) {
                case 0 :
                    axis[0] = sqrt(rotmat[n*0+0]-rotmat[n*1+1]-rotmat[n*2+2]+1)/2;
                    s = 1.0/(2*axis[0]);
                    axis[1] = s*rotmat[n*0+1];
                    axis[2] = s*rotmat[n*0+2];
                    break;
                case 1 :
                    axis[1] = sqrt(rotmat[n*1+1]-rotmat[n*0+0]-rotmat[n*2+2]+1)/2;
                    s = 1.0/(2*axis[1]);
                    axis[0] = s*rotmat[n*0+1];
                    axis[2] = s*rotmat[n*1+2];
                    break;
                case 2 :
                    axis[2] = sqrt(rotmat[n*2+2]-rotmat[n*0+0]-rotmat[n*1+1]+1)/2;
                    s = 1.0/(2*axis[2]);
                    axis[0] = s*rotmat[n*0+2];
                    axis[1] = s*rotmat[n*1+2];
                    break;
            }
        }
    }
    else {
        axis[0] = 1.0;
        axis[1] = 0.0;
        axis[2] = 0.0;
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

void aa_shift_vectors(const int n, const double* v, 
        const double axis[3], const double angle, const bool forward,
        double* vshift){

    double shiftmat[9]; /* 3 x 3 matrix as 1-D array */

    aa_shiftmat(axis, angle, forward, shiftmat);
    shift_vectors (n, v, shiftmat, forward, vshift);
}

/******************************************************************************/

void aa_shift_tensor2(const double* A, const double axis[3],
        const double angle, const bool forward, double* Ashift){

    double shiftmat[9]; /* 3 x 3 matrix as 1-D array */

    aa_shiftmat(axis, angle, forward, shiftmat);
    shift_tensor2(A, shiftmat, forward, Ashift);
}

/******************************************************************************/

void aa_shift_tensor3(const double* A, const double axis[3],
        const double angle, const bool forward, double* Ashift){

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    aa_shiftmat(axis, angle, forward, shiftmat);
    shift_tensor3(A, shiftmat, forward, Ashift);
}

/******************************************************************************/

void aa_to_dcm(const double axis[3], const double angle,
        double* dcm){

    aa_shiftmat(axis, angle, true, dcm);
}

/******************************************************************************/

void aa_to_euler(const double axis[3], const double angle, 
        enum EulangSeq seq, const bool world, double euler[3]){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    aa_rotmat(axis, angle, rotmat);
    factor_rotmat(rotmat, seq, world, euler);
}

/******************************************************************************/

void aa_to_quat(const double axis[3], const double angle, double q[4]){

    double half_angle = 0.5*angle;
    double sin_half_angle = sin(half_angle);
    double nrm;

    q[0] = cos(half_angle);
    q[1] = sin_half_angle*axis[0];
    q[2] = sin_half_angle*axis[1];
    q[3] = sin_half_angle*axis[2];

    /* Normalize the quaternion */
    nrm = norm(4, q);
    for (int i=0; i<4; ++i){
        q[i] /= nrm;
    }
}

/******************************************************************************/
                         /* DIRECTION COSINE MATRIX */
/******************************************************************************/

bool mat_is_dcm(const double* mat){

    return mat_is_rotmat(mat);
}

/******************************************************************************/

void dcm_from_axes(const double* A, const double* B, double* dcm) {

    double A_transpose[9];  /* 3 x 3 matrix as 1-D array */

    transpose(3, 3, A, A_transpose);
    mm_mult(3, 3, B, 3, A_transpose, dcm);
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
    shift_vectors (n, v, shiftmat, forward, vshift);
}

/******************************************************************************/

void dcm_shift_tensor2(const double* A, const double* dcm, const bool forward,
        double* Ashift){

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    dcm_shiftmat(dcm, forward, shiftmat);
    shift_tensor2(A, shiftmat, forward, Ashift);
}

/******************************************************************************/

void dcm_shift_tensor3(const double* A, const double* dcm, const bool forward,
        double* Ashift){

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    dcm_shiftmat(dcm, forward, shiftmat);
    shift_tensor3(A, shiftmat, forward, Ashift);
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
        double euler[3]){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    dcm_rotmat(dcm, rotmat);
    factor_rotmat(rotmat, seq, world, euler);
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
    shift_vectors (n, v, shiftmat, forward, vshift);
}

/******************************************************************************/

void euler_shift_tensor2(const double* A, const double euler[3],
        enum EulangSeq seq, const bool world, const bool forward,
        double* Ashift) {

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    euler_shiftmat(euler, seq, world, forward, shiftmat);
    shift_tensor2(A, shiftmat, forward, Ashift);
}

/******************************************************************************/

void euler_shift_tensor3(const double* A, const double euler[3],
        enum EulangSeq seq, const bool world, const bool forward,
        double* Ashift) {

    double shiftmat[9];  /* 3 x 3 matrix as 1-D array */

    euler_shiftmat(euler, seq, world, forward, shiftmat);
    shift_tensor3(A, shiftmat, forward, Ashift);
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
        enum EulangSeq to_seq, const bool to_world, double to_euler[3]){

    double rotmat[9];  /* 3 x 3 matrix as 1-D array */

    euler_rotmat(euler, seq, world, rotmat);
    factor_rotmat(rotmat, to_seq, to_world, to_euler);
}

/******************************************************************************/
