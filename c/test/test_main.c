#include "utils_math.h"
#include "rotlib.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//******************************************************************************

void test_quat_to_aa_roundtrip(void){
    double angle;
    double qf[4], qb[4], axis[3];

    printf("Testing quaternion -> axis angle (roundtrip)\n");
    printf("--------------------------------------------\n");
    //Generate a random quaternion
    quat_rand(qf);
    printf(" Quaternion %f %f %f %f\n", qf[0], qf[1], qf[2], qf[3]);
    //Convert qf to axis angle
    quat_to_aa(qf, axis, &angle);
    printf(" Axis %f %f %f  Angle %f\n", axis[0], axis[1], axis[2], angle);
    //Convert axis angle back to quaternion
    aa_to_quat(axis, angle, qb);
    printf(" Quaternion from axis angle %f %f %f %f\n",
            qb[0], qb[1], qb[2], qb[3]);
    if ( allclose(4, qf, qb, 1e-8, 1e-14) ){
        printf("Passed\n");
    } else {
        printf("Failed\n");
    }
}

//******************************************************************************

void test_quat_to_dcm_roundtrip(void){
    double qf[4], qb[4], dcm[9];

    printf("\n");
    printf("Testing quaternion -> dcm (roundtrip)\n");
    printf("-------------------------------------\n");
    //Generate a random quaternion
    quat_rand(qf);
    //qf[0] = 1.0; qf[1] = 2.0; qf[2] = 3.0; qf[3] = 4.0;
    //unitize(4, qf);
    printf(" Quaternion %f %f %f %f\n", qf[0], qf[1], qf[2], qf[3]);
    //Convert qf to dcm
    quat_to_dcm(qf, dcm);
    printf( " Direction cosine matrix\n");
    for (int i=0; i<3; ++i){
        for (int j=0; j<3; ++j){
            printf(" %f ", dcm[3*i+j]);
        }
        printf("\n");
    }
    //Convert dcm back to quaternion
    dcm_to_quat(dcm, qb);
    printf(" Quaternion from dcm %f %f %f %f\n",
            qb[0], qb[1], qb[2], qb[3]);
    if ( allclose(4, qf, qb, 1e-8, 1e-14) ){
        printf("Passed\n");
    } else {
        printf("Failed\n");
    }
}

//******************************************************************************

void test_quat_to_euler_roundtrip(void){
    double qf[4], qb[4], euler[3];

    bool world[] = {true, false};
    enum EulangSeq seq[] = {ES_XYZ, ES_XZY, ES_YXZ, ES_YZX, ES_ZXY, ES_ZYX};

    printf("\n");
    printf("Testing quaternion -> euler (roundtrip)\n");
    printf("--------------------------------------------\n");
    //Generate a random quaternion
    quat_rand(qf);
    //qf[0] = 1.0; qf[1] = 2.0; qf[2] = 3.0; qf[3] = 4.0;
    //unitize(4, qf);
    printf(" Quaternion %f %f %f %f\n", qf[0], qf[1], qf[2], qf[3]);
    //Convert qf to euler
    for (int i=0; i<6; ++i){
        for (int j=1; j<2; ++j){
            quat_to_euler(qf, seq[i], world[j], euler);
            printf( " Euler %d %d %f %f %f\n", seq[i], world[j], 
                    euler[0], euler[1], euler[2]);
            //Convert euler back to quaternion
            euler_to_quat(euler, seq[i], world[j], qb);
            printf(" Quaternion from euler %f %f %f %f\n",
                    qf[0], qf[1], qf[2], qf[3]);

            if ( allclose(4, qf, qb, 1e-8, 1e-14) ){
                printf("Passed\n");
            } else {
                printf("Failed\n");
            }
        }
    }
}

//******************************************************************************

void test_aa_to_dcm_roundtrip(void){
    double axisf[3], axisb[3], dcm[9];
    double anglef, angleb;

    printf("\n");
    printf("Testing axis angle -> dcm (roundtrip)\n");
    printf("-------------------------------------\n");
    //Generate a random axis angle
    aa_rand(axisf, &anglef);
    printf(" Axis %f %f %f. Angle %f\n", axisf[0], axisf[1], axisf[2], anglef);
    //Convert axis angle to dcm
    aa_to_dcm(axisf, anglef, dcm);
    printf( " Direction cosine matrix\n");
    for (int i=0; i<3; ++i){
        for (int j=0; j<3; ++j){
            printf(" %f ", dcm[3*i+j]);
        }
        printf("\n");
    }
    //Convert dcm back to axis angle
    dcm_to_aa(dcm, axisb, &angleb);
    printf(" Axis %f %f %f. Angle %f\n", axisb[0], axisb[1], axisb[2], angleb);
    if ( allclose(3, axisf, axisb, 1e-8, 1e-14) && 
         isclose(anglef, angleb, 1e-8, 1e-14) ){
        printf("Passed\n");
    } else {
        printf("Failed\n");
    }
}

//******************************************************************************

void test_aa_to_euler_roundtrip(void){
    double axisf[3], axisb[3], euler[3];
    double anglef, angleb;

    bool world[] = {true, false};
    enum EulangSeq seq[] = {ES_XYZ, ES_XZY, ES_YXZ, ES_YZX, ES_ZXY, ES_ZYX};

    printf("\n");
    printf("Testing axis angle -> euler (roundtrip)\n");
    printf("--------------------------------------------\n");
    //Generate a random axis angle
    aa_rand(axisf, &anglef);
    printf(" Axis %f %f %f. Angle %f\n", axisf[0], axisf[1], axisf[2], anglef);
    //Convert qf to euler
    for (int i=0; i<6; ++i){
        for (int j=1; j<2; ++j){
            aa_to_euler(axisf, anglef, seq[i], world[j], euler);
            printf( " Euler %d %d %f %f %f\n", seq[i], world[j], 
                    euler[0], euler[1], euler[2]);
            //Convert euler back to axis angle
            euler_to_aa(euler, seq[i], world[j], axisb, &angleb);
            printf(" Axis %f %f %f. Angle %f\n", 
                    axisb[0], axisb[1], axisb[2], angleb);

            if ( allclose(3, axisf, axisb, 1e-8, 1e-14) && 
                 isclose(anglef, angleb, 1e-8, 1e-14) ){
                printf("Passed\n");
            } else {
                printf("Failed\n");
            }
            printf("\n");
        }
    }
}

//******************************************************************************

void test_dcm_to_euler_roundtrip(void){
    double dcmf[9], dcmb[9], euler[3], axis[3];
    double angle;

    bool world[] = {true, false};
    enum EulangSeq seq[] = {ES_XYZ, ES_XZY, ES_YXZ, ES_YZX, ES_ZXY, ES_ZYX};

    printf("\n");
    printf("Testing dcm -> euler (roundtrip)\n");
    printf("--------------------------------------------\n");
    //Generate a random axis angle and convert to dcm
    aa_rand(axis, &angle);
    aa_to_dcm(axis, angle, dcmf);
    printf( " Direction cosine matrix\n");
    for (int i=0; i<3; ++i){
        for (int j=0; j<3; ++j){
            printf(" %f ", dcmf[3*i+j]);
        }
        printf("\n");
    }
    //Convert dcm to euler
    for (int i=0; i<6; ++i){
        for (int j=1; j<2; ++j){
            dcm_to_euler(dcmf, seq[i], world[j], euler);
            printf( " Euler %d %d %f %f %f\n", seq[i], world[j], 
                    euler[0], euler[1], euler[2]);
            //Convert euler back to dcm
            euler_to_dcm(euler, seq[i], world[j], dcmb);
            printf( " Direction cosine matrix\n");
            for (int i=0; i<3; ++i){
                for (int j=0; j<3; ++j){
                    printf(" %f ", dcmb[3*i+j]);
                }
                printf("\n");
            }
            if ( allclose(9, dcmf, dcmb, 1e-8, 1e-14) ){
                printf("Passed\n");
            } else {
                printf("Failed\n");
            }
            printf("\n");
        }
    }
}

//******************************************************************************

void test_euler_to_euler_roundtrip(void){
    double eulerf[3], eulerb[3], euler[3], axis[3];
    double angle;

    bool world[] = {true, false};
    enum EulangSeq seq[] = {ES_XYZ, ES_XZY, ES_YXZ, ES_YZX, ES_ZXY, ES_ZYX};

    printf("\n");
    printf("Testing euler -> euler (roundtrip)\n");
    printf("--------------------------------------------\n");
    //Generate a random axis angle and convert to euler
    aa_rand(axis, &angle);
    
    printf(" Axis %f %f %f  Angle %f\n", axis[0], axis[1], axis[2], angle);
    for (int i=0; i<5; ++i){
        for (int j=0; j<1; ++j){
            printf("\n");
            aa_to_euler(axis, angle, seq[i], world[j], eulerf);
            printf( " Euler %d %d %f %f %f\n", seq[i], world[j], 
                            eulerf[0], eulerf[1], eulerf[2]);
            //Convert one euler sequence to another
            for (int ii=i; ii <6; ++ii){
                for (int jj=j; jj <2; ++jj){
                    euler_to_euler(eulerf, seq[i], world[j], 
                        seq[ii], world[jj], euler);
                    printf( " Euler %d %d %f %f %f\n", seq[ii], world[jj], 
                            euler[0], euler[1], euler[2]);
                    //Convert the euler sequence back
                    euler_to_euler(euler, seq[ii], world[jj], 
                        seq[i], world[j], eulerb);
                    printf( " Euler %d %d %f %f %f\n", seq[i], world[j], 
                            eulerb[0], eulerb[1], eulerb[2]);

                    if ( allclose(3, eulerf, eulerb, 1e-8, 1e-14) ){
                        printf("Passed\n");
                    } else {
                        printf("Failed\n");
                    }
                    printf("\n");
                }
            }
        }
    }
}

//******************************************************************************

void test_angvel_to_qdot_roundtrip(void){
    /*
     * Inconvenient to do qdot -> angvel round trip as then for a given q, qdot
     * has to be chosen such that it is orthogonal to q. qdot does not have to
     * be a unit vector.
    */

    double q[4], qdot[4], omegaf[3], omegab[3];

    printf("Testing angvel -> qdot (roundtrip)\n");
    printf("----------------------------------\n");
    //Generate a random quaternion
    quat_rand(q);

    //Generate a random angular velocity between -1 and 1
    omegaf[0] = 2.0*rand()/RAND_MAX - 1.0;
    omegaf[1] = 2.0*rand()/RAND_MAX - 1.0;
    omegaf[2] = 2.0*rand()/RAND_MAX - 1.0;
    printf(" Omega %f %f %f\n", omegaf[0], omegaf[1], omegaf[2]);

    //Convert omegaf to qdot
    quat_deriv_from_angvel(q, omegaf, qdot);
    printf(" qdot %f %f %f %f\n", qdot[0], qdot[1], qdot[2], qdot[3]);

    //Convert qdot back to angular velocity
    quat_deriv_to_angvel(q, qdot, omegab);
    printf(" Omega <- qdot %f %f %f\n", omegab[0], omegab[1], omegab[2]);

    if ( allclose(3, omegaf, omegab, 1e-8, 1e-14) ){
        printf("Passed\n");
    } else {
        printf("Failed\n");
    }
}

/******************************************************************************/

int main(void){
    //Seed the random number generator
    srand( (unsigned int) time((time_t *) NULL) );

    test_quat_to_aa_roundtrip();
    test_quat_to_dcm_roundtrip();
    test_quat_to_euler_roundtrip();
    test_aa_to_dcm_roundtrip();
    test_aa_to_euler_roundtrip();
    test_dcm_to_euler_roundtrip();
    test_euler_to_euler_roundtrip();
    test_angvel_to_qdot_roundtrip();
    return 0;
}
