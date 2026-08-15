#!/usr/bin/env python

from pathlib import Path
import sys
import numpy as np

srcdir = Path('rotlib').resolve()
print(srcdir)
sys.path.insert(0, str(srcdir))

import rotlib

rng = np.random.default_rng()

def test_quat_to_aa_roundtrip():
    print("Testing quaternion -> axis angle (roundtrip)")
    print("--------------------------------------------")
    #Generate a random quaternion
    qf = rotlib.quat_rand(rng)
    print(" Quaternion " + ' '.join([f"{x:g}" for x in qf]) )
    #Convert qf to axis angle
    axis, angle = rotlib.quat_to_aa(qf)
    print(f" Axis " + ' '.join([f"{x:g}" for x in axis]) 
          + f"  Angle {angle:g}")
    #Convert axis angle back to quaternion
    qb = rotlib.aa_to_quat(axis, angle)
    print(f" Quaternion <- axis angle " 
          + ' '.join([f"{x:g}" for x in qb]))
    assert np.allclose(qf, qb, 1e-8, 1e-14)
    

def test_quat_to_dcm_roundtrip():
    print("Testing quaternion -> dcm (roundtrip)")
    print("-------------------------------------")
    #Generate a random quaternion
    qf = rotlib.quat_rand(rng)
    print(" Quaternion " + ' '.join([f"{x:g}" for x in qf]) )
    #Convert qf to dcm
    dcm = rotlib.quat_to_dcm(qf)
    print(" Direction cosine matrix")
    buf = ''
    for i in range(3):
        buf += f"  {dcm[i,0]:g} {dcm[i,1]:g} {dcm[i,2]:g}\n"
    print(buf[:-1])
    #Convert dcm back to quaternion
    qb = rotlib.dcm_to_quat(dcm)
    print(" Quaternion <- dcm " + ' '.join([f"{x:g}" for x in qb]) )
    assert np.allclose(qf, qb, 1e-8, 1e-14)


def test_quat_to_euler_roundtrip():
    seq = ['XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'] 
    world = [True, False]

    print("Testing quaternion -> euler (roundtrip)")
    print("---------------------------------------")
    #Generate a random quaternion
    qf = rotlib.quat_rand(rng)
    print(" Quaternion " + ' '.join([f"{x:g}" for x in qf]) )
    #Convert qf to euler
    for seqi in seq:
        for worldj in world:
            print("\n")
            euler = rotlib.quat_to_euler(qf, seqi, worldj)
            print(f" Euler {seqi} world={worldj}"
                   f" {euler[0]:g} {euler[1]:g} {euler[2]:g}")
            #Convert euler back to quaternion
            qb = rotlib.euler_to_quat(euler, seqi, worldj)
            print(" Quaternion <- euler " + ' '.join([f"{x:g}" for x in qb]))
            assert np.allclose(qf, qb, 1e-8, 1e-14)


def test_aa_to_dcm_roundtrip():
    print("Testing axis angle -> dcm (roundtrip)")
    print("-------------------------------------")
    #Generate a random axis angle
    axisf, anglef = rotlib.aa_rand(rng)
    print(f" Axis " + ' '.join([f"{x:g}" for x in axisf]) 
          + f"  Angle {anglef:g}")
    #Convert axis angle to dcm
    dcm = rotlib.aa_to_dcm(axisf, anglef)
    print(" Direction cosine matrix")
    buf = ''
    for i in range(3):
        buf += f"  {dcm[i,0]:g} {dcm[i,1]:g} {dcm[i,2]:g}\n"
    print(buf[:-1])
    #Convert dcm back to axis angle
    axisb, angleb = rotlib.dcm_to_aa(dcm)
    print(f" Axis " + ' '.join([f"{x:g}" for x in axisb]) 
          + f"  Angle {angleb:g}")
    assert ( np.allclose(axisf, axisb, 1e-8, 1e-14) and
             np.isclose(anglef, angleb, 1e-8, 1e-14) )


def test_aa_to_euler_roundtrip():
    seq = ['XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'] 
    world = [True, False]

    print("Testing axis angle -> euler (roundtrip)")
    print("---------------------------------------")
    #Generate a random axis angle
    axisf, anglef = rotlib.aa_rand(rng)
    print(f" Axis " + ' '.join([f"{x:g}" for x in axisf]) 
          + f"  Angle {anglef:g}")
    #Convert axis angle to euler
    for seqi in seq:
        for worldj in world:
            print("\n")
            euler = rotlib.aa_to_euler(axisf, anglef, seqi, worldj)
            print(f" Euler {seqi} world={worldj}"
                   f" {euler[0]:g} {euler[1]:g} {euler[2]:g}")
            #Convert euler back to axis angle
            axisb, angleb = rotlib.euler_to_aa(euler, seqi, worldj)
            print(f" Axis " + ' '.join([f"{x:g}" for x in axisb]) 
                  + f"  Angle {angleb:g}")
            assert ( np.allclose(axisf, axisb, 1e-8, 1e-14) and
                     np.isclose(anglef, angleb, 1e-8, 1e-14) )


def test_dcm_to_euler_roundtrip():
    seq = ['XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'] 
    world = [True, False]

    print("Testing dcm -> euler (roundtrip)")
    print("--------------------------------")
    #Generate a random dcm
    axis, angle = rotlib.aa_rand(rng)
    dcmf = rotlib.aa_to_dcm(axis, angle)
    print(" Direction cosine matrix")
    buf = ''
    for i in range(3):
        buf += f"  {dcmf[i,0]:g} {dcmf[i,1]:g} {dcmf[i,2]:g}\n"
    print(buf[:-1])
    #Convert dcm to euler
    for seqi in seq:
        for worldj in world:
            print("\n")
            euler = rotlib.dcm_to_euler(dcmf, seqi, worldj)
            print(f" Euler {seqi} world={worldj}"
                   f" {euler[0]:g} {euler[1]:g} {euler[2]:g}")
            #Convert euler back to dcm
            dcmb = rotlib.euler_to_dcm(euler, seqi, worldj)
            print(" Direction cosine matrix")
            buf = ''
            for i in range(3):
                buf += f"  {dcmb[i,0]:g} {dcmb[i,1]:g} {dcmb[i,2]:g}\n"
            print(buf[:-1])
            assert np.allclose(dcmf, dcmb, 1e-8, 1e-14)


def test_euler_to_euler_roundtrip():
    seq = ['XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'] 
    world = [True, False]

    print("Testing euler -> euler (roundtrip)")
    print("----------------------------------")
    #Generate a random axis angle and convert to euler
    axis, angle = rotlib.aa_rand(rng)

    for i in range(len(seq)-1):
        for j in range(len(world)-1):
            print("\n")
            eulerf = rotlib.aa_to_euler(axis, angle, seq[i], world[j])
            print(f" Euler {seq[i]} world={world[j]}"
                   f" {eulerf[0]:g} {eulerf[1]:g} {eulerf[2]:g}")
            #Convert one euler sequence to another
            for ii in range(i, len(seq)):
                for jj in range(len(world)):
                    print("\n")
                    euler = rotlib.euler_to_euler(eulerf, seq[i], world[j], 
                                                    seq[ii], world[jj])
                    print(f"   Euler {seq[ii]} world={world[jj]}"
                           f" {euler[0]:g} {euler[1]:g} {euler[2]:g}")
                    #Convert the euler sequence back 
                    eulerb = rotlib.euler_to_euler(euler, seq[ii], world[jj], 
                                                    seq[i], world[j])
                    print(f"   Euler {seq[i]} world={world[j]}"
                           f" {eulerb[0]:g} {eulerb[1]:g} {eulerb[2]:g}")
                    assert np.allclose(eulerf, eulerb, 1e-8, 1e-14)


def test_quat_from_any():
    q = rotlib.quat_rand(rng)
    axis, angle = rotlib.quat_to_aa(q)
    dcm = rotlib.quat_to_dcm(q)
    euler = rotlib.quat_to_euler(q, 'XYZ', True)

    q_from_any = rotlib.quat_from_any({'repr':'quat', 'quat': q})
    assert np.allclose(q, q_from_any, 1e-8, 1e-14)

    q_from_any = rotlib.quat_from_any({'repr':'axis_angle', 'axis': axis,
                                       'angle': angle})
    assert np.allclose(q, q_from_any, 1e-8, 1e-14)

    q_from_any = rotlib.quat_from_any({'repr':'dcm', 'dcm': dcm})
    assert np.allclose(q, q_from_any, 1e-8, 1e-14)

    q_from_any = rotlib.quat_from_any(
            {'repr':'euler', 'euler': euler, 'seq': 'XYZ', 'world': True})
    assert np.allclose(q, q_from_any, 1e-8, 1e-14)


def test_aa_from_any():
    q = rotlib.quat_rand(rng)
    axis, angle = rotlib.quat_to_aa(q)
    dcm = rotlib.quat_to_dcm(q)
    euler = rotlib.quat_to_euler(q, 'XYZ', True)

    axisfa, anglefa = rotlib.aa_from_any({'repr':'axis_angle', 'axis': axis,
                                       'angle': angle})
    assert np.allclose(axisfa, axis, 1e-8, 1e-14)
    assert np.isclose(anglefa, angle, 1e-8, 1e-14)

    axisfa, anglefa = rotlib.aa_from_any({'repr':'quat', 'quat': q})
    assert np.allclose(axisfa, axis, 1e-8, 1e-14)
    assert np.isclose(anglefa, angle, 1e-8, 1e-14)

    axisfa, anglefa = rotlib.aa_from_any({'repr':'dcm', 'dcm': dcm})
    assert np.allclose(axisfa, axis, 1e-8, 1e-14)
    assert np.isclose(anglefa, angle, 1e-8, 1e-14)

    axisfa, anglefa = rotlib.aa_from_any(
            {'repr':'euler', 'euler': euler, 'seq': 'XYZ', 'world': True})
    assert np.allclose(axisfa, axis, 1e-8, 1e-14)
    assert np.isclose(anglefa, angle, 1e-8, 1e-14)


def test_angvel_to_qdot_roundtrip():
    #Inconvenient to do qdot -> angvel round trip as then for a given q, qdot
    #has to be chosen such that it is orthogonal to q. qdot does not have to be
    #a unit vector.

    print("Testing angvel -> qdot (roundtrip)")
    print("----------------------------------")

    q = rotlib.quat_rand(rng)
    omegaf = rng.uniform(-1.0, 1.0, (3,))
    print(" omegaf " + ' '.join([f"{x:g}" for x in omegaf]) )

    qdot = rotlib.quat_deriv_from_angvel(q, omegaf)
    print(" qdot " + ' '.join([f"{x:g}" for x in qdot]) )

    omegab = rotlib.quat_deriv_to_angvel(q, qdot)
    print(" omegab " + ' '.join([f"{x:g}" for x in omegab]) )

    assert np.allclose(omegaf, omegab, 1e-8, 1e-14)


def test_align_vector_1():
    #Generate a random vector
    v = rng.random((3,))
    #Generate a random axis-angle
    axis, angle = rotlib.aa_rand(rng)
    #Rotate `v` and the x-axis by axis angle
    x = np.array([1, 0, 0], dtype=np.float64)
    vx = np.vstack((v,x))
    vx_rot = rotlib.aa_rotate_vectors(vx, axis, angle)
    vrot = vx_rot[0,:]
    xrot = vx_rot[1,:]
    #Angle between v and x
    angle_vx = np.acos(np.dot(v,x)/np.linalg.norm(v))
    #Align `v` such that x is along xrot
    valigned = rotlib.align(v, x, xrot)
    #Angle between valigned and xrot
    angle_va_xrot = np.acos(np.dot(valigned,xrot)/np.linalg.norm(valigned))

    #Angle must be preserved on alignment
    assert np.isclose(angle_vx, angle_va_xrot, 1e-8, 1e-14)

    #Check self alignment
    xaligned = rotlib.align(x, x, xrot)
    assert np.allclose(xrot, xaligned, 1e-8, 1e-14)


def test_align_vector_2():
    #Generate a random vector
    v = rng.random((3,))
    #Generate a random axis-angle
    axis, angle = rotlib.aa_rand(rng)
    #Rotate `v` and the x and y axes by axis angle
    x = np.array([1, 0, 0], dtype=np.float64)
    y = np.array([0, 1, 0], dtype=np.float64)
    vxy = np.vstack((v,x,y))
    vxy_rot = rotlib.aa_rotate_vectors(vxy, axis, angle)
    vrot = vxy_rot[0,:]
    xrot = vxy_rot[1,:]
    yrot = vxy_rot[2,:]
    #Align `v` such that x is along xrot and y along yrot
    valigned = rotlib.align(v, vxy[1:,:], vxy_rot[1:,:])
    assert np.allclose(vrot, valigned, 1e-8, 1e-14)

    #Check self alignment
    xaligned = rotlib.align(x, vxy[1:,:], vxy_rot[1:,:])
    yaligned = rotlib.align(y, vxy[1:,:], vxy_rot[1:,:])
    assert np.allclose(xrot, xaligned, 1e-8, 1e-14)
    assert np.allclose(yrot, yaligned, 1e-8, 1e-14)


def test_align_vector_3():
    #Generate a random vector
    v = rng.random((4,3))
    #Generate a random axis-angle
    axis, angle = rotlib.aa_rand(rng)
    #coordinate axes
    axes = np.identity(3, dtype=np.float64)
    #Rotate `v` and the coordinate axes by axis angle
    vrot = rotlib.aa_rotate_vectors(v, axis, angle)
    axes_rot = rotlib.aa_rotate_vectors(axes, axis, angle)
    #Align `v` such that axes is along axes_rot
    valigned = rotlib.align(v, axes, axes_rot)
    assert np.allclose(vrot, valigned, 1e-8, 1e-14)

    #Check self alignment
    axes_aligned = rotlib.align(axes, axes, axes_rot)
    assert np.allclose(axes_rot, axes_aligned, 1e-8, 1e-14)
