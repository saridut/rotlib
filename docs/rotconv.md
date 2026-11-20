# Conventions

## Active and passive rotations

In the following, all functions dealing with rotations perform *active*
rotations following the right-hand rule. E.g., `aa_rotate_vectors` actively
rotates vectors about an axis through an angle in a given frame.

*Passive* rotation operations are referred to as *shifting*. They are performed
using a *shifter* matrix. Conceptually, these can be thought of as transforming
(or shifting) the components of a fixed vector (or tensor) from one frame to
another. E.g. `aa_shift_vectors` transforms the components of a set of
vectors from frame A to frame B, where frame B is obtained by rotating frame A
about an axis through an angle and the set of vectors as well as the axis are
resolved in frame A. The transformation from A to B is referred to as the
*forward* direction, whereas that from B to A is called the *reverse*
direction. Note that the axis is resolved in frame A irrespective of the
direction.

## Quaternions

Quaternions are represented by a vector of four elements, where the first
element represents the scalar part and the latter three represent the
vector part. Thus, assuming 0-based indexing, `q = q[0] + q[1]*i + q[2]*j + q[3]*k`.

A unit (or normalized) quaternion, used for representing orientations (or
rotations) in 3-space is a quaternion whose magnitude is unity. It is related
to a 3-D rotation about an axis `u` through an angle `theta` via the
following relations (assuming 0-based indexing):

```text
q[0] = cos(theta/2)
q[1] = u[0] * sin(theta/2)
q[2] = u[1] * sin(theta/2)
q[3] = u[2] * sin(theta/2)
```

In the canonical form, the first non-zero element of a unit quaternion is
positive.

Multiplication of quaternions and rotations by unit quaternions follow
Hamilton's convention. See, e.g.,

* Kuipers, J. B., *Quaternions and Rotation Sequences*, 1999, Princeton
  University Press, Princeton, NJ
* Hanson, A. J., *Visualizing Quaternions*, 2006, Morgan Kaufmann Publishers,
  San Francisco, CA

## Axis-angle

The axis-angle representation consists of a unit vector pointing along the axis
of rotation and an angle (in radian) of rotation. The rotation follows the
right-hand rule for positive angles. In the canonical form the angle is
restricted to lie in [0, pi].

## Direction cosine matrix (DCM)

The DCM gives the orientation of a frame B with respect to another frame A, and
transforms vector components from A to B. Its entries are the dot products of
B's basis vectors with A's basis vectors, with both sets expressed in the same
orthonormal frame F. Let `a_i` be the basis vectors of frame A and `b_i` those
of frame B resolved in a fixed frame F. The DCM is `D_ij = b_i . a_j`. In
matrix notation, `D = B . A^T` for row-wise storage of the basis vectors and
`D = B^T . A` assuming column-wise storage. If A and F coincide, then
`D = B` (row-based) or `D = B^T` (column-based).


## Euler angles

Euler angles are represented as three-element vectors, along with a rotation
sequence and a boolean flag `world` specifying the frame of rotation. The
elements represent rotations about X, Y, and Z, respectively. The order of
rotations is specified by a separate argument `seq`.

E.g.  `ea = [phi, theta, psi], seq = 'YZX', world=True` indicates a rotation
by `theta` along 'Y', followed by another rotation by `psi` along 'Z', and
finally a third rotation about 'X' by `phi`, where the rotations are with
respect to a fixed frame called the *world* frame. If `world=False`, the
rotations are performed with respect to the moving *body* frame.

The following table describes ranges of returned principal representations,
they do not restrict angles accepted by the construction functions.

| Sequence | phi | theta | psi |
| --- | --- | --- | --- |
| XYZ, ZYX | [-pi, pi] | [-pi/2, pi/2] | [-pi, pi] |
| XZY, YZX | [-pi, pi] | [-pi, pi] | [-pi/2, pi/2] |
| ZXY, YXZ | [-pi/2, pi/2] | [-pi, pi] | [-pi, pi] |

Functions that extract Euler angles accept a non-negative tolerance parameter
`tol`. In factorizing a rotation matrix near a singularity, if the absolute
value of the cosine of the middle angle is less than or equal to `tol`, one of
the outer angles is set to zero as per the following table.

| Sequence | Angle set to zero |
| --- | --- |
| XYZ, XZY | phi |
| YXZ, YZX | theta |
| ZXY, ZYX | psi |

Note that applying the above convention near, but not exactly at a singularity,
can introduce reconstruction error of the order of `tol`; it is not a
guaranteed reconstruction-error bound.

*Reference:* D. Eberly,
[Euler angle formulas](http://www.geometrictools.com/Documentation/EulerAngles.pdf).
