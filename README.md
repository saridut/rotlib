# rotlib

Rotlib converts between quaternions, axis-angle representations, direction
cosine matrices, and Euler angles. It also rotates vectors and transforms
vector and tensor components between frames.

The Python, C, and Fortran implementations are maintained separately and follow
shared [conventions](docs/rotconv.md). See the
[API documentation](https://saridut.github.io/rotlib) for function contracts.

## Installation

### Python

Requires Python 3.12+ and NumPy 2.0+. In a virtual environment, install from the
`main` branch using Git:

```bash
python -m pip install --upgrade "rotlib @ git+https://github.com/saridut/rotlib.git@main#subdirectory=python"
```

### C and Fortran

Clone the repository and build from its root:

```bash
git clone https://github.com/saridut/rotlib.git
cd rotlib
make -C c/build lib
make -C fortran/build lib
```

| Language | Compiler | Source directory | Library |
| --- | --- | --- | --- |
| C | GCC | `c/src` | `c/lib/librotlib.a` |
| Fortran | gfortran | `fortran/src` | `fortran/lib/librotlib.a` |

Fortran module files are generated in `fortran/build`.

## Development

Clone the whole repository to work across all three implementations. Install
Python tooling in a local virtual environment with the `dev` extra:

```bash
git clone https://github.com/saridut/rotlib.git
cd rotlib
python3 -m venv python/.venv
source python/.venv/bin/activate
python -m pip install -e "./python[dev]"
```

The `dev` extra includes Python test, lint, and documentation tools. Install
native compilers, Doxygen, FORD, and native test frameworks separately.

Run the Python checks from the repository root:

```bash
python -m pytest python/tests
python -m ruff check python
```

Use `python -m tox -c python/pyproject.toml` for the configured Python-version
matrix and lint checks; the corresponding interpreters must be installed.

The existing C and Fortran test drivers still target older interfaces.
CppUTest/pFUnit integration and shared cases are planned.

## Building documentation

Run from the repository root:

| API | Tool | Command | HTML output |
| --- | --- | --- | --- |
| Python | Sphinx | `make -C python/docs doc` | `docs/python_api` |
| C | Doxygen | `make -C c/build doc` | `docs/c_api` |
| Fortran | FORD | `make -C fortran/build doc` | `docs/fortran_api` |

For Python, run `make -C python/docs rst` initially and when the module
structure changes to regenerate API stubs. Each HTML output directory contains
an `index.html` for local preview.

## Examples

Convert Euler angles to a unit quaternion and back. Angles are in radians,
stored in X, Y, Z order; these examples use the fixed-axis `XYZ` sequence.
Extracted Euler angles use principal ranges and may differ from the input
while representing the same rotation.

### Python

```python
import numpy as np
import rotlib

ea = np.array([0.2, 0.5, 1.0])
q = rotlib.euler_to_quat(ea, seq="XYZ", world=True)
ea_back = rotlib.quat_to_euler(q, seq="XYZ", world=True, tol=1e-14)
```

### C

```c
#include "rotlib.h"

int main(void) {
    const double ea[3] = {0.2, 0.5, 1.0};
    double q[4], ea_back[3];

    euler_to_quat(ea, ES_XYZ, true, q);
    quat_to_euler(q, ES_XYZ, true, 1e-14, ea_back);
    return 0;
}
```

### Fortran

```fortran
program example
    use constants_m, only: rp
    use rotlib_m, only: euler_to_quat, quat_to_euler
    implicit none
    real(rp) :: ea(3), q(4), ea_back(3)

    ea = [0.2_rp, 0.5_rp, 1.0_rp]
    call euler_to_quat(ea, 'XYZ', .true., q)
    call quat_to_euler(q, 'XYZ', .true., 1e-14_rp, ea_back)
end program example
```

Python and C store vector batches and frame bases as rows; Fortran stores them
as columns. Quaternions are scalar-first in all three implementations.

Further Python examples are available in `python/tests`.
