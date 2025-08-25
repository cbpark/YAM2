# YAM2: Yet Another Library for the M₂ Variables

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Pipeline](https://gitlab.com/cbpark/yam2/badges/master/pipeline.svg)](https://gitlab.com/cbpark/yam2/-/pipelines/)

**YAM2** is a free and open-source C++ library for calculating the family of kinematic variables M<sub>2</sub>, widely used in collider phenomenology with invisible particles. It provides multiple numerical optimization backends (e.g., SQP and augmented Lagrangian) for solving constrained minimization problems efficiently. YAM2 is suitable for analyses at both **hadron** and **lepton** colliders.

---

## Table of Contents

- [Installation](#installation)
  - [Dependencies](#dependencies)
  - [Build with Make](#build-with-make)
  - [Shared Library](#shared-library)
  - [Build with CMake](#build-with-cmake)
- [Usage](#usage)
  - [Basic Example](#basic-example)
  - [Input Format](#input-format)
  - [M2Cons Example](#m2cons-example)
- [Citation](#citation)
- [References](#references)

---

## Installation

### Dependencies

- A C++17 compliant compiler (GCC ≥ 7, Clang ≥ 5, or Xcode ≥ 10 on macOS).
- [NLopt](https://github.com/stevengj/nlopt) (≥ 2.6.2).

You can install NLopt from source or via package managers:

``` bash
# Arch Linux / Manjaro
sudo pacman -S nlopt

# CentOS / Fedora
sudo dnf install NLopt-devel

# Debian / Ubuntu
sudo apt-get install libnlopt-cxx-dev

# openSUSE
sudo zypper install nlopt

# macOS (Homebrew)
brew install nlopt
```

For Ubuntu versions prior to 20.04 LTS, use `libnlopt-dev`.

---

### Build with Make

```
make
```

This generates a static library `lib/libYAM2.a`.

If NLopt is installed in a non-standard location (e.g., `/usr/local`), you can specify its path as a prefix to the `make` command:

```
NLOPT=/usr/local make
```

To install headers and libraries (e.g., into `/usr/local`):

```
sudo DESTDIR=/usr/local make install
```

---

### Shared Library

For dynamic linking (e.g., in [ROOT](https://root.cern.ch/) macros):

```
make lib
```

This produces:

- Linux: `lib/libYAM2.so`
- macOS: `lib/libYAM2.dylib`

Example usage in ROOT:

``` c++
R__LOAD_LIBRARY(/usr/local/lib/libYAM2.so)
R__LOAD_LIBRARY(/usr/lib/libnlopt.so)

gSystem->AddIncludePath("/usr/local/include/YAM2");
```

---

### Build with CMake

```
mkdir -p build && cd build
cmake -Dnlopt_DIR=/usr/local -DCMAKE_INSTALL_PREFIX=/usr/local ..
make
sudo make install
```

If ROOT is installed but fails to configure due to missing dependencies (e.g., `Vdt`), install [vdt](https://github.com/dpiparo/vdt).

---

## Usage

The main interface is provided in [`yam2.h`](./include/YAM2/yam2.h):

``` c++
#include <YAM2/yam2.h>
```

### Basic Example

The type signature of the function for calculating M<sub>2CC</sub> is given by:

``` c++
std::optional<M2Solution> m2CCSQP(
    const std::optional<InputKinematics> &inp,
    double eps = EPS, int neval = NEVAL);
```

This function computes M<sub>2CC</sub> using the sequential quadratic programming (SQP) method.

The return type is `std::optional` of [`M2Solution`](./include/YAM2/yam2.h) to safely handle cases where the calculation does not succeed. You should always check if the `optional` contains a value before accessing it. If the computation fails due to invalid input or failure to converge to a minimum, the function returns an empty `std::optional`. Otherwise, it contains a valid `M2Solution` object.

Once the M<sub>2</sub> calculation succeeds, the result can be accessed using the `value()` method of `std::optional`.

``` c++
const auto input =
    yam2::mkInput({a1, a2}, {b1, b2}, ptmiss, yam2::Mass{m_invis});
const auto m2sol = yam2::m2CCSQP(input);

if (!m2sol) {
    std::cerr << "Failed.\n";
} else {
    std::cout << "M2CC = " << m2sol.value().m2() << '\n'
              << "Invisible particle 1 momentum (k1): "
              << m2sol.value().k1() << '\n'
              << "Invisible particle 1 momentum (k1): "
              << m2sol.value().k2() << '\n';
}
```

Available solvers include:
- `m2CCSQP`: M<sub>2CC</sub> with sequential quadratic programming
- `m2XCSQP`: M<sub>2XC</sub> with SQP
- `m2CCAugLagBFGS`: M<sub>2CC</sub> with augmented Lagrangian (BFGS update)

---

### Input Format

The function `mkInput` constructs the required kinematics:

``` c++
std::optional<InputKinematics> mkInput(
    const std::vector<FourMomentum> &as,
    const std::vector<FourMomentum> &bs,
    const TransverseMomentum &ptmiss,
    const Mass &minv);
```

- `as`, `bs`: visible particle four-momenta
- `ptmiss`: missing transverse momentum
- `minv`: invisible particle mass

The input should be validated before use:

``` c++
const auto input =
    yam2::mkInput({a1, a2}, {b1, b2}, ptmiss, yam2::Mass{m_invis});
if (!input) {
    std::cerr << "Invalid input.\n";
}
```

Optional arguments: tolerance `eps` (default `1e-3`) and max iterations `neval` (default `5000`).

---

### M2Cons Example

From version 2.0, YAM2 supports M<sub>2Cons</sub> [arXiv:1509.00298](https://arxiv.org/abs/1509.00298):

``` c++
const auto input = yam2::mkInput(a1, a2, ptmiss, m_invis, {}, mY);
const auto m2sol = yam2::m2Cons(input);
```

Here, `mY` is the on-shell mass of the parent particle in the _antler_ decay topology. At lepton colliders, `sqrt_s` (the center-of-mass energy) and `pz` (the total longitudinal momentum) can be used as additional constraints.

``` c++
const auto input = yam2::mkInput(a1, a2, ptmiss, m_invis, {}, sqrt_s, {pz});
const auto m2sol = yam2::m2Cons(input);
```

Here, `pz` is the longitudinal momentum of the total system. See [examples/m2cons.cc](./examples/m2cons.cc) and [YAM2-ditau](https://github.com/cbpark/YAM2-ditau).

---

## Citation

If you use YAM2 in your work, please cite:

``` bibtex
@article{Park:2020bsu,
    author = "Park, Chan Beom",
    title = "{YAM2: Yet another library for the M2 variables using sequential quadratic programming}",
    eprint = "2007.15537",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "CTPU-PTC-20-18",
    doi = "10.1016/j.cpc.2021.107967",
    journal = "Comput. Phys. Commun.",
    volume = "264",
    pages = "107967",
    year = "2021"
}
```

Once the YAM2 2.0 paper is published, its BibTeX entry will be added here. In the meantime, please cite the original YAM2 paper.

---

## References

- C.B. Park, *YAM2: Yet another library for the M₂ variables using sequential quadratic programming*, Comput. Phys. Commun. 264 (2021) 107967, [arXiv:2007.15537](https://arxiv.org/abs/2007.15537).
- W.S. Cho et al., *On-shell constrained M₂ variables with applications to mass measurements and topology disambiguation*, JHEP 08 (2014) 070, [arXiv:1401.1449](https://arxiv.org/abs/1401.1449).
- W.S. Cho et al., *OPTIMASS: A Package for the Minimization of Kinematic Mass Functions with Constraints*, JHEP 01 (2016) 026, [arXiv:1508.00589](https://arxiv.org/abs/1508.00589).
- J. Nocedal and S. Wright, *Numerical Optimization*, Springer, 2006.
