# YAM2: Yet another library for the M2 variables

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

YAM2 is free software library for calculating a set of kinematic variables M<sub>2</sub> using several numerical optimization algorithms.

## How to build

For building the source codes, it is required to have a C++ compiler, supporting the C++17 features, and [NLopt](https://github.com/stevengj/nlopt). For C++ compiler, check the pages given below.

* [C++ Standards Support in GCC](https://gcc.gnu.org/projects/cxx-status.html),
* [C++ Support in Clang](https://clang.llvm.org/cxx_status.html).

In theory, it will work with use GCC >= 7 or Clang >= 5.

Detailed instructions for installing NLopt by building the source code are given in [NLopt installation](https://nlopt.readthedocs.io/en/latest/NLopt_Installation/). For some Linux distributions, it can be installed by using the system package manager. For example,

```
# For Arch Linux/Manjaro:
sudo pacman -S nlopt

# For CentOS/Fedora:
sudo dnf install NLopt-devel

# For Debian/Ubuntu:
sudo apt-get install libnlopt-cxx-dev
```

For Ubuntu prior to 20.04 LTS (Focal Fossa), install `libnlopt-dev`. In CentOS, the [EPEL](https://fedoraproject.org/wiki/EPEL) repository must be installed and enabled. In macOS, one can install NLopt using [Homebrew](https://brew.sh/). We have tested our codes with the NLopt version >= 2.6.2.

The source code of the YAM2 library can be built by running `make`. If the path to NLopt is `/usr/local`, append the path to the `make` command:

```
NLOPT=/usr/local make
```

The command will build all the source codes, and then generate a static library file, `libYAM2.a`, in the `lib` directory. If the shared library is necessary, run `make lib`. See [`Makefile`](./Makefile) for the detail of the compilation flags and path settings.

The header and library files can also be installed to the other destination path outside build directory. If the path to be installed is `/usr/local`, run the command as follows.

```
DESTDIR=/usr/local make install
```

## How to use

The interfaces for using YAM2 are defined in the header file [`yam2.h`](./src/yam2.h). Users have to add the header to their analysis code through include directive.

``` c++
#include <yam2.h>
```

The type signature of the function for calculating M<sub>2CC</sub> can be seen in the following function declaration.

``` c++
std::optional<M2Solution> m2CCSQP(
    const std::optional<InputKinematics> &inp,
    double eps = EPS, int neval = NEVAL);
```

It will calculate M<sub>2CC</sub> using the sequential quadratic programming (SQP) method. For M<sub>2XC</sub>, the function to use is [`m2XCSQP`](./src/yam2.h). The function for calculating M<sub>2CC</sub> using the augmented Lagrangian method with the BFGS update is [`m2CCAugLagBFGS`](./src/yam2.h).

In the function declaration given above, one can see that the return type of the function is [`std::optional`](https://en.cppreference.com/w/cpp/utility/optional) of [`M2Solution`](./src/yam2.h). The class template `std::optional` causes a null value if the function has failed, or otherwise, it returns the contained value, that is, `M2Solution` in our case. The function fails if the input is incorrect or the function has eventually failed to find a minimum.

Once the calculation of the M<sub>2</sub> function is successful, the result can be extracted by the `value` method of `std::optional`.

``` c++
const auto m2sol = yam2::m2CCSQP(input.value());
if (!m2sol) {
    std::cerr << "Failed.\n";
} else {
    std::cout << "M2CC = " << m2sol.value().m2() << '\n'
              << "solution:\n"
              << "  k1: " << m2sol.value().k1() << '\n'
              << "  k2: " << m2sol.value().k2() << '\n';
}
```

As can be seen in the code snippet, the `M2Solution` class contains three methods: `m2` for the M<sub>2</sub> value, `k1` and `k2` for the M<sub>2</sub> solution to the invisible particle momenta. All the functions and classes are in the namespace of `yam2`.

### Input

There are three inputs to the functions for calculating M<sub>2</sub>. The first one is an instance of [`InputKinematics`](./src/input.h), which is for the particle momentum configuration of the given event. It can be constructed by using the [`mkInput`](./src/input.h) function,

``` c++
std::optional<InputKinematics> mkInput(
    const std::vector<FourMomentum> &as,
    const std::vector<FourMomentum> &bs,
    const TransverseMomentum &ptmiss, const Mass &minv);
```

Here `as` and `bs` correspond to the four-momenta of the visible particles for the decay process of

```
A1 + A2 --> a1 B1 + a2 B2 --> a1 b1 C1 + a2 b2 C2.
```

We stress that the order of the particle momenta should be set with care since it is not checked by the program: a<sub>i</sub> must be produced before having b<sub>i</sub> in the decay chain.

In addition to them, users have to insert the missing transverse momentum and the invisible particle mass into `ptmiss` and `minv`, respectively. See [`momentum.h`](./src/momentum.h) for the class definitions of `FourMomentum`, `TransverseMomentum`, and `Mass`.

The input momentum configuration should be validated before substituting it into the functions for calculating M<sub>2</sub>. An example code snippet using the `mkInput` is given below.

``` c++
const auto input =
    yam2::mkInput({a1, a2}, {b1, b2}, ptmiss, yam2::Mass{m_invis});
if (!input) {
    std::cerr << "Invalid input.\n";
}
const auto m2sol = yam2::m2CCSQP(input.value());
```

The other optional inputs to the [`m2CCSQP`](./src/yam2.h) function in the above are the tolerance `eps` and the maximal number of iterations `neval`. These will be set to the default values defined in [`yam2.h`](./src/yam2.h) unless users supply any input. In the current version of YAM2, their default values are `EPS` = 10<sup>-3</sup> and `NEVAL` = 5000. We recommend users to read the example analysis code enclosed with YAM2, [`examples/m2.cc`](./examples/m2.cc), before starting to write their analysis code for the M<sub>2</sub> variables. If you want to build the example code, run `make examples/m2`.

Supposing that the file name of the analysis code is `m2.cc` and the path to YAM2 is `/usr/local`, an example command for building an analysis code using YAM2 is as follows.

```
c++ -o m2.exe m2.cc -I/usr/local/include/YAM2 -L/usr/local/lib -lYAM2 -lnlopt
```

## References

* W.S. Cho et al, On-shell constrained M<sub>2</sub>​ variables with applications to mass measurements and topology disambiguation, JHEP 08 (2014) 070, [arXiv:1401.1449](https://arxiv.org/abs/1401.1449).

* W.S. Cho et al, OPTIMASS: A Package for the Minimization of Kinematic Mass Functions with Constraints, JHEP 01 (2016) 026, [arXiv:1508.00589](https://arxiv.org/abs/1508.00589).

* J. Nocedal and S. Wright, [Numerical Optimization](https://link.springer.com/book/10.1007/978-0-387-40065-5), Springer, 2006.
