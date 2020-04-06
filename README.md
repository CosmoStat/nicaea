# nicaea [![Build Status](https://travis-ci.org/CosmoStat/nicaea.svg?branch=master)](https://travis-ci.org/CosmoStat/nicaea) [![Documentation Status](https://readthedocs.org/projects/nicaea/badge/?version=latest)](http://nicaea.readthedocs.io/en/latest/?badge=latest)
                
Numerical routines to calculate cosmology and weak-lensing quantities.

NumerIcal Cosmology And lEnsing cAlculations

Documentation: http://nicaea.readthedocs.io

Web page: http://cosmostat.org/nicaea

## Information

### Authors:

  - Martin Kilbinger

  - Karim Benabed

  - Jean Coupon (HOD, halomodel)

  - Henry J. McCracken (HOD)

### Contributors:

  - Liping Fu (decomp_eb)

  - Catherine Heymans (intrinsic alignment)

  - Francois Lanusse (github, python support)
  
  - Shiming Gu (massive neutrinos)

### Version

2.7.4

### Download, compile, and run nicaea

#### Download the code

Clone the code from the github repository,

```bash
git clone https://github.com/CosmoStat/nicaea
```

A new directory `nicaea` is created automatically. Create a build directory, and configure the code as follows.
You can specify the installation path with the option to cmake `-DCMAKE_INSTALL_PREFIX=<PATH>`.

```bash
cd nicaea
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=<PATH>
```

On success, compile install the code:

```bash
make
[sudo] make install
```

The last command will copy the executable demo programs (e.g. lensingdemo)
to `PATH`/bin, the library libnicaea.a to `PATH`/lib, and the include
files to `PATH`/include/nicaea . If you did not specify the option `-DCMAKE_INSTALL_PREFIX=PATH` above, the default
location is used, in which case you might need to prepend `sudo` to the `make` command,
to get write access to the system-wide installation directory.

If the necessary libraries are found on the system, the python module
pynicaea is also installed.

The code can be tested with::
```sh
ctest -vv
```
To run the demo programs, go to `nicaea_2.7/par_files`.

## Reference

To reference nicaea, please use Kilbinger, Benabed et al. (2009), A&A, 497, 677 (https://arxiv.org/abs/0810.5129), in which something that resembles the first version of nicaea has been used.
