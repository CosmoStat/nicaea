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

### Version

2.7

### Download, compile, and run nicaea

#### Download the code

Clone the code from the github repository,

```bash
git clone https://github.com/CosmoStat/nicaea
```

A new directory `nicaea` is created automatically. Change into its build directory, and configure the code as follows:

```bash
cd nicaea/build
cmake ..
```

On success, compile and install the code:

```bash
make
make install
```

The last command will copy the executable demo programs (e.g. lensingdemo)
to `BASE`/bin, the library libnicaea.a to `BASE`/lib, and the include
files to `BASE`/include/nicaea . The default base directory is
`BASE`=nicaea_2.7 .

If the necessary libraries are found on the system, the python module
pynicaea is also installed.

The code can be tested with::
```sh
ctest -vv
```
To run the demo programs, go to `nicaea_2.7/par_files`.

