:numbered:

Download, compile, and run nicaea
=================================

Download the code
-----------------

Download the file nicaea_2.7.tgz from http://cosmostat.org/nicaea and un-tar
the archive. The packages fftw3 and gsl are required to compile and run nicaea.
You can install fftw3 from http://www.fftw.org, and gsl from
www.gnu.org/software/gsl.

Compile and install the code
----------------------------

Two options to compile nicaea exist. If nicaea is to be used as a library,
option 1 is recommended.

Option 1: using cmake, *recommended*::

	cd build
	cmake ..
	make && make install

The last command will copy the executable demo programs (e.g. lensingdemo)
to <BASE>/bin, the library libnicaea.a to <BASE>/lib, and the include
files to <BASE>/include/nicaea . The default base directory is
<BASE>=nicaea_2.7 .

If the necessary libraries are found on the system, the python module
pynicaea is also installed.

The code can be tested with::

	ctest -vv

To run the demo programs (see below), go to nicaea_2.7/par_files .

Option 2: using make.::

	cd Demo
	make

If fftw3 and gsl are not installed in a standard directory (e.g. /usr,
/usr/local), set the variables 'FFTW' and 'GSL' in the Makefile. The header
file fftw3.h is looked for in $(FFTW)/include and libfftw3.a in $(FFTW)/lib.
The gsl header files are looked for in $(GSL)/include, the libraries libgsl.a
and libgslcblas.a in $(GSL)/lib.

Various demo programs can be run in ./Demo, see below.

Run the demo programs
---------------------

The demo programs need parameter files in the working directory, which can be
found in par_files.

+------------------------+--------------+-----------------------------------------------------------------------+
| Program name           | Category     | Functionality                                                       	|
+========================+==============+=======================================================================+
| lensingdemo		 | Weak lensing | density- and lensing power spectrum, lensing second-order functions 	|
+------------------------+--------------+-----------------------------------------------------------------------+
| sn1demo 		 | SNIa         | Luminosity distance, distance module				      	|
+------------------------+--------------+-----------------------------------------------------------------------+
| halomodeldemo		 | Halo model   | Power spectrum						      	|
+------------------------+--------------+-----------------------------------------------------------------------+
| cmb_bao_demo 		 | CMB, BAO     | geometrical quantities, e.g. sound horizon, angular diameter distance	|
+------------------------+--------------+-----------------------------------------------------------------------+
| decomp_eb_demo	 | Weak lensing | E-/B-mode decomposition (generaalized ring statistic)			|
+------------------------+--------------+-----------------------------------------------------------------------+
| cosebi_demo		 | Weak lensing | E-/B-mode decomposition (COSEBIs)					|
+------------------------+--------------+-----------------------------------------------------------------------+
| third_order_demo	 | Weak lensing | Third-order aperture-mass moments					|
+------------------------+--------------+-----------------------------------------------------------------------+
