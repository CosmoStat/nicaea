Known bugs and shortcomings
===========================

- Some parameter combinations cause undefined behaviour of the
  program. These are (hopefully) intercepted and an error is created
  (see Sect. 5). E.g., for n_spec<0.7, f_NL (Peacock&Dodds) is not
  defined. For a closed Universe, the probed redshift can be larger
  than the maximum redshift.

- a=1.0 very rarely creates an error, use 0.99999... instead.

- The code is not well suited for Fisher matrix calculations. In particular
  for the inverse Fisher matrix, numerical derivatives have to be very
  accurate, and the interpolations between tabulated values (linear and
  spline) in nicaea introduce numerical noise that can render the Fisher
  matrix numerically singular (Wolz et al. 2012).

- Dark-energy models, in particular with varying w(z), are not recommended
  for the non_linear models smith03, and smith03_de. Instead, use the
  revised halofit model with smith03_revised.

In case of problems please don't hesitate to contact me at
martin.kilbinger@cea.fr . Questions and comments are welcome!


Changes compared to the Rob Smith's original halofit
====================================================

Parts of the program 'cosmo.c' is based on Rob Smiths' halofit (Smith et al.
2003). The code for determining the non-linear power spectrum has been improved
and made more efficient. The main changes are listed below. The code also
includes the non-linear fitting formulae of Peacock & Dodds (1996).

- Tabulation of the linear and non-linear power spectrum, constants
  are calculated only once.
- Integration cutoff for determination of non-linear scale knl
  flexible, as function of smoothing scale rmid; using Romberg
  integration.
- Bisection to find knl is iterative: if the bisection gets stuck at one
  end of the bisecting interval, the interval is shifted accordingly and
  a new bisection is started. If knl is larger than knlstern (I chose
  10^6 h/Mpc), the bisection is canceled and the linear power spectrum
  is used.
- Slope and curvature are calculated only once, after knl is fixed.
- The Eisenstein&Hu (1998) fit for the transfer function is used
  instead of Bond&Efstathiou (1984).
- The exact linear growth factor is used instead of the CPT92 fitting
  formula. Dark energy models are incorporated.
