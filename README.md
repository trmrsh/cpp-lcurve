This is a modified version of lcurve for the python wrapper [pyclurve](https://github.com/lidihei/pylcurve) 

- major modifications:
- - trm/lcurve.h add Pparam(double value, double range, double dstep, bool vary, bool defined) and
    Constructor from arguments added by lijiao
    Model(//Binary and stars
           double q_value, double q_range, double q_dstep, bool q_vary, bool q_defined,...
- - src/Makefile.am add lib_LTLIBRARIES = libpylcurve.la ...
- - src/trm_lcurve.cc  add Lcurve::Model::Model(//Binary and stars...
- - add libpylcurve.cc and pylight_curve_comp.cc

'lcurve' is for modelling of light curves.

Installation order: subs --> colly, binary --> roche --> lcurve

The file called 'Lcurve' that is generated in this directory defines
aliases (csh-style) and prints a help message when sourced.

Installation can be painful I'm afraid, but I don't have the time to make
a better job of it I am afraid. For details, please see:

 https://cygnus.astro.warwick.ac.uk/phsaap/software


This version allows python to call function in DLLs or shared libraries ([ctypes.CDLL](https://docs.python.org/3/library/ctypes.html)). 
The python wrapper of lcurve is [pylcurve](https://github.com/lidihei/pylcurve), which is written by Jiao Li (lijiao@bao.ac.cn).

