/*

!!begin
!!title    Computes the light-curve of a sphere and a Roche distorted star
!!author   Jiao Li
!!created  08 Oct  2022
!!descr    Computes the light-curve of a sphere and a Roche distorted star
!!index    lroche
!!lib      libpylcurve
!!css      style.css
!!class    Model
!!head1    Computes the light-curve of a sphere and a Roche distorted star

!!emph{lroche} computes the light curve equivalent to a model of a sphere and
a Roche-distorted star to model a white dwarf/main-sequence binary such as NN
Ser and can optionally include a disc and bright-spot as well. It includes
reprocessing and eclipses. The reprocessing is computed using the simple
addition of fluxes method, assuming that the irradiating star can be treated
as a point, and not including any 'back heating'. Phase 0 is defined as the
point when star 1 is furthest from the observer.  The star grids can de adapted 
somewhat to be focussed on eclipses and long exposures trapezoidally sub-divided 
to model smearing. Other physics included in !!emph{lroche}: Doppler beaming, 
gravitational lensing, Roemer time delays, asynchronous rotation of the stellar components.
Physics not included: eccentric orbits!

!!emph{lroche} is also designed as a test routine and so can add gaussian
noise to the result. It can either read a data file and compute at exactly the
right times or it can compute the light curve over a regularly spaced
sequence. This file is the main documentation for the data and model file
formats.

!!head2 Invocation

lroche model data [time1 time2 ntime expose] noise seed nfile [output] (device) [(roff) scale [ssfac]/[sstar1 sstar2 sdisc spot]]

!!head2 Arguments

!!table
!!arg{model}{File of parameters specifying the parameters. See below for a full description of these.}
!!arg{data}{Data file as a template to compute the fit. 'none' to indicate that
you do not have one and will define times instead. In the case of real data, the fit
will be scaled automatically to minimise chi**2. Data files requires six columns which
are: mid-exposure times, exposure times (same units as times), fluxes, uncertainties on fluxes, weight factors and 
finally integers to represent the number of subdivisions to model finite exposure effects. The weight factors
allow you to change the weighting of a point to the overall chi**2 without changing its uncertainty.}
!!arg{time1}{If data = 'none', first time to compute.}
!!arg{time2}{If data = 'none', last time to compute.}
!!arg{ntime}{If data = 'none', number of times to compute.}
!!arg{expose}{If data = 'none', length of exposure}
!!arg{noise}{If data = 'none', amount of noise to add to the results, RMS. In the case of real data, 
this is used a multiplier of the real error bars.}
!!arg{seed}{Seed integer}
!!arg{nfile}{Number of files to store, each with its own noise; only the first will be plotted if
nfile > 1. 0 is possible in which case there will be no output.}
!!arg{output}{File to save the results in the form of rows each with
time, exposure time, flux and uncertainty. If nfile > 1, then the files 
will have a number added, as in data001}
!!arg{device}{Plot device to use; 'none' or 'null' to ignore}
!!arg{roff}{Offset to add to the residuals when plotting a fit to data}
!!arg{scale}{Autoscale or not. If true, then the scale factors will be determined by minimisation of chi**2 using
linear scaling factors.}
!!arg{ssfac}{If scale=false, and a single global scaling factor is set, this is its value. Depending on the units of the fluxes, ssfac has the following interpretation. If ssfac is set =
(a/d)^2 (a = semi-major axis, d = distance), the result returned will be the flux density at Earth (fnu) in SI units (W/m^2/Hz). A milliJansky = 10^-29 W/m^2/Hz, so if you are fitting mJy
fluxes, the value of ssfac returned = 10^29 (a/d)^2. If you express a in solar radii and d in parsecs, and fit mJy fluxes, then ssfac = 5.0875 x 10^13 (a/d)^2.}
!!arg{sstar1}{If scale=false, scale factor for star 1}
!!arg{sstar2}{If scale=false, scale factor for star 2}
!!arg{sdisc}{If scale=false, scale factor for disc}
!!arg{sspot}{If scale=false, scale factor for spot}

!!table

!!head2 Parameter file

The model parameters come in two types, physical and computational. 'Physical' in this case are ones which have an initial
value, a range of plausible variation and a step size for derivative computation. 'Computational' are ones which cannot 
be varied and just have a value. !!emph{lroche} ignores the variable/non-variable distinction which are for the fitting
routines such as !!ref{levmarq.html}{levmarq}. The parameter file consists of a series of lines such as:
<pre>
q      =  0.12  0.01 0.0001 0 1
iangle =  82    2    0.01   1 1
r1     =  0.17  0.05 0.001  1 1
.
.
.
etc

delta_phase = 1.e-7 .  .  .  etc </pre> etc. The above would mean that q is
not to be varied (0 as the penultimate value), but iangle and r1 are. For
!!ref{simplex.html}{simplex}, !!ref{genetic.html}{genetic},
!!ref{simann.html}{simann} and !!ref{powell.html}{powell} the second parameter
specifies the range over which to vary the respective parameter (i.e. 0.01, 2
and 0.05 in this case), while the third parameter is used by
!!ref{levmarq.html}{levmarq} to compute numerical derivatives using finite
differences. Be careful to set this small enough to give accurate derivatives
but not so small that roundoff will be problematic. It should, at a minimum, be
smaller than the uncertainty in any parameter. The final integer is a relatively recent addition which specifies whether a given parameter is used at all and is designed for the star spot parameters.

Here is the full list of parameters.

!!head2 Physical parameters

!!head3 Binary and stars

!!table
!!arg{q}{Mass ratio, q = M2/M1}
!!arg{iangle}{Inclination angle, degrees}
!!arg{r1}{Radius of star 1, scaled by the binary separation}
!!arg{r2}{Radius of star 2, scaled by the binary separation. The radius is measured along the line of centres towards star 1. Set = -1 and hold fixed for
Roche lobe filling stars.}
!!arg{cphi3}{Third contact phase (star 1 starting to emerge from eclipse). This is an alternative way to specify the radii, based on a 
spherical approximation fot the two stars, i.e. unless the stars are spherical, it is not quite the true third contact. The radii will 
be computed from the contact phases according to the two equations r2+r1 = sqrt(1 - sin^2 i cos^2 (2*pi*cphi4)) and 
r2-r1 = sqrt(1 - sin^2 i cos^2 (2*pi*cphi3)). The radii returned are precise, just the interpretation as contact phases that is not 
precise. cphi3 and cphi4 need the boolean use_radii set to 0 to enabled. The reason for using them is to help with MCMC iterations as 
they prevent the nasty curved correlation between r1, r2 and i. This can save a huge amount of CPU time.}
!!arg{cphi4}{Fourth contact phase, star 1 fully emerged from eclipse. See cphi3 for details.}
!!arg{spin1}{This is the ratio of the spin frequency of star 1 to the orbital frequency. In this case a modified form of the Roche potential is used for star 1}
!!arg{spin2}{This is the ratio of the spin frequency of star 2 to the orbital frequency. In this case a modified form of the Roche potential is used for star 2}
!!arg{t1}{Temperature of star 1, Kelvin. This is really a substitute for surface brightness which is set assuming a black-body given this parameter. If it was not for irradiation that would be exactly what this is, a one-to-one replacement for surface brightness. Irradiation however introduces bolometric luminosities effectively and breaks the direct link. Some would then argue that one must use model atmospheres except at the moment irradiated model atmosphere are in their infancy.}
!!arg{t2}{Temperature of star 2, Kelvin. Set < 0 in order that it does not get scaled when using the iscale parameter.}
!!arg{ldc1_1, etc}{Limb darkening for stars is quite hard to specify precisely. Here we adopt a 4 coefficient 
approach which can either represent a straighforward polynimal expanion of the form I(mu) = 1 - \sum_i a_i (1-mu)^i, 
or rather better in some cases Claret's 4-coefficient formula I(mu) = 1 - \sum_i a_i (1 - mu^(i/2)) (i=1 to 4). 
You specify these by supplying the 4 coefficients for each star (which for form's sake are potentially variable but 
you would probably be unwise to let them be free) and later on a parameter to say whether it is the polynomial or 
Claret's representation. The polynomial allows one to use linear and quadratic limb darkening amongst others by setting 
the upper coefficients = 0. ldc1_1 is the first coefficient of star 1, ldc1_2 is the second, etc, while ldc2_1 
is the first coefficient for star 2 etc. See limb1, limb2, mucrit1, mucrit2 below.}
!!arg{velocity_scale}{Velocity scale,  sum of unprojected orbital speeds, used for accounting for 
Doppler beaming and gravitational lensing. On its own this makes little difference to the light curve, so you should
not usually let it be free, but you might want to if you have independent K1 or K2 information which you can apply as part of a prior.}
!!arg{beam_factor1}{The factor to use for Doppler beaming from star 1. This corresponds to the factor (3-alpha)
that multiplies -v_r/c in the standard beaming formula where alpha is related to the spectral shape. Use of this parameter
requires the velocity_scale to be set.}
!!arg{beam_factor2}{The factor to use for Doppler beaming from star 2. This corresponds to the factor (3-alpha)
that multiplies -v_r/c in the standard beaming formula where alpha is related to the spectral shape. Use of this parameter
requires the velocity_scale to be set.}
!!table

!!head3 General

!!table
!!arg{t0}{Zero point of ephemeris, marking time of mid-eclipse (or in general
superior conjunction) of star 1, same units as times.}
!!arg{period}{Orbital period, same units as times.}
!!arg{pdot}{Quadratic coefficient of ephemeris, same units as times}
!!arg{deltat}{Time shift between the primary and secondary eclipses to allow for small eccentricities and Roemer delays in the orbit. The 
sign is defined such that deltat > 0 implies that the secondary eclipse suffers a delay compared to the primary compared to precisely 0.5 
difference. deltat < 0 implies the secondary eclipse comes a little earlier than expected. Assuming that the "primary eclipse" is the eclipse 
of star 1, then, using the same sign convention, the Roemer delay is given by = P*(K1-K2)/(Pi*c) where P is the orbital period, K1 and K2 are
the usual projected radial velocity semi-amplitudes Pi = 3.14159.., and c = speed of light. See Kaplan (2010) for more details. The delay
is implemented by adjusting the orbital phase according to phi' = phi + (deltat/2/P)*(cos(2*Pi*phi)-1), i.e. there is no change at primary eclipse
but a delay of -deltat/P by the secondary eclipse.}
!!arg{gravity_dark}{Gravity darkening coefficient. Only matters for the Roche distorted case, but is prompted for always.
There are two alternatives for this. In the standard old method, the temperatures on the stars are set equal to t2*(g/gr)**gdark 
where g is the gravity at a given point and gr is the gravity at the point furthest from the primary (the 'backside' of the secondary). 
For a convectuive atmosphere, 0.08 is the usual value while 0.25 is the number for a radiative atmosphere. This is translated into
intensity using a blackbody approx. If you want to bypass the BB approx and invoke a direct relation flux ~ (g/gr)**gdark relation
you should set gdark_bolom (see below) to 0 (false.)}
!!arg{absorb}{The fraction of the irradiating flux from star 1 absorbed by star 2}
!!arg{slope, quad, cube}{Fudge factors to help cope with any trends in the data
as a result of e.g. airmass effects. The fit is multiplied by 
(1+x*(slope+x*(quad+x*cube))) where x is the time scaled so that it varies 
from -1 to 1 from start to end of the data. One should expect these number
to have absolute value << 1.}  
!!arg{third}{Third light contribution. Simply adds to whatever flux is
calculated and will be subject to auto-scaling like other flux. It only 
applies if global scaling rather than individual component scaling is used.
Third light is assumed strictly constant and is not subject to the slope, 
quad, cube parameters.}
!!table

!!head3 Spots

One spot allowed on each star (with some expectation that the
number may be increased if need be in the future):

!!table
!!arg{stsp11_long}{Longitude (degrees) of spot 1 on star 1, relative to
meridian defined by line of centres}
!!arg{stsp11_lat}{Latitude (degrees) of spot 1 on star 1}
!!arg{stsp11_lat}{FWHM (degrees) of spot 1 on star 1, as seen from its centre
of mass. Spot has gaussian distribution of temperature.}
!!arg{stsp11_tcen}{Central temp (K) of spot 1 on star 1}
!!arg{stsp21_long}{Longitude (degrees) of spot 1 on star 2, relative to
meridian defined by line of centres}
!!arg{stsp21_lat}{Latitude (degrees) of spot 1 on star 2}
!!arg{stsp21_lat}{FWHM (degrees) of spot 1 on star 2, as seen from its centre
of mass. Spot has gaussian distribution of temperature.}
!!arg{stsp21_tcen}{Central temp (K) of spot 1 on star 2}
!!table

!!head3 Disc

!!table
!!arg{rdisc1}{Inner radius of azimuthally symmetric disc. Set = -1 to set it equal to r1 (it should not be allowed to vary in this case)}
!!arg{rdisc2}{Outer radius of azimuthally symmetric disc. Set = -1 and hold fixed to clamp this to equal the bright spot radius.}
!!arg{height_disc}{Half height of disc at radius = 1. The height varies as a power law of radius}
!!arg{beta_disc}{Exponent of power law in radius of disc. Should be >= 1 to make concave disc; convex will not eclipse
properly.}
!!arg{temp_disc}{Temperature of outer part of disc. This is little more than a flux normalisation parameter but 
it is easier to think in terms of temperature}
!!arg{texp_disc}{Exponent of surface brightness (NB: not temperature) over disc}
!!arg{lin_limb_disc}{Linear limb darkening coefficient of the disc}
!!arg{quad_limb_disc}{Quadratic limb darkening coefficient of the disc}
!!arg{temp_edge}{Temperature at perpendicular edge of disc. Irradiation from the secondary is allowed so you should think of a bright rim at primary eclipse. Limb darkeining parameters of the 
disc are applied}
!!arg{absorb_edge}{Amount of secondary flux absorbed and reprocessed. This effect should lead 
to a sinusoidal variation with flux maximum at orbital phase 0.5. It was introduced to model
a possible accreting sdO/WD system discovered by Thomas Kupfer}

!!table

!!head3 Bright-spot

!!table
!!arg{radius_spot}{Distance from accretor of bright-spot (units of binary separation).}
!!arg{length_spot}{Length scale of spot (units of binary separation).}
!!arg{height_spot}{Height of spot (units of binary separation). This is only a normalisation constant.}
!!arg{expon_spot}{Spot is modeled as x**n*exp(-(x/l)**m). This parameter specifies the exponent 'n'}
!!arg{epow_spot}{This is the exponent m in the above expression}
!!arg{angle_spot}{This is the angle made by the line of elements of the spot measured in the direction of binary motion relative to 
the rim of the disc so that the "standard" value should be 0.}
!!arg{yaw_spot}{Allows the spot elements effectively to beam their light away from the perpendicular to the line of elements.
Measured as an angle in the same sense as angle_spot. 0 means standard perpendicular beaming.}
!!arg{temp_spot}{Normalises the surface brightness of the spot.}
!!arg{tilt_spot}{Allows spot to be other than perpendicular to the disc. 90 = perpendicular. If less than 90 then the
spot is visible for more than half a cycle.}
!!arg{cfrac_spot}{The fraction of the spot taken to be equally visible at all phases, i.e. pointing upwards.}
!!table

!!head2 Computational parameters
!!table
!!arg{delta_phase}{Accuracy in phase of eclipse computations. This determines the
accuracy of any Roche computations. Example: 1.e-7}
!!arg{nlat1f}{The number of latitudes for star 1's fine grid. This is used around the phase of primary eclipse (i.e. the
eclipse of star 1}
!!arg{nlat1c}{The number of latitudes for star 1's coarse grid. This is used away from primary eclipse.}
!!arg{nlat2f}{The number of latitudes for star 2's fine grid. This is used around the phase of secondary eclipse.}
!!arg{nlat2c}{The number of latitudes for star 2's coarse grid. This is used away from secondary eclipse.}
!!arg{npole}{True to set North pole of grid to the genuine stellar NP rather than substellar points. This is probably a good idea
when modelling well detached binaries, especially with extreme radius ratios because then it allows one to concentrate points
over a band of latitudes using the next two parameters}
!!arg{nlatfill}{Extra number of points to insert per normal latitude strip along the path of star 1 as it transits star 2. This is
designed to help tough extreme radius ratio cases. Take care to look at the resulting grid with visualise as the exact latitude range
chosen is a little approximate. This is only enabled if npole since only then do the latitude strips more-or-less line up with the movement
of the star.}
!!arg{nlngfill}{Extra number of points to insert per normal longitude strip along the path of star 1 as it transits star 2. This is
designed to help tough extreme radius ratio cases. Take care to look at the resulting grid with visualise as the exact latitude range
chosen is a little approximate.}

!!arg{lfudge}{The fine-grid latitude strip is computed assuming both stars are spherical. To allow for departures from this, this parameter
allows one to increase the latitude limits both up and down by an amount specified in degrees. Use the program !!ref{visualise.html}{visualise}
to judge how large this should be. However, one typically would like to avoid lfudge > 30*r1/r2 as that could more than double the width of the strip.}

!!arg{phase1}{this defines when star 1's fine grid is used abs(phase) < phase1. Thus phase1 = 0.05 will restrict
the fine grid use to phase 0.95 to 0.05.}
!!arg{llo, lhi}{These are experimental. They allow the user to fix the latitude limits of the fine strip which might be useful in preventing chi**2 variations
caused by variable grids. The values need to reflect the likely range of inclinations and can only really be set by trial and error using visualise. They are in degrees
following the usual convention for latitude on Earth. Set llo high and lhi low to stop them having any effect.}
!!arg{phase1}{this defines when star 1's fine grid is used abs(phase) < phase1. Thus phase1 = 0.05 will restrict
the fine grid use to phase 0.95 to 0.05.}
!!arg{phase2}{this defines when star 2's fine grid is used phase2 until 1-phase2. Thus phase2 = 0.45 will restrict
the fine grid use to phase 0.55 to 0.55.}
!!arg{nrad}{The number of radial strips over the disc}
!!arg{wavelength}{Wavelength (nm)}
!!arg{roche1}{Account for Roche distortion of star 1 or not}
!!arg{roche2}{Account for Roche distortion of star 2 or not}
!!arg{eclipse1}{Account for the eclipse of star 1 or not}
!!arg{eclipse2}{Account for the eclipse of star 2 or not}
!!arg{glens1}{Account for gravitational lensing by star 1. If you use this roche1 must be = 0 and the velocity_scale}
!!arg{use_radii}{If set = 1, the parameters r1 and r2 will be used to set the radii directly. If not, the third and fourth contact phases,
cphi3 and cphi4, will be used instead (see description for cphi3 for details).}
!!arg{tperiod}{The true orbital period in days. This is required, along with velocity_scale, if gravitational lensing is 
being applied to calculate proper dimensions in the system.}
!!arg{gdark_bolom}{True if the gravity darkening coefficient represents the bolometric value where T is proportional to
gravity to the power set by the coefficient. This is translated to flux variations using the black-body approximation.
If False, it represents a filter-integrated value 'y' coefficient such that the flux depends upon the gravity to the power 'y'.
This is itself an approximation and ideally should replaced by a proper function of gravity, but is probably good enough for
most purposes. Please see gravity_dark.}
!!arg{mucrit1}{Critical value of mu on star 1 below which intensity is assumed to be zero. This is to allow one to represent
Claret and Hauschildt's (2004) results where I(mu) drops steeply for mu < 0.08 or so. WARNING: this option is dangerous. I would normally
advise setting it = 0 unless you really know what you are doing as it leads to discontinuities.}
!!arg{mucrit2}{Critical value of mu on star 2 below which intensity is assumed to be zero. See comments on mucrit1 for more.}
!!arg{limb1}{String, either 'Poly' or 'Claret' determining the type of limb darkening law. See comments on ldc1_1 above.}
!!arg{limb2}{String, either 'Poly' or 'Claret' determining the type of limb darkening law. See comments on ldc1_1 above.}
!!arg{mirror}{Add any light not reprocessed in as if star reflected it or not as a crude approximation to the
effet of gray scattering}
!!arg{add_disc}{Add a disc or not}
!!arg{opaque}{Make disc opaque or not}
!!arg{iscale}{Individually scale the separate components or not. If set the
each component, star 1, star 2, disc and bright spot will be individually
scaled to minimise chi**2. Otherwise a single overall factor will be computed.
NB If you set this parameter then all temperature parameters (white dwarf,
secondary, disc and bright spot) must be held fixed otherwise near-total
degeneracy will result. The only reason it is not total is because of
reflection effect from irradiation of the secondary by the white dwarf, but
this is often very feeble and will not help, so, you have been warned.
Scaling should in general lead to faster convergence than not scaling.
You may have some cases where you do not want to include any secondary
star component. You can do this by setting t2 < 0. Note that if this is set
true, then the third light parameter will be ignored.}
!!table

!!end
* This file is modified by lijiao
*/

#include <climits>
#include <cfloat>
#include <cstdlib>
#include <iostream>
#include "cpgplot.h"
#include "trm/subs.h"
#include "trm/plot.h"
#include "trm/vec3.h"
#include "trm/input.h"
#include "trm/format.h"
#include "trm/roche.h"
#include "trm/lcurve.h"

int Lcurve::Fobj::neval = 0;
double Lcurve::Fobj::chisq_min;
Subs::Buffer1D<double> Lcurve::Fobj::scale_min;

// Main program
extern "C"{
    void pylcurve(double *times, double *exposes, int *ndivs, int Tsize,
                  double *calc, double *lcstar1, double *lcdisc,
                  double *lcedge, double *lcspot, double *lcstar2, double *wdwarflogrv,
                  //Binary and stars 
                  double q_value, double q_range, double q_dstep, bool q_vary, bool q_defined,
                  double iangle_value, double iangle_range, double iangle_dstep, bool iangle_vary, bool iangle_defined, 
                  double r1_value, double r1_range, double r1_dstep, bool r1_vary, bool r1_defined, 
                  double r2_value, double r2_range, double r2_dstep, bool r2_vary, bool r2_defined, 
                  double cphi3_value, double cphi3_range, double cphi3_dstep, bool cphi3_vary, bool cphi3_defined, 
                  double cphi4_value, double cphi4_range, double cphi4_dstep, bool cphi4_vary, bool cphi4_defined, 
                  double spin1_value, double spin1_range, double spin1_dstep, bool spin1_vary, bool spin1_defined, 
                  double spin2_value, double spin2_range, double spin2_dstep, bool spin2_vary, bool spin2_defined, 
                  double t1_value, double t1_range, double t1_dstep, bool t1_vary, bool t1_defined, 
                  double t2_value, double t2_range, double t2_dstep, bool t2_vary, bool t2_defined, 
                  double ldc1_1_value, double ldc1_1_range, double ldc1_1_dstep, bool ldc1_1_vary, bool ldc1_1_defined, 
                  double ldc1_2_value, double ldc1_2_range, double ldc1_2_dstep, bool ldc1_2_vary, bool ldc1_2_defined, 
                  double ldc1_3_value, double ldc1_3_range, double ldc1_3_dstep, bool ldc1_3_vary, bool ldc1_3_defined, 
                  double ldc1_4_value, double ldc1_4_range, double ldc1_4_dstep, bool ldc1_4_vary, bool ldc1_4_defined, 
                  double ldc2_1_value, double ldc2_1_range, double ldc2_1_dstep, bool ldc2_1_vary, bool ldc2_1_defined, 
                  double ldc2_2_value, double ldc2_2_range, double ldc2_2_dstep, bool ldc2_2_vary, bool ldc2_2_defined, 
                  double ldc2_3_value, double ldc2_3_range, double ldc2_3_dstep, bool ldc2_3_vary, bool ldc2_3_defined, 
                  double ldc2_4_value, double ldc2_4_range, double ldc2_4_dstep, bool ldc2_4_vary, bool ldc2_4_defined, 
                  double velocity_scale_value, double velocity_scale_range, double velocity_scale_dstep, bool velocity_scale_vary, bool velocity_scale_defined, 
                  double beam_factor1_value, double beam_factor1_range, double beam_factor1_dstep, bool beam_factor1_vary, bool beam_factor1_defined, 
                  double beam_factor2_value, double beam_factor2_range, double beam_factor2_dstep, bool beam_factor2_vary, bool beam_factor2_defined, 
                  //General
                  double t0_value, double t0_range, double t0_dstep, bool t0_vary, bool t0_defined, 
                  double period_value, double period_range, double period_dstep, bool period_vary, bool period_defined, 
                  double pdot_value, double pdot_range, double pdot_dstep, bool pdot_vary, bool pdot_defined, 
                  double deltat_value, double deltat_range, double deltat_dstep, bool deltat_vary, bool deltat_defined, 
                  double gravity_dark1_value, double gravity_dark1_range, double gravity_dark1_dstep, bool gravity_dark1_vary, bool gravity_dark1_defined, 
                  double gravity_dark2_value, double gravity_dark2_range, double gravity_dark2_dstep, bool gravity_dark2_vary, bool gravity_dark2_defined, 
                  double absorb_value, double absorb_range, double absorb_dstep, bool absorb_vary, bool absorb_defined, 
                  double slope_value, double slope_range, double slope_dstep, bool slope_vary, bool slope_defined, 
                  double quad_value, double quad_range, double quad_dstep, bool quad_vary, bool quad_defined, 
                  double cube_value, double cube_range, double cube_dstep, bool cube_vary, bool cube_defined, 
                  double third_value, double third_range, double third_dstep, bool third_vary, bool third_defined, 
                  // Star spots
                  double stsp11_long_value, double stsp11_long_range, double stsp11_long_dstep, bool stsp11_long_vary, bool stsp11_long_defined, 
                  double stsp11_lat_value, double stsp11_lat_range, double stsp11_lat_dstep, bool stsp11_lat_vary, bool stsp11_lat_defined, 
                  double stsp11_fwhm_value, double stsp11_fwhm_range, double stsp11_fwhm_dstep, bool stsp11_fwhm_vary, bool stsp11_fwhm_defined, 
                  double stsp11_tcen_value, double stsp11_tcen_range, double stsp11_tcen_dstep, bool stsp11_tcen_vary, bool stsp11_tcen_defined, 
                  double stsp12_long_value, double stsp12_long_range, double stsp12_long_dstep, bool stsp12_long_vary, bool stsp12_long_defined, 
                  double stsp12_lat_value, double stsp12_lat_range, double stsp12_lat_dstep, bool stsp12_lat_vary, bool stsp12_lat_defined, 
                  double stsp12_fwhm_value, double stsp12_fwhm_range, double stsp12_fwhm_dstep, bool stsp12_fwhm_vary, bool stsp12_fwhm_defined, 
                  double stsp12_tcen_value, double stsp12_tcen_range, double stsp12_tcen_dstep, bool stsp12_tcen_vary, bool stsp12_tcen_defined, 
                  double stsp13_long_value, double stsp13_long_range, double stsp13_long_dstep, bool stsp13_long_vary, bool stsp13_long_defined, 
                  double stsp13_lat_value, double stsp13_lat_range, double stsp13_lat_dstep, bool stsp13_lat_vary, bool stsp13_lat_defined, 
                  double stsp13_fwhm_value, double stsp13_fwhm_range, double stsp13_fwhm_dstep, bool stsp13_fwhm_vary, bool stsp13_fwhm_defined, 
                  double stsp13_tcen_value, double stsp13_tcen_range, double stsp13_tcen_dstep, bool stsp13_tcen_vary, bool stsp13_tcen_defined, 
                  double stsp21_long_value, double stsp21_long_range, double stsp21_long_dstep, bool stsp21_long_vary, bool stsp21_long_defined, 
                  double stsp21_lat_value, double stsp21_lat_range, double stsp21_lat_dstep, bool stsp21_lat_vary, bool stsp21_lat_defined, 
                  double stsp21_fwhm_value, double stsp21_fwhm_range, double stsp21_fwhm_dstep, bool stsp21_fwhm_vary, bool stsp21_fwhm_defined, 
                  double stsp21_tcen_value, double stsp21_tcen_range, double stsp21_tcen_dstep, bool stsp21_tcen_vary, bool stsp21_tcen_defined, 
                  double stsp22_long_value, double stsp22_long_range, double stsp22_long_dstep, bool stsp22_long_vary, bool stsp22_long_defined, 
                  double stsp22_lat_value, double stsp22_lat_range, double stsp22_lat_dstep, bool stsp22_lat_vary, bool stsp22_lat_defined, 
                  double stsp22_fwhm_value, double stsp22_fwhm_range, double stsp22_fwhm_dstep, bool stsp22_fwhm_vary, bool stsp22_fwhm_defined, 
                  double stsp22_tcen_value, double stsp22_tcen_range, double stsp22_tcen_dstep, bool stsp22_tcen_vary, bool stsp22_tcen_defined,
                  double uesp_long1_value, double uesp_long1_range, double uesp_long1_dstep, bool uesp_long1_vary, bool uesp_long1_defined,
                  double uesp_long2_value, double uesp_long2_range, double uesp_long2_dstep, bool uesp_long2_vary, bool uesp_long2_defined,                            
                  double uesp_lathw_value, double uesp_lathw_range, double uesp_lathw_dstep, bool uesp_lathw_vary, bool uesp_lathw_defined,
                  double uesp_taper_value, double uesp_taper_range, double uesp_taper_dstep, bool uesp_taper_vary, bool uesp_taper_defined,
                  double uesp_temp_value, double uesp_temp_range, double uesp_temp_dstep, bool uesp_temp_vary, bool uesp_temp_defined, 
                  // disc
                  double rdisc1_value, double rdisc1_range, double rdisc1_dstep, bool rdisc1_vary, bool rdisc1_defined,
                  double rdisc2_value, double rdisc2_range, double rdisc2_dstep, bool rdisc2_vary, bool rdisc2_defined, 
                  double height_disc_value, double height_disc_range, double height_disc_dstep, bool height_disc_vary, bool height_disc_defined, 
                  double beta_disc_value, double beta_disc_range, double beta_disc_dstep, bool beta_disc_vary, bool beta_disc_defined, 
                  double temp_disc_value, double temp_disc_range, double temp_disc_dstep, bool temp_disc_vary, bool temp_disc_defined, 
                  double texp_disc_value, double texp_disc_range, double texp_disc_dstep, bool texp_disc_vary, bool texp_disc_defined, 
                  double lin_limb_disc_value, double lin_limb_disc_range, double lin_limb_disc_dstep, bool lin_limb_disc_vary, bool lin_limb_disc_defined, 
                  double quad_limb_disc_value, double quad_limb_disc_range, double quad_limb_disc_dstep, bool quad_limb_disc_vary, bool quad_limb_disc_defined, 
                  double temp_edge_value, double temp_edge_range, double temp_edge_dstep, bool temp_edge_vary, bool temp_edge_defined,
                  double absorb_edge_value, double absorb_edge_range, double absorb_edge_dstep, bool absorb_edge_vary, bool absorb_edge_defined, 
                  //Bright-spot
                  double radius_spot_value, double radius_spot_range, double radius_spot_dstep, bool radius_spot_vary, bool radius_spot_defined,
                  double length_spot_value, double length_spot_range, double length_spot_dstep, bool length_spot_vary, bool length_spot_defined, 
                  double height_spot_value, double height_spot_range, double height_spot_dstep, bool height_spot_vary, bool height_spot_defined, 
                  double expon_spot_value, double expon_spot_range, double expon_spot_dstep, bool expon_spot_vary, bool expon_spot_defined, 
                  double epow_spot_value, double epow_spot_range, double epow_spot_dstep, bool epow_spot_vary, bool epow_spot_defined, 
                  double angle_spot_value, double angle_spot_range, double angle_spot_dstep, bool angle_spot_vary, bool angle_spot_defined, 
                  double yaw_spot_value, double yaw_spot_range, double yaw_spot_dstep, bool yaw_spot_vary, bool yaw_spot_defined, 
                  double temp_spot_value, double temp_spot_range, double temp_spot_dstep, bool temp_spot_vary, bool temp_spot_defined, 
                  double tilt_spot_value, double tilt_spot_range, double tilt_spot_dstep, bool tilt_spot_vary, bool tilt_spot_defined, 
                  double cfrac_spot_value, double cfrac_spot_range, double cfrac_spot_dstep, bool cfrac_spot_vary, bool cfrac_spot_defined,
                  // Computational parameters
                  double delta_phase, int nlat1f, int nlat2f, int nlat1c, int nlat2c, bool npole, 
                  int nlatfill, int nlngfill, double lfudge, double llo, double lhi, double phase1, double phase2, int nrad, double wavelength,
                  bool roche1, bool roche2, bool eclipse1, bool eclipse2, bool glens1, bool use_radii,
                  double tperiod, double gdark_bolom1, double gdark_bolom2, double mucrit1, double mucrit2, 
                  const char* pslimb1, const char* pslimb2, bool mirror, bool add_disc, bool opaque, bool add_spot, int nspot, bool iscale, bool info,
                  int parallel_threshold
                   ){
        try{
            //added by lijiao  
            Lcurve::Model model( // Binary and stars
                          q_value,  q_range,  q_dstep,  q_vary,  q_defined,
                          iangle_value,  iangle_range,  iangle_dstep,  iangle_vary,  iangle_defined, 
                          r1_value,  r1_range,  r1_dstep,  r1_vary,  r1_defined, 
                          r2_value,  r2_range,  r2_dstep,  r2_vary,  r2_defined, 
                          cphi3_value,  cphi3_range,  cphi3_dstep,  cphi3_vary,  cphi3_defined, 
                          cphi4_value,  cphi4_range,  cphi4_dstep,  cphi4_vary,  cphi4_defined, 
                          spin1_value,  spin1_range,  spin1_dstep,  spin1_vary,  spin1_defined, 
                          spin2_value,  spin2_range,  spin2_dstep,  spin2_vary,  spin2_defined, 
                          t1_value,  t1_range,  t1_dstep,  t1_vary,  t1_defined, 
                          t2_value,  t2_range,  t2_dstep,  t2_vary,  t2_defined, 
                          ldc1_1_value,  ldc1_1_range,  ldc1_1_dstep,  ldc1_1_vary,  ldc1_1_defined, 
                          ldc1_2_value,  ldc1_2_range,  ldc1_2_dstep,  ldc1_2_vary,  ldc1_2_defined, 
                          ldc1_3_value,  ldc1_3_range,  ldc1_3_dstep,  ldc1_3_vary,  ldc1_3_defined, 
                          ldc1_4_value,  ldc1_4_range,  ldc1_4_dstep,  ldc1_4_vary,  ldc1_4_defined, 
                          ldc2_1_value,  ldc2_1_range,  ldc2_1_dstep,  ldc2_1_vary,  ldc2_1_defined, 
                          ldc2_2_value,  ldc2_2_range,  ldc2_2_dstep,  ldc2_2_vary,  ldc2_2_defined, 
                          ldc2_3_value,  ldc2_3_range,  ldc2_3_dstep,  ldc2_3_vary,  ldc2_3_defined, 
                          ldc2_4_value,  ldc2_4_range,  ldc2_4_dstep,  ldc2_4_vary,  ldc2_4_defined, 
                          velocity_scale_value,  velocity_scale_range,  velocity_scale_dstep,  velocity_scale_vary,  velocity_scale_defined, 
                          beam_factor1_value,  beam_factor1_range,  beam_factor1_dstep,  beam_factor1_vary,  beam_factor1_defined, 
                          beam_factor2_value,  beam_factor2_range,  beam_factor2_dstep,  beam_factor2_vary,  beam_factor2_defined, 
                         //General
                          t0_value,  t0_range,  t0_dstep,  t0_vary,  t0_defined, 
                          period_value,  period_range,  period_dstep,  period_vary,  period_defined, 
                          pdot_value,  pdot_range,  pdot_dstep,  pdot_vary,  pdot_defined, 
                          deltat_value,  deltat_range,  deltat_dstep,  deltat_vary,  deltat_defined, 
                          gravity_dark1_value,  gravity_dark1_range,  gravity_dark1_dstep,  gravity_dark1_vary,  gravity_dark1_defined, 
                          gravity_dark2_value,  gravity_dark2_range,  gravity_dark2_dstep,  gravity_dark2_vary,  gravity_dark2_defined, 
                          absorb_value,  absorb_range,  absorb_dstep,  absorb_vary,  absorb_defined, 
                          slope_value,  slope_range,  slope_dstep,  slope_vary,  slope_defined, 
                          quad_value,  quad_range,  quad_dstep,  quad_vary,  quad_defined, 
                          cube_value,  cube_range,  cube_dstep,  cube_vary,  cube_defined, 
                          third_value,  third_range,  third_dstep,  third_vary,  third_defined, 
                         // Star spots
                          stsp11_long_value,  stsp11_long_range,  stsp11_long_dstep,  stsp11_long_vary,  stsp11_long_defined, 
                          stsp11_lat_value,  stsp11_lat_range,  stsp11_lat_dstep,  stsp11_lat_vary,  stsp11_lat_defined, 
                          stsp11_fwhm_value,  stsp11_fwhm_range,  stsp11_fwhm_dstep,  stsp11_fwhm_vary,  stsp11_fwhm_defined, 
                          stsp11_tcen_value,  stsp11_tcen_range,  stsp11_tcen_dstep,  stsp11_tcen_vary,  stsp11_tcen_defined, 
                          stsp12_long_value,  stsp12_long_range,  stsp12_long_dstep,  stsp12_long_vary,  stsp12_long_defined, 
                          stsp12_lat_value,  stsp12_lat_range,  stsp12_lat_dstep,  stsp12_lat_vary,  stsp12_lat_defined, 
                          stsp12_fwhm_value,  stsp12_fwhm_range,  stsp12_fwhm_dstep,  stsp12_fwhm_vary,  stsp12_fwhm_defined, 
                          stsp12_tcen_value,  stsp12_tcen_range,  stsp12_tcen_dstep,  stsp12_tcen_vary,  stsp12_tcen_defined, 
                          stsp13_long_value,  stsp13_long_range,  stsp13_long_dstep,  stsp13_long_vary,  stsp13_long_defined, 
                          stsp13_lat_value,  stsp13_lat_range,  stsp13_lat_dstep,  stsp13_lat_vary,  stsp13_lat_defined, 
                          stsp13_fwhm_value,  stsp13_fwhm_range,  stsp13_fwhm_dstep,  stsp13_fwhm_vary,  stsp13_fwhm_defined, 
                          stsp13_tcen_value,  stsp13_tcen_range,  stsp13_tcen_dstep,  stsp13_tcen_vary,  stsp13_tcen_defined, 
                          stsp21_long_value,  stsp21_long_range,  stsp21_long_dstep,  stsp21_long_vary,  stsp21_long_defined, 
                          stsp21_lat_value,  stsp21_lat_range,  stsp21_lat_dstep,  stsp21_lat_vary,  stsp21_lat_defined, 
                          stsp21_fwhm_value,  stsp21_fwhm_range,  stsp21_fwhm_dstep,  stsp21_fwhm_vary,  stsp21_fwhm_defined, 
                          stsp21_tcen_value,  stsp21_tcen_range,  stsp21_tcen_dstep,  stsp21_tcen_vary,  stsp21_tcen_defined, 
                          stsp22_long_value,  stsp22_long_range,  stsp22_long_dstep,  stsp22_long_vary,  stsp22_long_defined, 
                          stsp22_lat_value,  stsp22_lat_range,  stsp22_lat_dstep,  stsp22_lat_vary,  stsp22_lat_defined, 
                          stsp22_fwhm_value,  stsp22_fwhm_range,  stsp22_fwhm_dstep,  stsp22_fwhm_vary,  stsp22_fwhm_defined, 
                          stsp22_tcen_value,  stsp22_tcen_range,  stsp22_tcen_dstep,  stsp22_tcen_vary,  stsp22_tcen_defined,
                          uesp_long1_value, uesp_long1_range, uesp_long1_dstep, uesp_long1_vary, uesp_long1_defined,
                          uesp_long2_value, uesp_long2_range, uesp_long2_dstep, uesp_long2_vary, uesp_long2_defined,                            
                          uesp_lathw_value, uesp_lathw_range, uesp_lathw_dstep, uesp_lathw_vary, uesp_lathw_defined,
                          uesp_taper_value, uesp_taper_range, uesp_taper_dstep, uesp_taper_vary, uesp_taper_defined,
                          uesp_temp_value, uesp_temp_range, uesp_temp_dstep, uesp_temp_vary, uesp_temp_defined, 
                         // disc
                          rdisc1_value,  rdisc1_range,  rdisc1_dstep,  rdisc1_vary,  rdisc1_defined,
                          rdisc2_value,  rdisc2_range,  rdisc2_dstep,  rdisc2_vary,  rdisc2_defined, 
                          height_disc_value,  height_disc_range,  height_disc_dstep,  height_disc_vary,  height_disc_defined, 
                          beta_disc_value,  beta_disc_range,  beta_disc_dstep,  beta_disc_vary,  beta_disc_defined, 
                          temp_disc_value,  temp_disc_range,  temp_disc_dstep,  temp_disc_vary,  temp_disc_defined, 
                          texp_disc_value,  texp_disc_range,  texp_disc_dstep,  texp_disc_vary,  texp_disc_defined, 
                          lin_limb_disc_value,  lin_limb_disc_range,  lin_limb_disc_dstep,  lin_limb_disc_vary,  lin_limb_disc_defined, 
                          quad_limb_disc_value,  quad_limb_disc_range,  quad_limb_disc_dstep,  quad_limb_disc_vary,  quad_limb_disc_defined, 
                          temp_edge_value,  temp_edge_range,  temp_edge_dstep,  temp_edge_vary,  temp_edge_defined,
                          absorb_edge_value,  absorb_edge_range,  absorb_edge_dstep,  absorb_edge_vary,  absorb_edge_defined, 
                         //Bright-spot
                          radius_spot_value,  radius_spot_range,  radius_spot_dstep,  radius_spot_vary,  radius_spot_defined,
                          length_spot_value,  length_spot_range,  length_spot_dstep,  length_spot_vary,  length_spot_defined, 
                          height_spot_value,  height_spot_range,  height_spot_dstep,  height_spot_vary,  height_spot_defined, 
                          expon_spot_value,  expon_spot_range,  expon_spot_dstep,  expon_spot_vary,  expon_spot_defined, 
                          epow_spot_value,  epow_spot_range,  epow_spot_dstep,  epow_spot_vary,  epow_spot_defined, 
                          angle_spot_value,  angle_spot_range,  angle_spot_dstep,  angle_spot_vary,  angle_spot_defined, 
                          yaw_spot_value,  yaw_spot_range,  yaw_spot_dstep,  yaw_spot_vary,  yaw_spot_defined, 
                          temp_spot_value,  temp_spot_range,  temp_spot_dstep,  temp_spot_vary,  temp_spot_defined, 
                          tilt_spot_value,  tilt_spot_range,  tilt_spot_dstep,  tilt_spot_vary,  tilt_spot_defined, 
                          cfrac_spot_value,  cfrac_spot_range,  cfrac_spot_dstep,  cfrac_spot_vary,  cfrac_spot_defined,
                         // Computational parameters
                          delta_phase,  nlat1f,  nlat2f,  nlat1c,  nlat2c,  npole, 
                          nlatfill,  nlngfill,  lfudge,  llo,  lhi,  phase1,  phase2,  nrad,  wavelength,
                          roche1,  roche2,  eclipse1,  eclipse2,  glens1,  use_radii,
                          tperiod,  gdark_bolom1,  gdark_bolom2,  mucrit1,  mucrit2, 
                          pslimb1,  pslimb2,  mirror,  add_disc,  opaque,  add_spot,  nspot,  iscale
                        );
             
    
    
            // Compute light curve
            double wdwarf, chisq, wnok, logg1, logg2, rv1, rv2;
            Lcurve::pylight_curve_comp(model, times, exposes, ndivs, Tsize, info, calc,
                                       lcstar1, lcdisc, lcedge, lcspot, lcstar2,
                                       wdwarf, logg1, logg2, rv1, rv2, parallel_threshold);
    
            wdwarflogrv[0] = wdwarf;
            wdwarflogrv[1] = logg1;
            wdwarflogrv[2] = logg2;
            wdwarflogrv[3] = rv1;
            wdwarflogrv[4] = rv2;
            if(info){
                std::cout << "White dwarf's contribution = " << wdwarf
                          << std::endl;
                std::cout << "log10(g1 [cgs]) = " << logg1 << std::endl;
                std::cout << "log10(g2 [cgs]) = " << logg2 << std::endl;
                std::cout << "Vol-averaged r1 = " << rv1 << std::endl;
                std::cout << "Vol-averaged r2 = " << rv2 << std::endl;
                }
    
        } //try
        catch(const Roche::Roche_Error& err){
            std::cerr << "Roche::Roche_Error exception thrown" << std::endl;
            std::cerr << "libpylcurve: " << err.what() << std::endl;
            exit(EXIT_FAILURE);
        }
        catch(const Lcurve::Lcurve_Error& err){
            std::cerr << "Lcurve::Lcurve_Error exception thrown" << std::endl;
            std::cerr << err << std::endl;
            exit(EXIT_FAILURE);
        }
        catch(const std::string& err){
            std::cerr << err << std::endl;
            exit(EXIT_FAILURE);
        }
        catch(...){
            std::cerr << "Unknown exception caught in lroche" << std::endl;
            exit(EXIT_FAILURE);
        }
    }//void pylcurve(


   void pylcurve_smodel(const char* psmodel, 
                      double *times, double *exposes, int *ndivs, int Tsize,
                      bool info,
                      double *calc, double *lcstar1, double *lcdisc,
                      double *lcedge, double *lcspot, double *lcstar2,
                      double *wdwarflogrv, int parallel_threshold){
        try{
            const std::string smodel = psmodel;
            Lcurve::Model model(smodel);
            //free(psmodel);
            // Compute light curve
            double wdwarf, logg1, logg2, rv1, rv2;
            Lcurve::pylight_curve_comp(model, times, exposes, ndivs, Tsize, info, calc,
                                       lcstar1, lcdisc, lcedge, lcspot, lcstar2,
                                       wdwarf, logg1, logg2, rv1, rv2, parallel_threshold);
            wdwarflogrv[0] = wdwarf;
            wdwarflogrv[1] = logg1;
            wdwarflogrv[2] = logg2;
            wdwarflogrv[3] = rv1;
            wdwarflogrv[4] = rv2;
            if(info){
                std::cout << "White dwarf's contribution = " << wdwarf
                          << std::endl;
                std::cout << "log10(g1 [cgs]) = " << logg1 << std::endl;
                std::cout << "log10(g2 [cgs]) = " << logg2 << std::endl;
                std::cout << "Vol-averaged r1 = " << rv1 << std::endl;
                std::cout << "Vol-averaged r2 = " << rv2 << std::endl;
                }
    
        } //try
        catch(const Roche::Roche_Error& err){
            std::cerr << "Roche::Roche_Error exception thrown" << std::endl;
            std::cerr << "libpylcurve: " << err.what() << std::endl;
            exit(EXIT_FAILURE);
        }
        catch(const Lcurve::Lcurve_Error& err){
            std::cerr << "Lcurve::Lcurve_Error exception thrown" << std::endl;
            std::cerr << err << std::endl;
            exit(EXIT_FAILURE);
        }
        catch(const std::string& err){
            std::cerr << err << std::endl;
            exit(EXIT_FAILURE);
        }
        catch(...){
            std::cerr << "Unknown exception caught in lroche" << std::endl;
            exit(EXIT_FAILURE);
        }
   }//pylcurve_file
} //extern "C"
