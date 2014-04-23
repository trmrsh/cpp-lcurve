/*

!!begin
!!title    Computes the light-curve of a sphere and a Roche distorted star
!!author   T.R.Marsh 
!!created  02 Sep  2003
!!revised  20 Sep 2009
!!descr    Computes the light-curve of a sphere and a Roche distorted star
!!index    lroche
!!root     lroche
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
!!arg{ssfac}{If scale=false, and a single global scaling factor is set, this is its value.}
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
q      =  0.12  0.01 0.0001 0 
iangle =  82    2    0.01   1
r1     =  0.17  0.05 0.001  1
.
.
.
etc

delta_phase = 1.e-7
.
.
.
etc
</pre>
etc. The above would mean that q is not to be varied (0 at the end), but iangle and r1 are. For !!ref{simplex.html}{simplex}, 
!!ref{genetic.html}{genetic}, !!ref{simann.html}{simann} and !!ref{powell.html}{powell} the second parameter specifies
the range over which to vary the respective parameter (i.e. 0.01, 2 and 0.05 in this case), while the third parameter
is used by !!ref{levmarq.html}{levmarq} to compute numerical derivatives using finite differences. Be careful to set this small 
enough to give accurate derivatives but not so small that roundoff will be problematic. It should at a minimum be smaller than the
uncertainty in any parameter.

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
!!arg{beam_factor1}{The factor to use for Doppler beaming from star 2. This corresponds to the factor (3-alpha)
that multiplies -v_r/c in the standard beaming formula where alpha is related to the spectral shape. Use of this parameter
requires the velocity_scale to be set.}
!!table

!!head3 General

!!table
!!arg{t0}{Zero point of ephemeris, marking time of mid-eclipse (or in general superior conjunction) of star 1}
!!arg{period}{Orbital period, same units as time.}
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

!!arg{beta_disc}{Exponent of power law in radius of disc. Should be >= 1 to make concave disc; convex will not eclipse
properly.}
!!arg{temp_disc}{Temperature of outer part of disc. This is no more than a flux normalisation parameter but it s easier to think in
terms of temperature}
!!arg{texp_disc}{Exponent of surface brightness over disc}
!!arg{lin_limb_disc}{Linear limb darkening coefficient of the disc}
!!arg{quad_limb_disc}{Quadratic limb darkening coefficient of the disc}
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
!!arg{iscale}{Individually scale the separate components or not. If set the each component, star 1, star 2, disc and
bright spot will be individually scaled to minimise chi**2. Otherwise a single overall factor will be computed. NB
If you set this parameter then all temperature parameters (white dwarf, secondary, disc and bright spot) must be held
fixed otherwise near-total degeneracy will result. The only reason it is not total is because of reflection effect from
irradiation of the secondary by the white dwarf, but this is often very feeble and will not help, so, you have been warned.
Scaling should in general lead to faster convergence than not scaling. You may have some cases where you do not want to
include any secondary star component. You can do this by setting t2 < 0.}
!!table

!!end

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
int main(int argc, char* argv[]){

    try{

        // Construct Input object
        Subs::Input input(argc, argv, Lcurve::LCURVE_ENV, Lcurve::LCURVE_DIR);

        // Sign-in input variables
        input.sign_in("model",    Subs::Input::GLOBAL, Subs::Input::PROMPT);
        input.sign_in("data",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
        input.sign_in("time1",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("time2",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("ntime",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("expose",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("ndivide",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("noise",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("seed",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("nfile",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("output",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("device",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("roff",     Subs::Input::LOCAL,  Subs::Input::NOPROMPT);
        input.sign_in("scale",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("sstar1",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("sstar2",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("sdisc",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("sspot",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("ssfac",    Subs::Input::LOCAL,  Subs::Input::PROMPT);

        std::string smodel;
        input.get_value("model", smodel, "model",
                        "model file of parameter values");
        Lcurve::Model model(smodel);

        std::string sdata;
        input.get_value("data", sdata, "data",
                        "data file ('none' if you want to specify regularly-spaced times");
        bool no_file = (Subs::toupper(sdata) == "NONE");
        Lcurve::Data data, copy;
        if(!no_file){
            data.rasc(sdata);
            if(data.size() == 0)
                throw Lcurve::Lcurve_Error("No data read from file.");
            copy = data;
        }

        double time1, time2, expose, noise;
        int ntime, ndivide;
        if(no_file){
            input.get_value("time1", time1, 0.,
                            -DBL_MAX, DBL_MAX, "first time to compute");
            input.get_value("time2", time2, 1.,
                            -DBL_MAX, DBL_MAX, "last time to compute");
            input.get_value("ntime", ntime, 100, 2,
                            1000000, "number of times to compute");
            input.get_value("expose", expose, 0., 0.,
                            DBL_MAX, "exposure time (same units as ephemeris)");
            input.get_value("ndivide", ndivide, 1, 1,
                            INT_MAX, "ndivide factor to use");
            input.get_value("noise", noise, 0., 0.,
                            DBL_MAX, "RMS noise to add (after scaling)");
        }else{
            input.get_value("noise", noise, 0., 0.,
                            DBL_MAX, "RMS noise multiplier");
        }

        if(no_file){
            // Build fake data
            Lcurve::Datum datum = {0., expose, 0., noise, 1., ndivide};
            for(int i=0; i<ntime; i++){
                datum.time   = time1 + (time2-time1)*i/(ntime-1);
                data.push_back(datum);
            }
        }

        Subs::INT4 seed;
        input.get_value("seed", seed, 57565, INT_MIN, INT_MAX,
                        "random number seed");
        if(seed > 0) seed = -seed;
        int nfile;
        input.get_value("nfile", nfile, 1, 0, 20000,
                        "number of files to generate");
        std::string sout;
        if(nfile > 0)
            input.get_value("output", sout, "data", "file/root to save data");
        std::string device;
        input.get_value("device", device, "/xs", "plot device");
        double roff;
        input.get_value("roff", roff, 0., -DBL_MAX, DBL_MAX, "off to add to residuals when plotting a fit to data");

        bool scale = false;
        if(!no_file)
            input.get_value("scale", scale, true, "autoscale?");
        Subs::Buffer1D<double> sfac(4);
        if(!scale){
            if(model.iscale){
                input.get_value("sstar1", sfac[0], 1.,
                                -DBL_MAX, DBL_MAX, "star 1 scale factor");
                input.get_value("sstar2", sfac[1], 1.,
                                -DBL_MAX, DBL_MAX, "star 2 scale factor");
                input.get_value("sdisc",  sfac[2], 1.,
                                -DBL_MAX, DBL_MAX, "disc scale factor");
                input.get_value("sspot",  sfac[3], 1.,
                                -DBL_MAX, DBL_MAX, "spot scale factor");
            }else{
                input.get_value("ssfac", sfac[0], 1.,
                                -DBL_MAX, DBL_MAX, "global scale factor");
            }
        }
        input.save();

        // Construct function object
        Lcurve::Fobj func(model, data);

        // Compute light curve
        Subs::Array1D<double> fit;
        double wdwarf, chisq, wnok;
        Lcurve::light_curve_comp(model, data, scale, true, sfac,
                                 fit, wdwarf, chisq, wnok);

        Subs::Format form(12);

        // Save scale factors
        if(scale){
            std::cout << "Weighted chi**2 = " << form(chisq)
                      << ", wnok = " << form(wnok) << std::endl;
            if(model.iscale){
                input.set_default("sstar1", sfac[0]);
                input.set_default("sstar2", sfac[1]);
                input.set_default("sdisc",  sfac[2]);
                input.set_default("sspot",  sfac[3]);
                std::cout << "Scale factors = " << form(sfac[0]) << ", "
                          << form(sfac[1]) << ", " << form(sfac[2])
                          << ", " << form(sfac[3]) << std::endl;
            }else{
                input.set_default("ssfac", sfac[0]);
                std::cout << "Scale factor = " << form(sfac[0]) << std::endl;
            }
        }
        std::cout << "White dwarf's contribution = " << form(wdwarf)
                  << std::endl;

        if(!no_file){
            // Scale error bars
            for(size_t i=0; i<data.size(); i++)
                data[i].ferr *= noise;
        }

        int ndig = int(log10(double(nfile)+0.5))+1;

        // Add noise
        for(size_t i=0; i<data.size(); i++)
            data[i].flux = fit[i] + data[i].ferr*Subs::gauss2(seed);

        // make plot
        if((nfile == 0 || nfile == 1) && device != "none" && device != "null"){

            Subs::Plot plot(device);
            if(device.find("/xs") != std::string::npos){
                cpgscr(0,1,1,1);
                cpgscr(1,0,0,0);
            }
            cpgscr(2,0.7,0,0);
            cpgscr(3,0,0.6,0);
            cpgscr(4,0,0,0.5);
            cpgscr(5,0.7,0.7,0.7);
            double x1, x2;
            float y1, y2;
            x1 = x2 = data[0].time;
            y1 = data[0].flux - data[0].ferr;
            y2 = data[0].flux + data[0].ferr;
            for(size_t i=1; i<data.size(); i++){
                x1 = data[i].time > x1 ? x1 : data[i].time;
                x2 = data[i].time < x2 ? x2 : data[i].time;
                y1 = data[i].flux - data[i].ferr > y1 ? y1 :
                    data[i].flux - data[i].ferr;
                y2 = data[i].flux + data[i].ferr < y2 ? y2 :
                    data[i].flux + data[i].ferr;
                if(!no_file){
                    y1 = copy[i].flux - copy[i].ferr > y1 ? y1 :
                        copy[i].flux - copy[i].ferr;
                    y2 = copy[i].flux + copy[i].ferr < y2 ? y2 :
                        copy[i].flux + copy[i].ferr;
                    y1 = roff + copy[i].flux - data[i].flux - copy[i].ferr > y1 ? y1 : roff + copy[i].flux - data[i].flux - copy[i].ferr;
                    y2 = roff + copy[i].flux - data[i].flux + copy[i].ferr < y2 ? y2 : roff + copy[i].flux - data[i].flux + copy[i].ferr;
                }
            }
            double con = 0.;
            if(x2-x1 < 0.01*Subs::abs((x1+x2)/2.)){
                con = x1;
                x1 -= con;
                x2 -= con;
            }

            double range = x2-x1;
            x1   -= range/10.;
            x2   += range/10.;
            range = y2-y1;
            y1   -= range/10;
            y2   += range/10;

            cpgsch(1.5);
            cpgscf(2);
            cpgslw(4);
            cpgenv(float(x1), float(x2), y1, y2, 0, 0);
            std::string xlab = string("T - ") + Subs::str(con);
            cpglab(xlab.c_str(), " ", " ");

            if(!no_file){
                cpgsch(0.8);
                for(size_t i=0; i<copy.size(); i++){
                    cpgsci(5);
                    cpgslw(1);
                    cpgmove(copy[i].time-con, copy[i].flux - copy[i].ferr);
                    cpgdraw(copy[i].time-con, copy[i].flux + copy[i].ferr);

                    cpgsci(3);
                    cpgslw(3);
                    cpgpt1(copy[i].time-con, copy[i].flux, 17);

                    cpgsci(5);
                    cpgslw(1);
                    cpgmove(copy[i].time-con,
                            roff + copy[i].flux - data[i].flux - copy[i].ferr);
                    cpgdraw(copy[i].time-con,
                            roff + copy[i].flux - data[i].flux + copy[i].ferr);

                    cpgslw(3);
                    cpgsci(3);
                    cpgpt1(copy[i].time-con,
                           roff + copy[i].flux - data[i].flux, 17);
                }
            }
            if(noise == 0.0){
                cpgsci(1);
                cpgslw(5);
                cpgmove(data[0].time-con, data[0].flux);
                for(size_t i=1; i<data.size(); i++)
                    cpgdraw(data[i].time-con, data[i].flux);
            }else{
                for(size_t i=0; i<data.size(); i++){
                    cpgsci(2);
                    cpgmove(data[i].time-con, data[i].flux - data[i].ferr);
                    cpgdraw(data[i].time-con, data[i].flux + data[i].ferr);
                    cpgsci(3);
                    cpgpt1(data[i].time-con, data[i].flux, 1);
                }
            }
        }

        // Write out to disk
        if(nfile == 1){
            // note: noise will already have been added
            data.wrasc(sout);
            std::cout << "Written data to " << sout << std::endl;
        }else{
            for(int n=0; n<nfile; n++){
                // Add noise
                for(size_t i=0; i<data.size(); i++)
                    data[i].flux = fit[i] + data[i].ferr*Subs::gauss2(seed);

                std::string outfile = sout + Subs::str(n+1, ndig);
                data.wrasc(outfile);
                std::cout << "Written data to " << outfile << std::endl;
            }
        }
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
}
