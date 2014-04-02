#include <cstdlib>
#include <iostream>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/vec3.h"
#include "trm/lcurve.h"

/** set_star_continuum computes the continuum face-on brightness*area of each
 * element of the two stars assuming a black-body relation . The actual
 * contribution to the light-curve is the brightness times the area times (1 -
 * lin_limb*(1-mu) - quad_limb*(1-mu)**2), following Wade & Rucinski (1985)
 * where lin_limb and quad_limb are the linear and quadratic limb darkening
 * coefficients and mu is the cosine of the face angle.  This leads to a
 * correction factor of 1/(1-a/3-b/6) to the brightness to get the right
 * effective temperature.
 *
 * Irradiation is computed by adding the contribution of the other star to the
 * one being irradiated; the process is not iterated to compute
 * 'back-heating'. The finite size of the source is included only
 * approximately.
 *
 * \param mdl   Model defining parameters
 * \param star1        grid of elements over star 1, modified on exit
 * \param star2        grid of elements over star 2, modified on exit.
 */

void Lcurve::set_star_continuum(const Model& mdl,
                                Subs::Buffer1D<Lcurve::Point>& star1,
                                Subs::Buffer1D<Lcurve::Point>& star2){

    double r1, r2;
    mdl.get_r1r2(r1, r2);

    double rl1 = Roche::xl1(mdl.q);
    if(r1 < 0) r1 = rl1;

    double rl2 = 1.-rl1;
    if(r2 < 0) r2 = rl2;

    double mu, r, temp;
    Subs::Vec3 vec;
    const Subs::Vec3 cofm2(1.,0.,0.);

    // First compute irradiation of star 1 by star 2. 'geom' is the
    // geometrical factor by which the radiation from star 2 is reduced
    // by the time it hits the element in terms of flux per unit area
    // compared to its level as it leaves star 2
    double geom;

    // modify the gravity darkening coefficient to allow for two possibilities
    // the gravity darkening is implemented by modifying the temperature and
    // then calculating the flux in a BB approx. The 'bolometric' method does
    // this directly; the 'filter integrated' method modifies the exponent
    // used to give the desired behaviour of flux with gravity.
    const double GDCBOL1 = mdl.gdark_bolom1 ? mdl.gravity_dark1 :
        mdl.gravity_dark1 / Subs::dlpdlt(mdl.wavelength, mdl.t1);

    int nelem1 = star1.size();
    for(int i=0; i<nelem1; i++){
        vec = cofm2 - star1[i].posn;
        r   = vec.length();
        mu  = Subs::dot(star1[i].dirn,vec)/r;
        if(mu >= r2){

            // Full tilt irradiation
            geom = Subs::sqr(r2/r)*mu;
            temp = pow(pow(mdl.t1*pow(double(star1[i].gravity),GDCBOL1),4)
                       + mdl.absorb*pow(mdl.t2,4)*geom, 0.25);
        }else if(mu > -r2){

            // 'sunset' case
            double x0 = -mu/r2;
            // The following factor is a weighted version of 'mu' as the
            // secondary sets as far as this element is concerned.  When x0 =
            // -1 it equals r2 = mu. when x0 = 0 it equals 2*r2/(3*Pi) as
            // opposed to zero which it would be in the point source case.
            geom = Subs::sqr(r2/r)*r2*(sqrt(1.-x0*x0)*(2+x0*x0)/3 -
                                       x0*(Constants::PI/2-asin(x0)))/
                Constants::PI;
            temp = pow(pow(mdl.t1*pow(double(star1[i].gravity),GDCBOL1),4)
                       + mdl.absorb*pow(mdl.t2,4)*geom, 0.25);

        }else{

            // No irradiation
            geom = 0.;
            temp = mdl.t1*pow(double(star1[i].gravity),GDCBOL1);

        }

        // At this stage also add in a directly reflected part too
        star1[i].flux  = star1[i].area*Subs::planck(mdl.wavelength, temp);

        if(mdl.mirror) star1[i].flux  +=
                           star1[i].area*geom*
                           Subs::planck(mdl.wavelength,
                                        Subs::abs(double(mdl.t2)));

    }

    const Subs::Vec3 cofm1(0.,0.,0.);

    // See comments on GDCBOL1
    const double GDCBOL2 = mdl.gdark_bolom2 ? mdl.gravity_dark2 :
        mdl.gravity_dark2 / Subs::dlpdlt(mdl.wavelength,
                                         Subs::abs(double(mdl.t2)));

    int nelem2 = star2.size();
    for(int i=0; i<nelem2; i++){
        vec = cofm1 - star2[i].posn;
        r   = vec.length();
        mu  = Subs::dot(star2[i].dirn,vec)/r;

        if(mu >= r1){

            geom = Subs::sqr(r1/r)*mu;
            temp = pow(pow(Subs::abs(double(mdl.t2))*
                           pow(double(star2[i].gravity),GDCBOL2),4)
                       + mdl.absorb*pow(mdl.t1,4.)*geom, 0.25);

        }else if(mu > -r1){

            double x0 = -mu/r1;
            geom = Subs::sqr(r1/r)*r1*(sqrt(1.-x0*x0)*(2+x0*x0)/3
                                       - x0*(Constants::PI/2-asin(x0)))/
                Constants::PI;
            temp = pow(pow(Subs::abs(double(mdl.t2))*
                           pow(double(star2[i].gravity),GDCBOL2),4)
                       + mdl.absorb*pow(mdl.t1,4.)*geom, 0.25);

        }else{

            geom = 0.;
            temp = Subs::abs(double(mdl.t2))*pow(double(star2[i].gravity),
                                                 GDCBOL2);
        }

        star2[i].flux   = star2[i].area*Subs::planck(mdl.wavelength, temp);

        if(mdl.mirror) star2[i].flux  +=
                           star2[i].area*geom*
                           Subs::planck(mdl.wavelength, mdl.t1);
    }
}

