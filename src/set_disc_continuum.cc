#include <cstdlib>
#include <iostream>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/vec3.h"
#include "trm/lcurve.h"

/** set_disc_continuum computes the face-on brightness of each element of the
 * disc assuming a power law with radius.
 *
 * \param rdisc      radius of disc where the next parameter is defined
 * \param tdisc      the temperature at radius rdisc in the disc. This is
 *                   used to set the surface brightness scale.
 * \param texp       the exponent controlling the scaling of brightness with
 *                   radius. Negative to increase towards the centre of the
 *                   disc.
 * \param wave       wavelength of interest, nm
 * \param disc       grid of elements over disc
 */

void Lcurve::set_disc_continuum(double rdisc, double tdisc, double texp,
				double wave,
				Subs::Buffer1D<Lcurve::Point>& disc){

    // Reference surface brightness
    const double BRIGHT = Subs::planck(wave, tdisc);
    
    for(int i=0; i<disc.size(); i++){
	double r = disc[i].posn.length();
	disc[i].flux = BRIGHT*pow(r/rdisc, texp)*disc[i].area;
    }
}    

/** set_edge_continuum computes the face-on brightness of each element of the
 * edge of the disc assuming it intrinsically has temperature tedge
 * plus whatever it gets from irradiation. This to allow for unsual
 * sdO+WD accreting system found by Thomas Kupfer. Donor approximated
 * as point source of luminosity defined by its radius and temperature.
 * 
 * \param tedge      temperature at outer edge of disc
 * \param r2         radius of donor
 * \param t2         temperature of donor
 * \param absorb     amount of irradiation flux absorbed and reprocessed
 * \param wave       wavelength of interest, nm
 * \param edge       grid of elements over disc
 */

void Lcurve::set_edge_continuum(double tedge, double r2, double t2,
				double absorb, double wave,
				Subs::Buffer1D<Lcurve::Point>& edge){

    Subs::Vec3 vec;
    const Subs::Vec3 cofm2(1.,0.,0.);
    double temp, geom, mu, r;

    for(int i=0; i<edge.size(); i++){
        vec = cofm2 - edge[i].posn;
        r   = vec.length();
        mu  = Subs::dot(edge[i].dirn, vec)/r;

        if(mu >= r2){

            // Full tilt irradiation
            geom = Subs::sqr(r2/r)*mu;
            temp = pow(pow(tedge,4)  + absorb*pow(t2,4)*geom, 0.25);

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
            temp = pow(pow(tedge,4) + absorb*pow(t2,4)*geom, 0.25);

        }else{

            // No irradiation
            temp = tedge;

        }
        edge[i].flux = edge[i].area*Subs::planck(wave, temp);
    }
}    

