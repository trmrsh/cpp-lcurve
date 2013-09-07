#include <cstdlib>
#include <iostream>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/vec3.h"
#include "trm/lcurve.h"

/** set_disc_continuum computes the face-on brightness of each element of the disc assuming a power
 * law with radius.
 *
 * \param rdisc      radius of disc where the next parameter is defined
 * \param tdisc      the temperature at radius rdisc in the disc. This is used to set the surface brightness scale.
 * \param texp       the exponent controlling the scaling of brightness with radius. Negative to increase towards the centre of the disc.
 * \param wave       wavelength of interest, nm
 * \param disc       grid of elements over disc
 */

void Lcurve::set_disc_continuum(double rdisc, double tdisc, double texp, double wave, Subs::Buffer1D<Lcurve::Point>& disc){

    // Reference surface brightness
    const double BRIGHT = Subs::planck(wave, tdisc);
    
    for(int i=0; i<disc.size(); i++){
	double r = disc[i].posn.length();
	disc[i].flux = BRIGHT*pow(r/rdisc, texp)*disc[i].area;
    }
}    

