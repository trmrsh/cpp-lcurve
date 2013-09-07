#include <cstdlib>
#include <iostream>
#include "trm/subs.h"
#include "trm/vec3.h"
#include "trm/lcurve.h"

/** set_star_emission computes the emission line face-on brightness of elements over the 
 * secondary star assuming proportionality to intensity of flux from primary.
 * \param limb2    limb darkening coefficient, star 2
 * \param hbyr     H/R ratio to screen equator of secondary
 * \param star2    grid of elements over star 2
 * \param bright2  brightness array of emission line contributions for star 2
 */

void Lcurve::set_star_emission(double limb2, double hbyr, const Subs::Buffer1D<Lcurve::Point>& star2, Subs::Buffer1D<float>& bright2){
  
    double mu, r, hrtest;
    Subs::Vec3 vec;
    const Subs::Vec3 cofm1(0.,0.,0.);
    bright2.resize(star2.size());

    for(int i=0; i<star2.size(); i++){
	vec = cofm1 - star2[i].posn;
	r   = vec.length();
	mu  = Subs::dot(star2[i].dirn,vec)/r;  
	hrtest = Subs::abs(vec.z())/sqrt(Subs::sqr(vec.x())+Subs::sqr(vec.y()));
	if(mu > 0. && hrtest >= hbyr)
	    bright2[i] = mu/Subs::sqr(r)/(1.-limb2/3.);
	else
	    bright2[i] = 0.;
    }
}    

