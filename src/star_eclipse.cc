#include <cstdlib>
#include "trm_roche.h"
#include "trm_lcurve.h"

/** Covenience routine which wraps up the code to compute a star eclipse allowing for roche-dirstortion or not.
 * \param q mass ratio
 * \param r radius of star
 * \param spin spin/orbital frequency factor
 * \param ffac roche filling factor
 * \param iangle inclination angle
 * \param posn position of point
 * \param delta accuracy in phase
 * \param roche account for roche distortion or not
 * \param eclipses set of ingress/egress pairs 
 */
void Lcurve::star_eclipse(double q, double r, double spin, double ffac, double iangle, const Subs::Vec3& posn, double delta, bool roche, 
			  Roche::STAR star, Lcurve::Point::etype& eclipses){
    double ri = Subs::deg2rad(iangle);
    const double cosi  = cos(ri);
    const double sini  = sin(ri);
    double ingress, egress, lam1, lam2;
    Subs::Vec3 cofm = (star == Roche::PRIMARY) ? Subs::Vec3(0.,0.,0.) : Subs::Vec3(1.,0.,0.);

    if((roche  && Roche::ingress_egress(q, star, spin, ffac, iangle, delta, posn, ingress, egress)) || 
       (!roche && Roche::sphere_eclipse(cosi, sini, posn, cofm, r, ingress, egress, lam1, lam2)))
	eclipses.push_back(std::make_pair(ingress, egress));
}
