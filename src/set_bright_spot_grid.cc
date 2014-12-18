#include <cstdlib>
#include <iostream>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/vec3.h"
#include "trm/roche.h"
#include "trm/lcurve.h"

/**
 * set_bright_spot_grid sets up the elements needed to define the bright spot
 * at the edge of the disc. This is modelled as two straight lines of elements
 * coincident in position, but with one oriented parallel to the plane of the
 * disc and the other tilted by an arbitrary angle. This is to allow the
 * elements to make a hump while allowing flexibility in the step heights of
 * the bright-spot ingress and egress features. The brightness of the spot
 * along the lines is modelled as a function of distance x from the start of
 * each line as (x/l)**p1*exp(-(x/l)**p2) where l is a scale length and p1 and
 * p2 are power law exponents. The maximum of this occurs at x =
 * l*(p1/p2)**(1/p2), and is set to be the intersection of the gas stream with
 * the defined spot radius from the white dwarf.
 *
 * \param mdl       Model object defining parameters
 * \param spot      array representing the spot
 * \exception Exceptions are thrown if the specified radii over-fill the Roche lobes.
 */

void Lcurve::set_bright_spot_grid(const Model& mdl, Subs::Buffer1D<Lcurve::Point>& spot){

    double r1, r2;
    mdl.get_r1r2(r1, r2);

    double rl1 = Roche::xl11(mdl.q, mdl.spin1);
    if(r1 < 0)
        r1 = rl1;
    else if(r1 > rl1)
        throw Lcurve_Error("set_bright_spot_grid: the primary star is larger than its Roche lobe!");

    double rl2 = 1.-Roche::xl12(mdl.q, mdl.spin2);
    if(r2 < 0)
        r2 = rl2;
    else if(r2 > rl2)
        throw Lcurve_Error("set_bright_spot_grid: the secondary star is larger than its Roche lobe!");

    // Calculate a reference radius and potential for the two stars
    double rref1, pref1, ffac1 = r1/rl1;
    Roche::ref_sphere(mdl.q, Roche::PRIMARY, mdl.spin1, ffac1, rref1, pref1);

    double rref2, pref2, ffac2 = r2/rl2;
    Roche::ref_sphere(mdl.q, Roche::SECONDARY, mdl.spin2, ffac2, rref2, pref2);

    // Locate the spot position
    Subs::Vec3 bspot, v;
    Roche::strinit(mdl.q, bspot, v);
    Roche::stradv(mdl.q, bspot, v, mdl.radius_spot, 1.e-10, 1.e-3);

    // Now measure bright-spot angle relative to tangent to disc edge so we need
    // to add 90 + angle of bright-spot to the input value.
    double theta = atan2(bspot.y(), bspot.x()) + Constants::TWOPI/4. + Subs::deg2rad(mdl.angle_spot);

    double alpha = Subs::deg2rad(mdl.yaw_spot);
    double tilt  = Subs::deg2rad(mdl.tilt_spot);

    // The direction of the line of elements is set by angle_spot, but the
    // beaming direction adds in yaw_spot as well.
    Subs::Vec3 posn, dirn, bvec(cos(theta), sin(theta), 0);
    Subs::Vec3 pvec(0,0,1), tvec(sin(tilt)*sin(theta+alpha), -sin(tilt)*cos(theta+alpha), cos(tilt));

    // Length of bright spot in scale lengths
    const double BMAX = pow(mdl.expon_spot/mdl.epow_spot, 1./mdl.epow_spot);
    const double SFAC = 20. + BMAX;
    spot.resize(2*mdl.nspot);
    Lcurve::Point::etype eclipses;

    // This is where the spot height gets in
    const double AREA   = SFAC*mdl.length_spot*mdl.height_spot/(mdl.nspot-1);
    const double BRIGHT = Subs::planck(mdl.wavelength, mdl.temp_spot);

    for(int i=0; i<mdl.nspot; i++){

        // Position is adjusted to locate the impact point at the peak temperature point
        double dist = SFAC*i/(mdl.nspot-1);
        posn = bspot + mdl.length_spot*(dist-BMAX)*bvec;

        eclipses.clear();
        if(mdl.eclipse1) Lcurve::star_eclipse(mdl.q, r1, mdl.spin1, ffac1, mdl.iangle, posn, mdl.delta_phase, mdl.roche1, Roche::PRIMARY,   eclipses);
        if(mdl.eclipse2) Lcurve::star_eclipse(mdl.q, r2, mdl.spin2, ffac2, mdl.iangle, posn, mdl.delta_phase, mdl.roche2, Roche::SECONDARY, eclipses);

        // Factor here is adjusted to equal 1 at its peak
        double bright = BRIGHT*pow(dist/BMAX, mdl.expon_spot)*exp(mdl.expon_spot/mdl.epow_spot-pow(dist, mdl.epow_spot));

        // the tilted strip
        spot[i]            = Lcurve::Point(posn, tvec, AREA, 1., eclipses);
        spot[i].flux       = bright*(1-mdl.cfrac_spot)*spot[i].area;

        // the parallel strip
        spot[mdl.nspot+i]      = Lcurve::Point(posn, pvec, AREA, 1., eclipses);
        spot[mdl.nspot+i].flux = bright*mdl.cfrac_spot*spot[i].area;

    }
}

