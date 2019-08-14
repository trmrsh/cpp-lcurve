#include <cstdlib>
#include <iostream>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/vec3.h"
#include "trm/roche.h"
#include "trm/lcurve.h"

/**
 * comp_radius1 computes a volume-averaged scaled radius for star 1.
 *
 * \param star1 -- grid for star1
 * \return the value of volume-averaged R1/a
 */

double Lcurve::comp_radius1(const Subs::Buffer1D<Lcurve::Point>& star1){

    Subs::Vec3 vec, cofm1(0.,0.,0.);
    double sumsa=0., sumvol=0., r, rcosa;

    // Star 1.
    for(int i=0; i<star1.size(); i++){
        const Point& pt = star1[i];
        vec = pt.posn-cofm1;
        r = vec.length();
        rcosa = Subs::dot(pt.dirn, vec);

        // sum solid angle and 3x volume of all elements
        sumsa += pt.area*rcosa/std::pow(r,3);
        sumvol += pt.area*rcosa;
    }
    return std::pow(sumvol/sumsa,1./3.);
}

/**
 * comp_radius2 computes a volume-averaged scaled radius for star 2.
 *
 * \param star2 -- grid for star2
 * \return the value of volume-averaged R2/a
 */

double Lcurve::comp_radius2(const Subs::Buffer1D<Lcurve::Point>& star2){

    Subs::Vec3 vec, cofm2(1.,0.,0.);
    double sumsa=0., sumvol=0., r, rcosa;

    // Star 2.
    for(int i=0; i<star2.size(); i++){
        const Point& pt = star2[i];
        vec = pt.posn-cofm2;
        r = vec.length();
        rcosa = Subs::dot(pt.dirn, vec);

        // sum solid angle and 3x volume of all elements
        sumsa += pt.area*rcosa/std::pow(r,3);
        sumvol += pt.area*rcosa;
    }
    return std::pow(sumvol/sumsa,1./3.);
}
