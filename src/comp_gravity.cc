#include <cstdlib>
#include <iostream>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/vec3.h"
#include "trm/roche.h"
#include "trm/lcurve.h"

/**
 * comp_gravity2 computes a flux-weighted gravity for star2 as an approximate
 * estimate of its gravity. Returned as logg in cgs units, directly comparable
 * to usual surface gravity. Just uses fine grid.
 *
 * \param mdl -- the current model
 * \param star2 -- grid for star2
 * \return the value of logg
 */

double Lcurve::comp_gravity2(const Model& mdl,
			     const Subs::Buffer1D<Lcurve::Point>& star2){

    // Calculate the unit scaling factor to get CGS gravity
    double gm = std::pow(1000.*mdl.velocity_scale,3)*mdl.tperiod*Constants::DAY/Constants::TWOPI;
    double a = std::pow(gm/Subs::sqr(Constants::TWOPI/Constants::DAY/mdl.tperiod),1./3.);
    double gscale = 100*gm/Subs::sqr(a), gref;

    // radius stuff
    double r1, r2;
    mdl.get_r1r2(r1, r2);
    double rl2 = 1.- Roche::xl12(mdl.q, mdl.spin2);
    if(mdl.r2 < 0)
	r2 = rl2;


    // calculate the reference gravity in CGS units
    if(mdl.roche2){

	const double ACC = mdl.delta_phase/10.;
	double rref2, pref2, rad, ffac2=r2/rl2;
	Roche::ref_sphere(mdl.q, Roche::SECONDARY, mdl.spin2, ffac2, rref2, pref2);

	Subs::Vec3 dirn, posn, dvec;
        dirn.set(1.,0.,0.);
        Roche::face(mdl.q, Roche::SECONDARY, mdl.spin2, dirn, rref2, pref2, ACC, posn, dvec, rad, gref);
	gref *= gscale;

    }else{
        gref = gscale*mdl.q/(1+mdl.q)/Subs::sqr(r2);

    }

    double sumfg=0., sumf=0;

    // Star 2.
    for(int i=0; i<star2.size(); i++){
	const Point& pt = star2[i];
	// flux has built-in area factor
	sumfg += pt.flux*pt.gravity;
	sumf += pt.flux;
    }
    if((gref > 0) & (sumfg > 0) & (sumf > 0))
	return log10(gref*sumfg/sumf);
    else
	return 0.;
}
