#include <cstdlib>
#include <iostream>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/vec3.h"
#include "trm/roche.h"
#include "trm/lcurve.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

/**
 * set_star_grid sets up the elements needed to define the primary(secondary)
 * star in a co-rotating binary system. For each element set_star_grid computes
 * the position within the binary, the direction perpendicular to the element,
 * the area and gravity of the element, * (optionally) whether it is eclipsed by
 * the secondary(primary) star, and if so, what are the ingress * and egress
 * phases of the eclipse. This is optional to allow a saving of time in cases
 * where you are nor interested * in whether the elements are eclipsed.
 *
 * All numbers assume a separation = 1. The primary star is centred at (0,0,0),
 * the secondary at (1,0,0). The y-axis points in the direction of motion of the
 * secondary star.
 *
 * The grid elements are arranged in a series of constant latitude strips with
 * the 'North Pole' of the star taken as pointing towards its companion. Each
 * latitude strip has a constant latitude decrement and is split into a number
 * of equal longitude elements.  Measuring from the North pole, strip 'n' run
 * from Pi*n/nlat to Pi*(n+1)/nlat, with n running from 0 to nlat-1. At each
 * latitude a number of longitude strips equal to max(8,2*nlat*sin(theta))
 * will be selected where 'theta' is the angle from the North pole at the
 * centre of the latitude strip. This gives roughly square elements. By
 * defining the grid in terms of an integer we ensure that during optimsation,
 * the number of grid elements is fixed, although their sizes and exact
 * positions etc may vary as parameters such as mass ratio are varied.
 *
 * \param mdl        the Model with all parameters
 * \param which_star which star, primary or secondary.
 * \param fine       True for fine grid
 * \param star       array representing the primary(secondary) star (returned)
 * \exception Exceptions are thrown if the specified radii over-fill the Roche lobes.
 */

Subs::xy<double,double> envelope(double iangle, double lambda, double r1);

void Lcurve::set_star_grid(const Model& mdl, Roche::STAR which_star, bool fine,
                           Subs::Buffer1D<Lcurve::Point>& star){

    double r1, r2;
    mdl.get_r1r2(r1, r2);

    bool eclipse = which_star == Roche::PRIMARY ? mdl.eclipse1 : mdl.eclipse2;
    int nlat     = which_star == Roche::PRIMARY ? (fine ? mdl.nlat1f : mdl.nlat1c) :
        (fine ? mdl.nlat2f : mdl.nlat2c);
    int nlatfill = which_star == Roche::PRIMARY ? 0 : (fine ? mdl.nlatfill : 0);
    int nlngfill = which_star == Roche::PRIMARY ? 0 : (fine ? mdl.nlngfill : 0);
    double latlo = which_star == Roche::PRIMARY ? 0 : mdl.llo;
    double lathi = which_star == Roche::PRIMARY ? 0 : mdl.lhi;

    // Roche lobe stuff
    double rl1 = Roche::xl11(mdl.q, mdl.spin1);
    if(r1 < 0)
        r1 = rl1;
    else if(r1 > rl1)
        throw Lcurve_Error("set_star_grid: the primary star is larger than its Roche lobe!");

    double rl2 = 1.-Roche::xl12(mdl.q, mdl.spin2);
    if(r2 < 0)
        r2 = rl2;
    else if(r2 > rl2)
        throw Lcurve_Error("set_star_grid: the secondary star is larger than its Roche lobe!");

    if(mdl.glens1 && which_star == Roche::SECONDARY && mdl.roche1 && mdl.eclipse2)
        throw Lcurve_Error("set_star_grid: cannot have gravitational lensing, eclipse and Roche lobe geometry at the same time");

    // Gravitational lensing accounted for here by shrinking the lensing star
    // by a fixed amount when it is not the one having its grid computed. This
    // enables a proper computation of the eclipse of the star having its grid
    // computed by its companion. The proper correction to the radius is: r
    // --> r - 4GMd/(c^2 r) where d is the distance from the limb of the other
    // star to the point in question Here it is assumed that d = a - R/2, the
    // orbital separation minus half the radius of the star being eclipse,
    // which is approximate but should not be too bad.

    if(mdl.glens1 && which_star == Roche::SECONDARY){

        // Compute G(M1+M2), SI, and the separation a, SI.
        double gm = std::pow(1000.*mdl.velocity_scale,3)*mdl.tperiod*Constants::DAY/Constants::TWOPI;
        double a  = std::pow(gm/Subs::sqr(Constants::TWOPI/Constants::DAY/mdl.tperiod),1./3.);
        double rlens1    = 4.*gm/(1.+mdl.q)/a/Subs::sqr(Constants::C);

        r1 -= rlens1/(r1/(1.-r2/2.));
        if(r1 < 0) throw Lcurve_Error("set_star_grid: gravitational lensing correction more than the current program can cope with.");
    }

    // Calculate a reference radius and potential for the two stars
    double rref1, pref1, ffac1 = r1/rl1;
    Roche::ref_sphere(mdl.q, Roche::PRIMARY, mdl.spin1, ffac1, rref1, pref1);

    double rref2, pref2, ffac2 = r2/rl2;
    Roche::ref_sphere(mdl.q, Roche::SECONDARY, mdl.spin2, ffac2, rref2, pref2);

    // Compute latitude range over which extra points will be added. Only enabled
    // when setting the secondary grid and when the grid North pole is the genuine
    // North pole and when r2 > r1
    bool infill = mdl.npole && (which_star == Roche::SECONDARY) && (nlatfill > 0 || nlngfill > 0) && r2 > r1;
    double thelo = 0, thehi = 0;
    if(infill){

        double rangle = Subs::deg2rad(mdl.iangle);
        double cosi   = cos(rangle);
        double ratio  = (cosi+r1)/r2;

        // Lower latitude value comes from assuming star 1 is seen
        // at the limb of star 2
        if(ratio >= 1.){
            thehi = Constants::PI/2. + rangle;
        }else{

            // binary chop to discover where lower envelope of star 1 crosses edge of star 2
            double llo = 0., lhi = Constants::PI/2.;
            Subs::xy<double,double> xy;
            while(lhi > llo + 1.e-7){
                double lmid = (llo+lhi)/2.;
                xy = envelope(rangle, lmid, r1);
                if(xy.x*xy.x+xy.y*xy.y < r2*r2){
                    llo = lmid;
                }else{
                    lhi = lmid;
                }
            }
            double sini  = sin(rangle);

            // Final value, apply the lower latitude limit
            thehi = std::max(Constants::PI/2.-Subs::deg2rad(latlo),
                             std::min(Constants::PI/2. + rangle,
                                      acos(xy.y*sini/r2) + Subs::deg2rad(mdl.lfudge)));
        }

        // Upper latitude value comes from uppermost latitude covered when star 1 crosses
        // meridian of star 2.
        ratio = (cosi-r1)/r2;
        if(ratio >= 1.){
            infill = false;
        }else if(ratio <= -1.){
            thelo = 0.;
        }else{
            thelo = std::min(Constants::PI/2.-Subs::deg2rad(lathi),
                             std::max(0., Constants::PI/2. - acos(ratio) + rangle
                                      - Subs::deg2rad(mdl.lfudge)));
        }
    }

    // Compute number of faces needed
    const int NFACE = Lcurve::numface(nlat, infill, thelo, thehi, nlatfill, nlngfill);

    // Generate arrays over the star's face
    try{
        star.resize(NFACE);
    }catch(std::bad_alloc&){
        if(which_star == Roche::PRIMARY)
            throw Lcurve_Error("set_star_grid: failed to allocated enough memory for star 1");
        else
            throw Lcurve_Error("set_star_grid: failed to allocated enough memory for star 2");
    }

    const double ACC = mdl.delta_phase/10.;

    Subs::Vec3 dirn, posn, vec;

    // Compute reference gravity value, from the side of the star opposite from the L1 point
    // to ensure a non-zero value. Set to 1 if Roche distortion being ignored.
    double rad, gref;
    Subs::Vec3 dvec;

    if(which_star == Roche::PRIMARY && mdl.roche1){

        dirn.set(-1.,0.,0.);
        Roche::face(mdl.q, Roche::PRIMARY, mdl.spin1, dirn, rref1, pref1, ACC, posn, dvec, rad, gref);

    }else if(which_star == Roche::SECONDARY && mdl.roche2){

        dirn.set(1.,0.,0.);
        Roche::face(mdl.q, Roche::SECONDARY, mdl.spin2, dirn, rref2, pref2, ACC, posn, dvec, rad, gref);

    }else{

        gref = 1.;

    }

    // The grid starts at the North pole and ends at the South, proceeding in
    // a series of equi-latitudinal rings.  The polar axis is parallel to the
    // x-axis which points from one star to the other. The North pole is
    // defined to be the point closest to the other star (or the genuine North
    // Pole if npole is true). The angle theta is measured away from the North
    // pole, the angle phi is measured from the Y axis towards the Z axis.

    int nface = 0;
    double dtheta = Constants::PI/nlat;

    if(infill){
        add_faces(star, nface, 0., thelo, dtheta, 0, 0, mdl.npole, which_star, mdl.q, mdl.iangle,
                  r1, r2, rref1, rref2, mdl.roche1, mdl.roche2, mdl.spin1, mdl.spin2, eclipse, gref,
                  pref1, pref2, ffac1, ffac2, mdl.delta_phase);
        add_faces(star, nface, thelo, thehi, dtheta, nlatfill, nlngfill, mdl.npole, which_star, mdl.q,
                  mdl.iangle, r1, r2, rref1, rref2, mdl.roche1, mdl.roche2, mdl.spin1, mdl.spin2,
                  eclipse, gref, pref1, pref2, ffac1, ffac2, mdl.delta_phase);
        add_faces(star, nface, thehi, Constants::PI, dtheta, 0, 0, mdl.npole, which_star, mdl.q, mdl.iangle,
                  r1, r2, rref1, rref2, mdl.roche1, mdl.roche2, mdl.spin1, mdl.spin2, eclipse, gref,
                  pref1, pref2, ffac1, ffac2, mdl.delta_phase);
    }else{
        add_faces(star, nface, 0., Constants::PI, dtheta, 0, 0, mdl.npole, which_star, mdl.q, mdl.iangle,
                  r1, r2, rref1, rref2, mdl.roche1, mdl.roche2, mdl.spin1, mdl.spin2, eclipse, gref,
                  pref1, pref2, ffac1, ffac2, mdl.delta_phase);
    }
}

void Lcurve::add_faces(Subs::Buffer1D<Lcurve::Point>& star, int& nface, double tlo, double thi, double dtheta,
                       int nlatfill, int nlngfill, bool npole, Roche::STAR which_star, double q, double iangle,
                       double r1, double r2, double rref1, double rref2, bool roche1, bool roche2,
                       double spin1, double spin2, bool eclipse, double gref, double pref1, double pref2,
                       double ffac1, double ffac2, double delta){

    double ri = Subs::deg2rad(iangle);
    double cosi = cos(ri), sini = sin(ri);

    // Centre of masses of the stars
    const Subs::Vec3 cofm1(0.,0.,0.), cofm2(1.,0.,0.);

    // Can afford to be pretty careful on the location of faces as it is a fast computation
    const double ACC = delta/10.;

    bool infill = (nlatfill > 0) || (nlngfill > 0);
    int nlat, nlat1 = 0, nlat2 = 0;
    if(infill){
        // If infill is True we loop through the latitudes once on the side
        // facing the other star and once on the other side so that infilling
        // only occurs on one side. The infilled part will be included first
        nlat1 = int(std::ceil((1+nlatfill)*(thi-tlo)/dtheta));
        nlat2 = int(std::ceil((thi-tlo)/dtheta));
        nlat  = nlat1 + nlat2;
    }else{
        nlat = int(std::ceil((thi-tlo)/dtheta));
    }

    // vector of offsets needed for parallelisation as opposed
    // to counting up on the fly.
    std::vector<int> off(nlat);

    double sint;
    for(int nt=0; nt<nlat; nt++){
        off[nt] = nface;
        if(infill){
            if(nt < nlat1){
                sint   = sin(tlo + (thi-tlo)*(nt+0.5)/nlat1);
                nface += std::max(8,int(Constants::PI*sint*(1+nlngfill)/dtheta));
            }else{
                sint   = sin(tlo + (thi-tlo)*(nt-nlat1+0.5)/nlat2);
                nface += std::max(8,int(Constants::PI*sint/dtheta));
            }
        }else{
            sint   = sin(tlo + (thi-tlo)*(nt+0.5)/nlat);
            nface += std::max(16,int(2*Constants::PI*sint/dtheta));
        }
    }

    bool failed = false;
    std::string error_message;

#ifdef _OPENMP
    int mxth = std::min(16, omp_get_max_threads());
    omp_set_num_threads(mxth);
#pragma omp parallel for schedule(dynamic)
#endif

    for(int nt=0; nt<nlat; nt++){

        double phi1, phi2, theta, sint, cost;
        int nphi, nl;
        if(infill){
            if(nt < nlat1){
                theta = tlo + (thi-tlo)*(nt+0.5)/nlat1;
                sint  = sin(theta);
                cost  = cos(theta);
                nphi  = std::max(8,int(Constants::PI*sint*(1+nlngfill)/dtheta));
                nl    = nlat1;
                if(which_star == Roche::PRIMARY){
                    phi1 = -Constants::PI/2.;
                    phi2 = +Constants::PI/2.;
                }else{
                    phi1 = +Constants::PI/2.;
                    phi2 = +3.*Constants::PI/2.;
                }
            }else{
                theta = tlo + (thi-tlo)*(nt-nlat1+0.5)/nlat2;
                sint  = sin(theta);
                cost  = cos(theta);
                nphi = std::max(8,int(Constants::PI*sint/dtheta));
                nl   = nlat2;
                if(which_star == Roche::PRIMARY){
                    phi1 = +Constants::PI/2.;
                    phi2 = +3.*Constants::PI/2.;
                }else{
                    phi1 = -Constants::PI/2.;
                    phi2 = +Constants::PI/2.;
                }
            }
        }else{
            theta = tlo + (thi-tlo)*(nt+0.5)/nlat;
            sint  = sin(theta);
            cost  = cos(theta);
            nphi = std::max(16,int(2*Constants::PI*sint/dtheta));
            nl   = nlat;
            phi1 = 0.;
            phi2 = Constants::TWOPI;
        }

        for(int np=0; np<nphi; np++){

            try{
                double phi  = phi1 + (phi2-phi1)*(np+0.5)/nphi;
                double sinp = sin(phi);
                double cosp = cos(phi);

                Subs::Vec3 dirn, posn, dvec;
                if(npole){
                    dirn.set(sint*cosp, sint*sinp, cost);
                }else{
                    dirn.set(cost, sint*cosp, sint*sinp);
                }

                // Direction is now defined, so calculate radius and thus the
                // position according to whether we are accounting for Roche
                // geometry or not.
                double rad, gravity, lam1, lam2, area, ingress, egress;

                if(which_star == Roche::PRIMARY && roche1){
                    Roche::face(q, Roche::PRIMARY, spin1, dirn, rref1, pref1, ACC, posn, dvec, rad, gravity);
                }else if(which_star == Roche::SECONDARY && roche2){
                    Roche::face(q, Roche::SECONDARY, spin2, dirn, rref2, pref2, ACC, posn, dvec, rad, gravity);
                }else{

                    // Ignore Roche distortion
                    if(which_star == Roche::PRIMARY){
                        rad  = r1;
                        posn = cofm1 + rad*dirn;
                    }else{
                        rad  = r2;
                        posn = cofm2 + rad*dirn;
                    }
                    dvec    = dirn;
                    gravity = 1.;

                }

                // Area, accounting for angle of face
                area = ((phi2-phi1)/nphi*rad*sint)*((thi-tlo)/nl*rad)/dot(dirn, dvec);

                // Eclipse computation. We calculate whether a point is
                // eclipsed, and, if it is, its ingress and egress
                // phases. Account for spherical or Roche geometry of other
                // star. Since the stars are convex this calculation only
                // accounts for eclipse by the OTHER star.
                Lcurve::Point::etype eclipses;

                if(eclipse &&
                   ((which_star == Roche::PRIMARY &&
                     ((roche2  && Roche::ingress_egress(q, Roche::SECONDARY, spin2, ffac2, iangle, delta, posn, ingress, egress)) ||
                      (!roche2 && Roche::sphere_eclipse(cosi, sini, posn, cofm2, r2, ingress, egress, lam1, lam2)))) ||
                    (which_star == Roche::SECONDARY &&
                     ((roche1  && Roche::ingress_egress(q, Roche::PRIMARY, spin1, ffac1, iangle, delta, posn, ingress, egress)) ||
                      (!roche1 && Roche::sphere_eclipse(cosi, sini, posn, cofm1, r1, ingress, egress, lam1, lam2)))))){

                    eclipses.push_back(std::make_pair(ingress, egress));

                }

                star[off[nt]+np] = Lcurve::Point(posn, dvec, area, gravity/gref, eclipses);
            }
            catch(const Roche::Roche_Error& err){
                if(!failed){
                    failed = true;
                    error_message = std::string(err.what());
                }
            }
        }
    }
    // end of parallel loop
    if(failed)
        throw Lcurve::Lcurve_Error(error_message);
}

//! Returns x,y on lower envelope swept out by star 1 transiting star 2
/* See the document geometry to understand what is going on here
 */
Subs::xy<double,double> envelope(double rangle, double lambda, double r1){
    double sini = sin(rangle);
    double cosi = cos(rangle);
    double sinl = sin(lambda);
    double cosl = cos(lambda);
    double norm = std::sqrt(cosi*cosi+sini*sini*cosl*cosl);
    Subs::xy<double,double> temp;
    temp.x = sinl + r1*cosi*sinl/norm;
    temp.y = -cosi*cosl - r1*cosl/norm;
    return temp;
}
