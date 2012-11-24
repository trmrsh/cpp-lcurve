#include <cstdlib>
#include <iostream>
#include "trm_subs.h"
#include "trm_constants.h"
#include "trm_vec3.h"
#include "trm_roche.h"
#include "trm_lcurve.h"

/** 
 * comp_light computes a light curve point for a particular phase. It can
 * allow for finite exposures by trapezoidal integration.
 * \param iangle   orbital inclination
 * \param ldc1 limb darkening for star 1
 * \param ldc2 limb darkening for star 2
 * \param lin_limb_disc1  linear limb darkening for disc
 * \param quad_limb_disc1  quadratic limb darkening for disc
 * \param phase    orbital phase at centre of exposure
 * \param expose   length of exposure in terms of phase
 * \param ndiv     number of sub-divisions for exposure smearing using trapezoidal integration
 * \param q         mass ratio
 * \param beam_factor1 the 3-alpha factor for Doppler beaming for star1
 * \param beam_factor2 the 3-alpha factor for Doppler beaming for star2
 * \param spin1    spin/orbital ratio star 1 (makes a difference to the beaming)
 * \param spin2    spin/orbital ratio star 2 (makes a difference to the beaming)
 * \param vscale   the velocity scale V1+V2 (unprojected) for computation of Doppler beaming
 * \param glens1   account for gravitational lensing by star 1
 * \param rlens1   4x the gravitational radius (GM/c^2) of star 1, scaled by separation
 * \param gint     set of numbers for switching between coarse & fine grids
 * \param star1f   geometry/brightness array for star1, fine grid
 * \param star2f   geometry/brightness array for star2, fine grid
 * \param star1c   geometry/brightness array for star1, coarse grid
 * \param star2c   geometry/brightness array for star2, coarse grid
 * \param disc     geometry/brightness array for the disc
 * \param spot     geometry/brightness array for the bright spot
 * \return the light curve value desired.
 */

double Lcurve::comp_light(double iangle, const LDC& ldc1, const LDC& ldc2, double lin_limb_disc, double quad_limb_disc, 
			  double phase, double expose, int ndiv, double q, double beam_factor1, double beam_factor2,
			  double spin1, double spin2, float vscale, bool glens1, double rlens1, const Lcurve::Ginterp& gint,
			  const Subs::Buffer1D<Lcurve::Point>& star1f, const Subs::Buffer1D<Lcurve::Point>& star2f, 
			  const Subs::Buffer1D<Lcurve::Point>& star1c, const Subs::Buffer1D<Lcurve::Point>& star2c, 
			  const Subs::Buffer1D<Lcurve::Point>& disc, const Subs::Buffer1D<Lcurve::Point>& spot){

    const double XCOFM = q/(1.+q);
    const double cosi  = cos(Constants::TWOPI*iangle/360.);
    const double sini  = sin(Constants::TWOPI*iangle/360.);
    const double VFAC  = vscale/(Constants::C/1.e3);

    Subs::Vec3 earth, s;
    int ptype, i, nelem1, nelem2, nd;
    double sum=0., ssum, ssum2, mu, wgt, vr, phi, p, ph, mag, pd, d, phsq, rd;
    double vx, vy, vn, mud, ommu;

    for(nd=0; nd<ndiv; nd++){
      
        if(ndiv == 1){
	    phi = phase;
	    wgt = 1.;
	}else{
	    phi = phase + expose*(nd-double(ndiv-1)/2.)/(ndiv-1);
	    if(nd == 0 || nd == ndiv-1)
		wgt = 0.5;
	    else
		wgt = 1.;
	}
	
	earth = Roche::set_earth(cosi, sini, phi);
	
	ptype = gint.type(phi);
	const Subs::Buffer1D<Lcurve::Point>& star1 = ptype == 1 ? star1f : star1c;	
	const Subs::Buffer1D<Lcurve::Point>& star2 = ptype == 3 ? star2f : star2c;	
	nelem1 = star1.size();
	nelem2 = star2.size();

	ssum = 0.;

	// Star 1.
	for(i=0; i<nelem1; i++){
	    const Point& pt = star1[i];
	    if(pt.visible(phi)){
		mu = Subs::dot(earth, pt.dirn);
		if(ldc1.see(mu)){
		    if(beam_factor1 != 0.){
			vx  = -VFAC*spin1*pt.posn.y();
			vy  =  VFAC*(spin1*pt.posn.x()-XCOFM);
			vr  = -(earth.x()*vx + earth.y()*vy);
			vn  = pt.dirn.x()*vx + pt.dirn.y()*vy;
			mud = mu - mu*vr - vn;
			ssum += mu*pt.flux*(1.-beam_factor1*vr)*ldc1.imu(mud);
		    }else{
			ssum += mu*pt.flux*ldc1.imu(mu);
		    }
		}
	    }
	}
	ssum *= gint.scale1(phi);

	// Star 2.
	ssum2 = 0.;

	for(i=0; i<nelem2; i++){
	    const Point& pt = star2[i];
	    if(pt.visible(phi)){
		mu = Subs::dot(earth, pt.dirn);	    

		if(ldc2.see(mu)){
	      
		    // Account for magnifying effect of gravitational lensing here
		    mag = 1.;
		    if(glens1){
			// s = vector from centre of mass of star 1 to point of interest
			s = pt.posn;
		
			// d = distance along line of sight from star 1 to point in question 
			d = -Subs::dot(s, earth);
			if(d > 0.){
			    // p  = distance in plane of sky from cofm of star 1 to point in question
			    // pd = larger distance accounting for deflection by lensing. For p >> rlens1*d, q --> p,
			    // mag --> 1. Try to save time by avoiding square root if possible.
			    ph   = (p = (s+d*earth).length())/2.;
			    phsq = ph*ph;
			    rd   = rlens1*d;
			    if(phsq > 25.*rd){
				pd = p + rd/p;
			    }else{
				pd = ph + std::sqrt(phsq + rd);
			    }
			    mag = pd*pd/(pd-ph)/ph/4.;
			}
		    }
	      
		    if(beam_factor2 != 0.){
			vx  = -VFAC*spin2*pt.posn.y();
			vy  =  VFAC*(spin2*(pt.posn.x()-1.)+1.-XCOFM);
			vr  = -(earth.x()*vx + earth.y()*vy);
			vn  = pt.dirn.x()*vx + pt.dirn.y()*vy;
			mud = mu - mu*vr - vn;
			ssum2 += mu*mag*pt.flux*(1-beam_factor2*vr)*ldc2.imu(mud);
		    }else{
			ssum2 += mu*mag*pt.flux*ldc2.imu(mu);
		    }
		}
	    }
	}
    
	ssum += gint.scale2(phi)*ssum2;
	
	// Disc
	for(i=0; i<disc.size(); i++){
	    mu = Subs::dot(earth, disc[i].dirn);      
	    if(mu > 0. && disc[i].visible(phi)){
		ommu  = 1.-mu;
		ssum += mu*disc[i].flux*(1.-ommu*(lin_limb_disc+quad_limb_disc*ommu));
	    }
	}
	
	// Spot
	for(i=0; i<spot.size(); i++){
	    mu = Subs::dot(earth, spot[i].dirn);      
	    if(mu > 0. && spot[i].visible(phi))
		ssum += mu*spot[i].flux;
	}
	
	sum = sum + wgt*ssum;
    }
    
    return sum/std::max(1, ndiv-1);
    
}


/** 
 * comp_star1 computes the flux from star 1 for a particular phase. It can
 * allow for finite exposures by trapezoidal integration.
 * \param iangle   orbital inclination
 * \param ldc1 limb darkening for star 1
 * \param phase    orbital phase at centre of exposure
 * \param expose   length of exposure in terms of phase
 * \param ndiv     number of sub-divisions for exposure smearing using trapezoidal integration
 * \param q         mass ratio
 * \param beam_factor 3-alpha factor for Doppler beaming
 * \param vscale   the velocity scale V1+V2 (unprojected) for computation of Doppler beaming
 * \param gint     contains grid scaling factors
 * \param star1f   fine grid for star1
 * \param star1c   coarse grid for star1
 * \return the light curve value desired.
 */

double Lcurve::comp_star1(double iangle, const LDC& ldc1, double phase, double expose, int ndiv, 
			  double q, double beam_factor1, float vscale, const Lcurve::Ginterp& gint,
			  const Subs::Buffer1D<Lcurve::Point>& star1f, const Subs::Buffer1D<Lcurve::Point>& star1c){ 
  
    const double XCOFM = q/(1.+q);
    double ri = Subs::deg2rad(iangle);
    const double cosi  = cos(ri);
    const double sini  = sin(ri);
    const double VFAC  = vscale/(Constants::C/1.e3);  
    
    Subs::Vec3 earth;
    double phi, sum=0., ssum, mu, wgt, vr;
    double vx, vy, vn, mud;

    for(int nd=0; nd<ndiv; nd++){
	
	if(ndiv == 1){
	    phi = phase;
	    wgt = 1.;
	}else{
	    phi = phase + expose*(nd-double(ndiv-1)/2.)/(ndiv-1);
	    if(nd == 0 || nd == ndiv-1)
		wgt = 0.5;
	    else
		wgt = 1.;
	}

	earth = Roche::set_earth(cosi, sini, phi);

	// Define the grid to use
	const Subs::Buffer1D<Lcurve::Point>& star1 = gint.type(phi) == 1 ? star1f : star1c;	
	int nelem1 = star1.size();
	
	ssum = 0.;
	// Star 1.

	for(int i=0; i<nelem1; i++){
	    const Point& pt = star1[i];
	    if(pt.visible(phi)){
		mu = Subs::dot(earth, pt.dirn);
		if(ldc1.see(mu)){
		    if(beam_factor1 != 0.){
			vx  = -VFAC*pt.posn.y();
			vy  =  VFAC*(pt.posn.x()-XCOFM);
			vr  = -(earth.x()*vx + earth.y()*vy);
			vn  = pt.dirn.x()*vx + pt.dirn.y()*vy;
			mud = mu - mu*vr - vn;
			ssum += mu*pt.flux*(1.-beam_factor1*vr)*ldc1.imu(mud);
		    }else{
			vr  = VFAC*(earth.x()*pt.posn.y() - earth.y()*(pt.posn.x()-XCOFM));
			ssum += mu*pt.flux*(1.-beam_factor1*vr)*ldc1.imu(mu);
		    }
		}
	    }
	}
	
	sum += wgt*gint.scale1(phase)*ssum;
    }
    
    return sum/std::max(1, ndiv-1);
    
}


/** 
 * comp_star2 computes the flux from star 2 only. It can
 * allow for finite exposures by trapezoidal integration.
 * \param iangle   orbital inclination
 * \param ldc2   limb darkening for star 2
 * \param phase    orbital phase at centre of exposure
 * \param expose   length of exposure in terms of phase
 * \param ndiv     number of sub-divisions for exposure smearing using trapezoidal integration
 * \param q         mass ratio
 * \param beam_factor2 3-alpha Doppler beaming factor for star 2
 * \param vscale   the velocity scale V1+V2 (unprojected) for computation of Doppler beaming
 * \param glens1   account for gravitational lensing by star 1
 * \param rlens1   4x the gravitational radius (GM/c^2) of star 1, scaled by separation
 * \param gint     factors for switching grids
 * \param star2f   fine grid for star2
 * \param star2c   coarse grid for star2
 * \return the light curve value desired.
 */

double Lcurve::comp_star2(double iangle, const LDC& ldc2, double phase, double expose, int ndiv, 
			  double q, double beam_factor2, float vscale, bool glens1, double rlens1, const Lcurve::Ginterp& gint,
			  const Subs::Buffer1D<Lcurve::Point>& star2f, const Subs::Buffer1D<Lcurve::Point>& star2c){ 
    
    const double XCOFM = q/(1.+q);
    double ri = Subs::deg2rad(iangle);
    const double cosi  = cos(ri);
    const double sini  = sin(ri);
    const double VFAC  = vscale/(Constants::C/1.e3);  
    
    Subs::Vec3 earth, s;
    double phi, sum=0., ssum, mu, wgt, vr;
    double p, ph, mag, pd, d, phsq, rd, vx, vy, vn, mud;

    for(int nd=0; nd<ndiv; nd++){
	
	if(ndiv == 1){
	    phi = phase;
	    wgt = 1.;
	}else{
	    phi = phase + expose*(nd-double(ndiv-1)/2.)/(ndiv-1);
	    if(nd == 0 || nd == ndiv-1)
		wgt = 0.5;
	    else
		wgt = 1.;
	}

	earth = Roche::set_earth(cosi, sini, phi);
	
	// Define the grid to use
	const Subs::Buffer1D<Lcurve::Point>& star2 = gint.type(phi) == 3 ? star2f : star2c;	
	int nelem2 = star2.size();
	
	ssum = 0.;
	
	// Star 2.
	for(int i=0; i<nelem2; i++){
	    const Point& pt = star2[i];
	    if(pt.visible(phi)){
		mu = Subs::dot(earth, pt.dirn);	    
		if(ldc2.see(mu)){
		    
		    // Account for magnifying effect of gravitational lensing here
		    mag = 1.;
		    if(glens1){
			// s = vector from centre of mass of star 1 to point of interest
			s = pt.posn;
			
			// d = distance along line of sight from star 1 to point in question 
			d = -Subs::dot(s, earth);
			if(d > 0.){
			    // p  = distance in plane of sky from cofm of star 1 to point in question
			    // pd = larger distance accounting for lensing deflection. For p >> rlens1*d, pd --> p,
			    // mag --> 1. Try to save time by avoiding square root if possible.
			    ph   = (p = (s+d*earth).length())/2.;
			    phsq = ph*ph;
			    rd   = rlens1*d;
			    if(phsq > 25.*rd){
				pd = p + rd/p;
			    }else{
				pd = ph + std::sqrt(phsq + rd);
			    }
			    mag = pd*pd/(pd-ph)/ph/4.;
			}
		    }

		    if(beam_factor2 != 0.){
			vx  = -VFAC*pt.posn.y();
			vy  =  VFAC*(pt.posn.x()-XCOFM);
			vr  = -(earth.x()*vx + earth.y()*vy);
			vn  = pt.dirn.x()*vx + pt.dirn.y()*vy;
			mud = mu - mu*vr - vn;
			ssum += mu*mag*pt.flux*(1-beam_factor2*vr)*ldc2.imu(mud);
		    }else{
			vr  = VFAC*(earth.x()*pt.posn.y() - earth.y()*(pt.posn.x()-XCOFM));
			ssum += mu*mag*pt.flux*(1-beam_factor2*vr)*ldc2.imu(mu);
		    }		    
		}
	    }
	}
	
	sum += wgt*gint.scale2(phi)*ssum;
    }
    
    return sum/std::max(1, ndiv-1);
    
}

/** 
 * comp_disc computes the flux from the disc at a particular phase. It can
 * allow for finite exposures by trapezoidal integration.
 * \param iangle   orbital inclination
 * \param lin_limb_disc linear limb darkening coefficient for the disc
 * \param quad_limb_disc quadratic limb darkening for the disc
 * \param phase    orbital phase at centre of exposure
 * \param expose   length of exposure in terms of phase
 * \param ndiv     number of sub-divisions for exposure smearing using trapezoidal integration
 * \param q         mass ratio
 * \param vscale   the velocity scale V1+V2 (unprojected) for computation of Doppler beaming
 * \param disc     geometry/brightness array for the disc
 * \return the light curve value desired.
 */

double Lcurve::comp_disc(double iangle, double lin_limb_disc, double quad_limb_disc, double phase, double expose, 
			 int ndiv, double q, float vscale, const Subs::Buffer1D<Lcurve::Point>& disc){
    
    const Subs::Vec3 COFM(q/(1.+q),0.,0.), SPIN(0.,0.,1.);
    double ri = Subs::deg2rad(iangle);
    const double cosi  = cos(ri);
    const double sini  = sin(ri);
    
    Subs::Vec3 earth;
    double phi, sum=0., ssum, mu, wgt;
    for(int nd=0; nd<ndiv; nd++){
	
	if(ndiv == 1){
	    phi = phase;
	    wgt = 1.;
	}else{
	    phi = phase + expose*(nd-double(ndiv-1)/2.)/(ndiv-1);
	    if(nd == 0 || nd == ndiv-1)
		wgt = 0.5;
	    else
		wgt = 1.;
	}

	earth = Roche::set_earth(cosi, sini, phi);
	
	ssum = 0.;
	// Disc
	for(int i=0; i<disc.size(); i++){
	    mu = Subs::dot(earth, disc[i].dirn);      
	    if(mu > 0. && disc[i].visible(phi))
		ssum += mu*disc[i].flux*(1.-(1.-mu)*(lin_limb_disc+quad_limb_disc*(1.-mu)));
	}
	
	sum += wgt*ssum;
    }
    
    return sum/std::max(1, ndiv-1);
    
}

/** 
 * comp_spot computes a light curve point for a particular phase. It can
 * allow for finite exposures by trapezoidal integration.
 * \param iangle   orbital inclination
 * \param phase    orbital phase at centre of exposure
 * \param expose   length of exposure in terms of phase
 * \param ndiv     number of sub-divisions for exposure smearing using trapezoidal integration
 * \param q         mass ratio
 * \param vscale   the velocity scale V1+V2 (unprojected) for computation of Doppler beaming
 * \param spot     geometry/brightness array for the bright spot
 * \return the light curve value desired.
 */

double Lcurve::comp_spot(double iangle,double phase, double expose, int ndiv, double q, float vscale,
			 const Subs::Buffer1D<Lcurve::Point>& spot){
    
    const Subs::Vec3 COFM(q/(1.+q),0.,0.), SPIN(0.,0.,1.);
    double ri = Subs::deg2rad(iangle);
    const double cosi  = cos(ri);
    const double sini  = sin(ri);

    Subs::Vec3 earth;
    double phi, sum=0., ssum, mu, wgt;
    for(int nd=0; nd<ndiv; nd++){
	
	if(ndiv == 1){
	    phi = phase;
	    wgt = 1.;
	}else{
	    phi = phase + expose*(nd-double(ndiv-1)/2.)/(ndiv-1);
	    if(nd == 0 || nd == ndiv-1)
		wgt = 0.5;
	    else
		wgt = 1.;
	}
	earth = Roche::set_earth(cosi, sini, phi);
	
	ssum = 0.;
	// Spot
	for(int i=0; i<spot.size(); i++){
	    mu = Subs::dot(earth, spot[i].dirn);      
	    if(mu > 0. && spot[i].visible(phi))
		ssum += mu*spot[i].flux;
	}
	
	sum += wgt*ssum;
    }
    
    return sum/std::max(1, ndiv-1);
    
}
