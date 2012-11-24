/*

!!begin
!!title    Draws a Roche lobe, cvmovie style
!!author   T.R.Marsh 
!!created  20 October 2004
!!descr    Draws a Roche lobe, cvmovie style
!!index    picture
!!root     picture
!!css      style.css
!!class    Model
!!head1    picture -- draws a Roche lobe, cvmovie style

!!emph{picture} draws a Roche lobe, cvmovie style.

!!head2 Command invocation

picture q iangle phase r1 r2 rd roche delta ntheta1 nphi1 ntheta2 nphi2 rdisc1 rdisc2 height beta 
ndisct ndiscr nedge nplot spot device x1 x2 y1 y2 width reverse

!!head2 Arguments

!!table
!!arg{q}{Mass ratio, q = M2/M1}
!!arg{iangle}{Inclination angle, degrees}
!!arg{phase}{The orbital phase to compute}
!!arg{r1}{Radius of star 1, scaled by the binary separation}
!!arg{r2}{Radius of star 2, scaled by the binary separation. The radius is measured along the line of 
centres towards star 1. If larger than the L1 radius, it will be set equal to it.}
!!arg{rd}{Radius of the disc}
!!arg{roche}{true/false to account for Roche distortion of secondary or not}
!!arg{delta}{Accuracy in Roche lobe computations in terms of phase}
!!arg{ntheta1}{Number of lines of equal latitude on star 1}
!!arg{nphi1}{Number of lines of equal longitude on star 1}
!!arg{ntheta2}{Number of lines of equal latitude on star 2}
!!arg{nphi2}{Number of lines of equal longitude on star 2}
!!arg{rdisc1}{Inner radius of disc}
!!arg{rdisc2}{Outer radius of disc}
!!arg{height}{Height of disc at r=1}
!!arg{beta}{Flare exponent of disc}
!!arg{ndisct}{Number of lines of theta on disc}
!!arg{ndiscr}{Number of lines of radius on disc}
!!arg{nedge}{Number of vertical lines at edge of disc}
!!arg{nplot}{Number of points around equator of red star}
!!arg{spot}{Bright-spot size (0 to suppress)}
!!arg{device}{Plot device}
!!arg{x1}{left-hand X limit}
!!arg{x2}{right-hand X limit}
!!arg{y1}{lower Y limit}
!!arg{y2}{upper Y limit}
!!arg{width}{plot width in inches}
!!arg{reverse}{Reverse the asignment of black and white colour indices or not}
!!table

!!end

*/

#include <cstdlib>
#include <iostream>
#include "cpgplot.h"
#include "trm_subs.h"
#include "trm_plot.h"
#include "trm_vec3.h"
#include "trm_input.h"
#include "trm_roche.h"
#include "trm_lcurve.h"

// Main program
int main(int argc, char* argv[]){
  
  try{

    // Construct Input object
    Subs::Input input(argc, argv, Lcurve::LCURVE_ENV, Lcurve::LCURVE_DIR);

    // sign-in input variables
    input.sign_in("q",        Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("iangle",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("phase",    Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("r1",       Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("r2",       Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("rd",       Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("roche",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("delta",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("ntheta1",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nphi1",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("ntheta2",  Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nphi2",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("rdisc1",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("rdisc2",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("height",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("beta",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("ndisct",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("ndiscr",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nedge",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nplot",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("spot",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("device",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("x1",       Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("x2",       Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("y1",       Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("y2",       Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("width",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("reverse",  Subs::Input::LOCAL,  Subs::Input::PROMPT);

    double q;
    input.get_value("q", q, 0.2, 0.0001, 100., "mass ratio q (=M2/M1)");
    double iangle;
    input.get_value("iangle", iangle, 85., 0., 90., "orbital inclination angle (degrees)");
    double phase;
    input.get_value("phase", phase, 0., 0., 1., "orbital phase");
    double r1;
    input.get_value("r1", r1, 0.1, 1.e-8, 0.9999999999, "radius of star 1 (unit of binary separation)");
    double rl2 = 1.-Roche::xl1(q);
    double r2;
    input.get_value("r2", r2, std::min(0.1,rl2), -10., rl2, "radius of star 2 (unit of binary separation)");
    if(r2 <= 0) r2 = 0.99999999*rl2;
    double ffac = r2/rl2, rref, pref;
    double spin2 = 1.0;
    Roche::ref_sphere(q, Roche::SECONDARY, spin2, ffac, rref, pref);
    bool roche;
    input.get_value("roche", roche, true, "account for Roche distortion of secondary or not?");
    double delta;
    input.get_value("delta", delta, 1.e-6, 1.e-12, 1.e-2, "accuracy in phase for Roche lobe computations");
    int ntheta1;
    input.get_value("ntheta1", ntheta1, 10, 1, 1000, "number of latitude lines for star 1");
    int nphi1;
    input.get_value("nphi1", nphi1, 10, 1, 1000, "number of longitude lines for star 1");
    int ntheta2;
    input.get_value("ntheta2", ntheta2, 10, 1, 1000, "number of latitude lines for star 2");
    int nphi2;
    input.get_value("nphi2", nphi2, 10, 1, 1000, "number of longitude lines for star 2");

    double rdisc1;
    input.get_value("rdisc1", rdisc1, 0.1, r1, 1., "inner radius of disc (unit of binary separation)");
    double rdisc2;
    input.get_value("rdisc2", rdisc2, 0.4, r1, 1., "outer radius of disc (unit of binary separation)");
    double height;    input.get_value("height", height, 0.1, 0., 100., "height of disc at unit radius (units of binary separation)");
    double beta;
    input.get_value("beta", beta, 3., 0., 100., "flaring exponent of disc");
    int ndisct;
    input.get_value("ndisct", ndisct, 10, 1, 1000, "number of lines of theta on disc");
    int ndiscr;
    input.get_value("ndiscr", ndiscr, 10, 1, 1000, "number of lines of radius on disc");
    int nedge;
    input.get_value("nedge", nedge, 20, 1, 1000, "number of vertical lines at edge of disc");
    float spot;
    input.get_value("spot", spot, 1.f, 0.f, 20.f, "spot size (multiple of PGPLOT default)");
    int nplot;
    input.get_value("nplot", nplot, 100, 2, 100000, "number of points to plot around equator of red star");
    std::string device;
    input.get_value("device", device, "/xs", "plot device");
    float x1;
    input.get_value("x1", x1, -2.f, -100.f, 100.f, "left-hand X limit");
    float x2;
    input.get_value("x2", x2,  2.f, -100.f, 100.f, "right-hand X limit");
    float y1;
    input.get_value("y1", y1, -2.f, -100.f, 100.f, "lower Y limit");
    float y2;
    input.get_value("y2", y2,  2.f, -100.f, 100.f, "upper Y limit");
    float width;
    input.get_value("width", width,  4.f, 0.1f, 1000.f, "plot width, inches");
    bool reverse;
    input.get_value("reverse", reverse, false, "reverse black & white colour indices?");
    input.save();

    // Set standard vectors
    Subs::Vec3 dirn, rvec, dvec, earth;

    const Subs::Vec3 cofm1(0.,0.,0.), cofm2(1.,0.,0.), cofm(q/(1.+q),0.,0.);

    earth = Roche::set_earth(iangle, phase);

    // Compute sky basis vectors
    double cosp = cos(Constants::TWOPI*phase);
    double sinp = sin(Constants::TWOPI*phase);
    Subs::Vec3 xsky(sinp,cosp,0.);
    Subs::Vec3 ysky = Subs::cross(earth, xsky);

    // Plot
    Subs::Plot plot(device);

    float xs1, xs2, ys1, ys2;
    cpgqvsz(2,&xs1,&xs2,&ys1,&ys2);

    if(reverse){
      float bred, bgreen, bblue, fred, fgreen, fblue;
      cpgqcr(0, &bred, &bgreen, &bblue);
      cpgqcr(1, &fred, &fgreen, &fblue);
      cpgscr(1, bred, bgreen, bblue);
      cpgscr(0, fred, fgreen, fblue);
    }
    cpgpap(width,(y2-y1)/(x2-x1));
    cpgqvsz(2,&xs1,&xs2,&ys1,&ys2);
    cpgsvp(0.,1.,0.,1.);
    cpgwnad(x1, x2, y1, y2);
    cpgqvsz(2,&xs1,&xs2,&ys1,&ys2);
    cpgslw(2);

    double theta, sint, cost, phi, lam1, lam2, gravity, rad, gref;
    bool ready = false;
    const double acc = Constants::TWOPI*delta/10.;

    if(roche){
      // Compute reference gravity value.
      dirn.unitx();
      Roche::face(q, Roche::SECONDARY, spin2, dirn, rref, pref, acc, rvec, dvec, rad, gref);
    }else{
      gref = 1.;
    }

    cpgsci(2);

    // Loop over theta over star 2
    for(int nt=0; nt<ntheta2; nt++){

      cpgbbuf();

      theta = Constants::PI*(nt+0.5)/ntheta2;
      sint  = sin(theta);
      cost  = cos(theta);

      int nphi  = std::max(8,int(nplot*sint+0.5));

      for(int np=0; np<nphi; np++){

	phi  = Constants::TWOPI*np/(nphi-1);
	sinp = sin(phi);
	cosp = cos(phi);

	dirn.set(cost, sint*cosp, sint*sinp);

	if(roche){
	
	    Roche::face(q, Roche::SECONDARY, spin2, dirn, rref, pref, acc, rvec, dvec, rad, gravity);
	  
	}else{
	  
	  // Assume spherical
	  rvec    = r2*dirn + cofm2;
	  dvec    = dirn;
	  
	}
	
	// Eclipse & visibility computation, star1 assumed spherical.
	if(Subs::dot(earth, dvec) <= 0. || Roche::sphere_eclipse(earth, rvec, cofm1, r1, lam1, lam2) ||
	   Lcurve::disc_eclipse(iangle, phase, rdisc1, rdisc2, beta, height, rvec)){
	  ready = false;
	}else{
	  if(ready)
	    cpgdraw(Subs::dot(rvec, xsky), Subs::dot(rvec, ysky));
	  else
	    cpgmove(Subs::dot(rvec, xsky), Subs::dot(rvec, ysky));
	  ready = true;
	}
      }

      cpgebuf();

      ready = false;

    }

    // Loop over phi over star 2
    for(int np=0; np<nphi2; np++){

      cpgbbuf();

      phi   = Constants::PI*(np+0.5)/nphi2;
      sinp  = sin(phi);
      cosp  = cos(phi);

      int ntheta  = nplot;

      for(int nt=0; nt<ntheta; nt++){

	theta  = Constants::TWOPI*nt/(ntheta-1);
	sint   = sin(theta);
	cost   = cos(theta);

	dirn.set(cost, sint*cosp, sint*sinp);

	if(roche){
	
	    Roche::face(q, Roche::SECONDARY, spin2, dirn, rref, pref, acc, rvec, dvec, rad, gravity);
	  
	}else{
	  
	  // Assume spherical
	  rvec    = r2*dirn + cofm2;
	  dvec    = dirn;
	  
	}
	
	// Eclipse & visibility computation
	if(Subs::dot(earth, dvec) <= 0. || Roche::sphere_eclipse(earth, rvec, cofm1, r1, lam1, lam2) ||
	   Lcurve::disc_eclipse(iangle, phase, rdisc1, rdisc2, beta, height, rvec)){

	  ready = false;
	}else{
	  if(ready)
	    cpgdraw(Subs::dot(rvec, xsky), Subs::dot(rvec, ysky));
	  else
	    cpgmove(Subs::dot(rvec, xsky), Subs::dot(rvec, ysky));
	  ready = true;
	}
      }

      cpgebuf();

      ready = false;

    }

    cpgsci(4);

    // Loop over theta over star 1
    for(int nt=0; nt<ntheta1; nt++){

      cpgbbuf();

      theta = Constants::PI*(nt+0.5)/ntheta1;
      sint  = sin(theta);
      cost  = cos(theta);

      int nphi  = std::max(8,int(r1*nplot*sint/r2+0.5));

      for(int np=0; np<nphi; np++){

	phi  = Constants::TWOPI*np/(nphi-1);
	sinp = sin(phi);
	cosp = cos(phi);

	dirn.set(sint*cosp, sint*sinp, cost);

	// Assume spherical
	rvec    = r1*dirn + cofm1;
	dvec    = dirn;
	  
	// Eclipse & visibility computation
	if(Subs::dot(earth, dvec) <= 0. || 
	   (!roche && Roche::sphere_eclipse(earth, rvec, cofm2, rref, lam1, lam2)) || 
	   (roche  && Roche::fblink(q, Roche::SECONDARY, spin2, ffac, acc, earth, rvec)) ||
	   Lcurve::disc_eclipse(iangle, phase, rdisc1, rdisc2, beta, height, rvec)){
	  ready = false;
	}else{
	  if(ready)
	    cpgdraw(Subs::dot(rvec, xsky), Subs::dot(rvec, ysky));
	  else
	    cpgmove(Subs::dot(rvec, xsky), Subs::dot(rvec, ysky));
	  ready = true;
	}
      }

      cpgebuf();

      ready = false;

    }

    // Loop over phi over star 1
    for(int np=0; np<nphi1; np++){

      cpgbbuf();

      phi   = Constants::PI*(np+0.5)/nphi1;
      sinp  = sin(phi);
      cosp  = cos(phi);

      int ntheta  = int(r1*nplot/r2+0.5);

      for(int nt=0; nt<ntheta; nt++){

	theta  = Constants::TWOPI*nt/(ntheta-1);
	sint   = sin(theta);
	cost   = cos(theta);

	dirn.set(sint*cosp, sint*sinp, cost);

	// Assume spherical
	rvec    = r1*dirn + cofm1;
	dvec    = dirn;
	
	// Eclipse & visibility computation
	if(Subs::dot(earth, dvec) <= 0. || 
	   (!roche && Roche::sphere_eclipse(earth, rvec, cofm2, rref, lam1, lam2)) || 
	   (roche && Roche::fblink(q, Roche::SECONDARY, spin2, ffac, acc, earth, rvec)) ||
	   Lcurve::disc_eclipse(iangle, phase, rdisc1, rdisc2, beta, height, rvec)){
	  ready = false;
	}else{
	  if(ready)
	    cpgdraw(Subs::dot(rvec, xsky), Subs::dot(rvec, ysky));
	  else
	    cpgmove(Subs::dot(rvec, xsky), Subs::dot(rvec, ysky));
	  ready = true;
	}
      }

      cpgebuf();

      ready = false;


    }

    cpgsci(1);

    double rdisc;

    // Loop over theta for disc
    for(int nt=0; nt<ndisct; nt++){

      cpgbbuf();

      theta = Constants::TWOPI*(nt+0.5)/ndisct;

      int nrad  = std::max(2,int((rdisc2-rdisc1)*nplot/r2+0.5));

      for(int nr=0; nr<nrad; nr++){

	rdisc = rdisc1 + (rdisc2-rdisc1)*nr/(nrad-1);
	Lcurve::pos_disc(rdisc, theta, beta, height, rvec, dvec);
	  
	// Eclipse & visibility computation
	if(Subs::dot(earth, dvec) <= 0. || 
	   (!roche && Roche::sphere_eclipse(earth, rvec, cofm2, rref, lam1, lam2)) || 
	   (roche && Roche::fblink(q, Roche::SECONDARY, spin2, ffac, acc, earth, rvec)) ||
	   Roche::sphere_eclipse(earth, rvec, cofm1, r1, lam1, lam2) ||
	   Lcurve::disc_surface_eclipse(iangle, phase, rdisc1, rdisc2, beta, height, rvec)){
	  ready = false;
	}else{
	  if(ready)
	    cpgdraw(Subs::dot(rvec, xsky), Subs::dot(rvec, ysky));
	  else
	    cpgmove(Subs::dot(rvec, xsky), Subs::dot(rvec, ysky));
	  ready = true;
	}
      }

      cpgebuf();

      ready = false;

    }

    // Loop over radius of disc
    for(int nr=0; nr<ndiscr; nr++){

      cpgbbuf();

      rdisc   = rdisc1 + (rdisc2-rdisc1)*nr/(ndiscr-1);

      int ntheta  = int(rad*nplot/r2+0.5);

      for(int nt=0; nt<ntheta; nt++){

	theta  = Constants::TWOPI*nt/(ntheta-1);

	Lcurve::pos_disc(rdisc, theta, beta, height, rvec, dvec);
	  
	// Eclipse & visibility computation
	if(Subs::dot(earth, dvec) <= 0. || 
	   (!roche && Roche::sphere_eclipse(earth, rvec, cofm2, rref, lam1, lam2)) ||
	   (roche && Roche::fblink(q, Roche::SECONDARY, spin2, ffac, acc, earth, rvec)) ||
	   Roche::sphere_eclipse(earth, rvec, cofm1, r1, lam1, lam2) ||
	   Lcurve::disc_surface_eclipse(iangle, phase, rdisc1, rdisc2, beta, height, rvec)){
	  ready = false;
	}else{
	  if(ready)
	    cpgdraw(Subs::dot(rvec, xsky), Subs::dot(rvec, ysky));
	  else
	    cpgmove(Subs::dot(rvec, xsky), Subs::dot(rvec, ysky));
	  ready = true;
	}
      }

      cpgebuf();

      ready = false;

    }


    // Ring at upper edge of disc
    cpgbbuf();
    int nphi = int(rdisc2*nplot/r2+0.5);
    for(int np=0; np<nphi; np++){

      phi = Constants::TWOPI*np/(nphi-1);
      rvec.set(rdisc2*cos(phi), rdisc2*sin(phi), height*pow(rdisc2,beta));
 
      // Eclipse & visibility computation
      if(Subs::dot(earth, dvec) <= 0. || 
	 (!roche && Roche::sphere_eclipse(earth, rvec, cofm2, rref, lam1, lam2)) ||
	 (roche || Roche::fblink(q, Roche::SECONDARY, spin2, ffac, acc, earth, rvec)) ||
	 Roche::sphere_eclipse(earth, rvec, cofm1, r1, lam1, lam2)){ 
	ready = false;
      }else{
	if(ready)
	  cpgdraw(Subs::dot(rvec, xsky), Subs::dot(rvec, ysky));
	else
	  cpgmove(Subs::dot(rvec, xsky), Subs::dot(rvec, ysky));
	ready = true;
      }
    }

    cpgebuf();
    ready = false;

    // Ring at lower edge of disc
    cpgbbuf();
    nphi = int(rdisc2*nplot/r2+0.5);
    for(int np=0; np<nphi; np++){

      phi = Constants::TWOPI*np/(nphi-1);
      rvec.set(rdisc2*cos(phi), rdisc2*sin(phi), -height*pow(rdisc2,beta));
 
      // Eclipse & visibility computation
      if(Subs::dot(earth, dvec) <= 0. || 
	 (!roche && Roche::sphere_eclipse(earth, rvec, cofm2, rref, lam1, lam2)) ||
	 (roche || Roche::fblink(q, Roche::SECONDARY, spin2, ffac, acc, earth, rvec)) ||
	 Roche::sphere_eclipse(earth, rvec, cofm1, r1, lam1, lam2) ||
	 Lcurve::disc_eclipse(iangle, phase, rdisc1, rdisc2, beta, height, rvec)){
	ready = false;
      }else{
	if(ready)
	  cpgdraw(Subs::dot(rvec, xsky), Subs::dot(rvec, ysky));
	else
	  cpgmove(Subs::dot(rvec, xsky), Subs::dot(rvec, ysky));
	ready = true;
      }
    }

    cpgebuf();
    ready = false;


    // Ring at lower edge of disc
    cpgbbuf();
    double x, y;
    for(int ne=0; ne<nedge; ne++){

      phi  = Constants::TWOPI*ne/nedge;
      x    = rdisc2*cos(phi);
      y    = rdisc2*sin(phi);
      
      double h2 = height*pow(rdisc2,beta);
      int nvert = int(2*h2*nplot/r2+0.5);
      for(int nv=0; nv<nvert; nv++){
	
	rvec.set(x, y, -h2+2.*h2*nv/(nvert-1));
 
	// Eclipse & visibility computation
	if(Subs::dot(earth, dvec) <= 0. || 
	   (!roche && Roche::sphere_eclipse(earth, rvec, cofm2, rref, lam1, lam2)) ||
	   (roche && Roche::fblink(q, Roche::SECONDARY, spin2, ffac, acc, earth, rvec)) ||
	   Roche::sphere_eclipse(earth, rvec, cofm1, r1, lam1, lam2) ||
	   Lcurve::disc_eclipse(iangle, phase, rdisc1, rdisc2, beta, height, rvec)){
	  ready = false;
	}else{
	  if(ready)
	    cpgdraw(Subs::dot(rvec, xsky), Subs::dot(rvec, ysky));
	  else
	    cpgmove(Subs::dot(rvec, xsky), Subs::dot(rvec, ysky));
	  ready = true;
	}
      }
    
      cpgebuf();
      ready = false;
    }

    const int NSTREAM = std::max(2,int(nplot*(1-rl2)/r2+0.5));
    float xs[NSTREAM], ys[NSTREAM];
    Roche::streamr(q, rdisc2, xs, ys, NSTREAM);

    cpgslw(6);
    cpgsci(2);
    cpgbbuf();
    for(int ns=0; ns<NSTREAM; ns++){

      rvec.set(xs[ns], ys[ns], 0.);
 
      // Eclipse & visibility computation
      if((!roche && Roche::sphere_eclipse(earth, rvec, cofm2, rref, lam1, lam2)) ||
	 (roche && Roche::fblink(q, Roche::SECONDARY, spin2, ffac, acc, earth, rvec)) ||
	 Roche::sphere_eclipse(earth, rvec, cofm1, r1, lam1, lam2) ||
	 Lcurve::disc_eclipse(iangle, phase, rdisc1, rdisc2, beta, height, rvec)){
	ready = false;
      }else{
	if(ready){
	  cpgdraw(Subs::dot(rvec, xsky), Subs::dot(rvec, ysky));
	  if(ns == NSTREAM-1){
	    cpgsci(4);
	    cpgsch(spot);
	    cpgpt1(Subs::dot(rvec, xsky), Subs::dot(rvec, ysky), 15);
	  }
	}else{
	  cpgmove(Subs::dot(rvec, xsky), Subs::dot(rvec, ysky));
	}
	ready = true;
      }
    }
    
    cpgebuf();
    ready = false;
    
  }
  catch(const std::string& err){
    std::cerr << err << std::endl;
    exit(EXIT_FAILURE);
  }
}
