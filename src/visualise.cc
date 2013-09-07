/*

!!begin
!!title    Plots the appearance of Roche-distorted stars in a binary
!!author   T.R.Marsh 
!!created  17 Sep 2003
!!revised  29 May 2011
!!descr    Plots the appearance of Roche-distorted stars in a binary
!!index    visualise
!!root     visualise
!!css      style.css
!!class    Model
!!head1    visualise -- plots the appearance of Roche-distorted stars in a binary

!!emph{visualise} computes grids of points representing two stars and optionally a disc
and plots those visible at a particular phase projected onto the plane of the sky.
It can only handle stars less than or equal to their Roche lobes. This allows visualisation
and checking of the eclipse computations.

!!head2 Command invocation

visualise model nphase nphase (phase)/(phase1 phase2) device x1 x2 y1 y2

!!head2 Arguments

!!table
!!arg{model}{Parameter file defining the model, as used e.g. by !!ref{lroche.html}{lrcohe}. }
!!arg{nphase}{Number of orbital phases to display.}
!!arg{phase}{Orbital phase if nphase=1}
!!arg{phase1}{First orbital phase if nphase>1}
!!arg{phase2}{Last orbital phase if nphase>1}
!!arg{device}{Plot device}
!!arg{x1}{left-hand X limit}
!!arg{x2}{right-hand X limit}
!!arg{y1}{lower Y limit}
!!arg{y2}{upper Y limit}
!!arg{width}{width of plot in inches}
!!arg{reverse}{reverse black and white (or not)}
!!table

!!end

*/

#include <cstdlib>
#include <iostream>
#include "cpgplot.h"
#include "trm/subs.h"
#include "trm/plot.h"
#include "trm/vec3.h"
#include "trm/input.h"
#include "trm/roche.h"
#include "trm/lcurve.h"

int    Lcurve::Fobj::neval = 0;
double Lcurve::Fobj::chisq_min;
Subs::Buffer1D<double> Lcurve::Fobj::scale_min;

// Main program
int main(int argc, char* argv[]){

  // Defined at the end
  void plot_visible(const Subs::Buffer1D<Lcurve::Point>& object, const Subs::Vec3& earth, const Subs::Vec3& cofm, const Subs::Vec3& xsky, const Subs::Vec3& ysky, double phase);  
  
  try{
    
    // Construct Input object
    Subs::Input input(argc, argv, Lcurve::LCURVE_ENV, Lcurve::LCURVE_DIR);
    
    // sign-in input variables
    input.sign_in("model",    Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("nphase",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("phase",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("phase1",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("phase2",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("device",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("x1",       Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("x2",       Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("y1",       Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("y2",       Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("width",    Subs::Input::LOCAL,  Subs::Input::NOPROMPT);
    input.sign_in("reverse",  Subs::Input::LOCAL,  Subs::Input::NOPROMPT);
    
    std::string smodel;
    input.get_value("model", smodel, "model", "model file of parameter values");
    Lcurve::Model model(smodel);
    
    int nphase;
    input.get_value("nphase", nphase, 1, 1, 1000, "number of orbital phases");
    double phase, phase1, phase2;
    if(nphase == 1){
      input.get_value("phase", phase, 0., -10., 10., "orbital phase of interest");
    }else{
      input.get_value("phase1", phase1, 0., -10., 10., "first orbital phase of interest");
      input.get_value("phase2", phase2, 0., -10., 10., "last orbital phase of interest");
    }

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
    input.get_value("width", width,  8.f, 0.f, 100.f, "width of the plot (inches)");
    bool reverse;
    input.get_value("reverse", reverse, true, "reverse black and white?");

    input.save();

    double r1, r2, rdisc1=0., rdisc2=0.;
    model.get_r1r2(r1, r2);
    double rl1 = Roche::xl11(model.q,model.spin1);
    if(r1 <= 0) r1 = 0.99999999999*rl1;
    double rl2 = 1.-Roche::xl12(model.q,model.spin2);
    if(r2 <= 0) r2 = 0.99999999999*rl2;
    
    // Generate arrays over each star's face. 
    Subs::Buffer1D<Lcurve::Point> star1, star2, disc, outer_edge, inner_edge, bspot, stream;
    Lcurve::set_star_grid(model, Roche::PRIMARY, true, star1);
    Lcurve::set_star_grid(model, Roche::SECONDARY, true, star2);

    if(model.add_disc){

      rdisc1 = model.rdisc1 > 0. ? model.rdisc1 : r1;
      rdisc2 = model.rdisc2 > 0. ? model.rdisc2 : model.radius_spot;

      // note that the inner radius of the disc is set equal to that of the white dwarf if rdisc1 <= 0
      // while the outer disc is set equal to the spot radius
      Lcurve::set_disc_grid(model, disc);
      Lcurve::set_disc_edge(model, true, outer_edge);
      Lcurve::set_disc_edge(model, false, inner_edge);
	
      std::vector<std::pair<double,double> > eclipses;

      if(model.opaque){
	    
	// Apply eclipse by disc to star 1
	for(int i=0; i<star1.size(); i++){
	  eclipses =  Roche::disc_eclipse(model.iangle, rdisc1, rdisc2, model.beta_disc, model.height_disc, star1[i].posn);
	  for(size_t j=0; j<eclipses.size(); j++)
	    star1[i].eclipse.push_back(eclipses[j]);
	}
	    
	// Apply eclipse by disc to star 2
	for(int i=0; i<star2.size(); i++){
	  eclipses =  Roche::disc_eclipse(model.iangle, rdisc1, rdisc2, model.beta_disc, model.height_disc, star2[i].posn);
	  for(size_t j=0; j<eclipses.size(); j++)
	    star2[i].eclipse.push_back(eclipses[j]);  
	}
      }
    }

    if(model.add_spot){
      //      double wave = 5000, temp_spot = 10000, cfrac_spot = 0.5, height_spot = 0.01;
      //      Lcurve::set_bright_spot_grid(q, iangle, r1, r2, roche1, roche2, eclipse1, eclipse2, delta_phase, radius_spot, 
      //			   length_spot, height_spot, expon_spot, epow_spot, angle_spot, yaw_spot, temp_spot, tilt_spot, 
      //			   cfrac_spot, nspot, wave, bspot);
      Subs::Vec3 dir(1,0,0), posn, v;
	  
      double rl1 = Roche::xl1(model.q);
      
      // Calculate a reference radius and potential for the two stars
      double rref1, pref1, ffac1 = r1/rl1;
      Roche::ref_sphere(model.q, Roche::PRIMARY, model.spin1, ffac1, rref1, pref1);
      
      double rref2, pref2, ffac2 = r2/rl2;
      Roche::ref_sphere(model.q, Roche::SECONDARY, model.spin2, ffac2, rref2, pref2);
      
      dir.set(0,0,1);
      Roche::strinit(model.q, posn, v);
      
      Lcurve::Point::etype eclipses, edisc;
      Lcurve::star_eclipse(model.q, r1, model.spin1, ffac1, model.iangle, posn, model.delta_phase, model.roche1, Roche::PRIMARY,   eclipses);
      Lcurve::star_eclipse(model.q, r2, model.spin2, ffac2, model.iangle, posn, model.delta_phase, model.roche2, Roche::SECONDARY, eclipses);
      stream.push_back(Lcurve::Point(posn, dir, 0., 1., eclipses));
		       
      const int NSTREAM = int((rl1-model.radius_spot)/0.001);
      double radius;
      for(int i=0; i<NSTREAM; i++){
	radius = rl1 + (model.radius_spot-rl1)*(i+1)/NSTREAM;
	Roche::stradv(model.q, posn, v, radius, 1.e-10, 1.e-3);
	eclipses.clear();
	Lcurve::star_eclipse(model.q, r1, model.spin1, ffac1, model.iangle, posn, model.delta_phase, model.roche1, Roche::PRIMARY,   eclipses);
	Lcurve::star_eclipse(model.q, r2, model.spin2, ffac2, model.iangle, posn, model.delta_phase, model.roche2, Roche::SECONDARY, eclipses);
	if(model.add_disc){	       
	  edisc = Roche::disc_eclipse(model.iangle, rdisc1, rdisc2, model.beta_disc, model.height_disc, posn);
	  for(size_t j=0; j<edisc.size(); j++)
	    eclipses.push_back(edisc[j]);
	}
	stream.push_back(Lcurve::Point(posn, dir, 0., 1., eclipses));
      }	
    }
	
    // Plot
    Subs::Plot plot(device);

    cpgpap(width, (y2-y1)/(x2-x1));
    if(reverse){
      cpgscr(0,1,1,1);
      cpgscr(1,0,0,0);
      cpgscr(2,0.4,0,0);
      cpgscr(3,0,0.3,0);
      cpgscr(4,0,0,0.5);
    }
	
    // Make stars orbit around centre of mass of system
    const Subs::Vec3 cofm(model.q/(1.+model.q),0.,0.);
    Subs::Vec3 r, earth;  

    for(int np=0; np<nphase; np++){

      if(nphase > 1) phase = phase1 + (phase2-phase1)*np/double(nphase-1);

      cpgenv(x1, x2, y1, y2, 1, -2);
	  
      earth = Roche::set_earth(model.iangle, phase);
	  
      // Compute sky basis vectors
      double cosp = cos(Constants::TWOPI*phase);
      double sinp = sin(Constants::TWOPI*phase);
	  
      Subs::Vec3 xsky(sinp,cosp,0.);
      Subs::Vec3 ysky = Subs::cross(earth, xsky);

      cpgbbuf();
	  
      cpgsch(2);
	  
      // star 1
      cpgsci(4);
      plot_visible(star1, earth, cofm, xsky, ysky, phase);  
	  
      // star 2
      cpgsci(2);
      plot_visible(star2, earth, cofm, xsky, ysky, phase);  
	  
      if(model.add_disc){
	    
	// disc surface
	cpgsci(3);
	plot_visible(disc, earth, cofm, xsky, ysky, phase);  
	    
	// edges
	cpgsci(1);
	plot_visible(outer_edge, earth, cofm, xsky, ysky, phase);  
	plot_visible(inner_edge, earth, cofm, xsky, ysky, phase);  
	    
      }
	  
      if(model.add_spot){
	    
	cpgsci(2);
	plot_visible(stream, earth, cofm, xsky, ysky, phase);  
	    
	cpgsci(3);
	if(Subs::dot(earth, stream[stream.size()-1].posn) > 0. && stream[stream.size()-1].visible(phase)){
	  r = stream[stream.size()-1].posn - cofm;
	  cpgsch(4);
	  cpgslw(4);
	  cpgpt1(Subs::dot(r, xsky), Subs::dot(r, ysky), 18);
	}
      }
	  
      cpgebuf();
	  
    }
  }
  catch(const std::string& err){
    std::cerr << err << std::endl;
    exit(EXIT_FAILURE);
  }
}


// plots visible points
void plot_visible(const Subs::Buffer1D<Lcurve::Point>& object, const Subs::Vec3& earth, const Subs::Vec3& cofm, const Subs::Vec3& xsky, const Subs::Vec3& ysky, double phase){
  Subs::Vec3 r;
  for(int i=0; i<object.size(); i++){
    if(Subs::dot(earth, object[i].dirn) > 0. && object[i].visible(phase)){
      r = object[i].posn - cofm;
      cpgpt1(Subs::dot(r, xsky), Subs::dot(r, ysky), 1);
    }
  }
}
