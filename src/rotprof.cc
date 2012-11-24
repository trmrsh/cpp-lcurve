/*

!!begin
!!title    Computes rotational profile of one star in binary
!!author   T.R.Marsh 
!!created  20 Mar 2007
!!descr    Computes rotational profile of one star in binary
!!index    rotprof
!!root     rotprof
!!css      style.css
!!class    Model
!!head1    Computes rotational profile of one star in binary

!!emph{rotprof} computes the rotational profile from one of 
two stars in a binary

!!head2 Invocation

rotprof model time expose ndiv star spin vlo vhi nbin fwhm

!!head2 Arguments

!!table
!!arg{model}{File of parameters specifying the parameters. See !!ref{lroche.html}{lroche} for a full description.}
!!arg{time}{Time to compute.}
!!arg{expose}{Exposure time}
!!arg{ndiv}{Phase subdivision factor}
!!arg{star}{1 or 2 for primary or secondary}
!!arg{spin}{Spin frequency of star relative to orbit}
!!arg{vlo}{Lowest velocity of output (left-edge of lowest bin)}
!!arg{vhi}{Highest velocity of output (right-edge of highest bin)}
!!arg{nbin}{Number of bins}
!!arg{fwhm}{FWHM of blurring}
!!table

Note that the parameters vsini, vlo, vhi and fwhm are scaled in terms of K1+K2.
!!end

*/

#include <cfloat>
#include <cstdlib>
#include <iostream>
#include "trm_subs.h"
#include "trm_vec3.h"
#include "trm_input.h"
#include "trm_format.h"
#include "trm_roche.h"
#include "trm_lcurve.h"

int Lcurve::Fobj::neval = 0;
double Lcurve::Fobj::chisq_min;
Subs::Buffer1D<double> Lcurve::Fobj::scale_min;

// Main program
int main(int argc, char* argv[]){
  
  try{

    // Construct Input object
    Subs::Input input(argc, argv, Lcurve::LCURVE_ENV, Lcurve::LCURVE_DIR);

    // Sign-in input variables
    input.sign_in("model",    Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("time",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("expose",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("ndiv",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("star",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("vlo",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("vhi",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nbin",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("fwhm",     Subs::Input::LOCAL,  Subs::Input::PROMPT);

    std::string smodel;
    input.get_value("model", smodel, "model", "model file of parameter values");
    Lcurve::Model model(smodel);

    double time;
    input.get_value("time", time, 0., -DBL_MAX, DBL_MAX, "time to compute");
    double expose;
    input.get_value("expose", expose, 0., 0., DBL_MAX, "exposure time (same units as ephemeris)");
    int ndiv;
    input.get_value("ndiv", ndiv, 1, 1, 1000000, "phase division factor");
    int nstar;
    input.get_value("star", nstar, 1, 1, 2, "star number");
    double vlo;
    input.get_value("vlo", vlo, -10., -DBL_MAX, DBL_MAX, "lowest velocity to store");
    double vhi;
    input.get_value("vhi", vhi,  std::max(vlo, 10.), vlo, DBL_MAX, "highest velocity to store");
    int nbin;
    input.get_value("nbin", nbin, 1, 1, 1000000, "number of bins");
    double fwhm;
    input.get_value("fwhm", fwhm,  0.01, 0., 10000., "FWHM blurring");
    input.save();

    double r1, r2;
    model.get_r1r2(r1, r2);
    double rl1 = Roche::xl11(model.q,model.spin1);
    if(r1 <= 0) r1 = 0.99999999999*rl1;
    double rl2 = 1.-Roche::xl12(model.q,model.spin2);
    if(r2 <= 0) r2 = 0.99999999999*rl2;
    double spin    = (nstar == 1) ? model.spin1 : model.spin2;
    Lcurve::LDC ldc1( model.ldc1_1.value,  model.ldc1_2.value,  model.ldc1_3.value,  model.ldc1_4.value, model.mucrit1, model.limb1); 
    Lcurve::LDC ldc2( model.ldc2_1.value,  model.ldc2_2.value,  model.ldc2_3.value,  model.ldc2_4.value, model.mucrit2, model.limb2); 
    
    // Generate arrays over each star's face. Two sets are generated for star2 because it is wasteful
    // to compute a vast number of faces away from eclipse.
    Subs::Buffer1D<Lcurve::Point> star;
    
    Subs::Vec3 cstar;
    Lcurve::LDC ldc;
    if(nstar == 1){
      set_star_grid(model, Roche::PRIMARY, true, star);
	cstar.set(0,0,0);
	ldc = ldc1;
    }else if(nstar == 2){
      set_star_grid(model, Roche::SECONDARY, true, star);
	cstar.set(1,0,0);
	ldc = ldc2;
    }

    double xoff = cstar.x() - model.q/(1+model.q);

    // Set up accumulation buffer
    Subs::Array1D<double> fine(nbin);
    fine = 0.;
    Subs::Vec3 earth;

    for(int j=0; j<ndiv; j++){
	double phase;
	if(ndiv == 1){
	    phase = time;
	}else{
	    phase = time + expose*((j-1)/(ndiv-1)-0.5);
	}
	phase = (phase - model.t0)/model.period;

	earth = Roche::set_earth(model.iangle, phase);

	double cosp = cos(Constants::TWOPI*phase);
	double sinp = sin(Constants::TWOPI*phase);

	for(int i=0; i<star.size(); i++){
	    double mu = Subs::dot(earth, star[i].dirn);
      
	    if(ldc.see(mu) && star[i].visible(phase)){
		Subs::Vec3 off = star[i].posn - cstar;
		double v = xoff*sinp + spin*(off.y()*cosp + off.x()*sinp);
		int ind  = int(nbin*(v-vlo)/(vhi-vlo));
		if(ind >=0 && ind < nbin)
		    fine[ind] += star[i].area*ldc.imu(mu);
	    }
	}
    }

    Subs::Array1D<double> blurr(nbin);
    blurr = 0.;

    // Total width of 4*FWHM
    double vbin = (vhi-vlo)/nbin;
    const int NGAUSS = 2*int(2*fwhm/vbin)+1;
    const int NGMID  = NGAUSS/2;
    Subs::Array1D<double> gauss(NGAUSS);
    gauss = 0.;
    for(int i=0; i<NGAUSS; i++){
	gauss[NGMID+i] = exp(-pow(vbin*i/(fwhm/Constants::EFAC),2)/2);
    }
    for(int i=0; i<nbin; i++){
	for(int j=std::max(-NGMID,-i); j<std::min(NGMID+1, nbin-i); j++){
	    blurr[i+j] += gauss[NGMID+j]*fine[i];
	}
    }

    Subs::Format form(10);
    for(int i=0; i<nbin; i++){
	std::cout << form(vlo + vbin*(i+0.5)) << " " << form(blurr[i]) << std::endl;
    }

  }
  catch(const Lcurve::Lcurve_Error& err){
    std::cerr << "Lcurve::Lcurve_Error exception thrown" << std::endl;
    std::cerr << err << std::endl;
    exit(EXIT_FAILURE);
  }
  catch(const std::string& err){
    std::cerr << err << std::endl;
    exit(EXIT_FAILURE);
  }
}
