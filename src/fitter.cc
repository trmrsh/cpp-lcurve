/*

!!begin
!!title    Fits light-curve to a model with a sphere plus Roche-lobe distorted star
!!author   T.R.Marsh 
!!created  24 June 2005
!!revised  28 June 2005
!!descr    Fits light-curve to a model with a sphere plus Roche-lobe distorted star
!!index    fitter
!!root     fitter
!!css      style.css
!!class    Model
!!head1    Fits light-curve to a model with a sphere plus Roche-lobe distorted star

Using a file of supplied parameter values, fitter tries to optimise the fit between
a light curve and a model of a spherical and tidally distorted star. This is 
supposed to model a white dwarf/main-sequence binary such as NN Ser. It includes
reprocessing and eclipses. The reprocessing is computed using the simple addition of
fluxes method, assuming that the irradiating star can be treated as a point, and not 
including any 'back heating'. Phase 0 is defined as the point when star 1 is furthest from the observer.
Star 1 is the spherical one. Efforts are made to save time by using a coarser grid away from eclipse
and using no sub-division for these phases. This comes with some danger, mainly of discontinuities
at the cross over points. A rescaling is computed from the values of each grid at a cross-over point.
They should be close to 1. Setting all the element sizes to the same value will avoid any such problem, 
but will be slow.

The minimisation uses 'amoeba' which works with a group of N+1 points in N dimensional space and moves the points 
around to move downhill in chi**2.

!!table
!!arg{model}{A file of initial parameter values and an indication of whether they should be varied
or not for physical parameters and name = value for others. See !!ref{lroche}{lroche} for a full description of these.}
!!arg{data}{A light curve of times, exposure times, fluxes and uncertainties. # or blank to ignore,
negative errors ro mask.}
!!arg{ftol}{Once the N+1 corners are all within a fraction ftol of each other the routine will stop.}
!!arg{nmax}{The routine will also stop if there have this many model evaluations}
!!arg{output}{The file to store the final best model in}
!!table

!!end

*/

#include <cfloat>
#include <cstdlib>
#include <iostream>
#include "trm_subs.h"
#include "trm_format.h"
#include "trm_input.h"
#include "amoeba.h"
#include "trm_lcurve.h"

using Subs::operator+;

// Main program
int main(int argc, char* argv[]){
  
  try{

    // Construct Input object
    Subs::Input input(argc, argv, Lcurve::LCURVE_ENV, Lcurve::LCURVE_DIR);

    // Sign-in input variables
    input.sign_in("model",    Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("data",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("ftol",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nmax",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("output",   Subs::Input::GLOBAL, Subs::Input::PROMPT);

    std::string smodel;
    input.get_value("model", smodel, "model", "input model file");
    Lcurve::Model model(smodel);

    std::string sdata;
    input.get_value("data", sdata, "data", "data file (time, exposure, flux, error)");    
    Lcurve::Data data(sdata);

    double ftol;
    input.get_value("ftol", ftol, 1.e-5, DBL_MIN, 0.1, "fractional tolerance for convergence");    

    int nmax;
    input.get_value("nmax", nmax, 1000, 1, INT_MAX, "maximum number of model evaluations"); 

    std::string omodel;
    input.get_value("output", omodel, "model", "output model file");   

    // Construct function object
    Lcurve::Fobj func(model, data);

    // Load initial values and steps
    std::vector<double> corner, step;
    model.get_start_values(corner, step);
    std::vector<double> ncorner = corner;

    // Make corners for amoeba
    std::vector<std::pair<std::vector<double>, double> > params;
    int ndim = model.nvary();
    double chisq = func(corner);
    params.push_back(std::make_pair(corner, chisq));

    for(int i=1; i<ndim+1; i++){
      for(int j=0; j<ndim; j++)
	ncorner[j] = corner[j];
      ncorner[i-1] += step[i-1];
      chisq = func(ncorner);
      params.push_back(std::make_pair(ncorner, chisq));
    } 

    // Now run it
    int neval;
    Subs::amoeba(params, ftol, nmax, func, neval); 

    int nbest = 0;
    double cmin = params[0].second;
    for(int i=1; i<params.size(); i++){
      if(params[i].second < cmin){
	nbest = i;
	cmin = params[i].second;
      }
    }

    // Store results
    model.load_param(params[nbest].first);
    model.wrasc(omodel);

    Subs::Format form;
    std::cout << "Minimum chi**2 = " << form(func.chisq_min) << ", scale factor = " << form(func.scale_min) << ", number of model computations = " << func.neval << std::endl;
    std::cout << "Best model written to " << omodel << std::endl;

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
