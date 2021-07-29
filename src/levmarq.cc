/*

!!begin
!!title    Fits light-curve to a binary star model
!!author   T.R.Marsh 
!!created  18 May 2006
!!revised  12 Mar 2008
!!descr    Fits light-curve to a binary star model using Levenberg-Marquardt
!!index    levmarq
!!root     levmarq
!!css      style.css
!!class    Model
!!head1    Fits light-curve to a model with a binary star model using Levenberg-Marquardt

Using a file of supplied parameter values, !!emph{levmarq} tries to optimise
the fit between a light curve and a binary star model; see
!!ref{lroche.html}{lroche} for details of the model.  The minimisation uses
the Levenberg-Marquardt method which uses numerical derivatives. Centred
derivatives are used so that if N parameters are variable, one iteration can
require 2*N+1 full model calculations. The routine reports the value of the
lambda parameter that switches between inverse-Hessian (small lambda,
optimistic quadratic approximation to update the parameters) and steepest
descent (large lambda, needed when things become difficult, but slow in
general). The larger lambda becomes, the smaller the step size in the steepest
descent case.  Chi**2 values are reported whenever any decrease is
detected. The scale factor reported is the optimum multiplier from fit to data
and 'neval' is the number of model computations.  The routine skips the
derivative calculation in cases where chi**2 increases so the time between
steps during which lambda increases are shorter than those for which lambda
decreases.

It is always difficult to know when to stop in minimisation routines. This routine
stops for three reasons: (i) when a maximum iteration count is achieved, (ii) when
there has been a lack of progress measured by a failure to decrease in chi**2
three times in succession and (iii) when the value of lambda exceeds a user-defined
maximum.

!!table
!!arg{model}{A file of initial parameter values and an indication of whether they should be varied
or not for physical parameters and name = value for others. See !!ref{lroche.html}{lroche} for a full description of these.
The amounts of variation will be used as step sizes for the numerical derivatives so choose carefully. NR says that
for centred derivatives, the optimal choice of step size = (acc)**(1/3) times a characteristic scale of your
function (or more precisely times (f/f''')**(1/3)), where acc is the fractional accuracy of the machine.
If we take this to be 10**-15 for double precision, then this implies 10**-5 times the scale over which you think 
the function changes lots.}
!!arg{data}{A light curve of times, exposure times, fluxes and uncertainties. # or blank to ignore,
negative errors to mask.}
!!arg{nmax}{The routine will stop if there have been this many iterations.}
!!arg{delta}{The change in chi**2 to require. This criterion only applies to steps where progress is being made
(i.e. the alambda parameter is going down). If the last three such steps have all failed to decrease chi**2 by
at least delta, then iterations will be halted.}
!!arg{lmax}{The maximum value of lambda to allow. Should be much larger than 1, e.g. 2e10, because the 
value of lambda is multiplied or divided by 10 each step according to progress.}
!!arg{scale}{true/false to scale the erros to give a reduced chi**2 = 1 or not. Generally advise yes.}
!!arg{output}{The file to store the final best model in}
!!table

!!end

*/

#include <climits>
#include <cfloat>
#include <cstdlib>
#include <iostream>
#include "trm/subs.h"
#include "trm/array1d.h"
#include "trm/buffer2d.h"
#include "trm/format.h"
#include "trm/input.h"
#include "trm/lcurve.h"

// Class to store the data and calculate the Levenberg-Marquardt stuff

class Lmfunc {
    
public:
    
    // Total number of calls
    static int neval;
    
    // Minimum chisq
    static double chisq_min;
    
    // Scale factors for minimum chisq
    static Subs::Buffer1D<double> scale_min;
  
    // Constructor; stores the data needed to compute chi**2
    Lmfunc(const Lcurve::Model& model, const Lcurve::Data& data) : model(model), data(data) {}
  
    // The main work is done by this routine
    void lmcomp(Subs::Buffer2D<double>& alpha, Subs::Buffer1D<double>& beta, double& chisq);

    // Number of variable parameters
    int nvar() const {return model.nvary();}

    // Get the current parameters
    Subs::Array1D<double> get_param() const {
	return model.get_param();
    }

    // Set the current parameters
    void set_param(const Subs::Array1D<double>& param){
	model.set_param(param);
    }

private:
  
    Lcurve::Model model;
    const Lcurve::Data&  data;
  
};

void lmfit(Lmfunc& func, double& chisq, double& alambda, Subs::Buffer2D<double>& covar);

int    Lcurve::Fobj::neval = 0;
double Lcurve::Fobj::chisq_min;
Subs::Buffer1D<double> Lcurve::Fobj::scale_min;

int    Lmfunc::neval = 0;
double Lmfunc::chisq_min;
Subs::Buffer1D<double> Lmfunc::scale_min;

// Main program
int main(int argc, char* argv[]){

    try{

	// Construct Input object
	Subs::Input input(argc, argv, Lcurve::LCURVE_ENV, Lcurve::LCURVE_DIR);

	// Sign-in input variables
	input.sign_in("model",    Subs::Input::GLOBAL, Subs::Input::PROMPT);
	input.sign_in("data",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
	input.sign_in("nmax",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("delta",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("lmax",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("scale",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
	input.sign_in("output",   Subs::Input::GLOBAL, Subs::Input::PROMPT);

	std::string smodel;
	input.get_value("model", smodel, "model", "input model file");
	Lcurve::Model model(smodel);

	std::string sdata;
	input.get_value("data", sdata, "data", "data file (time, exposure, flux, error)");    
	Lcurve::Data data(sdata);

	int nmax;
	input.get_value("nmax", nmax, 10, 1, INT_MAX, "maximum number of iterations"); 

	double delta;
	input.get_value("delta", delta, 0.001, 0., DBL_MAX, "the minimum decrease in chi**2"); 

	double lmax;
	input.get_value("lmax", lmax, 1e12, 10., DBL_MAX, "the maximum value of lambda"); 

	bool scale;
	input.get_value("scale", scale, true, "scale errors to give reduced chisq = 1?"); 

	std::string omodel;
	input.get_value("output", omodel, "model", "output model file");   

	Lmfunc func(model, data);
	Lmfunc::neval = 0;

	double wdof = 0.;
	for(size_t i=0; i<data.size(); i++)
	    wdof += data[i].weight;
	if(wdof <= 0.)
	    throw Lcurve::Lcurve_Error("Total weight of input data <= 0!");

	Subs::Buffer2D<double> covar;
	double lambda = -1., lambda_old, chisq = 0, chisq_old = -1.;
	int ncount = 0, nfail = 0;
	const int NFMAX = 3;
	Subs::Format val(17), err(8);

	// Stop either when the maximum number of iterations has been reached, chi**2
	// has not decreased significantly for a while .
	while(ncount < nmax && (nfail < NFMAX || chisq > chisq_old) && lambda < lmax){
	    lambda_old = lambda;
	    chisq_old  = chisq;
	    lmfit(func, chisq, lambda, covar);
	    std::cout << "lambda: " << lambda_old << " ----> " << lambda << "; chisq: " << err(chisq_old) << " ----> " << err(chisq) << std::endl;

	    // We check whether any progress is being made. If yes, then we check whether
	    // it is enough. If not we add to a counter which is checked for above. 
	    if(ncount){
		if(chisq < chisq_old){
		    if(chisq_old - chisq < delta)
			nfail++;
		    else
			nfail = 0;
		}
	    }

	    // iteration counter.
	    ncount++;
	}

	// Report reason for halting
	if(ncount == nmax) std::cout << "Iterations halted on reaching maximum number = " << nmax << std::endl;
	if(nfail == NFMAX) std::cout << "Iterations halted after chi**2 decreased by less than " << delta << " " << NFMAX << " times in succession." << std::endl;
	if(lambda >= lmax) std::cout << "Iterations halted after lambda = " << lambda << " exceeded threshold = " << lmax << std::endl;
	  
	// set lambda = 0 to compute covariances
	lambda = 0.;
	lmfit(func, chisq, lambda, covar);

	Subs::Array1D<double> param = func.get_param();
	model.set_param(param);

	std::cout << "\nwdof= " << wdof << std::endl;

	std::cout << "\nParameter values:\n" << std::endl;
	for(int i=0; i<func.nvar(); i++){
	    if(scale)
		std::cout << model.get_name(i) << " = " << val(param[i]) << " +/- " << err(sqrt(Lmfunc::chisq_min/wdof*covar[i][i])) << std::endl; 
	    else
		std::cout << model.get_name(i) << " = " << val(param[i]) << " +/- " << err(sqrt(covar[i][i])) << std::endl; 
	}

	std::cout << "\nCorrelation coefficients:\n" << std::endl;
	for(int i=0; i<func.nvar(); i++){
	    for(int j=0; j<i; j++)
		std::cout << "r(" << model.get_name(i) << "," << model.get_name(j) << ") = " << err(covar[i][j]/sqrt(covar[i][i]*covar[j][j])) << std::endl;
	}
 
	// Store results
	model.wrasc(omodel);

	std::cout << "Minimum weighted chi**2 = " << err(Lmfunc::chisq_min) << ", scale factors = " 
		  << err(Lmfunc::scale_min[0]) << " " << err(Lmfunc::scale_min[1]) << " " << err(Lmfunc::scale_min[2]) << " " << err(Lmfunc::scale_min[3]) 
		  << ", number of model computations = " << Lmfunc::neval << std::endl;
	std::cout << "(Weighted) number of data points = " << wdof << ", chi**2/wdof = " << err(Lmfunc::chisq_min/wdof) << std::endl;
	if(scale)
	    std::cout << "Uncertainties were scaled so that effective chi**2/wdof = 1" << std::endl;
	else
	    std::cout << "Uncertainties are 'raw', i.e. they were not scaled so that effective chi**2/wdof = 1" << std::endl;

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


// This does the main work, equivalent to mrqcof in NR. Starting from the current model, 
// it computes the lightcurves and its derivatives with respect to each parameter and then 
// uses these to calculate the Levenberg-Marquardt matric and vector.
void Lmfunc::lmcomp(Subs::Buffer2D<double>& alpha, Subs::Buffer1D<double>& beta, double& chisq){
  
    bool start = (Lmfunc::neval == 0);

    // Retrieve the current variable parameters and step sizes
    Subs::Array1D<double> centre = model.get_param(), dstep = model.get_dstep();

    // Compute the fit
    Subs::Buffer1D<double> sfac(4), tsfac(4);
    Subs::Array1D<double> fit;
    double wdwarf, wnok, logg1, logg2, rv1, rv2;

    Lcurve::light_curve_comp(model, data, true, true, false, sfac, fit, wdwarf,
                             chisq, wnok, logg1, logg2, rv1, rv2);
    if(wnok == 0.0)
        throw Lcurve::Lcurve_Error("void Lmfunc::lmcomp: no good data!");
    Lmfunc::neval++;

    // Compute derivatives using finite differences centred on the
    // current point, displacing according to the step sizes.

    Subs::Array1D<double> buff(data.size()), tparam;
    std::vector<Subs::Array1D<double> > deriv(centre.size());
    double tchisq;

    for(int i=0; i<centre.size(); i++){
        tparam     = centre;
        tparam[i] += dstep[i];
        model.set_param(tparam);
        Lcurve::light_curve_comp(model, data, true, true, false, tsfac, deriv[i],
                                 wdwarf, tchisq, wnok, logg1, logg2, rv1, rv2);
        Lmfunc::neval++;

        // Next four lines are an attempt to reduce round-off error
        // in the finite difference.
        double temp = tparam[i];
        tparam[i] -= 2.*dstep[i];
        double h = temp - tparam[i];
        tparam[i] = temp - h;

        model.set_param(tparam);
        Lcurve::light_curve_comp(model, data, true, true, false, tsfac, buff,
                                 wdwarf, tchisq, wnok, logg1, logg2, rv1, rv2);
        Lmfunc::neval++;

        deriv[i] -= buff;
        deriv[i] /= h;

    }

    // Resave the stored parameters
    model.set_param(centre);

    // Now we can get on with the standard Levenberg-Marquardt stuff
    Subs::Buffer1D<double> dyda(model.nvary());
  
    // Initialise alpha and beta
    alpha = 0.;
    beta  = 0.;

    double wgt, dy, wt;
    for(int i=0; i<fit.size(); i++){
	if(data[i].weight > 0.0){
	    wgt = data[i].weight/Subs::sqr(data[i].ferr);
	    dy  = data[i].flux - fit[i];
	    
	    for(int l=0, j=0; l<model.nvary(); l++){
		wt = wgt*deriv[l][i];
		for(int m=0, k=0; m<=l; m++)
		    alpha[j][k++] += wt*deriv[m][i];
		beta[j++]  += wt*dy;
	    }
	}
    }

    for(int j=1; j<model.nvary(); j++)
	for(int k=0; k<j; k++) alpha[k][j] = alpha[j][k];

    if(start || chisq < Lmfunc::chisq_min){
	Lmfunc::scale_min = sfac;
	Lmfunc::chisq_min = chisq;
	Subs::Format form(12);
	std::cout << "Weighted chi**2 = " << form(chisq) << ", scale factors = " << 
	    form(sfac[0]) << ", " << form(sfac[1]) << ", " << form(sfac[2]) << ", " << form(sfac[3]) << ", neval = " << neval << std::endl;
    }

}  

/**
 * Implementation of mmrqmin for the problem in hand. 
 * \param param initial values of the variable parameters
 * \param func the function that does most of the work.
 * \param chisq the chii**2
 * \param lambda multiplier, < 0 to initialise
 * \param covar covariances returned when the routine is called with lambda = 0
 */

void lmfit(Lmfunc& func, double& chisq, double& lambda, Subs::Buffer2D<double>& covar){

    static Subs::Buffer1D<double> beta;
    static Subs::Array1D<double> da;
    static Subs::Buffer2D<double> oneda, alpha;
    static double ochisq;

    // Initialisation
    if(lambda < 0.){
	beta.resize(func.nvar());
	da.resize(func.nvar());
	oneda.resize(func.nvar(),1);
	alpha.resize(func.nvar(),func.nvar());
	covar.resize(func.nvar(),func.nvar());
	lambda=0.001;
	func.lmcomp(alpha, beta, chisq);
	ochisq = chisq;
    }

    // Alter linearised fitting matrix by augmenting diagonal elements
    covar = alpha;
    for (int j=0; j<func.nvar(); j++){
	covar[j][j] = alpha[j][j]*(1.0+lambda);
	oneda[j][0] = beta[j];
    }

    // Matrix solution to find a (hopefully) better solution
    Subs::gaussj(covar, oneda);
    for (int j=0; j<func.nvar(); j++) da[j] = oneda[j][0];

    // If lambda set = 0, this indicates convergence
    if(lambda == 0.) return;

    // Try out new solution
    Subs::Array1D<double> asave = func.get_param();

    func.set_param(asave + da);

    func.lmcomp(covar, da, chisq);

    if(chisq < ochisq){

	// Chi**2 reduced ==> success, so reduce lambda
	lambda  *= 0.1;
	ochisq   = chisq;
	alpha    = covar;
	beta     = da;

    }else{

	// Damn, Chi**2 increased ==> retrieve old parameters, back off on lambda
	lambda  *= 10.0;
	func.set_param(asave);
	chisq    = ochisq;
    }
}

