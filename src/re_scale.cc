#include "trm/subs.h"
#include "trm/lcurve.h"

/** This routine scales a given fit to obtain a minimum chi**2 (equivalent to assuming complete
 *  uncertainty in the distance and/or absolute scale of the system. It returns the optimum scaling 
 * factor and the re-scaled fit.
 * \param data  input data. Points with negative errors are ignored.
 * \param fit   input fit. On return this is scaled to the optimum value unless wnok = 0
 * \param chisq weighted chisq
 * \param wnok  weighted number of OK points
 * \return the optimum scale factor or 1 if nok = 0.
 */

double Lcurve::re_scale(const Lcurve::Data& data, Subs::Array1D<double>& fit, double& chisq, double& wnok){
    double wgt, sdy = 0., syy = 0.;
    wnok = 0.;
    for(int i=0; i<fit.size(); i++){
	if(data[i].weight > 0.){
	    wgt   = data[i].weight/Subs::sqr(data[i].ferr);
	    sdy  += wgt*data[i].flux*fit[i];
	    syy  += wgt*fit[i]*fit[i];
	    wnok += data[i].weight;
	}
    }
    double scale;
    if(wnok > 0.0){
	scale = sdy / syy;
	fit  *= scale;
	
	chisq = 0.;
	for(int i=0; i<fit.size(); i++)
	    if(data[i].weight > 0.0)
		chisq += data[i].weight*Subs::sqr((data[i].flux - fit[i])/data[i].ferr);
	
    }else{
	scale = 1.;
	chisq = 0.;
    }
    return scale;
}

