/*

!!begin
!!title    Computes the line profile from irradiated face of a Roche distorted star
!!author   T.R.Marsh 
!!created  21 Sep 2003
!!revised  29 May 2011
!!descr    Computes the line profile from irradiated face of a Roche distorted star
!!index    lprofile
!!root     lprofile
!!css      style.css
!!class    Model
!!head1    Computes the line profile from a Roche distorted star

!!emph{lprofile} computes the line profile from a Roche distorted
star, allowing for irradiation of one face.  The radiation from the
star is modelled as a slab of constant optical depth with a source
function that changes exponentially with depth. This allows one to
have a continuous transition from optically thin to thick and have
limb darkening or brightening. The computation is carried out for the
secondary star only.

!!table
!!arg{model}{File of parameters specifying the parameters. See !!ref{lroche.html}{lroche} for a full description of these.}
!!arg{tau}{Total vertical optical depth of emitting layer. Allows one to simulate optically thin (tau << 1) and 
optically thick (tau >> 1) cases.}
!!arg{efac}{Source function changes as exp(efac*tau). This must be less than 1. efac > 0 is equivalent to the usual 
limb darkened case. efac < 0 implies limb brightening.}
!!arg{aline}{Parameter controlling strength of absorption line which will have limb and gravity darkening as defined by 
the model file, i.e. it is assumed to be of constant equivalent width, which may not be accurate of course.}
!!arg{nphase}{number of phases}
!!arg{phase}{Phase to compute (if nphase == 1)}
!!arg{phase1}{First phase to compute (nphase > 1)}
!!arg{phase2}{Last phase to compute (nphase > 1)}
!!arg{store}{molly file to store results (nphase > 1)}
!!arg{expose}{Length of exposure in terms of orbital phase}
!!arg{ndiv}{Number of sub-divisions to account for finite exposures. The results will be trapezoidally
integrated over the length of the exposure.}
!!arg{v1}{left hand edge of profile array (same units as vscale)}
!!arg{v2}{right hand edge of profile array (same units as vscale)}
!!arg{nbin}{number of velocity bins}
!!arg{fwhm}{full-width half maximum of any blurring (same units as vscale)}
!!arg{nfine}{fine sub-division factor for blurring. Each contribution is added into the nearest pixel
in a fine array and then blurred.}
!!arg{align}{true to shift out the motion of the centre of mass of the secondary star. This allows the
profiles to be compared more easily although not necessarily directly with data.}
!!table

!!end

*/

#include <cstdlib>
#include <cfloat>
#include <iostream>
#include "trm/subs.h"
#include "trm/format.h"
#include "trm/vec3.h"
#include "trm/input.h"
#include "trm/roche.h"
#include "trm/colly.h"
#include "trm/lcurve.h"

int Lcurve::Fobj::neval = 0;
double Lcurve::Fobj::chisq_min;
Subs::Buffer1D<double> Lcurve::Fobj::scale_min;

// Main program
int main(int argc, char* argv[]){

    try{

        // Construct Input object
        Subs::Input input(argc, argv, Lcurve::LCURVE_ENV, Lcurve::LCURVE_DIR);

        // sign-in input variables
        input.sign_in("model",    Subs::Input::GLOBAL, Subs::Input::PROMPT);
        input.sign_in("tau",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("efac",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("aline",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("nphase",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("phase",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("phase1",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("phase2",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("store",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("expose",   Subs::Input::GLOBAL, Subs::Input::PROMPT);
        input.sign_in("ndiv",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
        input.sign_in("v1",       Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("v2",       Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("nbin",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("fwhm" ,    Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("nfine",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
        input.sign_in("align",    Subs::Input::LOCAL,  Subs::Input::PROMPT);

        std::string smodel;
        input.get_value("model", smodel, "model", "model file of parameter values");
        Lcurve::Model model(smodel);
        double tau;
        input.get_value("tau", tau, 1., 0., 1.e20, "total vertical optical depth");
        double befac;
        input.get_value("efac", befac, 0., -1000., 0.99999, "source function exponential factor");
        double aline;
        input.get_value("aline", aline, 0., -DBL_MAX, DBL_MAX, "absorption line strength factor");
        int nphase;
        input.get_value("nphase", nphase, 1, 1, 1000, "number of orbital phases");
        double phase, phase1, phase2;
        if(nphase == 1){
            input.get_value("phase", phase, 0., -10., 10., "orbital phase of interest");
        }else{
            input.get_value("phase1", phase1, 0., -10., 10., "first orbital phase of interest");
            input.get_value("phase2", phase2, 1., -10., 10., "last orbital phase of interest");
        }
        std::string store;
        input.get_value("store",  store, "profiles.mol", "molly file for storage of profiles");
        double expose;
        input.get_value("expose", expose, 0.001, 0., 0.5, "exposure length (orbital phase)");
        int ndiv;
        input.get_value("ndiv", ndiv, 1, 1, 1000, "number of sub-divisions for exposure smearing");
        double v1;
        input.get_value("v1", v1, -200., -10000., 10000., "velocity of left-edge of profile array (km/s)");
        double v2;
        input.get_value("v2", v2, 200., -10000., 10000., "velocity of right-edge of profile array (km/s)");
        int nbin;
        input.get_value("nbin", nbin, 100, 1, 10000, "number of velocity bins");
        double fwhm;
        input.get_value("fwhm", fwhm, 0.01, 1.e-5, 1000., "FWHM for blurring (km/s)");
        int nfine;
        input.get_value("nfine", nfine, 5, 1, 200, "sub-division factor for blurring");
        bool align;
        input.get_value("align", align, false, "align on mean motion of secondary star?");

        input.save();

        Lcurve::LDC ldc2 = model.get_ldc2();
        double r1, r2;
        model.get_r1r2(r1, r2);

        // Generate array over secondary star's face.
        Subs::Buffer1D<Lcurve::Point> star2;
        Lcurve::set_star_grid(model, Roche::SECONDARY, true, star2);

        // modify the gravity darkening coefficient to allow for two
        // possibilities the gravity darkening is implemented by
        // modifying the temperature and then calculating the flux in
        // a BB approx. The 'bolometric' method does this directly;
        // the 'filter integrated' method modifies the exponent used
        // to give the desired behaviour of flux with gravity.
        const double GDCBOL2 = model.gdark_bolom2 ? model.gravity_dark2 :
            model.gravity_dark2 / Subs::dlpdlt(model.wavelength, model.t2);

        // Calculate element surface brightnesses, applying gravity
        // darkening, i.e. assuming constant equivalent width
        Subs::Buffer1D<float> emm2, abs2(star2.size());
        Lcurve::set_star_emission(0., model.height_disc, star2, emm2);
        for(int i=0; i<star2.size(); i++){
            double temp = model.t2*pow(double(star2[i].gravity),GDCBOL2);
            abs2[i] = aline*Subs::planck(model.wavelength, temp);
        }

        if(model.add_disc && model.opaque){

            // account for eclipse by the disc
            std::vector<std::pair<double,double> > eclipses;

            // note that the inner radius of the disc is set equal to
            // that of the white dwarf if rdisc1 <= 0 while the outer
            // disc is set equal to the spot radius
            double rdisc1 = model.rdisc1 > 0. ? model.rdisc1 : r1;
            double rdisc2 = model.rdisc2 > 0. ? model.rdisc2 : model.radius_spot;

            // Apply eclipse by disc to star 2
            for(int i=0; i<star2.size(); i++){
                eclipses =  Roche::disc_eclipse(model.iangle, rdisc1, rdisc2,
                                                model.beta_disc, model.height_disc,
                                                star2[i].posn);
                for(size_t j=0; j<eclipses.size(); j++)
                    star2[i].eclipse.push_back(eclipses[j]);
            }
        }

        // Now compute line profile
        const int NFINE = nbin*nfine;
        Subs::Buffer1D<double> fine(NFINE);

        // blurr array stuff
        const int nblurr = int(3.*nfine*fwhm/((v2-v1)/nbin));
        const int nbtot  = 2*nblurr+1;
        float blurr[nbtot], sigma = fwhm/Constants::EFAC;
        float efac = Subs::sqr((v2-v1)/nbin/nfine/sigma)/2.;
        double sum = 0.;
        for(int k = -nblurr; k<= nblurr; k++)
            sum += (blurr[nblurr+k] = exp(-efac*k*k));
        for(int k=0; k< nbtot; k++)
            blurr[k] /= sum;

        // Make stars orbit around centre of mass of system
        const Subs::Vec3 cofm(model.q/(1.+model.q),0.,0.);
        Subs::Vec3 r, earth;
        float amount;
        double veloc, phi, wgt, mu, vk;
        int iadd;

        // Open molly file
        std::ofstream fout(store.c_str(), std::ios::out | std::ios::binary);
        if(!fout)
            throw Lcurve::Lcurve_Error("Failed to open " + store + ". Does it exist already?");
        Colly::molly spectrum;

        // Build up basic header
        Subs::Header head;
        std::vector<std::string> original;
        head.set("Object", new Subs::Hstring("lcurve profile"));
        head.set("q", new Subs::Hdouble(model.q));
        head.set("rangle", new Subs::Hdouble(model.iangle));
        head.set("r2", new Subs::Hdouble(r2));
        head.set("Xtra", new Subs::Hdirectory("Extra information on spectrum"));
        head.set("Xtra.FCODE", new Subs::Hint(4 ,"molly format code"));
        head.set("Xtra.UNITS", new Subs::Hstring("MILLIJANSKYS", "Flux units"));
        head.set("Xtra.NPIX",  new Subs::Hint(nbin,"Number of pixels"));
        head.set("Xtra.NARC",  new Subs::Hint(-2, "Number of arc coefficients"));
        double arc[2];
        arc[1] = nbin*log(1.+(v2-v1)/nbin/(Constants::C/1.e3));
        arc[0] = log(5000.) - arc[1]*(1.+nbin)/(2*nbin);
        head.set("Xtra.ARC",  new Subs::Hdvector(std::vector<double>(arc,arc+2),"Arc coefficients"));
        float *profile = new float[nbin];

        // Compute k1+k2
        double vscale = model.velocity_scale*sin(Subs::deg2rad(model.iangle));
        Subs::Format form(12);

        for(int np=0; np<nphase; np++){

            head.set("Record", new Subs::Hint(np+1));

            for(int i=0; i<NFINE; i++)
                fine[i] = 0.;

            if(nphase > 1)
                phase = phase1 + (phase2-phase1)*np/double(nphase-1);
            head.set("UTC",   new Subs::Hdouble(phase));
            head.set("Dwell", new Subs::Hfloat(expose));
            head.set("Orbital phase", new Subs::Hdouble(phase));

            double ommu = 1.-model.q/(1.+model.q);
            if(align){
                phi = phase - floor(phase);
                double sinp = sin(Constants::TWOPI*phi);
                vk = ommu*sinp;
            }else{
                vk = 0.;
            }

            double wmean  = 0., vwmean = 0.;
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
                phi = phi - floor(phi);
                earth = Roche::set_earth(model.iangle, phi);
                double cosp = cos(Constants::TWOPI*phi);
                double sinp = sin(Constants::TWOPI*phi);

                for(int i=0; i<star2.size(); i++){
                    mu = Subs::dot(earth, star2[i].dirn);

                    if(mu > 0. && star2[i].visible(phi)){
                        r = star2[i].posn - cofm;
                        amount =  wgt*mu*star2[i].area*(emm2[i]*(1.-exp((befac-1./mu)*tau))/(1.-befac*mu)-abs2[i]*ldc2.imu(mu));
                        veloc = vscale*(model.spin2*(cosp*r.y() + sinp*r.x()) + ommu*(1.-model.spin2)*sinp - vk);
                        wmean += amount;
                        vwmean += amount*veloc;
                        iadd = int(floor(NFINE*(veloc-v1)/(v2-v1)));
                        if(iadd >= 0 && iadd < NFINE)
                            fine[iadd] += amount;
                    }
                }
            }

            std::cout << form(phase) << " " << form(vwmean/wmean) << std::endl;

            // blurr
            for(int nb=0; nb<nbin; nb++){
                sum = 0.;
                int off = nfine*nb;
                for(int l=off; l<off+nfine; l++){
                    for(int k=0, j=l-nblurr; k<nbtot; k++, j++)
                        if(j >= 0 && j < NFINE) sum += blurr[k]*fine[j];
                }
                profile[nb]  = sum/std::max(1,ndiv-1)/((v2-v1)/nbin);
            }

            // write spectrum to molly file
            Colly::write_molly_head(fout, head, original, false);

            // write out data
            int nbyte = nbin*sizeof(float);
            fout.write((char*)&nbyte, sizeof(nbyte));
            fout.write((char*)profile, nbyte);
            fout.write((char*)&nbyte,sizeof(nbyte));
        }
        delete[] profile;
    }

    catch(const std::string& err){
        std::cerr << err << std::endl;
        exit(EXIT_FAILURE);
    }
}
