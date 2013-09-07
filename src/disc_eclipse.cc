#include "trm/lcurve.h"
#include "trm/roche.h"

/** 
 * Checks for eclipse by flared disc.
 * \param iangle inclination angle, degrees
 * \param phase orbital phase
 * \param rdisc1 inner disc  radius, units of separation
 * \param rdisc2 outer disc radius, units of separation
 * \param beta exponent of flaring, so that the height scales as r**beta. beta should be > 1
 * \param height disc height at unit radius in disc (even if it does not exist)
 * \param posn position of point in question.
 * \return 
 */
bool Lcurve::disc_eclipse(double iangle, double phase, double rdisc1, double rdisc2, double beta, double height, const Subs::Vec3& posn){

  double h1 = height*pow(rdisc1, beta);
  double h2 = height*pow(rdisc2, beta);

  Subs::Vec3 earth = Roche::set_earth(iangle, phase);

  if(earth.z() == 0.){
  
    // Special case of i=90
    if(posn.z() > h2 || posn.z() < h2){
      return false;
    }else{
      double fac = posn.sqr() - Subs::sqr(rdisc2);
      if(fac < 0.) return false;
      double bfac = Subs::dot(earth,posn);
      if(sqrt(fac) > bfac)
	return true;
      else
	return false;
    }

  }else{

    // Work out where line of sight crosses cylinder at edge of disc
    Subs::Vec3 flat_earth = earth;
    flat_earth.z() = 0;
    flat_earth.unit();
    Subs::Vec3 flat_posn = posn;
    flat_posn.z() = 0;

    double bquad  = Subs::dot(flat_earth, flat_posn);
    double fac    = Subs::sqr(bquad) - (flat_posn.sqr() - 0.9999*Subs::sqr(rdisc2));
    if(fac < 0.) return false;
    fac = sqrt(fac);
    if(fac < bquad) return false;

    double lam = -bquad - fac;
    if(lam > 0){
      double zh = posn.z() + lam*earth.z();
      if(zh > -h2 && zh < h2)
	return true;
    }

    lam = -bquad + fac;
    if(lam > 0){
      double zh = posn.z() + lam*earth.z();
      if(zh > -h2 && zh < h2)
	return true;
    }

    // Now work out where line of sight goes through planes at z = -h2 and +h1

    Subs::Vec3 cross;
    double lam1 = (-h2-posn.z())/earth.z();
    if(lam1 > 0.){
      cross = posn + lam1*earth;
      cross.z() = 0.;
      if(cross.sqr() < Subs::sqr(rdisc2))
	return true;
    }

    double lam2;
    if(Subs::sqr(posn.x()) + Subs::sqr(posn.y()) > 0.9999*Subs::sqr(rdisc2))
      lam2 = (+h2-posn.z())/earth.z();
    else
      lam2 = (+h1-posn.z())/earth.z();
    if(lam2 > 0.){
      cross = posn + lam2*earth;
      cross.z() = 0.;
      if(cross.sqr() < Subs::sqr(rdisc2))
	return true;
    }

  }
  return false;
}

/** 
 * Checks for eclipse by flared disc for points at its surface
 * \param iangle inclination angle
 * \param phase orbital phase
 * \param rdisc1 inner disc  radius
 * \param rdisc2 outer disc radius
 * \param beta exponent of flaring
 * \param height height at unit radius in disc (even if it does not exist)
 * \param posn position of point
 * \return 
 */
bool Lcurve::disc_surface_eclipse(double iangle, double phase, double rdisc1, double rdisc2, double beta, double height, const Subs::Vec3& posn){

  double h2 = height*pow(rdisc2, beta);

  Subs::Vec3 earth = Roche::set_earth(iangle, phase);

  if(earth.z() == 0.) return true;
    
  Subs::Vec3 cross;
    
  // Work out where line of sight goes through plane at z=h2
  double lam = (h2-posn.z())/earth.z();
  if(lam > 0.){
    cross = posn + lam*earth;
    cross.z() = 0.;
    if(cross.sqr() > Subs::sqr(rdisc2))
      return true;
  }
  return false;
}
