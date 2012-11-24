#include "trm_lcurve.h"

/** 
 * Computes position vector and direction cosine of point on upper surface of
 * disc
 * \param r radius
 * \param theta angle, radians measured from red star in diretcion of motion
 * \param beta exponent of flaring
 * \param height height at unit radius in disc (even if it does not exist)
 * \param posn position of point
 * \param dirn direction cosine of point
 */
void Lcurve::pos_disc(double r, double theta, double beta, double height, Subs::Vec3& posn, Subs::Vec3& dirn){

  static double old_theta = -1.e20, cost, sint;
  if(theta != old_theta){
    cost = cos(theta);
    sint = sin(theta);
    old_theta = theta;
  }

  double h = height*pow(r,beta);
  posn.set(r*cost,r*sint,h);
  double grad = beta*h/r;
  double z = 1./sqrt(1+grad*grad);
  double l = grad*z;
  dirn.set(-l*cost,-l*sint,z);
}
