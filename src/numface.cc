#include "trm/constants.h"
#include "trm/subs.h"
#include "trm/lcurve.h"

/** Calculates the number of elements needed given the number of latititude elements etc. See set_grid for
 * information on how the grid is defined.
 * \param nlat   basic number of latitude strips
 * \param infill whether to infill with a fine grid
 * \param thelo  theta limit at which fine grid starts
 * \param thehi  theta limit at which fine grid stops
 * \param nlatfill extra number of points between latitude strips
 * \param nlngfill extra number of points between longitude points
 */
int Lcurve::numface(int nlat, bool infill, double thelo, double thehi, int nlatfill, int nlngfill){

  int nface = 0;
  double theta, dphi, dtheta = Constants::PI/nlat;

  if(infill){
      int nl1 = int(std::ceil(thelo/dtheta));
      for(int i=0; i<nl1; i++){
	  theta  = thelo*(i+0.5)/nl1;
	  dphi   = dtheta/sin(theta);
	  nface += std::max(16,int(2*Constants::PI/dphi));
      }
      int nl2 = int(std::ceil((1+nlatfill)*(thehi-thelo)/dtheta));
      for(int i=0; i<nl2; i++){
	  theta  = thelo + (thehi-thelo)*(i+0.5)/nl2;
	  dphi   = dtheta/sin(theta)/(1+nlngfill);
	  nface += std::max(8,int(Constants::PI/dphi));
      }
      int nl3 = int(std::ceil((thehi-thelo)/dtheta));
      for(int i=0; i<nl3; i++){
	  theta  = thelo + (thehi-thelo)*(i+0.5)/nl3;
	  dphi   = dtheta/sin(theta);
	  nface += std::max(8,int(Constants::PI/dphi));
      }
      int nl4 = int(std::ceil((Constants::PI-thehi)/dtheta));
      for(int i=0; i<nl4; i++){
	  theta  = thehi + (Constants::PI-thehi)*(i+0.5)/nl4;
	  dphi   = dtheta/sin(theta);
	  nface += std::max(16,int(2*Constants::PI/dphi));
      }
  }else{
    nlat = int(std::ceil(Constants::PI/dtheta));
      for(int i=0; i<nlat; i++){
	  theta  = Constants::PI*(i+0.5)/nlat;
	  dphi   = dtheta/sin(theta);
	  nface += std::max(16,int(2*Constants::PI/dphi));
      }
  }
  return nface;
}
