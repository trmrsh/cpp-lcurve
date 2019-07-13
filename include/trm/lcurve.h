#ifndef TRM_LCURVE
#define TRM_LCURVE

#include <cmath>
#include <string>
#include "trm/subs.h"
#include "trm/array1d.h"
#include "trm/vec3.h"
#include "trm/roche.h"

//! Lcurve namespace. Stuff to do with light curve modelling.

namespace Lcurve {

  enum Star {STAR1, STAR2};

  //! Name of environment variable containing name of defaults directory
  const std::string LCURVE_ENV = "LCURVE_ENV";

  //! Default name of defaults directory
  const std::string LCURVE_DIR = ".lcurve";

  /** Lcurve::Lcurve_Error is the error class for the Roche programs.
   * It is inherited from the standard string class.
   */
  class Lcurve_Error : public std::string {
  public:

    //! Default constructor
    Lcurve_Error() : std::string() {}

    //! Constructor storing a message
    Lcurve_Error(const std::string& err) : std::string(err) {}
  };

  //! Structure defining a single element
  /** This defines the position, area, direction, gravity and brightness of an element
   * and also any phases during which it is eclipsed.
   */
  struct Point {

    typedef std::vector<std::pair<double,double> > etype;

    //! Default constructor
    Point() : posn(), dirn(), area(0.), gravity(1.), eclipse(), flux(0.) {}

    //! Constructor
    Point(const Subs::Vec3& posn_, const Subs::Vec3& dirn_, double area_, double gravity_, const etype& eclipse) :
      posn(posn_), dirn(dirn_), area(area_), gravity(gravity_), eclipse(eclipse), flux(0.) {}

    //! Position vector of element (units of binary separation)
    Subs::Vec3 posn;

    //! Outward facing direction of element (unit vector)
    Subs::Vec3 dirn;

    //! Area of element (units of binary separation**2)
    float area;

    //! Gravity of element
    float gravity;

    //! Ingress and egress phases of eclipses, if any
    etype eclipse;

    //! Brightness * area
    float flux;

    //! Computes whether a point is visible (else eclipsed)
    bool visible(double phase) const {
      double phi = phase - floor(phase);
      for(size_t i=0; i<eclipse.size(); i++){
	const std::pair<double,double> p = eclipse[i];
	if((phi >= p.first && phi <= p.second) || phi <= p.second-1.0)
	  return false;
      }
      return true;
    }

  };

  struct Ginterp {

    //! Start phase of coarse grid 0 -- 0.5
    double phase1;

    //! End phase of coarse grid 0 -- 0.5
    double phase2;

    //! Scale factor star 1 at phase1
    double scale11;

    //! Scale factor star 1 at 1-phase1
    double scale12;

    //! Scale factor star 2 at -phase2
    double scale21;

    //! Scale factor star 2 at phase2
    double scale22;

    //! Returns scale factor for star 1 at a given phase
    double scale1(double phase) const {
      // assume coarse grid outside -phase1 to + phase1
      double pnorm = phase - floor(phase);
      if(pnorm <= phase1 || pnorm >= 1.-phase1){
	return 1.;
      }else{
	return (scale11*(1.-phase1-pnorm)+scale12*(pnorm-phase1))/(1-2.*phase1);
      }
    }

    //! Returns scale factor for star 2 at a given phase
    double scale2(double phase) const {
      double pnorm = phase - floor(phase);
      if(pnorm >= phase2 && pnorm <= 1.-phase2){
	return 1.;
      }else if(pnorm < 0.5){
	return (scale22*(phase2-pnorm)+scale21*(pnorm+phase2))/(2.*phase2);
      }else{
	return (scale21*(1.+phase2-pnorm)+scale22*(pnorm-1.+phase2))/(2.*phase2);
      }
    }

    //! Returns integer type to represent situation we are in
    int type(double phase) const {
      double pnorm = phase - floor(phase);
      if(pnorm <= phase1 || pnorm >= 1.-phase1){
	// coarse grid for star 2, fine for star 1
	return 1;
      }else if((pnorm > phase1 && pnorm < phase2) || (pnorm > 1.-phase2 && pnorm < 1.-phase1)){
	// coarse grid for both stars
	return 2;
      }else{
	// coarse grid for star 1, fine for star 2
	return 3;
      }
    }
  };

  //! Dummy ASCII input operator to allow use of Buffer1D
  std::istream& operator>>(std::istream& s, Point& p);

  //! Dummy ASCII output operator to allow use of Buffer1D
  std::ostream& operator<<(std::ostream& s, const Point& p);

  //! Physical parameter structure
  /** Holds basic information for a physical parameter which are its value,
   * range for varying it, step size for derivative computation, whether it
   * is variable or not, and, optionally, whether it is defined or not
   */
  struct Pparam {

    //! Default constructor
  Pparam() : value(0), range(0), dstep(0), vary(false), defined(false) {}

    //! Constructor from a string
    /** Sets the values from a string of the form "0.12405 0.001"
     * giving the value and the step size. 0 is interpreted as holding
     * a parameter fixed.
     */
    Pparam(const std::string& entry){
      std::istringstream istr(entry);
      istr >> value >> range >> dstep >> vary;
      if(!istr)
        throw Lcurve::Lcurve_Error(std::string("Pparam: could not read entry = ") + entry);
      istr >> defined;
      if(!istr) defined = true;
    }

    //! Implicit conversion to a double by just returning the value
    operator double () const { return value; }

    //! The value of the parameter
    double value;

    //! The value of the range over which to vary it
    double range;

    //! The value of the step size for derivative computation
    double dstep;

    //! Whether the parameter varies
    bool vary;

    //! Whether the parameter has been defined
    bool defined;

  };


  //! ASCII input operator for a physical parameter
  std::ostream& operator<<(std::ostream& s, const Pparam& p);

  //! Lim darkening class
  class LDC {

  public:

    enum LDCtype {POLY, CLARET};

    //! Default. Sets all to zero.
    LDC() : ldc1(0.), ldc2(0.), ldc3(0.), ldc4(0.), mucrit(0.), ltype(POLY) {}

    //! Standard constructor
    LDC(double ldc1, double ldc2, double ldc3, double ldc4, double mucrit, LDCtype ltype) :
      ldc1(ldc1), ldc2(ldc2), ldc3(ldc3), ldc4(ldc4), mucrit(mucrit), ltype(ltype) {}

    //! Computes I(mu)
    double imu(double mu) const {
      if(mu <= 0){
	return 0.;
      }else{
	mu = std::min(mu, 1.);
	double ommu = 1.-mu, im = 1.;
	if(this->ltype == POLY){
	  im -= ommu*(this->ldc1+ommu*(this->ldc2+ommu*(this->ldc3+ommu*this->ldc4)));
	}else if(this->ltype == CLARET){
	  im -= this->ldc1+this->ldc2+this->ldc3+this->ldc4;
	  double msq = sqrt(mu);
	  im += msq*(this->ldc1+msq*(this->ldc2+msq*(this->ldc3+msq*this->ldc4)));
	}
	return im;
      }
    }

    //! To help applying mucrit
    bool see(double mu) const {return mu > this->mucrit;}

  private:

    double  ldc1;
    double  ldc2;
    double  ldc3;
    double  ldc4;
    double  mucrit;
    LDCtype ltype;

  };

  //! Model structure
  /** Defines the model to be used and which parameters are to be
   * varied in the fit.  The order of the paremeters here defines the
   * order in which they must occur in any parameter vector to be fed
   * to a generic minimisation routine such as amoeba.  Star 1 is
   * taken to be spherical. Star 2 can be tidally distorted.
   */
  struct Model {

    //! Constructor from a file
    Model(const std::string& file);

    //! Number of variable parameters
    int nvary() const;

    //! Uses a vector of numbers to update the variable parameter values
    void set_param(const Subs::Array1D<double>& vpar);

    //! Check legality of proposed set of parameters
    bool is_not_legal(const Subs::Array1D<double>& vpar) const;

    //! Returns the current variable parameter values
    Subs::Array1D<double> get_param() const;

    //! Returns the ranges of the variable parameters
    Subs::Array1D<double> get_range() const;

    //! Returns the step sizes of the variable parameters
    Subs::Array1D<double> get_dstep() const;

    //! Writes out model to an ASCII file
    void wrasc(const std::string& file) const;

    //! Returns the name of variable i 
    std::string get_name(int i) const;

    // Physical parameters (which can be varied)

    //! Mass ratio = M2 / M1
    Pparam q; 

    //! Inclination angle, degrees
    Pparam iangle; 

    //! Radius of star 1 / separation
    Pparam r1;

    //! Radius of star 1 / separation
    Pparam r2;

    //! Phase of third contact
    Pparam cphi3;

    //! Phase of fourth contact
    Pparam cphi4;

    //! Spin/orbital for star 1
    Pparam spin1;

    //! Spin/orbital for star 2
    Pparam spin2;

    //! Temperature in K star 1
    Pparam t1;

    //! Temperature in K star 2
    Pparam t2;

    //! Star 1, LDC #1
    Pparam ldc1_1;

    //! Star 1, LDC #2
    Pparam ldc1_2;

    //! Star 1, LDC #3
    Pparam ldc1_3;

    //! Star 1, LDC #4
    Pparam ldc1_4;

    //! Star 2, LDC #1
    Pparam ldc2_1;

    //! Star 2, LDC #2
    Pparam ldc2_2;

    //! Star 2, LDC #3
    Pparam ldc2_3;

    //! Star 2, LDC #4
    Pparam ldc2_4;

    //! Velocity scale  = v1+v2, unprojected
    Pparam velocity_scale;

    //! Doppler beaming 3-alpha that multiplies -v/c for star 1
    Pparam beam_factor1;

    //! Doppler beaming 3-alpha that multiplies -v/c for star 2
    Pparam beam_factor2;

    //! Zero point of ephemeris
    Pparam t0;

    //! Period
    Pparam period;

    //! Pdot term (third coeff of ephemeris)
    Pparam pdot;

    //! Roemer delay
    Pparam deltat;

    //! Gravity darkening coefficient 1
    Pparam gravity_dark1;

    //! Gravity darkening coefficient 2
    Pparam gravity_dark2;

    //! Absorbed fraction
    Pparam absorb;

    //! Slope (fudge factor really)
    Pparam slope;

    //! Quadratic term
    Pparam quad;

    //! Cubic term
    Pparam cube;

    //! Third light fraction
    Pparam third;

    //! Inner radius of disc
    Pparam rdisc1;

    //! Outer radius of disc
    Pparam rdisc2;

    //! Height of disc at r = 1
    Pparam height_disc;

    //! Radial exponent of disc height
    Pparam beta_disc;

    //! Temperature of disc
    Pparam temp_disc;

    //! Radial exponent of disc temperature
    Pparam texp_disc;

    //! Linear limb darkening of disc
    Pparam lin_limb_disc;

    //! Quadratic limb darkening of disc
    Pparam quad_limb_disc;

    //! Temperature of disc edge
    Pparam temp_edge;

    //! Reprocessing parameter at edge
    Pparam absorb_edge;

    //! Distance of bright spot from white dwarf
    Pparam radius_spot;

    //! Length scale of bright spot
    Pparam length_spot;

    //! Height of spot
    Pparam height_spot;

    //! Power law exponent for spot brightness
    Pparam expon_spot;

    //! Power law exponent inside exponent of spot brightness
    Pparam epow_spot;

    //! Angle of spot
    Pparam angle_spot;

    //! Extra angle of spot beaming
    Pparam yaw_spot;

    //! Peak spot temperature
    Pparam temp_spot;

    //! Tilt of spot
    Pparam tilt_spot;

    //! Constant fraction of spot
    Pparam cfrac_spot;

    //! longitude (in direction of orbit relative to other star) of spot
    Pparam stsp11_long;

    //! latitude of spot (degrees)
    Pparam stsp11_lat;

    //! fwhm of spot (degrees)
    Pparam stsp11_fwhm;

    //! central temperature of spot
    Pparam stsp11_tcen;

    //! longitude (in direction of orbit relative to other star) of spot
    Pparam stsp21_long;

    //! latitude of spot (degrees)
    Pparam stsp21_lat;

    //! fwhm of spot (degrees)
    Pparam stsp21_fwhm;

    //! central temperature of spot
    Pparam stsp21_tcen;

    // Computational ones

    //! Accuracy in phase for Roche lobe computations
    double delta_phase;

    //! Number of latitude strips on star 1, fine grid
    int nlat1f;

    //! Number of latitude strips on star 2, fine grid
    int nlat2f;

    //! Number of latitude strips on star 1, coarse grid
    int nlat1c;

    //! Number of latitude strips on star 2, coarse grid
    int nlat2c;

    //! Use genuine North pole for grid north pole rather than substellar points
    bool npole;

    //! Extra latitude points to finely sample star 2 along track of star 1
    int nlatfill;

    //! Extra longitude points to finely sample star 2 along track of star 1
    int nlngfill;

    //! Extra amount to add to fine latitude strip, degrees
    double lfudge;

    //! Lowest latitude for fine strip
    double llo;

    //! Highest latitude for fine strip
    double lhi;

    //! Start phase (0 to 0.5) of coarse grid
    double phase1;

    //! End phase (0 to 0.5) of coarse grid
    double phase2;

    //! Wavelength
    double wavelength;

    //! Account for Roche distortion of primary
    bool roche1;

    //! Account for Roche distortion of secondary
    bool roche2;

    //! Account for eclipse of star 1 by star 2
    bool eclipse1;

    //! Account for eclipse of star 2 by star 1
    bool eclipse2;

    //! Account for gravitational lensing by star 1
    bool glens1;

    //! Use radii (rather than contact phases)
    bool use_radii;

    //! True period in days (in case other is actually phase)
    double tperiod;

    //! Gravity darkening bolometric or filter integreted star 1
    bool gdark_bolom1;

    //! Gravity darkening bolometric or filter integreted star 2
    bool gdark_bolom2;

    //! Critical mu value for star1
    double mucrit1;

    //! Critical mu value for star2
    double mucrit2;

    //! Type of limb darkening star 1
    LDC::LDCtype limb1;

    //! Type of limb darkening star 2
    LDC::LDCtype limb2;

    //! Assume any flux not absorbed is truly reflected (or just ignore it)
    bool mirror;

    //! Add a disc or not
    bool add_disc;

    //! Number of radial strips for disc
    int nrad;

    //! Is disc opaque or not?
    bool opaque;

    //! Add a spot or not
    bool add_spot;

    //! Number of elements in spot
    int nspot;

    //! Individual scaling of components or not
    bool iscale;

    //! Returns relative radii, accounting for method of parameterising them
    void get_r1r2(double& rr1, double& rr2) const {
      if(this->use_radii){
	rr1 = this->r1;
	rr2 = this->r2;
      }else{
	double sini  = std::sin(Subs::deg2rad(this->iangle)); 
	double r2pr1 = std::sqrt(1.-std::pow(sini*std::cos(2.*M_PI*this->cphi4.value),2));
	double r2mr1 = std::sqrt(1.-std::pow(sini*std::cos(2.*M_PI*this->cphi3.value),2));
	rr1 = (r2pr1 - r2mr1)/2.;
	rr2 = (r2pr1 + r2mr1)/2.;
      }
      return;
    }
    
    //! Returns limb darkening struct for star 1
    LDC get_ldc1() const {
      return LDC(this->ldc1_1, this->ldc1_2, this->ldc1_3, this->ldc1_4, this->mucrit1, this->limb1);
    }

    //! Returns limb darkening struct for star 2
    LDC get_ldc2() const {
      return LDC(this->ldc2_1, this->ldc2_2, this->ldc2_3, this->ldc2_4, this->mucrit2, this->limb2);
    }

  };

  //! ASCII output operator for a Model
  std::ostream& operator<<(std::ostream& s, const Model& model);

  //! Holds all the data for a single point of a light curve
  struct Datum {

    //! The time
    double time;

    //! The exposure length in the same units as the time
    double expose;

    //! The flux
    double flux;

    //! The uncertainty on the flux in the same units
    double ferr;

    //! Weight factor for calculating goodness of fit
    double weight;

    //! Factor to split up data points to allow for finite exposures
    int ndiv;

  };

  //! ASCII input of a Datum (expects time expose flux ferr weight ndiv)
  std::istream& operator>>(std::istream& s, Datum& datum);

  //! ASCII output of a Datum (expects time expose flux ferr weight ndiv)
  std::ostream& operator<<(std::ostream& s, const Datum& datum);

  // Holds a light curve
  class Data : public std::vector<Datum> {

  public:

    //! Default constructor
    Data() : std::vector<Datum>() {}

    //! Constructor with pre-defined size
    Data(int n) : std::vector<Datum>(n) {}

    //! Constructor from a file
    Data(const std::string& file);

    //! Writes to an ASCII file
    void wrasc(const std::string& file) const;

    //! Reads from an ASCII file
    void rasc(const std::string& file);

  };

  //! Function object to compute chi**2

  /** This is inherited from the abstract class Subs::Afunc to define the basic
   * functionality needed for it to be used inside minimisation routines like amoeba
   */

  class Fobj : public Subs::Afunc {

  public:

    //! Total number of calls to Lcurve::chisq
    static int neval;

    //! Minimum chi**2 ever encountered
    static double chisq_min;

    //! Scale factors for minimum chi**2
    static Subs::Buffer1D<double> scale_min;

    //! Constructor; stores the data needed to compute chi**2
    Fobj(const Model& model, const Data& data) : model(model), data(data) {}

    //! Function call operator (overload)
    double operator()(const Subs::Array1D<double>& vpar);

    //! Computes chi-squared
    double chisq();

    //! Return fit values
    double operator[](int n) const;

  private:

    Model model;
    const Data&  data;
    Subs::Array1D<double> fit;

  };

  //! Computes position and direction on a disc
  void pos_disc(double r, double theta, double beta, double height, Subs::Vec3& posn, Subs::Vec3& dirn);

  //! Computes the number of faces given the number of latitude strips to be used
  int numface(int nlat, bool infill, double thelo, double thehi, int nlatfill, int nlngfill);

  //! Computes elements over the primary or secondary star
  void set_star_grid(const Model& mdl, Roche::STAR which_star, bool fine, Subs::Buffer1D<Lcurve::Point>& star);

  void add_faces(Subs::Buffer1D<Lcurve::Point>& star, int& nface, double tlo, double thi, double dtheta, int nlatfill, int nlngfill, 
		 bool npole, Roche::STAR which_star, double q, double iangle, double r1, double r2, double rref1, double rref2, 
		 bool roche1, bool roche2, double spin1, double spin2, bool eclipse, double gref, double pref1, double pref2, 
		 double ffac1, double ffac2, double delta);

  //! Computes elements over a disc
  void set_disc_grid(const Model& mdl, Subs::Buffer1D<Lcurve::Point>& disc);

  //! Bright spot elelemts
  void set_bright_spot_grid(const Model& mdl, Subs::Buffer1D<Lcurve::Point>& spot);
  
  //! Computes elements at rim of disc
  void set_disc_edge(const Model& mdl, bool outer,
		     Subs::Buffer1D<Lcurve::Point>& edge,
		     bool visual=true);

  //! Sets the continuum element contributions
  void set_star_continuum(const Model& mdl,
			  Subs::Buffer1D<Lcurve::Point>& star1,
			  Subs::Buffer1D<Lcurve::Point>& star2);

  //! Sets the disc continuum brightness.
  void set_disc_continuum(double rdisc, double tdisc, double texp,
			  double wave, Subs::Buffer1D<Lcurve::Point>& disc);

  //! Sets the disc edge continuum brightness.
  void set_edge_continuum(double tedge, double r2, double t2,
			  double absorb, double wave,
			  Subs::Buffer1D<Lcurve::Point>& edge);

  //! Sets the emission line brightness
  void set_star_emission(double limb2, double hbyr,
			 const Subs::Buffer1D<Lcurve::Point>& star2,
			 Subs::Buffer1D<float>& bright2);

  //! Computes the flux at a given phase
  double comp_light(double iangle, const LDC& ldc1, const LDC& ldc2,
		    double lin_limb_disc, double quad_limb_disc, 
		    double phase, double expose, int ndiv, double q,
		    double beam_factor1, double beam_factor2, 
		    double spin1, double spin2, float vscale, bool glens1,
		    double rlens1, const Ginterp& gint, 
		    const Subs::Buffer1D<Lcurve::Point>& star1f,
		    const Subs::Buffer1D<Lcurve::Point>& star2f,
		    const Subs::Buffer1D<Lcurve::Point>& star1c,
		    const Subs::Buffer1D<Lcurve::Point>& star2c,
		    const Subs::Buffer1D<Lcurve::Point>& disc,
		    const Subs::Buffer1D<Lcurve::Point>& edge,
		    const Subs::Buffer1D<Lcurve::Point>& spot);

  //! Computes flux from star 1 only
  double comp_star1(double iangle, const LDC& ldc1, double phase,
		    double expose, int ndiv, double q, double beam_factor1,
		    float vscale,  const Ginterp& gint, 
		    const Subs::Buffer1D<Lcurve::Point>& star1f,
		    const Subs::Buffer1D<Lcurve::Point>& star1c);

  //! Computes flux from star 2 only
  double comp_star2(double iangle, const LDC& ldc2, double phase,
		    double expose, int ndiv, double q, double beam_factor2,
		    float vscale, bool glens1, double rlens1,
		    const Ginterp& gint,
		    const Subs::Buffer1D<Lcurve::Point>& star2f,
		    const Subs::Buffer1D<Lcurve::Point>& star2c);

  //! Computes flux from disc
  double comp_disc(double iangle, double lin_limb_disc, double quad_limb_disc,
		   double phase, double expose, int ndiv, double q,
		   float vscale, const Subs::Buffer1D<Lcurve::Point>& disc);

  //! Computes flux from disc edge
  double comp_edge(double iangle, double lin_limb_disc, double quad_limb_disc,
		   double phase, double expose, int ndiv, double q,
		   float vscale, const Subs::Buffer1D<Lcurve::Point>& edge);

  //! Compute flux from spot
  double comp_spot(double iangle,double phase, double expose, int ndiv,
		   double q, float vscale, 
		   const Subs::Buffer1D<Lcurve::Point>& spot);

  //! Convenience routine 
  void star_eclipse(double q, double r, double spin, double ffac,
		    double iangle, const Subs::Vec3& posn, double delta, 
		    bool roche, Roche::STAR star, Point::etype& eclipses);

  //! Computes eclipse by a flared disc
  bool disc_eclipse(double iangle, double phase, double rdisc1,
		    double rdisc2, double beta, double height,
		    const Subs::Vec3& posn);

  //! Computes eclipse by a flared disc for points 
  bool disc_surface_eclipse(double iangle, double phase, double rdisc1,
			    double rdisc2, double beta, double height,
			    const Subs::Vec3& posn);
  
  //! Computes an entire light curve corresponding to a given data set.
  void light_curve_comp(const Lcurve::Model& model, const Lcurve::Data& data,
			bool scale, bool info, Subs::Buffer1D<double>& sfac, 
			Subs::Array1D<double>& calc, double& wdwarf,
			double& chisq, double& wnok);

  //! Re-scales a fit to minimise chi**2
  double re_scale(const Lcurve::Data& data, Subs::Array1D<double>& fit,
		  double& chisq, double& wnok);

};

#endif
