#include <map>
#include "trm/array1d.h"
#include "trm/lcurve.h"
#include "trm/roche.h"
#include "trm/format.h"

/** Outputs the value and step size of a physical parameter
 */
std::ostream& Lcurve::operator<<(std::ostream& s, const Pparam& p){
    static Subs::Format vform(17), sform(8);
    s << vform(p.value) << " " << sform(p.range) << " " <<
        sform(p.dstep) << " " << p.vary << " " << p.defined;
    return s;
}

/* Constructor from a parameter file
 * Reads in a file which defines both the physical and computational parameters
 * \param file the file of parameter values
 */
Lcurve::Model::Model(const std::string& file) {

    // List of parameter names to expect, both physical and computational used
    // for checking that we actually have read in all the one expected
    std::map<std::string,bool> names;
    typedef std::map<std::string,bool>::iterator NIT;

    // Physical 
    names["q"]              = false;
    names["iangle"]         = false;
    names["r1"]             = false;
    names["r2"]             = false;
    names["cphi3"]          = false;
    names["cphi4"]          = false;
    names["spin1"]          = false;
    names["spin2"]          = false;
    names["t1"]             = false;
    names["t2"]             = false;
    names["ldc1_1"]         = false;
    names["ldc1_2"]         = false;
    names["ldc1_3"]         = false;
    names["ldc1_4"]         = false;
    names["ldc2_1"]         = false;
    names["ldc2_2"]         = false;
    names["ldc2_3"]         = false;
    names["ldc2_4"]         = false;
    names["velocity_scale"] = false;
    names["beam_factor1"]   = false;
    names["beam_factor2"]   = false;
    names["t0"]             = false;
    names["period"]         = false;
    names["pdot"]           = false;
    names["deltat"]         = false;
    names["gravity_dark1"]  = false;
    names["gravity_dark2"]  = false;
    names["absorb"]         = false;
    names["slope"]          = false;
    names["quad"]           = false;
    names["cube"]           = false;
    names["third"]          = false;
    names["rdisc1"]         = false;
    names["rdisc2"]         = false;
    names["height_disc"]    = false;
    names["beta_disc"]      = false;
    names["temp_disc"]      = false;
    names["texp_disc"]      = false;
    names["lin_limb_disc"]  = false;
    names["quad_limb_disc"] = false;
    names["temp_edge"]      = false;
    names["absorb_edge"]    = false;
    names["radius_spot"]    = false;
    names["length_spot"]    = false;
    names["height_spot"]    = false;
    names["expon_spot"]     = false;
    names["epow_spot"]      = false;
    names["angle_spot"]     = false;
    names["yaw_spot"]       = false;
    names["temp_spot"]      = false;
    names["tilt_spot"]      = false;
    names["cfrac_spot"]     = false;

    // Star spots:
    //
    // 11 = star 1, spot 1
    // 21 = star 2, spot 1
    //
    names["stsp11_long"] = false;
    names["stsp11_lat"] = false;
    names["stsp11_fwhm"] = false;
    names["stsp11_tcen"] = false;

    names["stsp12_long"] = false;
    names["stsp12_lat"] = false;
    names["stsp12_fwhm"] = false;
    names["stsp12_tcen"] = false;

    names["stsp13_long"] = false;
    names["stsp13_lat"] = false;
    names["stsp13_fwhm"] = false;
    names["stsp13_tcen"] = false;

    names["stsp21_long"] = false;
    names["stsp21_lat"] = false;
    names["stsp21_fwhm"] = false;
    names["stsp21_tcen"] = false;

    names["stsp22_long"] = false;
    names["stsp22_lat"] = false;
    names["stsp22_fwhm"] = false;
    names["stsp22_tcen"] = false;
    
    names["uesp_long1"] = false;
    names["uesp_long2"] = false;
    names["uesp_lathw"] = false;
    names["uesp_taper"] = false;
    names["uesp_temp"] = false;

    // Computational
    names["delta_phase"] = false;
    names["nlat1f"] = false;
    names["nlat2f"] = false;
    names["nlat1c"] = false;
    names["nlat2c"] = false;
    names["npole"] = false;
    names["nlatfill"] = false;
    names["nlngfill"] = false;
    names["lfudge"] = false;
    names["llo"] = false;
    names["lhi"] = false;
    names["phase1"] = false;
    names["phase2"] = false;
    names["wavelength"] = false;
    names["roche1"] = false;
    names["roche2"] = false;
    names["eclipse1"] = false;
    names["eclipse2"] = false;
    names["glens1"] = false;
    names["use_radii"] = false;
    names["tperiod"] = false;
    names["gdark_bolom1"] = false;
    names["gdark_bolom2"] = false;
    names["mucrit1"] = false;
    names["mucrit2"] = false;
    names["limb1"] = false;
    names["limb2"] = false;
    names["mirror"] = false;
    names["add_disc"] = false;
    names["nrad"] = false;
    names["opaque"] = false;
    names["add_spot"] = false;
    names["nspot"] = false;
    names["iscale"] = false;

    // Parameter value pairs
    std::map<std::string, std::string> pv;

    // Read in the parameter values
    std::ifstream fin(file.c_str());
    if(!fin) throw Lcurve_Error("Failed to open " + file +
                                " for starting parameters.");

    const int MAX_LINE = 5000;
    int n = 0;
    char ch;
    std::string param, value;
    while(fin && !fin.eof()){
        n++;
        ch = fin.peek();
        if(fin.eof()){
            std::cerr << "End of file reached." << std::endl;
        }else if(ch == '#' || ch == ' ' || ch == '\t' || ch == '\n'){
            fin.ignore(MAX_LINE, '\n'); // skip to next line
        }else{

            if(fin >> param){

                NIT nit = names.find(param);
                if(nit == names.end())
                    throw Lcurve_Error("Parameter: " + param + " from line "
                                       + Subs::str(n) + " was not recognised.");
                nit->second = true;

                while(fin.get(ch) && ch != '='); // skips up to = sign
                if(!fin)
                    throw Lcurve_Error("Line " + Subs::str(n) +
                                       " starting: " + param + " is invalid");

                getline(fin,value);

                // Work out end of value string defined by first hash
                std::string::size_type ntemp, nhash = 0;
                while((ntemp = value.find('#',nhash)) != std::string::npos &&
                      ntemp > 0 && value[ntemp-1] == '\\'){
                    nhash = ntemp + 1;
                }
                if(ntemp != std::string::npos) value.erase(ntemp);

                // Chop off any space at start and end
                std::string::size_type n1 = value.find_first_not_of(" \t");
                std::string::size_type n2 = value.find_last_not_of(" \t");

                if(n1 != std::string::npos)
                    pv[param] = value.substr(n1,n2-n1+1);
                else
                    throw Lcurve_Error("Line " + Subs::str(n) + " starting: "
                                       + param + " has no values");
            }else if(fin.eof()){
                std::cerr << "End of parameter file reached." << std::endl;
            }else{
                throw Lcurve_Error("Parameter input failure on line "
                                   + Subs::str(n));
            }
        }
    }
    fin.close();

    std::cerr << n << " lines read from " << file << std::endl;

    // Check that everything has been initialised, except for the star spots
    // that do not have to be defined and some recently added parameters
    bool ok = true;
    for(NIT nit=names.begin(); nit!=names.end(); nit++){
        if(!nit->second){
            if(nit->first == "pdot" || nit->first == "third" ||
	       nit->first == "temp_edge" || nit->first == "absorb_edge"){
                std::cerr << nit->first
                          << " was not defined; will be set = 0" << std::endl;
                pv[nit->first] = "0 1.e-10 1.e-10 0 0";
            }else if(nit->first.substr(0,4) != "stsp" && nit->first.substr(0,4) != "uesp"){
                std::cerr << nit->first << " was not defined " << std::endl;
                ok = false;
            }
        }
    }
    if(!ok) throw Lcurve_Error("One or more parameters were not initialised");

    // Initialise physical parameters
    q                = Pparam(pv["q"]);
    iangle           = Pparam(pv["iangle"]);
    r1               = Pparam(pv["r1"]);
    r2               = Pparam(pv["r2"]);
    cphi3            = Pparam(pv["cphi3"]);
    cphi4            = Pparam(pv["cphi4"]);
    spin1            = Pparam(pv["spin1"]);
    spin2            = Pparam(pv["spin2"]);
    t1               = Pparam(pv["t1"]);
    t2               = Pparam(pv["t2"]);
    ldc1_1           = Pparam(pv["ldc1_1"]);
    ldc1_2           = Pparam(pv["ldc1_2"]);
    ldc1_3           = Pparam(pv["ldc1_3"]);
    ldc1_4           = Pparam(pv["ldc1_4"]);
    ldc2_1           = Pparam(pv["ldc2_1"]);
    ldc2_2           = Pparam(pv["ldc2_2"]);
    ldc2_3           = Pparam(pv["ldc2_3"]);
    ldc2_4           = Pparam(pv["ldc2_4"]);
    velocity_scale   = Pparam(pv["velocity_scale"]);
    beam_factor1     = Pparam(pv["beam_factor1"]);
    beam_factor2     = Pparam(pv["beam_factor2"]);
    t0               = Pparam(pv["t0"]);
    period           = Pparam(pv["period"]);
    pdot             = Pparam(pv["pdot"]);
    deltat           = Pparam(pv["deltat"]);
    gravity_dark1    = Pparam(pv["gravity_dark1"]);
    gravity_dark2    = Pparam(pv["gravity_dark2"]);
    absorb           = Pparam(pv["absorb"]);
    slope            = Pparam(pv["slope"]);
    quad             = Pparam(pv["quad"]);
    cube             = Pparam(pv["cube"]);
    third            = Pparam(pv["third"]);
    rdisc1           = Pparam(pv["rdisc1"]);
    rdisc2           = Pparam(pv["rdisc2"]);
    height_disc      = Pparam(pv["height_disc"]);
    beta_disc        = Pparam(pv["beta_disc"]);
    temp_disc        = Pparam(pv["temp_disc"]);
    texp_disc        = Pparam(pv["texp_disc"]);
    lin_limb_disc    = Pparam(pv["lin_limb_disc"]);
    quad_limb_disc   = Pparam(pv["quad_limb_disc"]);
    temp_edge        = Pparam(pv["temp_edge"]);
    absorb_edge      = Pparam(pv["absorb_edge"]);
    radius_spot      = Pparam(pv["radius_spot"]);
    length_spot      = Pparam(pv["length_spot"]);
    height_spot      = Pparam(pv["height_spot"]);
    expon_spot       = Pparam(pv["expon_spot"]);
    epow_spot        = Pparam(pv["epow_spot"]);
    angle_spot       = Pparam(pv["angle_spot"]);
    yaw_spot         = Pparam(pv["yaw_spot"]);
    temp_spot        = Pparam(pv["temp_spot"]);
    tilt_spot        = Pparam(pv["tilt_spot"]);
    cfrac_spot       = Pparam(pv["cfrac_spot"]);

    // star-spot parameters need not have been defined, but
    // if one of a group has, then all of them should be
    if(names["stsp11_long"] || names["stsp11_lat"] ||
       names["stsp11_fwhm"] || names["stsp11_tcen"]){

        if(!(names["stsp11_long"] && names["stsp11_lat"] &&
             names["stsp11_fwhm"] && names["stsp11_tcen"]))
            throw Lcurve_Error("One or more of the star spot 11 parameters were not initialised");
        stsp11_long = Pparam(pv["stsp11_long"]);
        stsp11_lat  = Pparam(pv["stsp11_lat"]);
        stsp11_fwhm = Pparam(pv["stsp11_fwhm"]);
        stsp11_tcen = Pparam(pv["stsp11_tcen"]);
    }

    if(names["stsp12_long"] || names["stsp12_lat"] ||
       names["stsp12_fwhm"] || names["stsp12_tcen"]){

        if(!(names["stsp12_long"] && names["stsp12_lat"] &&
             names["stsp12_fwhm"] && names["stsp12_tcen"]))
            throw Lcurve_Error("One or more of the star spot 12 parameters were not initialised");
        stsp12_long = Pparam(pv["stsp12_long"]);
        stsp12_lat  = Pparam(pv["stsp12_lat"]);
        stsp12_fwhm = Pparam(pv["stsp12_fwhm"]);
        stsp12_tcen = Pparam(pv["stsp12_tcen"]);
    }

    if(names["stsp13_long"] || names["stsp13_lat"] ||
       names["stsp13_fwhm"] || names["stsp13_tcen"]){

        if(!(names["stsp13_long"] && names["stsp13_lat"] &&
             names["stsp13_fwhm"] && names["stsp13_tcen"]))
            throw Lcurve_Error("One or more of the star spot 13 parameters were not initialised");
        stsp13_long = Pparam(pv["stsp13_long"]);
        stsp13_lat  = Pparam(pv["stsp13_lat"]);
        stsp13_fwhm = Pparam(pv["stsp13_fwhm"]);
        stsp13_tcen = Pparam(pv["stsp13_tcen"]);
    }

    if(names["stsp21_long"] || names["stsp21_lat"] ||
       names["stsp21_fwhm"] || names["stsp21_tcen"]){
        if(!(names["stsp21_long"] && names["stsp21_lat"] &&
             names["stsp21_fwhm"] && names["stsp21_tcen"]))
            throw Lcurve_Error("One or more of the star spot 21 parameters were not initialised");
        stsp21_long = Pparam(pv["stsp21_long"]);
        stsp21_lat  = Pparam(pv["stsp21_lat"]);
        stsp21_fwhm = Pparam(pv["stsp21_fwhm"]);
        stsp21_tcen = Pparam(pv["stsp21_tcen"]);
    }

    if(names["stsp22_long"] || names["stsp22_lat"] ||
       names["stsp22_fwhm"] || names["stsp22_tcen"]){

        if(!(names["stsp22_long"] && names["stsp22_lat"] &&
             names["stsp22_fwhm"] && names["stsp22_tcen"]))
            throw Lcurve_Error("One or more of the star spot 22 parameters were not initialised");
        stsp22_long = Pparam(pv["stsp22_long"]);
        stsp22_lat  = Pparam(pv["stsp22_lat"]);
        stsp22_fwhm = Pparam(pv["stsp22_fwhm"]);
        stsp22_tcen = Pparam(pv["stsp22_tcen"]);
    }
	
    if(names["uesp_long1"] || names["uesp_long2"] ||
       names["uesp_lathw"] || names["uesp_taper"] || names["uesp_temp"]){
      if(!(names["uesp_long1"] && names["uesp_long2"] &&
	   names["uesp_lathw"] && names["uesp_taper"] && names["uesp_temp"]))
	throw Lcurve_Error("One or more of uniform equatorial spot parameters were not initialised");
      uesp_long1 = Pparam(pv["uesp_long1"]);
      uesp_long2 = Pparam(pv["uesp_long2"]);
      uesp_lathw = Pparam(pv["uesp_lathw"]);
      uesp_taper = Pparam(pv["uesp_taper"]);
      uesp_temp = Pparam(pv["uesp_temp"]);
    }

    delta_phase = Subs::string_to_double(pv["delta_phase"]);
    nlat1f = Subs::string_to_int(pv["nlat1f"]);
    nlat2f = Subs::string_to_int(pv["nlat2f"]);
    nlat1c = Subs::string_to_int(pv["nlat1c"]);
    nlat2c = Subs::string_to_int(pv["nlat2c"]);
    npole = Subs::string_to_bool(pv["npole"]);

    nlatfill = Subs::string_to_int(pv["nlatfill"]);
    nlngfill = Subs::string_to_int(pv["nlngfill"]);
    lfudge = Subs::string_to_double(pv["lfudge"]);
    llo = Subs::string_to_double(pv["llo"]);
    lhi = Subs::string_to_double(pv["lhi"]);
    phase1 = Subs::string_to_double(pv["phase1"]);
    phase2 = Subs::string_to_double(pv["phase2"]);
    wavelength = Subs::string_to_double(pv["wavelength"]);

    roche1 = Subs::string_to_bool(pv["roche1"]);
    roche2 = Subs::string_to_bool(pv["roche2"]);
    eclipse1 = Subs::string_to_bool(pv["eclipse1"]);
    eclipse2 = Subs::string_to_bool(pv["eclipse2"]);
    glens1 = Subs::string_to_bool(pv["glens1"]);
    if(glens1 && roche1)
      throw Lcurve_Error("For reasons of simplicity, glens1 = 1 and roche1 = 1 are not simultaneously allowed");
    use_radii        = Subs::string_to_bool(pv["use_radii"]);
    tperiod          = Subs::string_to_double(pv["tperiod"]);
    gdark_bolom1     = Subs::string_to_bool(pv["gdark_bolom1"]);
    gdark_bolom2     = Subs::string_to_bool(pv["gdark_bolom2"]);
    mucrit1          = Subs::string_to_double(pv["mucrit1"]);
    mucrit2          = Subs::string_to_double(pv["mucrit2"]);

    if(pv["limb1"] == "Poly"){
        limb1 = LDC::POLY;
    }else if(pv["limb1"] == "Claret"){
        limb1 = LDC::CLARET;
    }else{
        throw Lcurve_Error("Could not recognize the value of limb1 = " +
                           pv["limb1"] + "; should be 'Poly' or 'Claret'");
    }
    if(pv["limb2"] == "Poly"){
        limb2 = LDC::POLY;
    }else if(pv["limb2"] == "Claret"){
        limb2 = LDC::CLARET;
    }else{
        throw Lcurve_Error("Could not recognize the value of limb2 = " +
                           pv["limb2"] + "; should be 'Poly' or 'Claret'");
    }

    mirror           = Subs::string_to_bool(pv["mirror"]);
    add_disc         = Subs::string_to_bool(pv["add_disc"]);
    nrad             = Subs::string_to_int(pv["nrad"]);
    opaque           = Subs::string_to_bool(pv["opaque"]);
    add_spot         = Subs::string_to_bool(pv["add_spot"]);
    nspot            = Subs::string_to_int(pv["nspot"]);
    iscale           = Subs::string_to_bool(pv["iscale"]);

}




/* Constructor from parameters added by lijiao
 */
Lcurve::Model::Model(//Binary and stars
                     double q_value, double q_range, double q_dstep, bool q_vary, bool q_defined,
                     double iangle_value, double iangle_range, double iangle_dstep, bool iangle_vary, bool iangle_defined, 
                     double r1_value, double r1_range, double r1_dstep, bool r1_vary, bool r1_defined, 
                     double r2_value, double r2_range, double r2_dstep, bool r2_vary, bool r2_defined, 
                     double cphi3_value, double cphi3_range, double cphi3_dstep, bool cphi3_vary, bool cphi3_defined, 
                     double cphi4_value, double cphi4_range, double cphi4_dstep, bool cphi4_vary, bool cphi4_defined, 
                     double spin1_value, double spin1_range, double spin1_dstep, bool spin1_vary, bool spin1_defined, 
                     double spin2_value, double spin2_range, double spin2_dstep, bool spin2_vary, bool spin2_defined, 
                     double t1_value, double t1_range, double t1_dstep, bool t1_vary, bool t1_defined, 
                     double t2_value, double t2_range, double t2_dstep, bool t2_vary, bool t2_defined, 
                     double ldc1_1_value, double ldc1_1_range, double ldc1_1_dstep, bool ldc1_1_vary, bool ldc1_1_defined, 
                     double ldc1_2_value, double ldc1_2_range, double ldc1_2_dstep, bool ldc1_2_vary, bool ldc1_2_defined, 
                     double ldc1_3_value, double ldc1_3_range, double ldc1_3_dstep, bool ldc1_3_vary, bool ldc1_3_defined, 
                     double ldc1_4_value, double ldc1_4_range, double ldc1_4_dstep, bool ldc1_4_vary, bool ldc1_4_defined, 
                     double ldc2_1_value, double ldc2_1_range, double ldc2_1_dstep, bool ldc2_1_vary, bool ldc2_1_defined, 
                     double ldc2_2_value, double ldc2_2_range, double ldc2_2_dstep, bool ldc2_2_vary, bool ldc2_2_defined, 
                     double ldc2_3_value, double ldc2_3_range, double ldc2_3_dstep, bool ldc2_3_vary, bool ldc2_3_defined, 
                     double ldc2_4_value, double ldc2_4_range, double ldc2_4_dstep, bool ldc2_4_vary, bool ldc2_4_defined, 
                     double velocity_scale_value, double velocity_scale_range, double velocity_scale_dstep, bool velocity_scale_vary, bool velocity_scale_defined, 
                     double beam_factor1_value, double beam_factor1_range, double beam_factor1_dstep, bool beam_factor1_vary, bool beam_factor1_defined, 
                     double beam_factor2_value, double beam_factor2_range, double beam_factor2_dstep, bool beam_factor2_vary, bool beam_factor2_defined, 
                     //General
                     double t0_value, double t0_range, double t0_dstep, bool t0_vary, bool t0_defined, 
                     double period_value, double period_range, double period_dstep, bool period_vary, bool period_defined, 
                     double pdot_value, double pdot_range, double pdot_dstep, bool pdot_vary, bool pdot_defined, 
                     double deltat_value, double deltat_range, double deltat_dstep, bool deltat_vary, bool deltat_defined, 
                     double gravity_dark1_value, double gravity_dark1_range, double gravity_dark1_dstep, bool gravity_dark1_vary, bool gravity_dark1_defined, 
                     double gravity_dark2_value, double gravity_dark2_range, double gravity_dark2_dstep, bool gravity_dark2_vary, bool gravity_dark2_defined, 
                     double absorb_value, double absorb_range, double absorb_dstep, bool absorb_vary, bool absorb_defined, 
                     double slope_value, double slope_range, double slope_dstep, bool slope_vary, bool slope_defined, 
                     double quad_value, double quad_range, double quad_dstep, bool quad_vary, bool quad_defined, 
                     double cube_value, double cube_range, double cube_dstep, bool cube_vary, bool cube_defined, 
                     double third_value, double third_range, double third_dstep, bool third_vary, bool third_defined, 
                     // Star spots
                     double stsp11_long_value, double stsp11_long_range, double stsp11_long_dstep, bool stsp11_long_vary, bool stsp11_long_defined, 
                     double stsp11_lat_value, double stsp11_lat_range, double stsp11_lat_dstep, bool stsp11_lat_vary, bool stsp11_lat_defined, 
                     double stsp11_fwhm_value, double stsp11_fwhm_range, double stsp11_fwhm_dstep, bool stsp11_fwhm_vary, bool stsp11_fwhm_defined, 
                     double stsp11_tcen_value, double stsp11_tcen_range, double stsp11_tcen_dstep, bool stsp11_tcen_vary, bool stsp11_tcen_defined, 
                     double stsp12_long_value, double stsp12_long_range, double stsp12_long_dstep, bool stsp12_long_vary, bool stsp12_long_defined, 
                     double stsp12_lat_value, double stsp12_lat_range, double stsp12_lat_dstep, bool stsp12_lat_vary, bool stsp12_lat_defined, 
                     double stsp12_fwhm_value, double stsp12_fwhm_range, double stsp12_fwhm_dstep, bool stsp12_fwhm_vary, bool stsp12_fwhm_defined, 
                     double stsp12_tcen_value, double stsp12_tcen_range, double stsp12_tcen_dstep, bool stsp12_tcen_vary, bool stsp12_tcen_defined, 
                     double stsp13_long_value, double stsp13_long_range, double stsp13_long_dstep, bool stsp13_long_vary, bool stsp13_long_defined, 
                     double stsp13_lat_value, double stsp13_lat_range, double stsp13_lat_dstep, bool stsp13_lat_vary, bool stsp13_lat_defined, 
                     double stsp13_fwhm_value, double stsp13_fwhm_range, double stsp13_fwhm_dstep, bool stsp13_fwhm_vary, bool stsp13_fwhm_defined, 
                     double stsp13_tcen_value, double stsp13_tcen_range, double stsp13_tcen_dstep, bool stsp13_tcen_vary, bool stsp13_tcen_defined, 
                     double stsp21_long_value, double stsp21_long_range, double stsp21_long_dstep, bool stsp21_long_vary, bool stsp21_long_defined, 
                     double stsp21_lat_value, double stsp21_lat_range, double stsp21_lat_dstep, bool stsp21_lat_vary, bool stsp21_lat_defined, 
                     double stsp21_fwhm_value, double stsp21_fwhm_range, double stsp21_fwhm_dstep, bool stsp21_fwhm_vary, bool stsp21_fwhm_defined, 
                     double stsp21_tcen_value, double stsp21_tcen_range, double stsp21_tcen_dstep, bool stsp21_tcen_vary, bool stsp21_tcen_defined, 
                     double stsp22_long_value, double stsp22_long_range, double stsp22_long_dstep, bool stsp22_long_vary, bool stsp22_long_defined, 
                     double stsp22_lat_value, double stsp22_lat_range, double stsp22_lat_dstep, bool stsp22_lat_vary, bool stsp22_lat_defined, 
                     double stsp22_fwhm_value, double stsp22_fwhm_range, double stsp22_fwhm_dstep, bool stsp22_fwhm_vary, bool stsp22_fwhm_defined, 
                     double stsp22_tcen_value, double stsp22_tcen_range, double stsp22_tcen_dstep, bool stsp22_tcen_vary, bool stsp22_tcen_defined, 
                     double uesp_long1_value, double uesp_long1_range, double uesp_long1_dstep, bool uesp_long1_vary, bool uesp_long1_defined,
                     double uesp_long2_value, double uesp_long2_range, double uesp_long2_dstep, bool uesp_long2_vary, bool uesp_long2_defined,
                     double uesp_lathw_value, double uesp_lathw_range, double uesp_lathw_dstep, bool uesp_lathw_vary, bool uesp_lathw_defined,
                     double uesp_taper_value, double uesp_taper_range, double uesp_taper_dstep, bool uesp_taper_vary, bool uesp_taper_defined,
                     double uesp_temp_value, double uesp_temp_range, double uesp_temp_dstep, bool uesp_temp_vary, bool uesp_temp_defined,
                     // disc
                     double rdisc1_value, double rdisc1_range, double rdisc1_dstep, bool rdisc1_vary, bool rdisc1_defined, //disc
                     double rdisc2_value, double rdisc2_range, double rdisc2_dstep, bool rdisc2_vary, bool rdisc2_defined, 
                     double height_disc_value, double height_disc_range, double height_disc_dstep, bool height_disc_vary, bool height_disc_defined, 
                     double beta_disc_value, double beta_disc_range, double beta_disc_dstep, bool beta_disc_vary, bool beta_disc_defined, 
                     double temp_disc_value, double temp_disc_range, double temp_disc_dstep, bool temp_disc_vary, bool temp_disc_defined, 
                     double texp_disc_value, double texp_disc_range, double texp_disc_dstep, bool texp_disc_vary, bool texp_disc_defined, 
                     double lin_limb_disc_value, double lin_limb_disc_range, double lin_limb_disc_dstep, bool lin_limb_disc_vary, bool lin_limb_disc_defined, 
                     double quad_limb_disc_value, double quad_limb_disc_range, double quad_limb_disc_dstep, bool quad_limb_disc_vary, bool quad_limb_disc_defined, 
                     double temp_edge_value, double temp_edge_range, double temp_edge_dstep, bool temp_edge_vary, bool temp_edge_defined,
                     double absorb_edge_value, double absorb_edge_range, double absorb_edge_dstep, bool absorb_edge_vary, bool absorb_edge_defined, 
                     //Bright-spot
                     double radius_spot_value, double radius_spot_range, double radius_spot_dstep, bool radius_spot_vary, bool radius_spot_defined, //Bright-spot
                     double length_spot_value, double length_spot_range, double length_spot_dstep, bool length_spot_vary, bool length_spot_defined, 
                     double height_spot_value, double height_spot_range, double height_spot_dstep, bool height_spot_vary, bool height_spot_defined, 
                     double expon_spot_value, double expon_spot_range, double expon_spot_dstep, bool expon_spot_vary, bool expon_spot_defined, 
                     double epow_spot_value, double epow_spot_range, double epow_spot_dstep, bool epow_spot_vary, bool epow_spot_defined, 
                     double angle_spot_value, double angle_spot_range, double angle_spot_dstep, bool angle_spot_vary, bool angle_spot_defined, 
                     double yaw_spot_value, double yaw_spot_range, double yaw_spot_dstep, bool yaw_spot_vary, bool yaw_spot_defined, 
                     double temp_spot_value, double temp_spot_range, double temp_spot_dstep, bool temp_spot_vary, bool temp_spot_defined, 
                     double tilt_spot_value, double tilt_spot_range, double tilt_spot_dstep, bool tilt_spot_vary, bool tilt_spot_defined, 
                     double cfrac_spot_value, double cfrac_spot_range, double cfrac_spot_dstep, bool cfrac_spot_vary, bool cfrac_spot_defined,
                     // Computational parameters
                     double delta_phase, int nlat1f, int nlat2f, int nlat1c, int nlat2c, bool npole, 
                     int nlatfill, int nlngfill, double lfudge, double llo, double lhi, double phase1, double phase2, int nrad, double wavelength,
                     bool roche1, bool roche2, bool eclipse1, bool eclipse2, bool glens1, bool use_radii,
                     double tperiod, bool gdark_bolom1, bool gdark_bolom2, double mucrit1, double mucrit2, 
                     const char* pslimb1, const char* pslimb2, bool mirror, bool add_disc, bool opaque, bool add_spot, int nspot, bool iscale
                    ) {


    // Star spots:
    //
    // 11 = star 1, spot 1
    // 21 = star 2, spot 1
    //

    // Computational

    // Parameter value pairs
    // Initialise physical parameters
    q                = Pparam(q_value, q_range, q_dstep, q_vary, q_defined);
    // std::cout << "lijiao trm_lroche.cc:  q2): " << q.value << std::endl; 
    iangle           = Pparam(iangle_value, iangle_range, iangle_dstep, iangle_vary, iangle_defined);
    r1               = Pparam(r1_value, r1_range, r1_dstep, r1_vary, r1_defined);
    r2               = Pparam(r2_value, r2_range, r2_dstep, r2_vary, r2_defined);
    cphi3            = Pparam(cphi3_value, cphi3_range, cphi3_dstep, cphi3_vary, cphi3_defined);
    cphi4            = Pparam(cphi4_value, cphi4_range, cphi4_dstep, cphi4_vary, cphi4_defined);
    spin1            = Pparam(spin1_value, spin1_range, spin1_dstep, spin1_vary, spin1_defined);
    spin2            = Pparam(spin2_value, spin2_range, spin2_dstep, spin2_vary, spin2_defined);
    t1               = Pparam(t1_value, t1_range, t1_dstep, t1_vary, t1_defined);
    t2               = Pparam(t2_value, t2_range, t2_dstep, t2_vary, t2_defined);
    ldc1_1           = Pparam(ldc1_1_value, ldc1_1_range, ldc1_1_dstep, ldc1_1_vary, ldc1_1_defined);
    ldc1_2           = Pparam(ldc1_2_value, ldc1_2_range, ldc1_2_dstep, ldc1_2_vary, ldc1_2_defined);
    ldc1_3           = Pparam(ldc1_3_value, ldc1_3_range, ldc1_3_dstep, ldc1_3_vary, ldc1_3_defined);
    ldc1_4           = Pparam(ldc1_4_value, ldc1_4_range, ldc1_4_dstep, ldc1_4_vary, ldc1_4_defined);
    ldc2_1           = Pparam(ldc2_1_value, ldc2_1_range, ldc2_1_dstep, ldc2_1_vary, ldc2_1_defined);
    ldc2_2           = Pparam(ldc2_2_value, ldc2_2_range, ldc2_2_dstep, ldc2_2_vary, ldc2_2_defined);
    ldc2_3           = Pparam(ldc2_3_value, ldc2_3_range, ldc2_3_dstep, ldc2_3_vary, ldc2_3_defined);
    ldc2_4           = Pparam(ldc2_4_value, ldc2_4_range, ldc2_4_dstep, ldc2_4_vary, ldc2_4_defined);
    velocity_scale   = Pparam(velocity_scale_value, velocity_scale_range, velocity_scale_dstep, velocity_scale_vary, velocity_scale_defined);
    beam_factor1     = Pparam(beam_factor1_value, beam_factor1_range, beam_factor1_dstep, beam_factor1_vary, beam_factor1_defined);
    beam_factor2     = Pparam(beam_factor2_value, beam_factor2_range, beam_factor2_dstep, beam_factor2_vary, beam_factor2_defined);
    t0               = Pparam(t0_value, t0_range, t0_dstep, t0_vary, t0_defined);
    period           = Pparam(period_value, period_range, period_dstep, period_vary, period_defined);
    pdot             = Pparam(pdot_value, pdot_range, pdot_dstep, pdot_vary, pdot_defined);
    deltat           = Pparam(deltat_value, deltat_range, deltat_dstep, deltat_vary, deltat_defined);
    gravity_dark1    = Pparam(gravity_dark1_value, gravity_dark1_range, gravity_dark1_dstep, gravity_dark1_vary, gravity_dark1_defined);
    gravity_dark2    = Pparam(gravity_dark2_value, gravity_dark2_range, gravity_dark2_dstep, gravity_dark2_vary, gravity_dark2_defined);
    absorb           = Pparam(absorb_value, absorb_range, absorb_dstep, absorb_vary, absorb_defined);
    slope            = Pparam(slope_value, slope_range, slope_dstep, slope_vary, slope_defined);
    quad             = Pparam(quad_value, quad_range, quad_dstep, quad_vary, quad_defined);
    cube             = Pparam(cube_value, cube_range, cube_dstep, cube_vary, cube_defined);
    third            = Pparam(third_value, third_range, third_dstep, third_vary, third_defined);
    rdisc1           = Pparam(rdisc1_value, rdisc1_range, rdisc1_dstep, rdisc1_vary, rdisc1_defined);
    rdisc2           = Pparam(rdisc2_value, rdisc2_range, rdisc2_dstep, rdisc2_vary, rdisc2_defined);
    height_disc      = Pparam(height_disc_value, height_disc_range, height_disc_dstep, height_disc_vary, height_disc_defined);
    beta_disc        = Pparam(beta_disc_value, beta_disc_range, beta_disc_dstep, beta_disc_vary, beta_disc_defined);
    temp_disc        = Pparam(temp_disc_value, temp_disc_range, temp_disc_dstep, temp_disc_vary, temp_disc_defined);
    texp_disc        = Pparam(texp_disc_value, texp_disc_range, texp_disc_dstep, texp_disc_vary, texp_disc_defined);
    lin_limb_disc    = Pparam(lin_limb_disc_value, lin_limb_disc_range, lin_limb_disc_dstep, lin_limb_disc_vary, lin_limb_disc_defined);
    quad_limb_disc   = Pparam(quad_limb_disc_value, quad_limb_disc_range, quad_limb_disc_dstep, quad_limb_disc_vary, quad_limb_disc_defined);
    temp_edge        = Pparam(temp_edge_value, temp_edge_range, temp_edge_dstep, temp_edge_vary, temp_edge_defined);
    absorb_edge      = Pparam(absorb_edge_value, absorb_edge_range, absorb_edge_dstep, absorb_edge_vary, absorb_edge_defined);
    radius_spot      = Pparam(radius_spot_value, radius_spot_range, radius_spot_dstep, radius_spot_vary, radius_spot_defined);
    length_spot      = Pparam(length_spot_value, length_spot_range, length_spot_dstep, length_spot_vary, length_spot_defined);
    height_spot      = Pparam(height_spot_value, height_spot_range, height_spot_dstep, height_spot_vary, height_spot_defined);
    expon_spot       = Pparam(expon_spot_value, expon_spot_range, expon_spot_dstep, expon_spot_vary, expon_spot_defined);
    epow_spot        = Pparam(epow_spot_value, epow_spot_range, epow_spot_dstep, epow_spot_vary, epow_spot_defined);
    angle_spot       = Pparam(angle_spot_value, angle_spot_range, angle_spot_dstep, angle_spot_vary, angle_spot_defined);
    yaw_spot         = Pparam(yaw_spot_value, yaw_spot_range, yaw_spot_dstep, yaw_spot_vary, yaw_spot_defined);
    temp_spot        = Pparam(temp_spot_value, temp_spot_range, temp_spot_dstep, temp_spot_vary, temp_spot_defined);
    tilt_spot        = Pparam(tilt_spot_value, tilt_spot_range, tilt_spot_dstep, tilt_spot_vary, tilt_spot_defined);
    cfrac_spot       = Pparam(cfrac_spot_value, cfrac_spot_range, cfrac_spot_dstep, cfrac_spot_vary, cfrac_spot_defined);

    // star-spot parameters need not have been defined, but
    // if one of a group has, then all of them should be
    if(stsp11_long_defined || stsp11_lat_defined ||
        stsp11_fwhm_defined || stsp11_tcen_defined){
        if(!(stsp11_long_defined && stsp11_lat_defined &&
             stsp11_fwhm_defined && stsp11_tcen_defined))
            throw Lcurve_Error("One or more of the star spot 11 parameters were not initialised");
        stsp11_long      = Pparam(stsp11_long_value, stsp11_long_range, stsp11_long_dstep, stsp11_long_vary, stsp11_long_defined);
        stsp11_lat       = Pparam(stsp11_lat_value, stsp11_lat_range, stsp11_lat_dstep, stsp11_lat_vary, stsp11_lat_defined);
        stsp11_fwhm      = Pparam(stsp11_fwhm_value, stsp11_fwhm_range, stsp11_fwhm_dstep, stsp11_fwhm_vary, stsp11_fwhm_defined);
        stsp11_tcen      = Pparam(stsp11_tcen_value, stsp11_tcen_range, stsp11_tcen_dstep, stsp11_tcen_vary, stsp11_tcen_defined);
       }

    if(stsp12_long_defined || stsp12_lat_defined ||
        stsp12_fwhm_defined || stsp12_tcen_defined){
        if(!(stsp12_long_defined && stsp12_lat_defined &&
             stsp12_fwhm_defined && stsp12_tcen_defined))
            throw Lcurve_Error("One or more of the star spot 12 parameters were not initialised");
        stsp12_long      = Pparam(stsp12_long_value, stsp12_long_range, stsp12_long_dstep, stsp12_long_vary, stsp12_long_defined);
        stsp12_lat       = Pparam(stsp12_lat_value, stsp12_lat_range, stsp12_lat_dstep, stsp12_lat_vary, stsp12_lat_defined);
        stsp12_fwhm      = Pparam(stsp12_fwhm_value, stsp12_fwhm_range, stsp12_fwhm_dstep, stsp12_fwhm_vary, stsp12_fwhm_defined);
        stsp12_tcen      = Pparam(stsp12_tcen_value, stsp12_tcen_range, stsp12_tcen_dstep, stsp12_tcen_vary, stsp12_tcen_defined);
       }

    if(stsp13_long_defined || stsp13_lat_defined ||
        stsp13_fwhm_defined || stsp13_tcen_defined){
        if(!(stsp13_long_defined && stsp13_lat_defined &&
             stsp13_fwhm_defined && stsp13_tcen_defined))
            throw Lcurve_Error("One or more of the star spot 13 parameters were not initialised");
        stsp13_long      = Pparam(stsp13_long_value, stsp13_long_range, stsp13_long_dstep, stsp13_long_vary, stsp13_long_defined);
        stsp13_lat       = Pparam(stsp13_lat_value, stsp13_lat_range, stsp13_lat_dstep, stsp13_lat_vary, stsp13_lat_defined);
        stsp13_fwhm      = Pparam(stsp13_fwhm_value, stsp13_fwhm_range, stsp13_fwhm_dstep, stsp13_fwhm_vary, stsp13_fwhm_defined);
        stsp13_tcen      = Pparam(stsp13_tcen_value, stsp13_tcen_range, stsp13_tcen_dstep, stsp13_tcen_vary, stsp13_tcen_defined);
       }

    if(stsp21_long_defined || stsp21_lat_defined ||
        stsp21_fwhm_defined || stsp21_tcen_defined){
        if(!(stsp21_long_defined && stsp21_lat_defined &&
             stsp21_fwhm_defined && stsp21_tcen_defined))
            throw Lcurve_Error("One or more of the star spot 21 parameters were not initialised");
        stsp21_long      = Pparam(stsp21_long_value, stsp21_long_range, stsp21_long_dstep, stsp21_long_vary, stsp21_long_defined);
        stsp21_lat       = Pparam(stsp21_lat_value, stsp21_lat_range, stsp21_lat_dstep, stsp21_lat_vary, stsp21_lat_defined);
        stsp21_fwhm      = Pparam(stsp21_fwhm_value, stsp21_fwhm_range, stsp21_fwhm_dstep, stsp21_fwhm_vary, stsp21_fwhm_defined);
        stsp21_tcen      = Pparam(stsp21_tcen_value, stsp21_tcen_range, stsp21_tcen_dstep, stsp21_tcen_vary, stsp21_tcen_defined);
         }

    if(stsp22_long_defined || stsp22_lat_defined ||
        stsp22_fwhm_defined || stsp22_tcen_defined){
        if(!(stsp22_long_defined && stsp22_lat_defined &&
             stsp22_fwhm_defined && stsp22_tcen_defined))
            throw Lcurve_Error("One or more of the star spot 22 parameters were not initialised");
        stsp22_long      = Pparam(stsp22_long_value, stsp22_long_range, stsp22_long_dstep, stsp22_long_vary, stsp22_long_defined);
        stsp22_lat       = Pparam(stsp22_lat_value, stsp22_lat_range, stsp22_lat_dstep, stsp22_lat_vary, stsp22_lat_defined);
        stsp22_fwhm      = Pparam(stsp22_fwhm_value, stsp22_fwhm_range, stsp22_fwhm_dstep, stsp22_fwhm_vary, stsp22_fwhm_defined);
        stsp22_tcen      = Pparam(stsp22_tcen_value, stsp22_tcen_range, stsp22_tcen_dstep, stsp22_tcen_vary, stsp22_tcen_defined);
       }

    if(uesp_long1_defined || uesp_long2_defined ||
        uesp_lathw_defined || uesp_taper_defined || uesp_temp_defined){
        if(!(uesp_long1_defined && uesp_long2_defined &&
             uesp_lathw_defined && uesp_taper_defined && uesp_temp_defined))
         throw Lcurve_Error("One or more of uniform equatorial spot parameters were not initialised");
        uesp_long1 = Pparam(uesp_long1_value, uesp_long1_range, uesp_long1_dstep, uesp_long1_vary, uesp_long1_defined);
        uesp_long2 = Pparam(uesp_long2_value, uesp_long2_range, uesp_long2_dstep, uesp_long2_vary, uesp_long2_defined);
        uesp_lathw = Pparam(uesp_lathw_value, uesp_lathw_range, uesp_lathw_dstep, uesp_lathw_vary, uesp_lathw_defined);
        uesp_taper = Pparam(uesp_taper_value, uesp_taper_range, uesp_taper_dstep, uesp_taper_vary, uesp_taper_defined);
        uesp_temp = Pparam(uesp_temp_value, uesp_temp_range, uesp_temp_dstep, uesp_temp_vary, uesp_temp_defined);
      }

    this->delta_phase = delta_phase;
    this->nlat1f = nlat1f, this->nlat2f = nlat2f, this->nlat1c = nlat1c, this->nlat2c = nlat2c;
    this->npole = npole; this->nlatfill = nlatfill, this->nlngfill = nlngfill;
    this->lfudge = lfudge, this->llo = llo, this->lhi = lhi;
    this->phase1 = phase1, this->phase2 = phase2, this->nrad = nrad;
    this->wavelength = wavelength;
    this->roche1 = roche1, this->roche2 = roche2, this->eclipse1 = eclipse1, this->eclipse2 = eclipse2;
    this->glens1 = glens1, this->use_radii = use_radii;
    this->tperiod = tperiod, this->gdark_bolom1 = gdark_bolom1, this->gdark_bolom2 = gdark_bolom2;
    this->mucrit1 = mucrit1, this->mucrit2 = mucrit2;
    this->mirror = mirror, this->add_disc = add_disc, this->opaque = opaque, this->add_spot = add_spot, this->nspot = nspot;
    this->iscale = iscale;
    if(glens1 && roche1)
      throw Lcurve_Error("For reasons of simplicity, glens1 = 1 and roche1 = 1 are not simultaneously allowed");
    std::string slimb1 = pslimb1;
    std::string slimb2 = pslimb2;
    //if(strcmp(pslimb1, "Poly") == 0){
    if(slimb1 == "Poly"){
        limb1 = LDC::POLY;
    }else if(slimb1 ==  "Claret"){
        limb1 = LDC::CLARET;
    }else{
        throw Lcurve_Error("Could not recognize the value of slimb1 = " +
                           slimb1 + "; should be 'Poly' or 'Claret'");
    }
    if(slimb2 == "Poly"){
        limb2 = LDC::POLY;
    }else if(slimb2 == "Claret"){
        limb2 = LDC::CLARET;
    }else{
        throw Lcurve_Error("Could not recognize the value of slimb2 = " +
                           slimb2 + "; should be 'Poly' or 'Claret'");
    }

}



/** Works out the number of variable physical parameters */
int Lcurve::Model::nvary() const {
    int n = 0;

    if(q.vary) n++;
    if(iangle.vary) n++;
    if(use_radii){
        if(r1.vary) n++;
        if(r2.vary) n++;
    }else{
        if(cphi3.vary) n++;
        if(cphi4.vary) n++;
    }
    if(spin1.vary) n++;
    if(spin2.vary) n++;
    if(t1.vary) n++;
    if(t2.vary) n++;
    if(ldc1_1.vary) n++;
    if(ldc1_2.vary) n++;
    if(ldc1_3.vary) n++;
    if(ldc1_4.vary) n++;
    if(ldc2_1.vary) n++;
    if(ldc2_2.vary) n++;
    if(ldc2_3.vary) n++;
    if(ldc2_4.vary) n++;
    if(velocity_scale.vary) n++;
    if(beam_factor1.vary) n++;
    if(beam_factor2.vary) n++;

    if(t0.vary) n++;
    if(period.vary) n++;
    if(pdot.vary) n++;
    if(deltat.vary) n++;
    if(gravity_dark1.vary) n++;
    if(gravity_dark2.vary) n++;
    if(absorb.vary) n++;
    if(slope.vary) n++;
    if(quad.vary) n++;
    if(cube.vary) n++;
    if(third.vary) n++;

    if(add_disc){
        if(rdisc1.vary) n++;
        if(rdisc2.vary) n++;
        if(height_disc.vary) n++;
        if(beta_disc.vary) n++;
        if(temp_disc.vary) n++;
        if(texp_disc.vary) n++;
        if(lin_limb_disc.vary) n++;
        if(quad_limb_disc.vary) n++;
        if(temp_edge.vary) n++;
        if(absorb_edge.vary) n++;
    }

    if(add_spot){
        if(radius_spot.vary) n++;
        if(length_spot.vary) n++;
        if(height_spot.vary) n++;
        if(expon_spot.vary) n++;
        if(epow_spot.vary) n++;
        if(angle_spot.vary) n++;
        if(yaw_spot.vary) n++;
        if(temp_spot.vary) n++;
        if(tilt_spot.vary) n++;
        if(cfrac_spot.vary) n++;
    }

    if(stsp11_long.defined && stsp11_long.vary) n++;
    if(stsp11_lat.defined  && stsp11_lat.vary) n++;
    if(stsp11_fwhm.defined && stsp11_fwhm.vary) n++;
    if(stsp11_tcen.defined && stsp11_tcen.vary) n++;

    if(stsp12_long.defined && stsp12_long.vary) n++;
    if(stsp12_lat.defined  && stsp12_lat.vary) n++;
    if(stsp12_fwhm.defined && stsp12_fwhm.vary) n++;
    if(stsp12_tcen.defined && stsp12_tcen.vary) n++;

    if(stsp13_long.defined && stsp13_long.vary) n++;
    if(stsp13_lat.defined  && stsp13_lat.vary) n++;
    if(stsp13_fwhm.defined && stsp13_fwhm.vary) n++;
    if(stsp13_tcen.defined && stsp13_tcen.vary) n++;

    if(stsp21_long.defined && stsp21_long.vary) n++;
    if(stsp21_lat.defined  && stsp21_lat.vary) n++;
    if(stsp21_fwhm.defined && stsp21_fwhm.vary) n++;
    if(stsp21_tcen.defined && stsp21_tcen.vary) n++;

    if(stsp22_long.defined && stsp22_long.vary) n++;
    if(stsp22_lat.defined  && stsp22_lat.vary) n++;
    if(stsp22_fwhm.defined && stsp22_fwhm.vary) n++;
    if(stsp22_tcen.defined && stsp22_tcen.vary) n++;
    
    if(uesp_long1.defined && uesp_long1.vary) n++;
    if(uesp_long2.defined && uesp_long2.vary) n++;
    if(uesp_lathw.defined && uesp_lathw.vary) n++;
    if(uesp_taper.defined && uesp_taper.vary) n++;
    if(uesp_temp.defined && uesp_temp.vary) n++;
    
    return n;
}

void Lcurve::Model::set_param(const Subs::Array1D<double>& vpar) {
    if(nvary() != vpar.size())
        throw Lcurve_Error("Lcurve::Model::set_param: conflicting numbers of variable parameters");
    int n = 0;

    if(q.vary)         q.value                   = vpar[n++];
    if(iangle.vary)    iangle.value              = vpar[n++];
    if(use_radii){
        if(r1.vary)        r1.value                  = vpar[n++];
        if(r2.vary)        r2.value                  = vpar[n++];
    }else{
        if(cphi3.vary)     cphi3.value               = vpar[n++];
        if(cphi4.vary)     cphi4.value               = vpar[n++];
    }
    if(spin1.vary)     spin1.value               = vpar[n++];
    if(spin2.vary)     spin2.value               = vpar[n++];
    if(t1.vary)        t1.value                  = vpar[n++];
    if(t2.vary)        t2.value                  = vpar[n++];
    if(ldc1_1.vary)    ldc1_1.value              = vpar[n++];
    if(ldc1_2.vary)    ldc1_2.value              = vpar[n++];
    if(ldc1_3.vary)    ldc1_3.value              = vpar[n++];
    if(ldc1_4.vary)    ldc1_4.value              = vpar[n++];
    if(ldc2_1.vary)    ldc2_1.value              = vpar[n++];
    if(ldc2_2.vary)    ldc2_2.value              = vpar[n++];
    if(ldc2_3.vary)    ldc2_3.value              = vpar[n++];
    if(ldc2_4.vary)    ldc2_4.value              = vpar[n++];
    if(velocity_scale.vary) velocity_scale.value = vpar[n++];
    if(beam_factor1.vary) beam_factor1.value     = vpar[n++];
    if(beam_factor2.vary) beam_factor2.value     = vpar[n++];

    if(t0.vary)        t0.value                  = vpar[n++];
    if(period.vary)    period.value              = vpar[n++];
    if(pdot.vary)      pdot.value                = vpar[n++];
    if(deltat.vary)    deltat.value              = vpar[n++];
    if(gravity_dark1.vary) gravity_dark1.value   = vpar[n++];
    if(gravity_dark2.vary) gravity_dark2.value   = vpar[n++];
    if(absorb.vary)    absorb.value              = vpar[n++];
    if(slope.vary) slope.value                   = vpar[n++];
    if(quad.vary) quad.value                     = vpar[n++];
    if(cube.vary) cube.value                     = vpar[n++];
    if(third.vary) third.value                   = vpar[n++];

    if(add_disc){
        if(rdisc1.vary) rdisc1.value                 = vpar[n++];
        if(rdisc2.vary) rdisc2.value                 = vpar[n++];
        if(height_disc.vary) height_disc.value       = vpar[n++];
        if(beta_disc.vary) beta_disc.value           = vpar[n++];
        if(temp_disc.vary) temp_disc.value           = vpar[n++];
        if(texp_disc.vary) texp_disc.value           = vpar[n++];
        if(lin_limb_disc.vary) lin_limb_disc.value   = vpar[n++];
        if(quad_limb_disc.vary) quad_limb_disc.value = vpar[n++];
        if(temp_edge.vary) temp_edge.value           = vpar[n++];
        if(absorb_edge.vary) absorb_edge.value       = vpar[n++];
    }

    if(add_spot){
        if(radius_spot.vary) radius_spot.value       = vpar[n++];
        if(length_spot.vary) length_spot.value       = vpar[n++];
        if(height_spot.vary) height_spot.value       = vpar[n++];
        if(expon_spot.vary) expon_spot.value         = vpar[n++];
        if(epow_spot.vary) epow_spot.value           = vpar[n++];
        if(angle_spot.vary) angle_spot.value         = vpar[n++];
        if(yaw_spot.vary) yaw_spot.value             = vpar[n++];
        if(temp_spot.vary) temp_spot.value           = vpar[n++];
        if(tilt_spot.vary) tilt_spot.value           = vpar[n++];
        if(cfrac_spot.vary) cfrac_spot.value         = vpar[n++];
    }

    if(stsp11_long.defined && stsp11_long.vary) stsp11_long.value = vpar[n++];
    if(stsp11_lat.defined  && stsp11_lat.vary) stsp11_lat.value = vpar[n++];
    if(stsp11_fwhm.defined && stsp11_fwhm.vary) stsp11_fwhm.value = vpar[n++];
    if(stsp11_tcen.defined && stsp11_tcen.vary) stsp11_tcen.value = vpar[n++];

    if(stsp12_long.defined && stsp12_long.vary) stsp12_long.value = vpar[n++];
    if(stsp12_lat.defined  && stsp12_lat.vary) stsp12_lat.value = vpar[n++];
    if(stsp12_fwhm.defined && stsp12_fwhm.vary) stsp12_fwhm.value = vpar[n++];
    if(stsp12_tcen.defined && stsp12_tcen.vary) stsp12_tcen.value = vpar[n++];

    if(stsp13_long.defined && stsp13_long.vary) stsp13_long.value = vpar[n++];
    if(stsp13_lat.defined  && stsp13_lat.vary) stsp13_lat.value = vpar[n++];
    if(stsp13_fwhm.defined && stsp13_fwhm.vary) stsp13_fwhm.value = vpar[n++];
    if(stsp13_tcen.defined && stsp13_tcen.vary) stsp13_tcen.value = vpar[n++];

    if(stsp21_long.defined && stsp21_long.vary) stsp21_long.value = vpar[n++];
    if(stsp21_lat.defined  && stsp21_lat.vary) stsp21_lat.value = vpar[n++];
    if(stsp21_fwhm.defined && stsp21_fwhm.vary) stsp21_fwhm.value = vpar[n++];
    if(stsp21_tcen.defined && stsp21_tcen.vary) stsp21_tcen.value = vpar[n++];

    if(stsp22_long.defined && stsp22_long.vary) stsp22_long.value = vpar[n++];
    if(stsp22_lat.defined  && stsp22_lat.vary) stsp22_lat.value = vpar[n++];
    if(stsp22_fwhm.defined && stsp22_fwhm.vary) stsp22_fwhm.value = vpar[n++];
    if(stsp22_tcen.defined && stsp22_tcen.vary) stsp22_tcen.value = vpar[n++];

    if(uesp_long1.defined && uesp_long1.vary) uesp_long1.value = vpar[n++];
    if(uesp_long2.defined && uesp_long2.vary) uesp_long2.value = vpar[n++];
    if(uesp_lathw.defined && uesp_lathw.vary) uesp_lathw.value = vpar[n++];
    if(uesp_taper.defined && uesp_taper.vary) uesp_taper.value = vpar[n++];
    if(uesp_temp.defined && uesp_temp.vary) uesp_temp.value = vpar[n++];
}

std::string Lcurve::Model::get_name(int i) const {
    if(i >= nvary())
        throw Lcurve_Error("Lcurve::Model::get_name: parameter index = " +
                           Subs::str(i) + " is out of range 0 to " + Subs::str(nvary()));

    int n = -1;
    if(q.vary) n++;
    if(n == i) return "q";

    if(iangle.vary) n++;
    if(n == i) return "iangle";

    if(use_radii){
        if(r1.vary) n++;
        if(n == i) return "r1";

        if(r2.vary) n++;
        if(n == i) return "r2";
    }else{
        if(cphi3.vary) n++;
        if(n == i) return "cphi3";

        if(cphi4.vary) n++;
        if(n == i) return "cphi4";
    }

    if(spin1.vary) n++;
    if(n == i) return "spin1";

    if(spin2.vary) n++;
    if(n == i) return "spin2";

    if(t1.vary) n++;
    if(n == i) return "t1";

    if(t2.vary) n++;
    if(n == i) return "t2";

    if(ldc1_1.vary) n++;
    if(n == i) return "ldc1_1";

    if(ldc1_2.vary) n++;
    if(n == i) return "ldc1_2";

    if(ldc1_3.vary) n++;
    if(n == i) return "ldc1_3";

    if(ldc1_4.vary) n++;
    if(n == i) return "ldc1_4";

    if(ldc2_1.vary) n++;
    if(n == i) return "ldc2_1";

    if(ldc2_2.vary) n++;
    if(n == i) return "ldc2_2";

    if(ldc2_3.vary) n++;
    if(n == i) return "ldc2_3";

    if(ldc2_4.vary) n++;
    if(n == i) return "ldc2_4";

    if(velocity_scale.vary)   n++;
    if(n == i) return "velocity_scale";

    if(beam_factor1.vary)   n++;
    if(n == i) return "beam_factor1";

    if(beam_factor2.vary)   n++;
    if(n == i) return "beam_factor2";

    if(t0.vary)        n++;
    if(n == i) return "t0";

    if(period.vary)    n++;
    if(n == i) return "period";

    if(pdot.vary)    n++;
    if(n == i) return "pdot";

    if(deltat.vary)    n++;
    if(n == i) return "deltat";

    if(gravity_dark1.vary) n++;
    if(n == i) return "gravity_dark1";

    if(gravity_dark2.vary) n++;
    if(n == i) return "gravity_dark2";

    if(absorb.vary) n++;
    if(n == i) return "absorb";

    if(slope.vary) n++;
    if(n == i) return "slope";

    if(quad.vary) n++;
    if(n == i) return "quad";

    if(cube.vary) n++;
    if(n == i) return "cube";

    if(third.vary) n++;
    if(n == i) return "third";

    if(add_disc){
        if(rdisc1.vary) n++;
        if(n == i) return "rdisc1";

        if(rdisc2.vary) n++;
        if(n == i) return "rdisc2";

        if(height_disc.vary) n++;
        if(n == i) return "height_disc";

        if(beta_disc.vary) n++;
        if(n == i) return "beta_disc";

        if(temp_disc.vary) n++;
        if(n == i) return "temp_disc";

        if(texp_disc.vary) n++;
        if(n == i) return "texp_disc";

        if(lin_limb_disc.vary) n++;
        if(n == i) return "lin_limb_disc";

        if(quad_limb_disc.vary) n++;
        if(n == i) return "quad_limb_disc";

        if(temp_edge.vary) n++;
        if(n == i) return "temp_edge";

        if(absorb_edge.vary) n++;
        if(n == i) return "absorb_edge";
    }

    if(add_spot){
        if(radius_spot.vary) n++;
        if(n == i) return "radius_spot";

        if(length_spot.vary) n++;
        if(n == i) return "length_spot";

        if(height_spot.vary) n++;
        if(n == i) return "height_spot";

        if(expon_spot.vary)  n++;
        if(n == i) return "expon_spot";

        if(epow_spot.vary)  n++;
        if(n == i) return "epow_spot";

        if(angle_spot.vary)  n++;
        if(n == i) return "angle_spot";

        if(yaw_spot.vary)  n++;
        if(n == i) return "yaw_spot";

        if(temp_spot.vary)   n++;
        if(n == i) return "temp_spot";

        if(tilt_spot.vary)   n++;
        if(n == i) return "tilt_spot";

        if(cfrac_spot.vary)   n++;
        if(n == i) return "cfrac_spot";

    }

    if(stsp11_long.defined && stsp11_long.vary) n++;
    if(n == i) return "stsp11_long";
    if(stsp11_lat.defined  && stsp11_lat.vary) n++;
    if(n == i) return "stsp11_lat";
    if(stsp11_fwhm.defined && stsp11_fwhm.vary) n++;
    if(n == i) return "stsp11_fwhm";
    if(stsp11_tcen.defined && stsp11_tcen.vary) n++;
    if(n == i) return "stsp11_tcen";

    if(stsp12_long.defined && stsp12_long.vary) n++;
    if(n == i) return "stsp12_long";
    if(stsp12_lat.defined  && stsp12_lat.vary) n++;
    if(n == i) return "stsp12_lat";
    if(stsp12_fwhm.defined && stsp12_fwhm.vary) n++;
    if(n == i) return "stsp12_fwhm";
    if(stsp12_tcen.defined && stsp12_tcen.vary) n++;
    if(n == i) return "stsp12_tcen";

    if(stsp13_long.defined && stsp13_long.vary) n++;
    if(n == i) return "stsp13_long";
    if(stsp13_lat.defined  && stsp13_lat.vary) n++;
    if(n == i) return "stsp13_lat";
    if(stsp13_fwhm.defined && stsp13_fwhm.vary) n++;
    if(n == i) return "stsp13_fwhm";
    if(stsp13_tcen.defined && stsp13_tcen.vary) n++;
    if(n == i) return "stsp13_tcen";

    if(stsp21_long.defined && stsp21_long.vary) n++;
    if(n == i) return "stsp21_long";
    if(stsp21_lat.defined  && stsp21_lat.vary) n++;
    if(n == i) return "stsp21_lat";
    if(stsp21_fwhm.defined && stsp21_fwhm.vary) n++;
    if(n == i) return "stsp21_fwhm";
    if(stsp21_tcen.defined && stsp21_tcen.vary) n++;
    if(n == i) return "stsp21_tcen";

    if(stsp22_long.defined && stsp22_long.vary) n++;
    if(n == i) return "stsp22_long";
    if(stsp22_lat.defined  && stsp22_lat.vary) n++;
    if(n == i) return "stsp22_lat";
    if(stsp22_fwhm.defined && stsp22_fwhm.vary) n++;
    if(n == i) return "stsp22_fwhm";
    if(stsp22_tcen.defined && stsp22_tcen.vary) n++;
    if(n == i) return "stsp22_tcen";
    
    if(uesp_long1.defined && uesp_long1.vary) n++;
    if(n == i) return "uesp_long1";
    if(uesp_long2.defined && uesp_long2.vary) n++;
    if(n == i) return "uesp_long2";
    if(uesp_lathw.defined && uesp_lathw.vary) n++;
    if(n == i) return "uesp_lathw";
    if(uesp_taper.defined && uesp_taper.vary) n++;
    if(n == i) return "uesp_taper";
    if(uesp_temp.defined && uesp_temp.vary) n++;
    if(n == i) return "uesp_temp";

    return "UNKNOWN";
}

/** Check legality of a proposed set of variable parameters. It is assumed
 * that the constant parameters are OK
 * \param vpar the parameters to check
 * \return true if not OK
 */
bool Lcurve::Model::is_not_legal(const Subs::Array1D<double>& vpar) const {
    if(nvary() != int(vpar.size()))
        throw Lcurve_Error("are_legal: conflicting numbers of variable parameters");
    int n = 0;

    if(q.vary) {
        if(vpar[n] <= 0. || vpar[n] > 100.) return true;
        n++;
    }
    double xl11 = Roche::xl11(q.value, spin1.value);
    double xl12 = 1.-Roche::xl12(q.value, spin2.value);

    if(iangle.vary) {
        if(vpar[n] < 0. || vpar[n] > 90.) return true;
        n++;
    }

    if(use_radii){

        if(r1.vary) {
            if(vpar[n] <= 0. || vpar[n] >= xl11) return true;
            n++;
        }

        if(r2.vary) {
            if(vpar[n] <= 0. || vpar[n] > xl12) return true;
            n++;
        }

    }else{

        if(cphi3.vary) {
            if(vpar[n] <= 0. || vpar[n] > 0.25) return true;
            n++;
        }

        if(cphi4.vary) {
            if(vpar[n] <= 0. || vpar[n] > 0.25) return true;
            n++;
        }
    }

    if(spin1.vary) {
        if(vpar[n] <= 0. || vpar[n] > 1000.) return true;
        n++;
    }

    if(spin2.vary) {
        if(vpar[n] <= 0. || vpar[n] > 1000.) return true;
        n++;
    }

    if(t1.vary) {
        if(vpar[n] <= 500. || vpar[n] > 1.e6) return true;
        n++;
    }

    if(t2.vary) {
        if(vpar[n] <= 500. || vpar[n] > 1.e6) return true;
        n++;
    }

    if(ldc1_1.vary){
        if(vpar[n] < -1. || vpar[n] > 1.) return true;
        n++;
    }

    if(ldc1_2.vary){
        if(vpar[n] < -1. || vpar[n] > 1.) return true;
        n++;
    }

    if(ldc1_3.vary){
        if(vpar[n] < -1. || vpar[n] > 1.) return true;
        n++;
    }

    if(ldc1_4.vary){
        if(vpar[n] < -1. || vpar[n] > 1.) return true;
        n++;
    }

    if(ldc2_1.vary){
        if(vpar[n] < -1. || vpar[n] > 1.) return true;
        n++;
    }

    if(ldc2_2.vary){
        if(vpar[n] < -1. || vpar[n] > 1.) return true;
        n++;
    }

    if(ldc2_3.vary){
        if(vpar[n] < -1. || vpar[n] > 1.) return true;
        n++;
    }

    if(ldc2_4.vary){
        if(vpar[n] < -1. || vpar[n] > 1.) return true;
        n++;
    }

    if(velocity_scale.vary) n++;

    if(beam_factor1.vary) n++;

    if(beam_factor2.vary) n++;

    if(t0.vary) n++;

    if(period.vary) {
        if(vpar[n] <= 1.e-3 || vpar[n] > 100.) return true;
        n++;
    }

    if(pdot.vary) n++;

    if(deltat.vary) {
        if(vpar[n] <= -1. || vpar[n] > 1.) return true;
        n++;
    }

    if(gravity_dark1.vary){
        if(vpar[n] < -1. || vpar[n] > 1.) return true;
        n++;
    }

    if(gravity_dark2.vary){
        if(vpar[n] < -1. || vpar[n] > 1.) return true;
        n++;
    }

    if(absorb.vary){
        if(vpar[n] < 0. || vpar[n] > 10.) return true;
        n++;
    }

    if(slope.vary){
        if(vpar[n] < -2. || vpar[n] > 2.) return true;
        n++;
    }

    if(quad.vary){
        if(vpar[n] < -2. || vpar[n] > 2.) return true;
        n++;
    }

    if(cube.vary){
        if(vpar[n] < -2. || vpar[n] > 2.) return true;
        n++;
    }

    if(third.vary) n++;

    if(add_disc){

        if(rdisc1.vary){
            if(vpar[n] < 0. || vpar[n] > 1) return true;
            n++;
        }

        if(rdisc2.vary){
            if(vpar[n] < 0 || vpar[n] > 1) return true;
            n++;
        }

        if(height_disc.vary){
            if(vpar[n] < 0) return true;
            n++;
        }

        if(beta_disc.vary){
            if(vpar[n] < 1 || vpar[n] > 100) return true;
            n++;
        }

        if(temp_disc.vary){
            if(vpar[n] < 500 || vpar[n] > 1.e6) return true;
            n++;
        }

        if(texp_disc.vary){
            if(vpar[n] < -100 || vpar[n] > 100) return true;
            n++;
        }

        if(lin_limb_disc.vary){
            if(vpar[n] < 0. || vpar[n] > 1.) return true;
            n++;
        }

        if(quad_limb_disc.vary){
            if(vpar[n] < 0. || vpar[n] > 1.) return true;
            n++;
        }

        if(temp_edge.vary){
            if(vpar[n] <= 0 || vpar[n] > 1.e6) return true;
            n++;
        }

        if(absorb_edge.vary){
            if(vpar[n] < 0 || vpar[n] > 1000) return true;
            n++;
        }
    }

    if(add_spot){

        if(radius_spot.vary){
            if(vpar[n] < 0. || vpar[n] > 1.) return true;
            n++;
        }

        if(length_spot.vary){
            if(vpar[n] < 0. || vpar[n] > 1.) return true;
            n++;
        }

        if(height_spot.vary){
            if(vpar[n] < 0. || vpar[n] > 1.) return true;
            n++;
        }

        if(expon_spot.vary){
            if(vpar[n] <= 0. || vpar[n] > 100.) return true;
            n++;
        }

        if(epow_spot.vary){
            if(vpar[n] <= 0. || vpar[n] > 10.) return true;
            n++;
        }

        if(angle_spot.vary){
            if(vpar[n] < -1000. || vpar[n] > 1000.) return true;
            n++;
        }

        if(yaw_spot.vary){
            if(vpar[n] < -1000. || vpar[n] > 1000.) return true;
            n++;
        }

        if(temp_spot.vary){
            if(vpar[n] < 0. || vpar[n] > 1.e10) return true;
            n++;
        }

        if(tilt_spot.vary){
            if(vpar[n] < -1000. || vpar[n] > 1000.) return true;
            n++;
        }

        if(cfrac_spot.vary){
            if(vpar[n] < 0. || vpar[n] > 1.) return true;
            n++;
        }

    }

    if(stsp11_long.defined && stsp11_long.vary){
        if(vpar[n] < -400. || vpar[n] > 400.) return true;
        n++;
    }
    if(stsp11_lat.defined && stsp11_lat.vary){
        if(vpar[n] < -90. || vpar[n] > 90.) return true;
        n++;
    }
    if(stsp11_fwhm.defined  && stsp11_fwhm.vary){
        if(vpar[n] <= 0. || vpar[n] > 180.) return true;
        n++;
    }
    if(stsp11_tcen.defined  && stsp11_tcen.vary){
        if(vpar[n] <= 0.) return true;
        n++;
    }

    if(stsp12_long.defined && stsp12_long.vary){
        if(vpar[n] < -400. || vpar[n] > 400.) return true;
        n++;
    }
    if(stsp12_lat.defined && stsp12_lat.vary){
        if(vpar[n] < -90. || vpar[n] > 90.) return true;
        n++;
    }
    if(stsp12_fwhm.defined  && stsp12_fwhm.vary){
        if(vpar[n] <= 0. || vpar[n] > 180.) return true;
        n++;
    }
    if(stsp12_tcen.defined  && stsp12_tcen.vary){
        if(vpar[n] <= 0.) return true;
        n++;
    }

    if(stsp13_long.defined && stsp13_long.vary){
        if(vpar[n] < -400. || vpar[n] > 400.) return true;
        n++;
    }
    if(stsp13_lat.defined && stsp13_lat.vary){
        if(vpar[n] < -90. || vpar[n] > 90.) return true;
        n++;
    }
    if(stsp13_fwhm.defined  && stsp13_fwhm.vary){
        if(vpar[n] <= 0. || vpar[n] > 180.) return true;
        n++;
    }
    if(stsp13_tcen.defined  && stsp13_tcen.vary){
        if(vpar[n] <= 0.) return true;
        n++;
    }
    
    if(stsp21_long.defined && stsp21_long.vary){
        if(vpar[n] < -400. || vpar[n] > 400.) return true;
        n++;
    }
    if(stsp21_lat.defined && stsp21_lat.vary){
        if(vpar[n] < -90. || vpar[n] > 90.) return true;
        n++;
    }
    if(stsp21_fwhm.defined  && stsp21_fwhm.vary){
        if(vpar[n] <= 0. || vpar[n] > 180.) return true;
        n++;
    }
    if(stsp21_tcen.defined  && stsp21_tcen.vary){
        if(vpar[n] <= 0.) return true;
        n++;
    }

    if(stsp22_long.defined && stsp22_long.vary){
        if(vpar[n] < -400. || vpar[n] > 400.) return true;
        n++;
    }
    if(stsp22_lat.defined && stsp22_lat.vary){
        if(vpar[n] < -90. || vpar[n] > 90.) return true;
        n++;
    }
    if(stsp22_fwhm.defined  && stsp22_fwhm.vary){
        if(vpar[n] <= 0. || vpar[n] > 180.) return true;
        n++;
    }
    if(stsp22_tcen.defined  && stsp22_tcen.vary){
        if(vpar[n] <= 0.) return true;
        n++;
    }
    
    if(uesp_long1.defined && uesp_long1.vary) {
      if(vpar[n] < -360. || vpar[n] > 360.) return true;
      n++;
    }

    if(uesp_long2.defined && uesp_long2.vary) {
      if(vpar[n] < -360. || vpar[n] > 720.) return true;
      n++;
    }

    if(uesp_lathw.defined && uesp_lathw.vary) {
      if(vpar[n] <= 0. || vpar[n] > 80.) return true;
      n++;
    }
    
    if(uesp_taper.defined && uesp_taper.vary) {
      if(vpar[n] <= 0.) return true;
      n++;
    }
    
    if(uesp_temp.defined && uesp_temp.vary) {
      if(vpar[n] <= 0.) return true;
      n++;
    }
    
    return false;

}

/** This routine returns the current values of all variable
 * parameters.
 * \return an Array1D of the variable parameters only. They come
 * in a standard order
 * q, iangle, r1, r2, etc.
 */
Subs::Array1D<double> Lcurve::Model::get_param() const {

    Subs::Array1D<double> temp;

    if(q.vary) temp.push_back(q.value);
    if(iangle.vary) temp.push_back(iangle.value);
    if(use_radii){
        if(r1.vary)     temp.push_back(r1.value);
        if(r2.vary)     temp.push_back(r2.value);
    }else{
        if(cphi3.vary)  temp.push_back(cphi3.value);
        if(cphi4.vary)  temp.push_back(cphi4.value);
    }
    if(spin1.vary)  temp.push_back(spin1.value);
    if(spin2.vary)  temp.push_back(spin2.value);
    if(t1.vary)     temp.push_back(t1.value);
    if(t2.vary)     temp.push_back(t2.value);
    if(ldc1_1.vary) temp.push_back(ldc1_1.value);
    if(ldc1_2.vary) temp.push_back(ldc1_2.value);
    if(ldc1_3.vary) temp.push_back(ldc1_3.value);
    if(ldc1_4.vary) temp.push_back(ldc1_4.value);
    if(ldc2_1.vary) temp.push_back(ldc2_1.value);
    if(ldc2_2.vary) temp.push_back(ldc2_2.value);
    if(ldc2_3.vary) temp.push_back(ldc2_3.value);
    if(ldc2_4.vary) temp.push_back(ldc2_4.value);
    if(velocity_scale.vary) temp.push_back(velocity_scale.value);
    if(beam_factor1.vary) temp.push_back(beam_factor1.value);
    if(beam_factor2.vary) temp.push_back(beam_factor2.value);

    if(t0.vary) temp.push_back(t0.value);
    if(period.vary) temp.push_back(period.value);
    if(pdot.vary) temp.push_back(pdot.value);
    if(deltat.vary) temp.push_back(deltat.value);
    if(gravity_dark1.vary) temp.push_back(gravity_dark1.value);
    if(gravity_dark2.vary) temp.push_back(gravity_dark2.value);
    if(absorb.vary) temp.push_back(absorb.value);
    if(slope.vary) temp.push_back(slope.value);
    if(quad.vary) temp.push_back(quad.value);
    if(cube.vary) temp.push_back(cube.value);
    if(third.vary) temp.push_back(third.value);

    if(add_disc){
        if(rdisc1.vary) temp.push_back(rdisc1.value);
        if(rdisc2.vary) temp.push_back(rdisc2.value);
        if(height_disc.vary) temp.push_back(height_disc.value);
        if(beta_disc.vary) temp.push_back(beta_disc.value);
        if(temp_disc.vary) temp.push_back(temp_disc.value);
        if(texp_disc.vary) temp.push_back(texp_disc.value);
        if(lin_limb_disc.vary) temp.push_back(lin_limb_disc.value);
        if(quad_limb_disc.vary) temp.push_back(quad_limb_disc.value);
        if(temp_edge.vary) temp.push_back(temp_edge.value);
        if(absorb_edge.vary) temp.push_back(absorb_edge.value);
    }

    if(add_spot){
        if(radius_spot.vary) temp.push_back(radius_spot.value);
        if(length_spot.vary) temp.push_back(length_spot.value);
        if(height_spot.vary) temp.push_back(height_spot.value);
        if(expon_spot.vary)  temp.push_back(expon_spot.value);
        if(epow_spot.vary)   temp.push_back(epow_spot.value);
        if(angle_spot.vary)  temp.push_back(angle_spot.value);
        if(yaw_spot.vary)    temp.push_back(yaw_spot.value);
        if(temp_spot.vary)   temp.push_back(temp_spot.value);
        if(tilt_spot.vary)   temp.push_back(tilt_spot.value);
        if(cfrac_spot.vary)  temp.push_back(cfrac_spot.value);
    }

    if(stsp11_long.defined && stsp11_long.vary)
        temp.push_back(stsp11_long.value);
    if(stsp11_lat.defined && stsp11_lat.vary)
        temp.push_back(stsp11_lat.value);
    if(stsp11_fwhm.defined && stsp11_fwhm.vary)
        temp.push_back(stsp11_fwhm.value);
    if(stsp11_tcen.defined && stsp11_tcen.vary)
        temp.push_back(stsp11_tcen.value);

    if(stsp12_long.defined && stsp12_long.vary)
        temp.push_back(stsp12_long.value);
    if(stsp12_lat.defined && stsp12_lat.vary)
        temp.push_back(stsp12_lat.value);
    if(stsp12_fwhm.defined && stsp12_fwhm.vary)
        temp.push_back(stsp12_fwhm.value);
    if(stsp12_tcen.defined && stsp12_tcen.vary)
        temp.push_back(stsp12_tcen.value);

    if(stsp13_long.defined && stsp13_long.vary)
        temp.push_back(stsp13_long.value);
    if(stsp13_lat.defined && stsp13_lat.vary)
        temp.push_back(stsp13_lat.value);
    if(stsp13_fwhm.defined && stsp13_fwhm.vary)
        temp.push_back(stsp13_fwhm.value);
    if(stsp13_tcen.defined && stsp13_tcen.vary)
        temp.push_back(stsp13_tcen.value);
    
    if(stsp21_long.defined && stsp21_long.vary)
        temp.push_back(stsp21_long.value);
    if(stsp21_lat.defined && stsp21_lat.vary)
        temp.push_back(stsp21_lat.value);
    if(stsp21_fwhm.defined && stsp21_fwhm.vary)
        temp.push_back(stsp21_fwhm.value);
    if(stsp21_tcen.defined && stsp21_tcen.vary)
        temp.push_back(stsp21_tcen.value);

    if(stsp22_long.defined && stsp22_long.vary)
        temp.push_back(stsp22_long.value);
    if(stsp22_lat.defined && stsp22_lat.vary)
        temp.push_back(stsp22_lat.value);
    if(stsp22_fwhm.defined && stsp22_fwhm.vary)
        temp.push_back(stsp22_fwhm.value);
    if(stsp22_tcen.defined && stsp22_tcen.vary)
        temp.push_back(stsp22_tcen.value);
    
    if(uesp_long1.defined && uesp_long1.vary)
      temp.push_back(uesp_long1.value);
    if(uesp_long2.defined && uesp_long2.vary) 
      temp.push_back(uesp_long2.value);
    if(uesp_lathw.defined && uesp_lathw.vary) 
      temp.push_back(uesp_lathw.value);
    if(uesp_taper.defined && uesp_taper.vary) 
      temp.push_back(uesp_taper.value);
    if(uesp_temp.defined && uesp_temp.vary) 
      temp.push_back(uesp_temp.value);

    return temp;
}

/** This routine returns the values of all the ranges of the variable
 * parameters.
 * \return an Array1D of the range of the variable parameters only.
 * They come in a standard order
 * q, iangle, r1, r2, etc.
 */
Subs::Array1D<double> Lcurve::Model::get_range() const {

    Subs::Array1D<double> temp;

    if(q.vary) temp.push_back(q.range);
    if(iangle.vary) temp.push_back(iangle.range);
    if(use_radii){
        if(r1.vary) temp.push_back(r1.range);
        if(r2.vary) temp.push_back(r2.range);
    }else{
        if(cphi3.vary) temp.push_back(cphi3.range);
        if(cphi4.vary) temp.push_back(cphi4.range);
    }
    if(spin1.vary) temp.push_back(spin1.range);
    if(spin2.vary) temp.push_back(spin2.range);
    if(t1.vary) temp.push_back(t1.range);
    if(t2.vary) temp.push_back(t2.range);
    if(ldc1_1.vary) temp.push_back(ldc1_1.range);
    if(ldc1_2.vary) temp.push_back(ldc1_2.range);
    if(ldc1_3.vary) temp.push_back(ldc1_3.range);
    if(ldc1_4.vary) temp.push_back(ldc1_4.range);
    if(ldc2_1.vary) temp.push_back(ldc2_1.range);
    if(ldc2_2.vary) temp.push_back(ldc2_2.range);
    if(ldc2_3.vary) temp.push_back(ldc2_3.range);
    if(ldc2_4.vary) temp.push_back(ldc2_4.range);
    if(velocity_scale.vary) temp.push_back(velocity_scale.range);
    if(beam_factor1.vary) temp.push_back(beam_factor1.range);
    if(beam_factor2.vary) temp.push_back(beam_factor2.range);

    if(t0.vary) temp.push_back(t0.range);
    if(period.vary) temp.push_back(period.range);
    if(pdot.vary) temp.push_back(pdot.range);
    if(deltat.vary) temp.push_back(deltat.range);
    if(gravity_dark1.vary) temp.push_back(gravity_dark1.range);
    if(gravity_dark2.vary) temp.push_back(gravity_dark2.range);
    if(absorb.vary) temp.push_back(absorb.range);
    if(slope.vary) temp.push_back(slope.range);
    if(quad.vary) temp.push_back(quad.range);
    if(cube.vary) temp.push_back(cube.range);
    if(third.vary) temp.push_back(third.range);

    if(add_disc){
        if(rdisc1.vary) temp.push_back(rdisc1.range);
        if(rdisc2.vary) temp.push_back(rdisc2.range);
        if(height_disc.vary) temp.push_back(height_disc.range);
        if(beta_disc.vary) temp.push_back(beta_disc.range);
        if(temp_disc.vary) temp.push_back(temp_disc.range);
        if(texp_disc.vary) temp.push_back(texp_disc.range);
        if(lin_limb_disc.vary) temp.push_back(lin_limb_disc.range);
        if(quad_limb_disc.vary) temp.push_back(quad_limb_disc.range);
        if(temp_edge.vary) temp.push_back(temp_edge.range);
        if(absorb_edge.vary) temp.push_back(absorb_edge.range);
    }

    if(add_spot){
        if(radius_spot.vary) temp.push_back(radius_spot.range);
        if(length_spot.vary) temp.push_back(length_spot.range);
        if(height_spot.vary) temp.push_back(height_spot.range);
        if(expon_spot.vary)  temp.push_back(expon_spot.range);
        if(epow_spot.vary)   temp.push_back(epow_spot.range);
        if(angle_spot.vary)  temp.push_back(angle_spot.range);
        if(yaw_spot.vary)    temp.push_back(yaw_spot.range);
        if(temp_spot.vary)   temp.push_back(temp_spot.range);
        if(tilt_spot.vary)   temp.push_back(tilt_spot.range);
        if(cfrac_spot.vary)  temp.push_back(cfrac_spot.range);
    }

    if(stsp11_long.defined && stsp11_long.vary)
        temp.push_back(stsp11_long.range);
    if(stsp11_lat.defined && stsp11_lat.vary)
        temp.push_back(stsp11_lat.range);
    if(stsp11_fwhm.defined && stsp11_fwhm.vary)
        temp.push_back(stsp11_fwhm.range);
    if(stsp11_tcen.defined && stsp11_tcen.vary)
        temp.push_back(stsp11_tcen.range);

    if(stsp12_long.defined && stsp12_long.vary)
        temp.push_back(stsp12_long.range);
    if(stsp12_lat.defined && stsp12_lat.vary)
        temp.push_back(stsp12_lat.range);
    if(stsp12_fwhm.defined && stsp12_fwhm.vary)
        temp.push_back(stsp12_fwhm.range);
    if(stsp12_tcen.defined && stsp12_tcen.vary)
        temp.push_back(stsp12_tcen.range);

    if(stsp13_long.defined && stsp13_long.vary)
        temp.push_back(stsp13_long.range);
    if(stsp13_lat.defined && stsp13_lat.vary)
        temp.push_back(stsp13_lat.range);
    if(stsp13_fwhm.defined && stsp13_fwhm.vary)
        temp.push_back(stsp13_fwhm.range);
    if(stsp13_tcen.defined && stsp13_tcen.vary)
        temp.push_back(stsp13_tcen.range);

    if(stsp21_long.defined && stsp21_long.vary)
        temp.push_back(stsp21_long.range);
    if(stsp21_lat.defined && stsp21_lat.vary)
        temp.push_back(stsp21_lat.range);
    if(stsp21_fwhm.defined && stsp21_fwhm.vary)
        temp.push_back(stsp21_fwhm.range);
    if(stsp21_tcen.defined && stsp21_tcen.vary)
        temp.push_back(stsp21_tcen.range);

    if(stsp22_long.defined && stsp22_long.vary)
        temp.push_back(stsp22_long.range);
    if(stsp22_lat.defined && stsp22_lat.vary)
        temp.push_back(stsp22_lat.range);
    if(stsp22_fwhm.defined && stsp22_fwhm.vary)
        temp.push_back(stsp22_fwhm.range);
    if(stsp22_tcen.defined && stsp22_tcen.vary)
        temp.push_back(stsp22_tcen.range);
    
    if(uesp_long1.defined && uesp_long1.vary)
      temp.push_back(uesp_long1.range);
    if(uesp_long2.defined && uesp_long2.vary) 
      temp.push_back(uesp_long2.range);
    if(uesp_lathw.defined && uesp_lathw.vary) 
      temp.push_back(uesp_lathw.range);
    if(uesp_taper.defined && uesp_taper.vary) 
      temp.push_back(uesp_taper.range);
    if(uesp_temp.defined && uesp_temp.vary) 
      temp.push_back(uesp_temp.range);

    return temp;
}

/** This routine returns the step sizes for calculating derivatives 
 * of all variable
 * parameters.
 * \return an Array1D of the setp sizes for variable parameters only. They come in a standard order
 * q, iangle, r1, r2, etc.
 */
Subs::Array1D<double> Lcurve::Model::get_dstep() const {

    Subs::Array1D<double> temp;

    if(q.vary) temp.push_back(q.dstep);
    if(iangle.vary) temp.push_back(iangle.dstep);
    if(use_radii){
        if(r1.vary)     temp.push_back(r1.dstep);
        if(r2.vary)     temp.push_back(r2.dstep);
    }else{
        if(cphi3.vary) temp.push_back(cphi3.dstep);
        if(cphi4.vary) temp.push_back(cphi4.dstep);
    }
    if(spin1.vary) temp.push_back(spin1.dstep);
    if(spin2.vary) temp.push_back(spin2.dstep);
    if(t1.vary) temp.push_back(t1.dstep);
    if(t2.vary) temp.push_back(t2.dstep);
    if(ldc1_1.vary) temp.push_back(ldc1_1.dstep);
    if(ldc1_2.vary) temp.push_back(ldc1_2.dstep);
    if(ldc1_3.vary) temp.push_back(ldc1_3.dstep);
    if(ldc1_4.vary) temp.push_back(ldc1_4.dstep);
    if(ldc2_1.vary) temp.push_back(ldc2_1.dstep);
    if(ldc2_2.vary) temp.push_back(ldc2_2.dstep);
    if(ldc2_3.vary) temp.push_back(ldc2_3.dstep);
    if(ldc2_4.vary) temp.push_back(ldc2_4.dstep);
    if(velocity_scale.vary) temp.push_back(velocity_scale.dstep);
    if(beam_factor1.vary) temp.push_back(beam_factor1.dstep);
    if(beam_factor2.vary) temp.push_back(beam_factor2.dstep);

    if(t0.vary) temp.push_back(t0.dstep);
    if(period.vary) temp.push_back(period.dstep);
    if(pdot.vary) temp.push_back(pdot.dstep);
    if(deltat.vary) temp.push_back(deltat.dstep);
    if(gravity_dark1.vary) temp.push_back(gravity_dark1.dstep);
    if(gravity_dark2.vary) temp.push_back(gravity_dark2.dstep);
    if(absorb.vary) temp.push_back(absorb.dstep);
    if(slope.vary) temp.push_back(slope.dstep);
    if(quad.vary) temp.push_back(quad.dstep);
    if(cube.vary) temp.push_back(cube.dstep);
    if(third.vary) temp.push_back(third.dstep);

    if(add_disc){
        if(rdisc1.vary) temp.push_back(rdisc1.dstep);
        if(rdisc2.vary) temp.push_back(rdisc2.dstep);
        if(height_disc.vary) temp.push_back(height_disc.dstep);
        if(beta_disc.vary) temp.push_back(beta_disc.dstep);
        if(temp_disc.vary) temp.push_back(temp_disc.dstep);
        if(texp_disc.vary) temp.push_back(texp_disc.dstep);
        if(lin_limb_disc.vary) temp.push_back(lin_limb_disc.dstep);
        if(quad_limb_disc.vary) temp.push_back(quad_limb_disc.dstep);
        if(temp_edge.vary) temp.push_back(temp_edge.dstep);
        if(absorb_edge.vary) temp.push_back(absorb_edge.dstep);
    }

    if(add_spot){
        if(radius_spot.vary) temp.push_back(radius_spot.dstep);
        if(length_spot.vary) temp.push_back(length_spot.dstep);
        if(height_spot.vary) temp.push_back(height_spot.dstep);
        if(expon_spot.vary) temp.push_back(expon_spot.dstep);
        if(epow_spot.vary) temp.push_back(epow_spot.dstep);
        if(angle_spot.vary) temp.push_back(angle_spot.dstep);
        if(yaw_spot.vary) temp.push_back(yaw_spot.dstep);
        if(temp_spot.vary) temp.push_back(temp_spot.dstep);
        if(tilt_spot.vary) temp.push_back(tilt_spot.dstep);
        if(cfrac_spot.vary) temp.push_back(cfrac_spot.dstep);
    }

    if(stsp11_long.defined && stsp11_long.vary)
        temp.push_back(stsp11_long.dstep);
    if(stsp11_lat.defined && stsp11_lat.vary)
        temp.push_back(stsp11_lat.dstep);
    if(stsp11_fwhm.defined && stsp11_fwhm.vary)
        temp.push_back(stsp11_fwhm.dstep);
    if(stsp11_tcen.defined && stsp11_tcen.vary)
        temp.push_back(stsp11_tcen.dstep);

    if(stsp12_long.defined && stsp12_long.vary)
        temp.push_back(stsp12_long.dstep);
    if(stsp12_lat.defined && stsp12_lat.vary)
        temp.push_back(stsp12_lat.dstep);
    if(stsp12_fwhm.defined && stsp12_fwhm.vary)
        temp.push_back(stsp12_fwhm.dstep);
    if(stsp12_tcen.defined && stsp12_tcen.vary)
        temp.push_back(stsp12_tcen.dstep);

    if(stsp13_long.defined && stsp13_long.vary)
        temp.push_back(stsp13_long.dstep);
    if(stsp13_lat.defined && stsp13_lat.vary)
        temp.push_back(stsp13_lat.dstep);
    if(stsp13_fwhm.defined && stsp13_fwhm.vary)
        temp.push_back(stsp13_fwhm.dstep);
    if(stsp13_tcen.defined && stsp13_tcen.vary)
        temp.push_back(stsp13_tcen.dstep);
    
    if(stsp21_long.defined && stsp21_long.vary)
        temp.push_back(stsp21_long.dstep);
    if(stsp21_lat.defined && stsp21_lat.vary)
        temp.push_back(stsp21_lat.dstep);
    if(stsp21_fwhm.defined && stsp21_fwhm.vary)
        temp.push_back(stsp21_fwhm.dstep);
    if(stsp21_tcen.defined && stsp21_tcen.vary)
        temp.push_back(stsp21_tcen.dstep);

    if(stsp22_long.defined && stsp22_long.vary)
        temp.push_back(stsp22_long.dstep);
    if(stsp22_lat.defined && stsp22_lat.vary)
        temp.push_back(stsp22_lat.dstep);
    if(stsp22_fwhm.defined && stsp22_fwhm.vary)
        temp.push_back(stsp22_fwhm.dstep);
    if(stsp22_tcen.defined && stsp22_tcen.vary)
        temp.push_back(stsp22_tcen.dstep);
    
    if(uesp_long1.defined && uesp_long1.vary)
      temp.push_back(uesp_long1.dstep);
    if(uesp_long2.defined && uesp_long2.vary) 
      temp.push_back(uesp_long2.dstep);
    if(uesp_lathw.defined && uesp_lathw.vary) 
      temp.push_back(uesp_lathw.dstep);
    if(uesp_taper.defined && uesp_taper.vary) 
      temp.push_back(uesp_taper.dstep);
    if(uesp_temp.defined && uesp_temp.vary) 
      temp.push_back(uesp_temp.dstep);
    
    return temp;
}

/** Outputs the current values of the model parameters */
std::ostream& Lcurve::operator<<(std::ostream& s, const Model& model){
    s << "q              = " << model.q              << "\n";
    s << "iangle         = " << model.iangle         << "\n";
    s << "r1             = " << model.r1             << "\n";
    s << "r2             = " << model.r2             << "\n";
    s << "cphi3          = " << model.cphi3          << "\n";
    s << "cphi4          = " << model.cphi4          << "\n";
    s << "spin1          = " << model.spin1          << "\n";
    s << "spin2          = " << model.spin2          << "\n";
    s << "t1             = " << model.t1             << "\n";
    s << "t2             = " << model.t2             << "\n";
    s << "ldc1_1         = " << model.ldc1_1         << "\n";
    s << "ldc1_2         = " << model.ldc1_2         << "\n";
    s << "ldc1_3         = " << model.ldc1_3         << "\n";
    s << "ldc1_4         = " << model.ldc1_4         << "\n";
    s << "ldc2_1         = " << model.ldc2_1         << "\n";
    s << "ldc2_2         = " << model.ldc2_2         << "\n";
    s << "ldc2_3         = " << model.ldc2_3         << "\n";
    s << "ldc2_4         = " << model.ldc2_4         << "\n";
    s << "velocity_scale = " << model.velocity_scale << "\n";
    s << "beam_factor1   = " << model.beam_factor1    << "\n";
    s << "beam_factor2   = " << model.beam_factor2    << "\n\n";

    s << "t0             = " << model.t0             << "\n";
    s << "period         = " << model.period         << "\n";
    s << "pdot           = " << model.pdot           << "\n";
    s << "deltat         = " << model.deltat         << "\n";
    s << "gravity_dark1  = " << model.gravity_dark1  << "\n";
    s << "gravity_dark2  = " << model.gravity_dark2  << "\n";
    s << "absorb         = " << model.absorb         << "\n";
    s << "slope          = " << model.slope          << "\n\n";
    s << "quad           = " << model.quad           << "\n\n";
    s << "cube           = " << model.cube           << "\n\n";
    s << "third          = " << model.third          << "\n\n";

    s << "rdisc1         = " << model.rdisc1         << "\n";
    s << "rdisc2         = " << model.rdisc2         << "\n";
    s << "height_disc    = " << model.height_disc    << "\n";
    s << "beta_disc      = " << model.beta_disc      << "\n";
    s << "temp_disc      = " << model.temp_disc      << "\n";
    s << "texp_disc      = " << model.texp_disc      << "\n";
    s << "lin_limb_disc  = " << model.lin_limb_disc  << "\n";
    s << "quad_limb_disc = " << model.quad_limb_disc << "\n";
    s << "temp_edge      = " << model.temp_edge      << "\n";
    s << "absorb_edge    = " << model.absorb_edge    << "\n\n";

    s << "radius_spot    = " << model.radius_spot    << "\n";
    s << "length_spot    = " << model.length_spot    << "\n";
    s << "height_spot    = " << model.height_spot    << "\n";
    s << "expon_spot     = " << model.expon_spot     << "\n";
    s << "epow_spot      = " << model.epow_spot      << "\n";
    s << "angle_spot     = " << model.angle_spot     << "\n";
    s << "yaw_spot       = " << model.yaw_spot       << "\n";
    s << "temp_spot      = " << model.temp_spot      << "\n";
    s << "tilt_spot      = " << model.tilt_spot      << "\n";
    s << "cfrac_spot     = " << model.cfrac_spot     << "\n\n";

    s << "stsp11_long    = " << model.stsp11_long    << "\n";
    s << "stsp11_lat     = " << model.stsp11_lat     << "\n";
    s << "stsp11_fwhm    = " << model.stsp11_fwhm    << "\n";
    s << "stsp11_tcen    = " << model.stsp11_tcen    << "\n\n";

    s << "stsp12_long    = " << model.stsp12_long    << "\n";
    s << "stsp12_lat     = " << model.stsp12_lat     << "\n";
    s << "stsp12_fwhm    = " << model.stsp12_fwhm    << "\n";
    s << "stsp12_tcen    = " << model.stsp12_tcen    << "\n\n";

    s << "stsp13_long    = " << model.stsp13_long    << "\n";
    s << "stsp13_lat     = " << model.stsp13_lat     << "\n";
    s << "stsp13_fwhm    = " << model.stsp13_fwhm    << "\n";
    s << "stsp13_tcen    = " << model.stsp13_tcen    << "\n\n";
    
    s << "stsp21_long    = " << model.stsp21_long    << "\n";
    s << "stsp21_lat     = " << model.stsp21_lat     << "\n";
    s << "stsp21_fwhm    = " << model.stsp21_fwhm    << "\n";
    s << "stsp21_tcen    = " << model.stsp21_tcen    << "\n\n";

    s << "stsp22_long    = " << model.stsp22_long    << "\n";
    s << "stsp22_lat     = " << model.stsp22_lat     << "\n";
    s << "stsp22_fwhm    = " << model.stsp22_fwhm    << "\n";
    s << "stsp22_tcen    = " << model.stsp22_tcen    << "\n\n";

    s << "uesp_long1     = " << model.uesp_long1   << "\n";
    s << "uesp_long2     = " << model.uesp_long2   << "\n";
    s << "uesp_lathw     = " << model.uesp_lathw   << "\n";
    s << "uesp_taper     = " << model.uesp_taper   << "\n";
    s << "uesp_temp      = " << model.uesp_temp    << "\n\n";

    s << "delta_phase    = " << model.delta_phase    << "\n";
    s << "nlat1f         = " << model.nlat1f         << "\n";
    s << "nlat2f         = " << model.nlat2f         << "\n";
    s << "nlat1c         = " << model.nlat1c         << "\n";
    s << "nlat2c         = " << model.nlat2c         << "\n";
    s << "npole          = " << model.npole          << "\n";
    s << "nlatfill       = " << model.nlatfill       << "\n";
    s << "nlngfill       = " << model.nlngfill       << "\n";
    s << "lfudge         = " << model.lfudge         << "\n";
    s << "llo            = " << model.llo            << "\n";
    s << "lhi            = " << model.lhi            << "\n";
    s << "phase1         = " << model.phase1         << "\n";
    s << "phase2         = " << model.phase2         << "\n";

    s << "wavelength     = " << model.wavelength     << "\n";
    s << "roche1         = " << model.roche1         << "\n";
    s << "roche2         = " << model.roche2         << "\n";
    s << "eclipse1       = " << model.eclipse1       << "\n";
    s << "eclipse2       = " << model.eclipse2       << "\n";
    s << "glens1         = " << model.glens1         << "\n";
    s << "use_radii      = " << model.use_radii      << "\n";

    s << "tperiod        = " << model.tperiod        << "\n";
    s << "gdark_bolom1   = " << model.gdark_bolom1   << "\n";
    s << "gdark_bolom2   = " << model.gdark_bolom2   << "\n";
    s << "mucrit1        = " << model.mucrit1        << "\n";
    s << "mucrit2        = " << model.mucrit2        << "\n";
    if(model.limb1 == LDC::POLY){
        s << "limb1          = Poly\n";
    }else if(model.limb1 == LDC::CLARET){
        s << "limb1          = Claret\n";
    }
    if(model.limb2 == LDC::POLY){
        s << "limb2          = Poly\n";
    }else if(model.limb2 == LDC::CLARET){
        s << "limb2          = Claret\n";
    }
    s << "mirror         = " << model.mirror         << "\n";
    s << "add_disc       = " << model.add_disc       << "\n";
    s << "nrad           = " << model.nrad           << "\n";
    s << "opaque         = " << model.opaque         << "\n";
    s << "add_spot       = " << model.add_spot       << "\n";
    s << "nspot          = " << model.nspot          << "\n";
    s << "iscale         = " << model.iscale         << "\n";
    return s;
}

/** Writes model out to an ASCII file
 */
void Lcurve::Model::wrasc(const std::string& file) const {

    std::ofstream fout(file.c_str());
    if(!fout)
        throw Lcurve_Error("Lcurve::Model::wrasc: failed to open " +
                           file + " for output.");

    fout << *this;
    fout.close();
}

/** ASCII input of a Datum. This expects the entries time,
 * exposure length, flux, error on flux, weight factor and ndiv
 * separated by spaces.
 */
std::istream& Lcurve::operator>>(std::istream& s, Datum& datum){
    std::string buff;
    getline(s, buff);
    std::istringstream istr(buff);
    istr >> datum.time >> datum.expose >> datum.flux >> datum.ferr >>
        datum.weight >> datum.ndiv;
    if(!istr)
        throw Lcurve_Error("Lcurve::operator>>: failed to read t,e,f,fe,wgt,ndiv of a Datum");
    if(datum.ferr < 0.){
        datum.ferr -= datum.ferr;
        datum.weight = 0.;
    }
    return s;
}


/** ASCII output of a Datum.
 */
std::ostream& Lcurve::operator<<(std::ostream& s, const Datum& datum) {
    static Subs::Format tform(17), dform(10);
    s << tform(datum.time) << " " << datum.expose << " " << dform(datum.flux)
      << " " << datum.ferr << " " << datum.weight << " " << datum.ndiv;
    return s;
}

/** Reads in data from a file in the form of a series of lines
 * specifying each 'Datum'. Lines starting with # or blank are ignored
 * \param file the ASCII file to load.
 */
Lcurve::Data::Data(const std::string& file) : std::vector<Datum>() {
    this->rasc(file);
}

/** Reads in data from a file in the form of a series of lines
 * specifying each 'Datum'. Lines starting with # or blank are ignored
 * \param file the ASCII file to load.
 */
void Lcurve::Data::rasc(const std::string& file) {

    // Read in the parameter values
    std::ifstream fin(file.c_str());
    if(!fin)
        throw Lcurve_Error("Lcurve::Data::Data: failed to open " + file +
                           " for data.");

    this->clear();

    const int MAX_LINE = 5000;
    int n = 0;
    char ch;
    Datum datum;
    while(fin && !fin.eof()){
        n++;
        ch = fin.peek();
        if(fin.eof()){
            std::cout << "End of file reached." << std::endl;
        }else if(ch == '#' || ch == ' ' || ch == '\t' || ch == '\n'){
            fin.ignore(MAX_LINE, '\n'); // skip to next line
        }else{
            if(fin >> datum){
                this->push_back(datum);
            }else if(fin.eof()){
                std::cout << "End of data file reached." << std::endl;
            }else{
                throw Lcurve_Error("Data file input failure on line " +
                                   Subs::str(n));
            }
        }
    }
    fin.close();

    std::cout << this->size() << " lines of data read from " << file
              << "\n" << std::endl;

}

/** Writes to an ASCII file with a series of rows
 * each with a time, exposure, flux and uncertainty
 */
void Lcurve::Data::wrasc(const std::string& file) const {

    std::ofstream fout(file.c_str());
    if(!fout)
        throw Lcurve_Error("Lcurve::Data::wrasc: failed to open " +
                           file + " for output.");

    for(size_t i=0; i<this->size(); i++)
        fout << (*this)[i] << std::endl;

}

/** Operator of function object to allow it to be called as func(vpar).
 * \param vpar vector of numbers of the variable parameters in the order
 * in which the parameters occur in Model
 * \return Returns chi**2
 */
double Lcurve::Fobj::operator()(const Subs::Array1D<double>& vpar){
    // Load up the variable parameters
    model.set_param(vpar);
    return chisq();
}

/** Access the fit values
 */
double Lcurve::Fobj::operator[](int n) const {
    return fit[n];
}

/** Computes chisq for current model and stores fit
 * \return Returns chi**2
 */
double Lcurve::Fobj::chisq() {

    Subs::Buffer1D<double> sfac(4);
    double wdwarf, chisq, wnok, logg1, logg2, rv1, rv2;
    Lcurve::light_curve_comp(model, data, true, true, false, sfac, fit, wdwarf,
                             chisq, wnok, logg1, logg2, rv1, rv2);
    if(wnok == 0.0)
        throw Lcurve::Lcurve_Error("Lcurve::Fobj::chisq: no good data");
    Lcurve::Fobj::neval++;

    Subs::Format form(12);
    if(Lcurve::Fobj::neval == 1 || chisq < Lcurve::Fobj::chisq_min){
        Lcurve::Fobj::chisq_min = chisq;
        Lcurve::Fobj::scale_min = sfac;
        if(model.iscale){
            std::cout << "Weighted chi**2 = " << form(chisq)
                      << ", scale factors = ";
            std::cout << form(sfac[0]) << ", " <<  form(sfac[1])
                      << ", " <<  form(sfac[2]) << ", "
                      << form(sfac[3]) << ", ";
            std::cout << ", neval = " << neval << ", wnok = "
                      << form(wnok) << std::endl;
        }else{
            std::cout << "Weighted chi**2 = " << form(chisq)
                      << ", scale factor = " << form(sfac[0])
                      << ", neval = " << neval << ", wnok = "
                      << form(wnok) << std::endl;
        }
    }
    return chisq;
}

std::istream& Lcurve::operator>>(std::istream& s, Point& p) {
    throw Lcurve_Error("Attempt to use operator>>(std::istream& s, Point&) is an error");
}

//! Dummy ASCII output operator to allow use of Buffer1D
std::ostream& Lcurve::operator<<(std::ostream& s, const Point& p) {
    throw Lcurve_Error("Attempt to use operator<<(std::ostream& s, const Point&) is an error");
}
