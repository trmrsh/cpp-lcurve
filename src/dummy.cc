#include <map>
#include "trm_lcurve.h"
#include "trm_roche.h"
#include "trm_format.h"

using Subs::operator+;

std::istream& Lcurve::operator>>(std::istream& s, Point& p) {
  throw Lcurve_Error("Attempt to use operator>>(std::istream& s, Point&) is an error");
}

std::ostream& Lcurve::operator<<(std::ostream& s, const Point& p) {
  throw Lcurve_Error("Attempt to use operator<<(std::ostream& s, const Point&) is an error");
}
