#include <map>
#include "trm/lcurve.h"
#include "trm/roche.h"
#include "trm/format.h"

using Subs::operator+;

std::istream& Lcurve::operator>>(std::istream& s, Point& p) {
  throw Lcurve_Error("Attempt to use operator>>(std::istream& s, Point&) is an error");
}

std::ostream& Lcurve::operator<<(std::ostream& s, const Point& p) {
  throw Lcurve_Error("Attempt to use operator<<(std::ostream& s, const Point&) is an error");
}
