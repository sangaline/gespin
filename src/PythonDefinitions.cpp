//boost include files
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/tuple.hpp>

//gespin include files
#include "Nucleon.h"

using namespace boost::python;

/***************** Nucleon helpers *********************/
namespace NucleonHelpers {
    tuple Position(Nucleon& nucleon) {
        return make_tuple(nucleon.X(), nucleon.Y(), nucleon.Z());
    }
}
/********************************************************/

BOOST_PYTHON_MODULE(core)
{
/***************** Nucleon class ***********************/
    class_<Nucleon>("Nucleon", init<double, double, double>())
        .def(init<>())
        .def("X", &Nucleon::X)
        .def("Y", &Nucleon::Y)
        .def("Z", &Nucleon::Z)
        .def("position", &NucleonHelpers::Position)
    ;
/********************************************************/
}
