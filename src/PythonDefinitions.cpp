//boost include files
#include <boost/python.hpp>
#include <boost/python/module.hpp>

//gespin include files
#include "Nucleon.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(core)
{
/***************** Nucleon class ***********************/
    class_<Nucleon>("Nucleon", init<double, double, double>())
        .def(init<>())
        .def("X", &Nucleon::X)
        .def("Y", &Nucleon::Y)
        .def("Z", &Nucleon::Z)
    ;
/********************************************************/
}
