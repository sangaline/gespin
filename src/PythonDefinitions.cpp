//boost include files
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/extract.hpp>

//gespin include files
#include "Nucleon.h"

using namespace boost::python;

/***************** Nucleon helpers *********************/
namespace NucleonHelpers {
    tuple GetPosition(Nucleon& nucleon) {
        return make_tuple(nucleon.X(), nucleon.Y(), nucleon.Z());
    }
    void SetPosition(Nucleon& nucleon, tuple position) {
        nucleon.SetPosition(extract<float>(position[0]),
                            extract<float>(position[1]),
                            extract<float>(position[2]));
    }
}
/********************************************************/

BOOST_PYTHON_MODULE(core)
{
/***************** Nucleon class ***********************/
    class_<Nucleon>("Nucleon", init<double, double, double>())
        .def(init<>())
        .add_property("x", &Nucleon::X, &Nucleon::SetX)
        .add_property("y", &Nucleon::Y, &Nucleon::SetY)
        .add_property("z", &Nucleon::Z, &Nucleon::SetZ)
        .add_property("position", &NucleonHelpers::GetPosition, &NucleonHelpers::SetPosition)
    ;
/********************************************************/
}
