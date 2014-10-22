//boost include files
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/errors.hpp>

//gespin include files
#include "Nucleon.h"
#include "NucleonCollection.h"

using namespace boost::python;

//for general use in deep copying
template<typename T> const T DeepCopy(const T& v, dict d) { return T(v); }

/***************** Nucleon helpers **********************/
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
/***************** Nucleon identities ***********************/
    enum_<NucleonIdentity>("nucleon_identities")
        .value("unspecified", NucleonIdentity::unspecified)
        .value("proton", NucleonIdentity::proton)
        .value("neutron", NucleonIdentity::neutron)
        .value("antiproton", NucleonIdentity::antiproton)
        .value("antineutron", NucleonIdentity::antineutron)
    ;
/********************************************************/

/***************** Nucleon class ***********************/
    class_<Nucleon>("Nucleon", init<double, double, double>())
        .def(init<>())
        .def("__deepcopy__", &DeepCopy<Nucleon>)
        .add_property("x", &Nucleon::X, &Nucleon::SetX)
        .add_property("y", &Nucleon::Y, &Nucleon::SetY)
        .add_property("z", &Nucleon::Z, &Nucleon::SetZ)
        .add_property("r", &Nucleon::R, &Nucleon::SetR)
        .add_property("theta", &Nucleon::Theta, &Nucleon::SetTheta)
        .add_property("phi", &Nucleon::Phi, &Nucleon::SetPhi)
        .add_property("position", &NucleonHelpers::GetPosition, &NucleonHelpers::SetPosition)
        .add_property("radius", &Nucleon::Radius, &Nucleon::SetRadius)
        .add_property("identity", &Nucleon::Identity, &Nucleon::SetIdentity)
    ;
/********************************************************/
}
