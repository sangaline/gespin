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

/***************** NucleonCollection helpers ************/
namespace NucleonCollectionHelpers {
    //the return list is by proxy
    list GetNucleonList(NucleonCollection& nucleon_collection) {
        list l;
        for(unsigned int i = 0; i < nucleon_collection.NucleonCount(); i++) {
            l.append(Nucleon(nucleon_collection[i]));
        }
        return l;
    }

    //the setting is by proxy as well
    void SetNucleonList(NucleonCollection& nucleon_collection, list l) {
        nucleon_collection.Reset();
        for(unsigned int i = 0; i < len(l); i++) {
            nucleon_collection.AddNucleon(extract<Nucleon>(l[i]));
        }
    }

    //[] access is by reference, so nucleon_collection[0].x += 5 works
    Nucleon& GetNucleon(NucleonCollection& nucleon_collection, int index) {
        if(index >= int(nucleon_collection.NucleonCount())) {
            PyErr_SetString(PyExc_StopIteration, "No more data.");
            throw_error_already_set();
        }
        if(index < 0) {
            index = index % nucleon_collection.NucleonCount();
        }
        return nucleon_collection[index];
    }

    void SetNucleon(NucleonCollection& nucleon_collection, int index, Nucleon& nucleon) {
        //the nucleon assignment operator must preserve the parent
        nucleon_collection[index] = nucleon;
    }
};
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

/***************** Nucleon Collection *******************/
    class_<NucleonCollection>("NucleonCollection", init<>())
        .def(init<double, unsigned int, double>())
        .def("__len__", &NucleonCollection::NucleonCount)
        .def("__getitem__", &NucleonCollectionHelpers::GetNucleon, boost::python::return_internal_reference<>())
        .def("__setitem__", &NucleonCollectionHelpers::SetNucleon)
        .def("__deepcopy__", &DeepCopy<NucleonCollection>)
        .def("append", &NucleonCollection::AddNucleon)
        .def("reset", &NucleonCollection::Reset)
        .def("Nucleons", &NucleonCollectionHelpers::GetNucleonList)
        .def("SetNucleons", &NucleonCollectionHelpers::SetNucleonList)
        .add_property("likelihood", &NucleonCollection::Likelihood)
    ;
/********************************************************/
}
