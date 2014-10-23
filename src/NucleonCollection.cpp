#include "NucleonCollection.h"
#include "Nucleon.h"

#include "math.h"

NucleonCollection::NucleonCollection(double pairwise_max, unsigned int units, double length)
  : likelihood(1), nucleon_count(0), units(units), length(length), cube_length(length/units) {
    pairwise_units = ceil(pairwise_max/cube_length);
    cubes = new nucleon_array** [2*units];
    for(unsigned int i = 0; i < 2*units; i++) {
        cubes[i] = new nucleon_array* [2*units];
        for(unsigned int j = 0; j < 2*units; j++) {
            cubes[i][j] = new nucleon_array [2*units];
        }
    }
}

NucleonCollection::NucleonCollection(const NucleonCollection &nucleon_collection) {
    nucleon_count = 0;
    units = nucleon_collection.units;
    length = nucleon_collection.length;
    cube_length = nucleon_collection.cube_length;

    for(unsigned int i = 0; i < nucleon_collection.NucleonCount(); i++) {
        AddNucleon(nucleon_collection[i]);
    }
}

NucleonCollection::~NucleonCollection() {
    //this deletes all of the nucleons
    Reset();

    for(unsigned int i = 0; i < 2*units; i++) {
        for(unsigned int j = 0; j < 2*units; j++) {
            delete [] cubes[i][j];
        }
        delete [] cubes[i];
    }
    delete cubes;
}

unsigned int NucleonCollection::AddNucleon(const Nucleon& nucleon) {
    Nucleon *new_nucleon = new Nucleon(nucleon);
    ordered.push_back(new_nucleon);

    new_nucleon->cube = &FindCube(nucleon);
    new_nucleon->cube->push_back(new_nucleon);
    new_nucleon->parent = this;
    return ++nucleon_count;
}

NucleonCollection::nucleon_array& NucleonCollection::FindCube(double x, double y, double z, int &i, int &j, int &k) {
    remquo(x, cube_length, &i);
    remquo(y, cube_length, &j);
    remquo(z, cube_length, &k);

    static const int units_minus_one = units - 1;
    static const unsigned int units_times_two = units*2;

    i = (i + units_minus_one)%units_times_two;
    j = (j + units_minus_one)%units_times_two;
    k = (k + units_minus_one)%units_times_two;

    return cubes[i][j][k];
}

NucleonCollection::nucleon_array& NucleonCollection::FindCube(double x, double y, double z) const {
    int i, j, k;
    return FindCube(x, y, z, i, j, k);
}

NucleonCollection::nucleon_array& NucleonCollection::FindCube(const Nucleon &nucleon) const {
    return FindCube(nucleon.X(), nucleon.Y(), nucleon.Z());
}

void NucleonCollection::Reset(bool delete_nucleons) {
    for(nucleon_array::iterator it = ordered.begin(); it != ordered.end(); it++) {
        (*it)->cube->clear();
        if(delete_nucleons) {
            delete *it;
        }
    }
    ordered.clear();
    nucleon_count = 0;
}

void NucleonCollection::SetNucleonPosition(Nucleon *nucleon, double x, double y, double z) {
    nucleon_array new_cube = FindCube(x, y, z);

    if(&new_cube != nucleon->cube) {
        RemoveNucleonFromCube(nucleon);
        new_cube.push_back(nucleon);
        nucleon->cube = &new_cube;

        //keep it instide the region
        while(x > length) x -= length;
        while(x < -length) x += length;
        while(y > length) y -= length;
        while(y < -length) y += length;
        while(z > length) z -= length;
        while(z < -length) z += length;
    }

    nucleon->x = x;
    nucleon->y = y;
    nucleon->z = z;
}

void NucleonCollection::RemoveNucleonFromCube(Nucleon *nucleon) {
    for(nucleon_array::iterator it = nucleon->cube->begin(); it != nucleon->cube->end(); it++) {
        if(*it == nucleon) {
            nucleon->cube->erase(it);
            break;
        }
    }
}

void NucleonCollection::RemoveNucleonFromOrdered(Nucleon *nucleon) {
    for(nucleon_array::iterator it = ordered.begin(); it != ordered.end(); it++) {
        if(*it == nucleon) {
            ordered.erase(it);
            break;
        }
    }
}

void NucleonCollection::RemoveNucleon(Nucleon *nucleon) {
    RemoveNucleonFromCube(nucleon);
    RemoveNucleonFromOrdered(nucleon);
}
