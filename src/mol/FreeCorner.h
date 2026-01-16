//
// Created by jones on 1/15/2026.
//

#pragma once

#include <GraphMol/GraphMol.h>

namespace Gape {
    using namespace RDKit;

    class FreeCorner {

        Atom* atomX;
        Atom* atomA;
        Atom* atomB;
        Atom* atomC;
        Atom* atomD;
        Bond* bondAB;
        Bond* bondBX;
        Bond* bondXC;
        Bond* bondCD;

        std::vector<int> ringIndices;

        const static double CORNER_FLAP_TOL;
    public:
        static bool checkFreeCorder(const Atom *atomX);
    };
} // Gape
