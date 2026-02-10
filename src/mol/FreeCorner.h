//
// Created by jones on 1/15/2026.
//

#pragma once

#include <GraphMol/GraphMol.h>

#include "gape/RotatableBond.h"

namespace Gape {
    using namespace RDKit;

    class CornerRotation : RotatableBond {
        const Atom *rootAtom;
        const Atom *checkAtom;

    public:
        CornerRotation(const Bond *bond, const SuperpositionMolecule *superpositionMolecule,
                       const Atom *rootAtom, const Atom *checkAtom);
    };

    class FreeCorner {
        const Atom *atomX;
        const Atom *atomA;
        const Atom *atomB;
        const Atom *atomC;
        const Atom *atomD;
        const Bond *bondAB;
        const Bond *bondBX;
        const Bond *bondXC;
        const Bond *bondCD;
        const std::vector<int> ringIndices;
        const std::unique_ptr<CornerRotation> rotationCD;
        const std::unique_ptr<CornerRotation> rotationXC;
        const std::unique_ptr<CornerRotation> rotationAB;


        RDGeom::Point3D findOtherCorner(const Conformer &conformer) const;

    public:
        FreeCorner(const Atom *atom_x, const Atom *atom_a, const Atom *atom_b, const Atom *atom_c, const Atom *atom_d,
                   const Bond *bond_ab, const Bond *bond_bx, const Bond *bond_xc, const Bond *bond_cd,
                   const std::vector<int> &ring_indices, std::unique_ptr<CornerRotation> &rotation_cd,
                   std::unique_ptr<CornerRotation> &rotation_xc, std::unique_ptr<CornerRotation> &rotation_ab
                   )
            : atomX(atom_x),
              atomA(atom_a),
              atomB(atom_b),
              atomC(atom_c),
              atomD(atom_d),
              bondAB(bond_ab),
              bondBX(bond_bx),
              bondXC(bond_xc),
              bondCD(bond_cd),
              ringIndices(ring_indices) ,
              rotationCD(std::move(rotation_cd)),
              rotationXC(std::move(rotation_xc)),
              rotationAB(std::move(rotation_ab)) {
        }

        const static double CORNER_FLAP_TOL;

        static std::shared_ptr<FreeCorner>
        isFreeCorner(const SuperpositionMolecule *superpositionMolecule, Atom *atomX);
    };
} // Gape
