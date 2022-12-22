//
// Created by jones on 11/25/2022.
//

#ifndef GAPE_ROTATABLEBOND_H
#define GAPE_ROTATABLEBOND_H

#include <GraphMol/GraphMol.h>
#include "SuperpositionMolecule.h"
#include <ForceField/MMFF/Contribs.h>

namespace Gape {
    using namespace RDKit;

    class SuperpositionMolecule;

    class TorsionInfo {
    public:
        const unsigned int index0, index1, index2, index3;
        const ForceFields::MMFF::MMFFTor mmffTorsion;

        TorsionInfo(unsigned int idx0, unsigned int idx1, unsigned int idx2, unsigned int idx3,
                    ForceFields::MMFF::MMFFTor mmffTor) : index0(idx0), index1(idx1), index2(idx2), index3(idx3),
                                                          mmffTorsion(mmffTor) {}
    };

    class RotatableBond {
    private:
        const Atom *atom1, *atom2;
        SuperpositionMolecule *molecule;
        std::vector<const Atom *> atom1List, atom2List;
        RotatableBondType rotatableBondType;
        std::vector<TorsionInfo> torsions;

    public:
        RotatableBond(RotatableBondType rotatableBondType, const Bond * bond,
                      SuperpositionMolecule *superpositionMolecule);

        bool isSeparatedByBond(const Atom *a1, const Atom *a2) const;

        RotatableBondType getRotatableBondType() const { return rotatableBondType; }

        ~RotatableBond();
    };


}

#endif //GAPE_ROTATABLEBOND_H
