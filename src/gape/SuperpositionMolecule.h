//
// Created by gareth on 10/18/22.
//

#ifndef GAPE_SUPERPOSITIONMOLECULE_H
#define GAPE_SUPERPOSITIONMOLECULE_H

#include <GraphMol/GraphMol.h>
#include "Gape.h"

using namespace RDKit;
namespace ForceFields {
    class ForceField;
}

namespace Gape {

    enum RotatableBondType {
        None, Flip, Full
    };
    class SuperpositionMolecule {

    public:
        explicit SuperpositionMolecule(const ROMol &mol, const Gape &settings);

        virtual ~SuperpositionMolecule();

        SuperpositionMolecule(const SuperpositionMolecule &) = delete;

        SuperpositionMolecule &operator=(SuperpositionMolecule &) = delete;
        
        std::string ToMolBlock() const;

    private:
        RWMol mol;
        ForceFields::ForceField *forceField;
        const Gape &settings;

        void findFreelyRotatableBonds();

        bool isO2(const Atom &atom) const;
        bool isO3(const Atom &atom) const;
        bool isAmideBond(const Bond &bond) const;
        bool isNpl3Atom(const Atom &atom) const;
        static bool isTerminalBond(const Bond &bond);
        bool isArginineCarbon(const Atom &atom) const;
        static bool isSp2Carbon(const Atom &atom) ;
        bool atomIsInRing(const Atom &atom) const;
        bool isCOOHCarbon(const Atom &atom, Atom * &o2Atom, Atom *&o3Atom) const;
        RotatableBondType isRotatableBond(const Bond &bond, bool &canFlatten) const;

    };

} // namespace GAPE

#endif //GAPE_SUPERPOSITIONMOLECULE_H
