//
// Created by gareth on 10/18/22.
//

#ifndef GAPE_SUPERPOSITIONMOLECULE_H
#define GAPE_SUPERPOSITIONMOLECULE_H

#include <GraphMol/GraphMol.h>

using namespace RDKit;
namespace ForceFields {
    class ForceField;
}

namespace Gape {

    class SuperpositionMolecule {

    public:
        explicit SuperpositionMolecule(const ROMol &mol);

        virtual ~SuperpositionMolecule();

        SuperpositionMolecule(const SuperpositionMolecule &) = delete;

        SuperpositionMolecule &operator=(SuperpositionMolecule &) = delete;
        
        std::string ToMolBlock() const;

    private:
        RWMol mol;
        ForceFields::ForceField *forceField;

    };

} // namespace GAPE

#endif //GAPE_SUPERPOSITIONMOLECULE_H
