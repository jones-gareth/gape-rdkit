//
// Created by gareth on 10/18/22.
//

// #include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>

#include "SuperpositionMolecule.h"
#include <GraphMol/MolOps.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/FileParsers/FileParsers.h>

using namespace RDKit;

namespace Gape {

    SuperpositionMolecule::SuperpositionMolecule(const ROMol &inputMol) {
        mol = inputMol;
        // TODO solvate
        MolOps::addHs(mol);
        DGeomHelpers::EmbedParameters embedParameters;
        auto confId = EmbedMolecule(mol, embedParameters);
        MMFF::MMFFMolProperties mmffMolProperties(mol);
        if (mmffMolProperties.isValid()) {
            forceField = MMFF::constructForceField(mol);
            ForceFieldsHelper::OptimizeMolecule(*forceField);
        }
    }

    SuperpositionMolecule::~SuperpositionMolecule() {
        delete forceField;
    }

    std::string SuperpositionMolecule::ToMolBlock() const{
        return MolToMolBlock(mol);
    }

} // namespace GAPE
