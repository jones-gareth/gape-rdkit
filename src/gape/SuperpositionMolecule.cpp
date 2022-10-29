//
// Created by gareth on 10/18/22.
//

// #include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>

#include "SuperpositionMolecule.h"
#include <GraphMol/MolOps.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <ForceField/MMFF/Contribs.h>

using namespace RDKit;

namespace Gape {

    SuperpositionMolecule::SuperpositionMolecule(const ROMol &inputMol) {
        mol = inputMol;
        // TODO solvate
        MolOps::addHs(mol);
        DGeomHelpers::EmbedParameters embedParameters;
        auto confId = EmbedMolecule(mol, embedParameters);
        MMFF::MMFFMolProperties mmffMolProperties(mol);
        assert(mmffMolProperties.isValid());
        forceField = MMFF::constructForceField(mol);
        ForceFieldsHelper::OptimizeMolecule(*forceField);

        for (auto contrib: forceField->contribs()) {
           const ForceFields::ForceFieldContrib *ptr = contrib.get();
           const auto *torsion = dynamic_cast<const MMFF::TorsionAngleContrib *>(ptr);
           if (torsion != nullptr) {
               // ForceFields::MMFF::Utils::calcTorsionEnergy
               std::cerr << "Hello" << std::endl;
           }
        }
    }

    SuperpositionMolecule::~SuperpositionMolecule() {
        delete forceField;
    }

    std::string SuperpositionMolecule::ToMolBlock() const {
        return MolToMolBlock(mol);
    }

} // namespace GAPE
