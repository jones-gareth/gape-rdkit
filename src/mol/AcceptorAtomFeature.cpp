//
// Created by Gareth Jones on 3/23/2023.
//

#include "AcceptorAtomFeature.h"
#include "../util/Reporter.h"

namespace Gape {

    thread_local double AcceptorAtomFeature::hBondLen = 2.9;
    thread_local double AcceptorAtomFeature::chargeFactor = 2.0;
    thread_local double AcceptorAtomFeature::matchFactor = 1.0;
    thread_local bool AcceptorAtomFeature::scaleLonePairs = true;

    std::vector<std::shared_ptr<Feature>>
    AcceptorAtomFeature::findAcceptorAtoms(const SuperpositionMolecule *superpositionMolecule) {
        auto& mol = superpositionMolecule->getMol();
        auto& acceptors = superpositionMolecule->getAcceptors();
        std::vector<std::shared_ptr<Feature>> features;

        int featureSetNumber = 0;
        for (auto [atom, t]: acceptors) {
            auto acceptorAtom = mol.getAtomWithIdx(atom->getIdx());
            auto feature = std::make_shared<AcceptorAtomFeature>(featureSetNumber++, superpositionMolecule, acceptorAtom);
            features.push_back(std::static_pointer_cast<Feature>(feature));
        }
       return features;
    }

    AcceptorAtomFeature::AcceptorAtomFeature(const int featureSetNum, const SuperpositionMolecule *spMol,
                                             const Atom *featureAtom): AcceptorAtomFeature(featureSetNum) {
        molecule = spMol;
        atom = featureAtom;
        const auto& acceptors = molecule->getAcceptors();
        const auto acceptorIt = acceptors.find(atom);
        assert(acceptorIt != acceptors.end());
        hydrogenBondingType = acceptorIt->second.get();
        charged = false; //TODO add charge to hydrogen bonding types
        REPORT(Reporter::DEBUG) << "N Lone Pairs " << nLonePairs;
        acceptorAtom = std::make_unique<AcceptorAtom>(molecule, atom);
    }


    std::string AcceptorAtomFeature::pharmLabel() const
    {
        return  "" ;
    }

    std::string AcceptorAtomFeature::info() const
    {
	    return "";
    }




 } // Gape