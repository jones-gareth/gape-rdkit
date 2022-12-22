//
// Created by gareth on 12/1/22.
//

#ifndef GAPE_MOLUTIL_H
#define GAPE_MOLUTIL_H

#include <GraphMol/GraphMol.h>

using namespace RDKit;

namespace Gape {

    void findAtoms(const Atom *startAtom, const std::vector<const Atom *> &stopAtoms,
                   std::set<const Atom *> &foundAtoms);

    void
    findAtoms(const Atom *startAtom, const Atom *stopAtom, std::set<const Atom *> &foundAtoms);


} // namespace Gape
#endif //GAPE_MOLUTIL_H
