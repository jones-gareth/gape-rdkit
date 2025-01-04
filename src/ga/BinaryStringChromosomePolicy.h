//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef BINARYSTRINGCHROMOSOMEPOLICY_H
#define	BINARYSTRINGCHROMOSOMEPOLICY_H

#include "util/export.h"
#include "util/RandomUtil.h"

namespace Gape {
    class GA_EXPORT BinaryStringChromosomePolicy {
    public:
        explicit BinaryStringChromosomePolicy(RandomUtil& rng_);

        BinaryStringChromosomePolicy(const BinaryStringChromosomePolicy& orig) = delete;

        BinaryStringChromosomePolicy& operator=(const BinaryStringChromosomePolicy& other) = delete;

        ~BinaryStringChromosomePolicy() = default;

        bool mutate(int pos, bool currentValue) const;

        bool initialize(int pos) const;

        bool isAllowSwitch() const { return allowSwitch; }

        void setAllowSwitch(const bool allow) { allowSwitch = allow; }

    private:
        bool allowSwitch = false;
        RandomUtil& rng;
    };
}


#endif	/* BINARYSTRINGCHROMOSOMEPOLICY_H */
