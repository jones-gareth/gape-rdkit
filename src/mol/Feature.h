//
// Created by Gareth Jones on 1/19/23.
//

#ifndef GAPE_FEATURE_H
#define GAPE_FEATURE_H

namespace Gape {

    enum FeatureType {
        HydrophobicAtom, DonorInteractionPoint, AcceptorAtom, AromaticRing
        // TODO UserFeatures
    };

    class Feature {
    protected:
        const FeatureType featureType;
    };

} // Gape

#endif //GAPE_FEATURE_H
