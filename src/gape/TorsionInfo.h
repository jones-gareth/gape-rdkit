#pragma once

namespace Gape {
    using namespace RDKit;

    class TorsionInfo
    {
    public:
        const unsigned int index0, index1, index2, index3;
        const ForceFields::MMFF::MMFFTor mmffTorsion;
        double referenceAngle = .0;

        TorsionInfo(unsigned int idx0, unsigned int idx1, unsigned int idx2, unsigned int idx3,
                    ForceFields::MMFF::MMFFTor mmffTor) : index0(idx0), index1(idx1), index2(idx2), index3(idx3),
                                                          mmffTorsion(mmffTor)
        {
        }

        double torsionEnergy(const Conformer& conformer) const;

        bool operator==(const TorsionInfo &other) const;
    };
}