//
// Created by jones on 1/15/2026.
//

#include "FreeCorner.h"

#include "MolUtil.h"
#include "util/Reporter.h"

namespace Gape {
    const double FreeCorner::CORNER_FLAP_TOL = 0.15;

    std::shared_ptr<FreeCorner> FreeCorner::isFreeCorner(const SuperpositionMolecule *superpositionMolecule,
                                                         Atom *atomX) {
        const auto &mol = atomX->getOwningMol();
        const auto ringInfo = mol.getRingInfo();
        const auto atomRings = ringInfo->atomMembers(atomX->getIdx());
        if (atomRings.size() != 1) {
            return nullptr;
        }

        const auto ringIndices = ringInfo->atomRings()[atomRings[0]];
        if (ringIndices.size() < 4) {
            return nullptr;
        }

        const Atom *atomC = nullptr;
        const Atom *atomB = nullptr;
        const Bond *bondBX = nullptr;
        const Bond *bondXC = nullptr;

        for (int idx: ringIndices) {
            auto bond = mol.getBondBetweenAtoms(idx, atomX->getIdx());
            if (bond == nullptr) {
                continue;
            }
            if (bond->getBondType() != Bond::BondType::SINGLE) {
                return nullptr;
            }
            if (atomB == nullptr) {
                atomB = mol.getAtomWithIdx(idx);
                bondBX = bond;
            } else if (atomC == nullptr) {
                atomC = mol.getAtomWithIdx(idx);
                bondXC = bond;
            }
        }

        assert(atomB != nullptr);
        assert(atomC != nullptr);

        const Atom *atomA = nullptr;
        const Atom *atomD = nullptr;
        const Bond *bondAB = nullptr;
        const Bond *bondCD = nullptr;
        for (int idx: ringIndices) {
            if (idx == atomB->getIdx() || idx == atomC->getIdx() || idx == atomX->getIdx()) {
                continue;
            }
            auto bondToB = mol.getBondBetweenAtoms(idx, atomB->getIdx());
            if (bondToB != nullptr) {
                if (bondToB->getBondType() != Bond::BondType::SINGLE) {
                    return nullptr;
                }
                atomA = mol.getAtomWithIdx(idx);
                bondAB = bondToB;
                continue;
            }

            auto bondToC = mol.getBondBetweenAtoms(idx, atomC->getIdx());
            if (bondToC != nullptr) {
                if (bondToC->getBondType() != Bond::BondType::SINGLE) {
                    return nullptr;
                }
                atomD = mol.getAtomWithIdx(idx);
                bondCD = bondToC;
            }
        }

        assert(atomA != nullptr);
        assert(atomD != nullptr);

        auto rotationCD = std::make_unique<CornerRotation>(bondCD, superpositionMolecule, atomC, atomB);
        auto rotationXC = std::make_unique<CornerRotation>(bondXC, superpositionMolecule, atomX, atomB);
        auto rotationAB = std::make_unique<CornerRotation>(bondAB, superpositionMolecule, atomB, atomX);

        return std::make_shared<FreeCorner>(atomX, atomA, atomB, atomC, atomD, bondAB, bondBX, bondXC, bondCD,
                                            ringIndices, rotationCD, rotationXC, rotationAB);
    }

    CornerRotation::CornerRotation(const Bond *bond, const SuperpositionMolecule *superpositionMolecule,
                                   const Atom *rootAtom, const Atom *checkAtom) : RotatableBond(RotatableBondType::Full,
            bond, superpositionMolecule),
        rootAtom(rootAtom), checkAtom(checkAtom) {
        if (bond->getBeginAtom() == rootAtom) {
            atom1 = bond->getBeginAtom();
            atom2 = bond->getEndAtom();
        } else {
            atom1 = bond->getEndAtom();
            atom2 = bond->getBeginAtom();
        }

        std::set<const Atom *> atoms1, atoms2;
        const std::vector<const Atom *> stop1Atoms{atom2, checkAtom};
        const std::vector<const Atom *> stop2Atoms{atom1, checkAtom};
        findAtoms(atom1, stop1Atoms, atoms1);
        findAtoms(atom2, stop2Atoms, atoms1);
        atom1List.assign(atoms1.begin(), atoms1.end());
        atom2List.assign(atoms2.begin(), atoms2.end());
    }

    /**
     * quadraticRoots[] are solutions of ax^2 + bx + c. A check is made for
     * roots close to zero.
     *
     * @param a
     * @param b
     * @param c
     * @throws GaException
     */
    void solveQuadratic(const double a, const double b, const double c, double &root1, double &root2) {
        double sqrTerm = b * b - 4.0 * a * c;

        REPORT(Reporter::DEBUG) << "Solve Quad a " << a << " b " << b << " c " << c;
        REPORT(Reporter::DEBUG) << "Solve Qaaud b**2-4ac " << sqrTerm;

        if (sqrTerm < .0) {
            double check = b * b / (4 * a * c);
            if (check > 0.9999999 && check < 1.00000001)
                sqrTerm = .0;
            else
                throw std::runtime_error("solveQuadratic: complex root\n");
        }

        double x = sqrt(sqrTerm);
        root1 = (-b + x) / (2.0 * a);
        root2 = (-b - x) / (2.0 * a);
    }

    /**
     * Returns true if two distances are nearly the same. The parameter
     * CORNER_FLAP_TOL is used to check for this. This routine is used to
     * determine which of the two alternate points for the free corner
     * corresponds to the current corner.
     *
     * @param p1
     * @param p2
     * @return
     */
    bool nearlyEqual(double p1, double p2) {
        double diff = p1 - p2;
        if (diff > FreeCorner::CORNER_FLAP_TOL || -diff > FreeCorner::CORNER_FLAP_TOL)
            return false;
        return true;
    }

    RDGeom::Point3D FreeCorner::findOtherCorner(const Conformer &conformer) const {
        auto pointA = conformer.getAtomPos(atomA->getIdx());
        auto pointB = conformer.getAtomPos(atomB->getIdx());
        const auto &pointX = conformer.getAtomPos(atomX->getIdx());
        const auto &pointC = conformer.getAtomPos(atomC->getIdx());
        const auto &pointD = conformer.getAtomPos(atomD->getIdx());
        const auto midPoint = (pointC + pointD) / 2.0;

        int no = 0;
        bool fitted = false;

        double r1sqr = 0, r2sqr = 0, r3sqr = 0, x12 = 0, x13 = 0, y12 = 0, y13 = 0, y23 = 0, z12 = 0, z13 = 0, l1Sqr = 0
                , l2Sqr = 0, l3Sqr = 0, p12Sqr = 0, p13Sqr = 0, epsilon = 0;

        while (!fitted) {
            r1sqr = squareDistance(pointX, pointA);
            r2sqr = squareDistance(pointX, pointB);
            r3sqr = squareDistance(pointX, midPoint);

            REPORT(Reporter::DEBUG) << "r1sqr " << r1sqr;
            REPORT(Reporter::DEBUG) << "r2sqr " << r2sqr;
            REPORT(Reporter::DEBUG) << "r3sqr " << r3sqr;

            x12 = pointA[0] - pointB[0];
            x13 = pointA[0] - midPoint[0];
            y12 = pointA[1] - pointB[1];
            y13 = pointA[1] - midPoint[1];
            y23 = pointB[1] - midPoint[1];
            z12 = pointA[2] - pointB[2];
            z13 = pointA[2] - midPoint[2];

            REPORT(Reporter::DEBUG) << "x12 " << x12;
            REPORT(Reporter::DEBUG) << "x13 " << x13;
            REPORT(Reporter::DEBUG) << "y12 " << y12;
            REPORT(Reporter::DEBUG) << "y13 " << y13;
            REPORT(Reporter::DEBUG) << "y23 " << y23;
            REPORT(Reporter::DEBUG) << "z12 " << z12;
            REPORT(Reporter::DEBUG) << "z13 " << z13;

            l1Sqr = pointA[0] * pointA[0] + pointA[1] * pointA[1] + pointA[2] * pointA[2] - r1sqr;
            l2Sqr = pointB[0] * pointB[0] + pointB[1] * pointB[1] + pointB[2] * pointB[2] - r2sqr;
            l3Sqr = midPoint[0] * midPoint[0] + midPoint[1] * midPoint[1] + midPoint[2] * midPoint[2] - r3sqr;

            REPORT(Reporter::DEBUG) << "l1Sqr " << l1Sqr;
            REPORT(Reporter::DEBUG) << "l2Sqr " << l2Sqr;
            REPORT(Reporter::DEBUG) << "l3Sqr " << l3Sqr;

            p12Sqr = 0.5 * (l1Sqr - l2Sqr);
            p13Sqr = 0.5 * (l1Sqr - l3Sqr);

            REPORT(Reporter::DEBUG) << "p12Sqr " << p12Sqr;
            REPORT(Reporter::DEBUG) << "p13Sqr " << p13Sqr;
            REPORT(Reporter::DEBUG) << "pointA[0] " << pointA[0];
            REPORT(Reporter::DEBUG) << "pointB[0] " << pointB[0];
            REPORT(Reporter::DEBUG) << "midPoint[0] " << midPoint[0];
            double denom = -y23 * pointA[0] + y13 * pointB[0] - y12 * midPoint[0];
            epsilon = 1.0 / denom;
            REPORT(Reporter::DEBUG) << "denom " << denom;
            REPORT(Reporter::DEBUG) << "epsilon " << epsilon;

            if (abs(epsilon) < 1.0e-5) {
                pointA = conformer.getAtomPos(atomB->getIdx());
                pointB = conformer.getAtomPos(atomA->getIdx());
            } else {
                fitted = true;
            }
            no++;
            if (no > 2)
                throw std::runtime_error("Epsilon Test failed");
        }

        double alpha = epsilon * (y12 * p13Sqr - y13 * p12Sqr) - pointA[0];
        double beta = epsilon * (x13 * p12Sqr - x12 * p13Sqr) - pointA[1];
        double gamma = epsilon * (y13 * z12 - y12 * z13);
        double delta = epsilon * (x12 * z13 - x13 * z12);

        double c1 = 1.0 + gamma * gamma + delta * delta;
        double c2 = -2.0 * (pointA[2] - alpha * gamma - beta * delta);
        double c3 = pointA[2] * pointA[2] - r1sqr + alpha * alpha + beta * beta;

        double root1, root2;
        solveQuadratic(c1, c2, c3, root1, root2);

        RDGeom::Point3D altPoint;
        if (nearlyEqual(root1, pointX[2]))
            altPoint[2] = root2;
        else if (nearlyEqual(root2, pointX[2]))
            altPoint[2] = root1;
        else
            throw std::runtime_error("find other corner: didn't find original point");

        altPoint[0] = alpha + pointA[0] + gamma * altPoint[2];
        altPoint[1] = beta + pointA[1] + delta * altPoint[2];

        REPORT(Reporter::DEBUG) << "Other corner " << altPoint;
        return altPoint;
    }

    /**
     * Flips the free corner. See Payne and Glen for all the gory details
     */
    void FreeCorner::flipFreeCorner(SuperpositionCoordinates &superpositionCoordinates) const {
        REPORT(Reporter::DEBUG) << "Flipping corner";
        REPORT(Reporter::DEBUG) << "Atom A " << atomA->getIdx();
        REPORT(Reporter::DEBUG) << "Atom B " << atomB->getIdx();
        REPORT(Reporter::DEBUG) << "Atom X " << atomX->getIdx();
        REPORT(Reporter::DEBUG) << "Atom C " << atomC->getIdx();
        REPORT(Reporter::DEBUG) << "Atom D " << atomD->getIdx();

        auto conformer = superpositionCoordinates.getConformer();
        auto altPoint = findOtherCorner(conformer);

        for (int i = 0; i < 4; i++) {
            const Atom *otherAtom = nullptr;
            if (i == 0)
                otherAtom = atomA;
            else if (i == 1)
                otherAtom = atomB;
            else if (i == 2)
                otherAtom = atomC;
            else if (i == 3)
                otherAtom = atomD;

            const auto d1 = squareDistance(altPoint, conformer.getAtomPos(otherAtom->getIdx()));
            const auto d2 = squareDistance(conformer.getAtomPos(atomX->getIdx()),
                                           conformer.getAtomPos(otherAtom->getIdx()));
            const double diff = d1 - d2;
            REPORT(Reporter::DEBUG) << "Point distance check " << diff;
            if (abs(diff) > CORNER_FLAP_TOL) {
                REPORT(Reporter::DEBUG) << "corner not flipped (" << i << ")";
                return;
            }
        }

        auto yVec = conformer.getAtomPos(atomX->getIdx());
        auto zVec = conformer.getAtomPos(atomB->getIdx());

        double rot1 = torsionAngle(altPoint, conformer.getAtomPos(atomA->getIdx()),
                                   conformer.getAtomPos(atomB->getIdx()), yVec);
        REPORT(Reporter::DEBUG) << "rotationAB " << rot1;
        RDGeom::Transform3D rotation1Transform;
        rotationAB->rotateBond(rot1, superpositionCoordinates, rotation1Transform);
        rotation1Transform.TransformPoint(yVec);

        if (Reporter::getMinReportingLevel() <= Reporter::DEBUG) {
            const double sqrDistance = squareDistance(altPoint, yVec);
            if (sqrDistance > CORNER_FLAP_TOL * CORNER_FLAP_TOL)
                REPORT(Reporter::DEBUG) << "flipCorner: (1) point mismatch, sqrDistance " << sqrDistance;
            else
                REPORT(Reporter::DEBUG) << "flipCorner: (1) sqrDistance check " << sqrDistance;
        }

        double rot2 = torsionAngle(altPoint, conformer.getAtomPos(atomC->getIdx()),
                                   conformer.getAtomPos(atomD->getIdx()), conformer.getAtomPos(atomX->getIdx()));
        REPORT(Reporter::DEBUG) << "rotationCD " << rot2;
        RDGeom::Transform3D rotation2Transform;
        rotationCD->rotateBond(rot2, superpositionCoordinates, rotation2Transform);
        rotation2Transform.TransformPoint(zVec);

        if (Reporter::getMinReportingLevel() <= Reporter::DEBUG) {
            double sqrDistance = squareDistance(altPoint, conformer.getAtomPos(atomX->getIdx()));
            if (sqrDistance > CORNER_FLAP_TOL * CORNER_FLAP_TOL)
                REPORT(Reporter::DEBUG) << "flipCorner: (2) point mismatch, sqrDistance " << sqrDistance;
            else
                REPORT(Reporter::DEBUG) << "flipCorner: (2) sqrDistance check " << sqrDistance;
        }

        double rot3 = torsionAngle(zVec, conformer.getAtomPos(atomX->getIdx()),
                                   conformer.getAtomPos(atomC->getIdx()), conformer.getAtomPos(atomB->getIdx()));
        REPORT(Reporter::DEBUG) << "rotationXC " << rot3;
        RDGeom::Transform3D rotation3Transform;
        rotationXC->rotateBond(-rot3, superpositionCoordinates, rotation3Transform);
        rotation3Transform.TransformPoint(zVec);

        if (Reporter::getMinReportingLevel() <= Reporter::DEBUG) {
            double sqrDistance = squareDistance(conformer.getAtomPos(atomB->getIdx()), zVec);
            if (sqrDistance > CORNER_FLAP_TOL * CORNER_FLAP_TOL)
                REPORT(Reporter::DEBUG) << "flipCorner: (3) point mismatch, sqrDistance " << sqrDistance;
            else
                REPORT(Reporter::DEBUG) << "flipCorner: (3) sqrDistance check " << sqrDistance;
            // TODO - add atomZ and atomY to molecule to visualize corner flip
            // conformer.setAtomPos(atomZ->getIdx(), zVec);
            // conformer.setAtomPos(atomY->getIdx(), yVec);
        }
    }

    void FreeCorner::addTorsions(const SuperpositionMolecule &superpositionMolecule,
                                 std::vector<TorsionInfo> &torsions) const {
        const auto mmffMolProperties = superpositionMolecule.getMMFFMolProperties();
        for (int j = 0; j < 4; j++) {
            const Bond *bond = nullptr;
            if (j == 0)
                bond = bondAB;
            else if (j == 1)
                bond = bondBX;
            else if (j == 2)
                bond = bondXC;
            else if (j == 3)
                bond = bondCD;

            auto atom1 = bond->getBeginAtom();
            auto atom2 = bond->getEndAtom();
            const auto &mol = superpositionMolecule.getMol();

            for (const auto atom0: mol.atomNeighbors(atom1)) {
                if (atom0 != atom2) {
                    for (const auto atom3: mol.atomNeighbors(atom2)) {
                        if (atom3 != atom1) {
                            ForceFields::MMFF::MMFFTor mmffTorsion;
                            const auto index0 = atom0->getIdx();
                            const auto index1 = atom1->getIdx();
                            const auto index2 = atom2->getIdx();
                            const auto index3 = atom3->getIdx();
                            unsigned int torsionType;
                            if (mmffMolProperties->getMMFFTorsionParams(mol, index0, index1, index2, index3,
                                                                        torsionType,
                                                                        mmffTorsion)) {
                                TorsionInfo torsionInfo(index0, index1, index2, index3, mmffTorsion);

                                if (auto it = std::find(torsions.begin(), torsions.end(), torsionInfo);
                                    it == torsions.end()) {
                                    torsions.push_back(torsionInfo);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
} // Gape
