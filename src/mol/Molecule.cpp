/*
 * Molecule.cpp
 *
 *  Created on: Apr 28, 2014
 *      Author: gjones
 */

#include "Molecule.h"

#include <boost/format.hpp>
#include <cassert>
#include <fstream>
#include <functional>
#include <numeric>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <memory>

#include "../util/GzipReader.h"
#include "../util/Reporter.h"
#include "../util/FileSystemUtil.h"
#include "AtomType.h"
#include "BondType.h"
#include "MoleculeSssr.h"
#include "Ring.h"
#include "Huckel.h"

namespace GarethMol {

using namespace std;
using namespace GarethUtil;

using AtomTypeId = AtomType::AtomTypeId;

static bool readToMol2Section(istream & in, string section) {
    string line;
    while (getline(in, line)) {
        if (startsWith(line, section)) {
            return true;
        }
    }

    return false;
}

bool Molecule::loadMol2Molecule(istream & in) {
    string line;

    if (!readToMol2Section(in, "@<TRIPOS>MOLECULE")) {
        if (in.eof()) {
            return false;
        } else
            throw runtime_error(
                    "unable to find molecule section in Mol2 input!");
    }

    // first line is molecule name
    getline(in, line);
    name = removeTrailingLF(line);

    // second line no atoms and bonds
    int nAtoms = 0, nBonds = 0;
    getline(in, line);
    REPORT(Reporter::TRACE) << "got size line " << line;
    stringstream(line) >> nAtoms >> nBonds;

    // reserve for atoms plus hydrogen and LP

    atoms.reserve(nAtoms + 30);
    bonds.reserve(nBonds + 30);
    coords = CoordMatrix(4, nAtoms);

    REPORT(Reporter::DETAIL) << "Reading MOL2 molecule " << name << " n atoms "
            << nAtoms << " n bonds " << nBonds;

    // third and forth line molecule type and charge type - ignore for now
    //getline(in, line);
    //getline(in, line);

    if (!readToMol2Section(in, "@<TRIPOS>ATOM"))
        throw runtime_error("unable to find atom section in Mol2 input!");
    for (int no = 0; no < nAtoms; ++no) {
        int atomNo;
        string atomName, atomTypeName;
        CoordVector c;
        c(3) = 1.0;

        getline(in, line);
        stringstream(line) >> atomNo >> atomName >> c[0] >> c[1] >> c[2]
                >> atomTypeName;
        REPORT(Reporter::DEBUG) << "reading atom no " << atomNo << " type "
                << atomTypeName;
        atomNo--;
        assert(atomNo == no);
        auto & atomType = AtomType::typeFromName(atomTypeName);
        atoms.push_back(make_unique<Atom>(atomNo, atomType, atomName));
        this->coords.col(no) = c;
    }

    if (!readToMol2Section(in, "@<TRIPOS>BOND"))
        throw runtime_error("unable to find bond section in Mol2 input!");
    for (int no = 0; no < nBonds; ++no) {
        int bondNo, atom1No, atom2No;
        string bondTypeStr;

        getline(in, line);
        stringstream(line) >> bondNo >> atom1No >> atom2No >> bondTypeStr;
        REPORT(Reporter::DEBUG) << "reading bond no " << bondNo << " type "
                << bondTypeStr;
        --bondNo;
        --atom1No;
        --atom2No;
        assert(bondNo == no);
        auto & bondType = BondType::typeFromName(bondTypeStr);
        auto bond = make_unique<Bond>(bondNo, bondType, atoms.at(atom1No).get(),
                atoms.at(atom2No).get());
        bonds.push_back(move(bond));
    }

    hasCoordinates = true;
    return true;
}

void Molecule::writeMol2File(ostream & out, const bool & includeLonePairs,
        const string & comment) const {
    writeMol2File(out, getCoords(), includeLonePairs, comment);
}

void Molecule::writeMol2File(ostream & out, const CoordMatrix & otherCoords,
        const bool & includeLonePairs, const string & comment) const {
    int nAtoms = atoms.size();

    assert(otherCoords.cols() == nAtoms);

    int nBonds = bonds.size();
    REPORT(Reporter::NORMAL) << "Writing MOL2 molecule " << name << " n atoms "
            << nAtoms << " n bonds " << nBonds;

    out << "#       Name:                   " << name << endl;
    out << "#       Creating User Name:     " << getUserName() << endl;
    out << "#       Creation Time:          " << currentTime() << endl;
    out << "#       " << endl;
    out << "#       " << comment << endl;
    out << endl << endl;

    out << "@<TRIPOS>MOLECULE" << endl;
    out << name << endl;
    out << "   " << nAtoms << "    " << bonds.size() << "    0" << endl;
    out << "SMALL" << endl << "USER_CHARGES" << endl << endl << endl;

    map<int, int> atomNoMap = mapAtomsForOutput(includeLonePairs);

    out << "@<TRIPOS>ATOM" << endl;
    for (int no = 0; no < nAtoms; ++no) {
        auto & atom = atoms.at(no);
        const CoordVector & c = otherCoords.col(no);
        int atomNo = atomNoMap.at(atom->getAtomNo()) + 1;
        if (atomNo < 1)
            continue;
        out << "   " << atomNo << "   " << atom->getName() << "   " << c(0)
                << "   " << c(1) << "   " << c(2) << "   "
                << atom->getAtomType().getName() << "   1   MOLECULE";
        out << endl;
    }

    out << "@<TRIPOS>BOND" << endl;
    int no = 1;
    for (auto & bond : bonds) {
        int atom1No = atomNoMap.at(bond->getAtom1().getAtomNo()) + 1;
        int atom2No = atomNoMap.at(bond->getAtom2().getAtomNo()) + 1;
        if (atom1No < 1 || atom2No < 1)
            continue;
        const auto bondTypeId = bond->sybylBondType(*this);
        const auto & bondType = BondType::typeFromTypeId(bondTypeId);
        out << "   " << no << "    " << atom1No << "   " << atom2No << "   "
                << bondType.getName() << "   " << endl;
        ++no;
    }

}

bool Molecule::loadSdfMolecule(istream & in) {

    string line;
    getline(in, line);
    if (in.eof()) {
        return false;
    }
    removeTrailingLF(line);
    if (line.compare("$$$$") == 0) {
        getline(in, line);
        removeTrailingLF(line);
        if (in.eof()) {
            return false;
        }
    }
    name = line;
    getline(in, line);
    getline(in, line);
    getline(in, line);

    if (in.eof())
        return false;
    int nAtoms = 0, nBonds = 0;
    bool ok;
    nAtoms = convertString<int>(line, &ok, 0, 3);
    if (!ok) {
        return false;
    }
    nBonds = convertString<int>(line, &ok, 3, 6);
    if (!ok) {
        return false;
    }

    // reserve for atoms plus hydrogen and LP

    atoms.reserve(nAtoms + 30);
    bonds.reserve(nBonds + 30);
    coords = CoordMatrix(4, nAtoms);

    REPORT(Reporter::DEBUG) << "Reading SDF molecule " << name << " n atoms "
            << nAtoms << " n bonds " << nBonds;

    sdfChiral = false;

    if (line.length() > 15) {
        string chiralStr = line.substr(12, 3);
        stringstream ss(chiralStr);
        int chiral = 0;
        ss >> chiral;
        if (!ss.fail() && chiral == 1) {
            sdfChiral = true;
        }
    }

    for (int i = 0; i < nAtoms; i++) {
        getline(in, line);

        CoordVector c;
        c(3) = 1.0;

        c(0) = convertString<double>(line, &ok, 0, 10);
        c(1) = convertString<double>(line, &ok, 10, 10);
        c(2) = convertString<double>(line, &ok, 20, 10);
        string element = convertString<string>(line, &ok, 30, 4);

        auto & atomType = AtomType::typeFromName(element);
        auto atom = make_unique<Atom>(i, atomType, element);
        coords.col(i) = c;

        int massDifference = convertString<int>(line, &ok, 34, 3);
        atom->setSdfMassDifference(massDifference);
        int formalCharge = convertString<int>(line, &ok, 36, 3);
        if (formalCharge != 0) {
            atom->setFormalCharge(4 - formalCharge);
            atom->setFormalChargeSet(true);
        }
        atom->setFormalChargeSet(false);

        atom->setSdfStereo(Atom::SdfStereo::NONE);
        int stereo = convertString<int>(line, &ok, 39, 3);
        if (ok && stereo > 0) {
            if (stereo == 1)
                atom->setSdfStereo(Atom::SdfStereo::ODD);
            else if (stereo == 2)
                atom->setSdfStereo(Atom::SdfStereo::EVEN);
            else if (stereo == 3)
                atom->setSdfStereo(Atom::SdfStereo::EITHER);
        }

        int valence = convertString<int>(line, &ok, 48, 3);
        if (ok && valence != 0) {
            valence = 15;
            atom->setSdfValence(valence);
            atom->setSdfValenceSet(true);
        } else {
            atom->setSdfValence(0);
            atom->setSdfValenceSet(false);
        }

        atoms.push_back(move(atom));
    }
    hasCoordinates = true;

    for (int i = 0; i < nBonds; i++) {
        getline(in, line);
        bool ok;
        int atom1No = convertString<int>(line, &ok, 0, 3);
        int atom2No = convertString<int>(line, &ok, 3, 6);
        int type = convertString<int>(line, &ok, 6, 9);
        atom1No--;
        atom2No--;
        auto & bondType = BondType::sdfType(type);
        auto bond = make_unique<Bond>(i, bondType, atoms.at(atom1No).get(),
                atoms.at(atom2No).get());
        bonds.push_back(move(bond));
    }

    regex pat(R"(> *<(.*)>)");
    smatch matches;

    while (true) {
        getline(in, line);
        if (in.eof())
            return true;
        if (startsWith(line, "$$$$"))
            return true;
        trim(line);
        if (startsWith(line, ">")) {
            if (regex_search(line, matches, pat)) {
                string sdField = matches[1];
                vector<string> sdValues;
                while (!line.empty()) {
                    getline(in, line);
                    trim(line);
                    if (!line.empty()) {
                        sdValues.push_back(line);
                    }
                }
                if (!sdValues.empty()) {
                    REPORT(Reporter::TRACE) << "adding sd field " << sdField
                            << " values " << collectionToString(sdValues);
                    addSdfField(sdField, sdValues);
                }
            }
        }
    }

    return true;
}

void Molecule::writeSdfFile(ostream & out, const string & comment) const {
    writeSdfFile(out, getCoords(), comment);
}

void Molecule::writeSdfFile(ostream & out, const CoordMatrix & otherCoords,
        const string & comment) const {
    out << name << endl;
    out << " -GAPE---08220413213D" << endl << endl;

    // determine output atoms
    map<int, int> atomNoMap = mapAtomsForOutput();

    // count atoms and bonds to be output
    auto atomAcc = [](int s, map<int, int>::value_type pair) {
        if (pair.second == -1) return s; else return ++s;
    };
    int nAtoms = std::accumulate(atomNoMap.begin(), atomNoMap.end(), 0,
            atomAcc);
    auto bondAcc = [atomNoMap](int s, const unique_ptr<Bond> & bond) {
        int a1 = atomNoMap.at(bond->getAtom1().getAtomNo());
        int a2 = atomNoMap.at(bond->getAtom2().getAtomNo());
        if (a1 == -1 || a2 == -1)
        return s;
        return ++s;
    };
    int nBonds = std::accumulate(bonds.begin(), bonds.end(), 0, bondAcc);

    REPORT(Reporter::NORMAL) << "Writing SDF molecule " << name << " n atoms "
            << nAtoms << " n bonds " << nBonds;

    out << boost::format("%3d%3d") % nAtoms % nBonds;
    out << "  0  0";
    if (sdfChiral)
        out << "  1";
    else
        out << "  0";
    out << "  0  0  0  0  0999 V2000" << endl;

    // atom block
    for (size_t i = 0; i < atoms.size(); ++i) {
        auto & atom = atoms.at(i);
        int atomNo = atomNoMap.at(atom->getAtomNo());
        if (atomNo == -1)
            continue;

        const CoordVector & c = otherCoords.col(i);
        out
                << boost::format("%10.4f%10.4f%10.4f %-3s") % c(0) % c(1) % c(2)
                        % atom->getAtomType().sdType();
        if (atom->getSdfMassDifference() < 0) {
            out << boost::format("-%1d") % -atom->getSdfMassDifference();
        } else {
            out << boost::format(" %1d") % atom->getSdfMassDifference();
        }
        int charge = 0;
        if (atom->isFormalChargeSet()) {
            charge = 4 - atom->getFormalCharge();
        }
        int stereo = 0;
        if (atom->getSdfStereo() == Atom::SdfStereo::ODD)
            stereo = 1;
        else if (atom->getSdfStereo() == Atom::SdfStereo::EVEN)
            stereo = 2;
        else if (atom->getSdfStereo() == Atom::SdfStereo::EITHER)
            stereo = 3;
        int valence = atom->getSdfValence();
        if (valence == 0 && atom->isSdfValenceSet())
            valence = 15;
        assert(!atom->isSdfValenceSet() && valence != 15);
        out
                << boost::format("%3d%3d  0  0%3d  0  0  0  0  0  0") % charge
                        % stereo % valence << endl;
    }

    // bond block
    for (auto & bond : bonds) {
        int atom1no = atomNoMap.at(bond->getAtom1().getAtomNo());
        int atom2no = atomNoMap.at(bond->getAtom2().getAtomNo());
        if (atom1no == -1 || atom2no == -1)
            continue;
        ++atom1no;
        ++atom2no;
        int stereo = 0;
        if (bond->getStereo() == Bond::SdfStereo::UP)
            stereo = 1;
        if (bond->getStereo() == Bond::SdfStereo::DOWN)
            stereo = 4;
        if (bond->getStereo() == Bond::SdfStereo::EITHER)
            stereo = 6;
        if (bond->getStereo() == Bond::SdfStereo::CIS_TRANS)
            stereo = 3;
        out
                << boost::format("%3d%3d%3d%3d  0  0  0") % atom1no % atom2no
                        % bond->getBondType().sdType() % stereo << endl;
    }

    out << "M  END" << endl;

    if (!comment.empty()
            && !isKeyPresent<string, vector<string>>("comment", sdfFields)) {
        out << ">  <COMMENT>" << endl << comment << endl << endl;
    }
    for (auto iter = sdfFields.begin(); iter != sdfFields.end(); ++iter) {
        out << ">  <" << iter->first << ">" << " (" << iter->second.size()
                << ")" << endl;
        for (string value : iter->second) {
            out << value << endl;
        }
        out << endl;
    }

    out << "$$$$" << endl;
}

const map<int, int> Molecule::mapAtomsForOutput(
        const bool & includeLonePairs) const {
    map<int, int> atomNoMap;

    int no = 0;
    for (auto & atom : atoms) {
        if (!includeLonePairs
                && atom->getAtomType().getType() == AtomType::AtomTypeId::LP) {
            atomNoMap[atom->getAtomNo()] = -1;
        } else {
            atomNoMap[atom->getAtomNo()] = no;
            ++no;
        }
    }

    return atomNoMap;
}

// helper function for reading molecules
static vector<Molecule::MoleculePtr> readMolecules(istream & in,
        const Molecule::FileType fileType) {
    vector<Molecule::MoleculePtr> mols;
    while (true) {
        auto molecule = make_unique<Molecule>();
        if (!molecule->loadMolecule(fileType, in)) {
            break;
        }
        mols.push_back(move(molecule));
    }
    return mols;
}

bool Molecule::readMoleculeFromFile(const string & fileName) {
    const FileType fileType = Molecule::fileNameToType(fileName);
    auto func =
            [this, fileType] (istream &in) -> bool {return loadMolecule(fileType, in);};
    return readObjectFromFile<bool>(fileName, func);
}

unique_ptr<Molecule> Molecule::readAndInitializeMolecule(
        const string & fileName) {
    auto molecule = make_unique<Molecule>();
    auto ok = molecule->readMoleculeFromFile(fileName);
    assert(ok);
    molecule->initialize();
    return molecule;
}

vector<Molecule::MoleculePtr> Molecule::readMoleculesFromFile(
        const string & fileName) {

    auto fileType = fileNameToType(fileName);
    auto func = [fileType] (istream & in) {return readMolecules(in, fileType);};
    return readObjectFromFile<vector<Molecule::MoleculePtr>>(fileName, func);
}

void Molecule::writeToFile(const FileType & fileType, ostream & out,
        const bool & includeLonePairs, const string & comment) const {
    writeToFile(fileType, out, getCoords(), includeLonePairs, comment);
}

void Molecule::writeToFile(const FileType & fileType, ostream & out,
        const CoordMatrix & otherCoords, const bool & includeLonePairs,
        const string & comment) const {
    switch (fileType) {
    case FileType::MOL2:
        writeMol2File(out, otherCoords, includeLonePairs, comment);
        return;
    case FileType::SDF:
        writeSdfFile(out, otherCoords, comment);
        return;
    default:
        throw std::invalid_argument("Unknown molecule file type");
    }
}

// helper function for writing molecules
static void writeMolecules(ostream & out, vector<Molecule::MoleculePtr> & mols,
        const Molecule::FileType fileType) {
    for (auto & molecule : mols) {
        molecule->writeToFile(fileType, out);
    }
}

void Molecule::writeMoleculesToFile(string fileName,
        vector<MoleculePtr> & mols) {
    const auto fileType = Molecule::fileNameToType(fileName);

    auto func =
            [ & mols, fileType](ostream & out) {writeMolecules(out, mols, fileType);};
    writeObjectToFile(fileName, func);
}

Molecule::FileType Molecule::fileNameToType(const string & fileName) {
    string test(fileName);
    toLowerCase(test);
    const auto fileType =
            endsWith(test, "mol2") || endsWith(test, "mol2.gz") ?
                    FileType::MOL2 : FileType::SDF;
    return fileType;
}

bool Molecule::loadMolecule(FileType fileType, istream &in) {
    switch (fileType) {
    case FileType::MOL2:
        return loadMol2Molecule(in);
    case FileType::SDF:
        return loadSdfMolecule(in);
    default:
        throw std::invalid_argument("Unknown molecule file type");
    }
    return false;
}

/**
 * Test to see if two atoms have the same 2D atom type.
 *
 * @param thisAtom
 * @param otherAtom
 * @return
 */
static bool same2DAtomType(const Atom & thisAtom, const Atom & otherAtom) {
    auto thisType = thisAtom.getAtomType().getType();
    auto otherType = otherAtom.getAtomType().getType();

    if (thisType == otherType) {
        return true;
    }

    // for CCDC multi conformer molecules there can be trigonal NPL3 or pyrimidal
    // N3 typing for the same atom depending on 3D geometry- if we have many more of
    // these cases we may need to set atom types for each conformer while loading
    // multiconformer molecules (e.g. see CCDC multiconformers for 1ywr)
    if (thisType == AtomTypeId::N3 && otherType == AtomTypeId::NPL3) {
        return true;
    }
    if (thisType == AtomTypeId::NPL3 && otherType == AtomTypeId::N3) {
        return true;
    }

    return false;
}

bool Molecule::same2dMolecule(const Molecule & otherMol) const {
    auto nAtoms = getAtoms().size();
    if (otherMol.getAtoms().size() != nAtoms) {
        return false;
    }
    auto nBonds = getBonds().size();
    if (otherMol.getBonds().size() != nBonds) {
        return false;
    }

    for (size_t i = 0; i < nAtoms; ++i) {
        if (!same2DAtomType(getAtom(i), otherMol.getAtom(i))) {
            return false;
        }
    }

    for (size_t i = 0; i < nBonds; ++i) {
        const Bond & bond1 = getBond(i);
        const Bond & bond2 = otherMol.getBond(i);
        if (!bond1.sameBondType(bond2)) {
            return false;
        }
        if (!same2DAtomType(bond1.getAtom1(), bond2.getAtom1())) {
            return false;
        }
        if (!same2DAtomType(bond1.getAtom2(), bond2.getAtom2())) {
            return false;
        }
    }

    return true;
}

boost::optional<const Bond &> Molecule::isBonded(const Atom & atom1,
        const Atom & atom2) const {
    return isBonded(atom1.getAtomNo(), atom2.getAtomNo());
}

boost::optional<const Bond &> Molecule::isBonded(const int atom1No,
        const int atom2No) const {
    auto pred = [atom1No, atom2No] (const unique_ptr<Bond> & bond) -> bool {
        int testAtom1No = bond->getAtom1().getAtomNo();
        int testAtom2No = bond->getAtom2().getAtomNo();
        return (testAtom1No == atom1No && testAtom2No == atom2No) ||
        (testAtom2No == atom1No && testAtom1No == atom2No);
    };
    auto value = findFirstInUniquePtrList<Bond>(bonds, pred);
    if (value) {
        return *value.get();
    }
    return boost::none;
}

void Molecule::createAtomNeighbourhood() {
    for (auto & atom : atoms) {
        atom->buildNeighbourhood(*this);
    }
    if (bondTable == nullptr) {
        createBondTable();
    }

    updateNeighbourHoodWithRings();
}

void Molecule::initialize() {
    createBondTable();
    MoleculeSssr moleculeSssr(*this);
    moleculeSssr.doSssr();
    createAtomNeighbourhood();
    assignAtomTypes();
    Huckel huckel(*this);
    huckel.findAromaticRings();
}

void Molecule::createBondTable() {
    int nAtoms = this->nAtoms();
    bondTable = make_unique<Array2D<Bond *>>(nAtoms, nAtoms);
    for (auto & bond : bonds) {
        bondTable->get(bond->getAtom1().getAtomNo(),
                bond->getAtom2().getAtomNo()) = bond.get();
        bondTable->get(bond->getAtom2().getAtomNo(),
                bond->getAtom1().getAtomNo()) = bond.get();
    }

    bonds13 = make_unique<Array2D<bool>>(nAtoms, nAtoms);
    int nBonds = this->nBonds();
    angles.clear();
    angles.reserve(nBonds * 4);
    int angleNo = 0;

    for (auto i = 0; i < nBonds; i++) {

        auto & bond1 = getBond(i);
        if (!AtomType::isRealAtom(bond1.getAtom1().getAtomTypeId())) {
            continue;
        }
        if (!AtomType::isRealAtom(bond1.getAtom2().getAtomTypeId())) {
            continue;
        }

        for (auto j = i + 1; j < nBonds; j++) {

            auto & bond2 = getBond(j);
            if (!AtomType::isRealAtom(bond2.getAtom1().getAtomTypeId())) {
                continue;
            }
            if (!AtomType::isRealAtom(bond2.getAtom2().getAtomTypeId())) {
                continue;
            }

            auto atom1 = bond1.getAtom1().getAtomNo();
            auto atom2 = bond1.getAtom2().getAtomNo();
            auto atom3 = bond2.getAtom1().getAtomNo();
            auto atom4 = bond2.getAtom2().getAtomNo();
            unique_ptr<Angle> angle = nullptr;

            if (atom1 == atom3) {
                angle = make_unique<Angle>(angleNo++, atom2, atom1, atom4);
            } else if (atom1 == atom4) {
                angle = make_unique<Angle>(angleNo++, atom2, atom1, atom3);
            }
            else if (atom2 == atom3) {
                angle = make_unique<Angle>(angleNo++, atom1, atom2, atom4);
            } else if (atom2 == atom4) {
                angle = make_unique<Angle>(angleNo++, atom1, atom2, atom3);
            }

            if (angle != nullptr) {
                bonds13->get(angle->getAtom1(), angle->getAtom3()) = true;
                bonds13->get(angle->getAtom3(), angle->getAtom1()) = true;
                angles.push_back(move(angle));
            }
        }
    }
    assert(angleNo == static_cast<int>(angles.size()));

    bonds14 = make_unique<Array2D<bool>>(nAtoms, nAtoms);
    torsions.clear();
    torsions.reserve(nBonds * 16);
    int torsionNo = 0;

    for (auto i = 0; i < nBonds; i++) {
        auto & bond2 = getBond(i);

        auto atom2 = bond2.getAtom1().getAtomNo();
        auto atom3 = bond2.getAtom2().getAtomNo();

        for (auto j = 0; j < nBonds; j++) {
            if (i == j) {
                continue;
            }
            auto & bond1 = getBond(j);
            int atom1 = -1;
            if (bond1.getAtom1().getAtomNo() == atom2) {
                atom1 = bond1.getAtom2().getAtomNo();
            }
            if (bond1.getAtom2().getAtomNo() == atom2) {
                atom1 = bond1.getAtom1().getAtomNo();
            }
            if (atom1 != -1) {
                for (auto k = 0; k < nBonds; k++) {
                    if (k == i || k == j) {
                        continue;
                    }
                    auto & bond3 = getBond(k);
                    int atom4 = -1;
                    if (bond3.getAtom1().getAtomNo() == atom3) {
                        atom4 = bond3.getAtom2().getAtomNo();
                    }
                    if (bond3.getAtom2().getAtomNo() == atom3) {
                        atom4 = bond3.getAtom1().getAtomNo();
                    }

                    if (atom4 != -1) {
                        auto torsion = make_unique<Torsion>(torsionNo++, atom1, atom2, atom3,
                                atom4, bond2.getBondNo());
                        bonds14->get(atom1, atom4) = true;
                        bonds14->get(atom4, atom1) = true;
                        torsions.push_back(move(torsion));
                    }
                }

            }
        }

    }
    assert(torsionNo == static_cast<int>(torsions.size()));
}

Bond * Molecule::getBond(int atom1, int atom2) const {
    auto bond = bondTable->get(atom1, atom2);
    return bond;
}

void Molecule::updateNeighbourHoodWithRings() {
    for (auto const & atom : atoms) {
        atom->getNeighbourhoodToEdit( { }).clearRings();
    }

    for (auto const & ring : rings) {
        for (auto const & atom : ring->getAtoms()) {
            atom->getNeighbourhoodToEdit( { }).addRing(ring.get());
        }
    }

}

int Molecule::assignAtomTypes() {
    REPORT(Reporter::DEBUG) << "Assigning atom types for molecule "
            << getName();
    int typesChanged = 0;
    AtomEditKey atomEditKey;
    for (auto & atom : atoms) {
        auto checkType = atom->checkAtomType(*this);
        assert(checkType != AtomType::AtomTypeId::ATM_NONE);
        if (checkType != atom->getAtomTypeId()) {
            typesChanged++;
            const auto & newAtomType = AtomType::typeFromTypeId(checkType);
            REPORT(Reporter::DEBUG) << "Set atom type for atom "
                    << atom->info() << " to type " << newAtomType.getName();
            atom->setAtomType(atomEditKey, newAtomType);
        }
    }

    BondEditKey bondEditKey;
    for (auto & bond : bonds) {
        auto checkType = bond->checkBondType(*this);
        if (checkType != bond->getBondTypeId()) {
            const auto & newBondType = BondType::typeFromTypeId(checkType);
            REPORT(Reporter::DEBUG) << "Set bond type for bond "
                    << bond->info() << " to type " << newBondType.getName();
            bond->setBondType(bondEditKey, newBondType);
        }
    }
#ifndef NDEBUG
    for (auto & atom : atoms) {
        auto checkType = atom->checkAtomType(*this);
        if (checkType != atom->getAtomTypeId()) {
            const auto & newAtomType = AtomType::typeFromTypeId(checkType);
            REPORT(Reporter::FATAL) << "Atom " << atom->info()
                    << " assigned type " << atom->getAtomType().getName()
                    << " second check assignment is " << newAtomType.getName();
        }
        assert(checkType == atom->getAtomTypeId());
    }
#endif

    return typesChanged;
}

void Molecule::centerCoords() {
    coords = centerCoordMatrix(coords);
}

CoordVector Molecule::centroid() const {
    return GarethUtil::centroid(coords);
}

int Molecule::nRealAtoms() const {
    auto accum = [] (int count, const unique_ptr<Atom> & atom) -> int {
        auto inc = AtomType::isRealAtom(atom->getAtomTypeId()) ? 1 : 0;
        return count + inc;
    };
    return reduce<unique_ptr<Atom>, int>(atoms, 0, accum);
}

int Molecule::nHeavyAtoms() const {
    auto accum = [] (int count, const unique_ptr<Atom> & atom) -> int {
        auto inc = AtomType::isHeavy(atom->getAtomTypeId()) ? 1 : 0;
        return count + inc;
    };
    return reduce<unique_ptr<Atom>, int>(atoms, 0, accum);
}

int Molecule::nRealBonds() const {
    auto accum =
            [] (int count, const unique_ptr<Bond> & bond) -> int {
                auto realBond = AtomType::isRealAtom(bond->getAtom1().getAtomTypeId()) &&
                AtomType::isRealAtom(bond->getAtom2().getAtomTypeId());
                auto inc = realBond ? 1 : 0;
                return count + inc;
            };
    return reduce<unique_ptr<Bond>, int>(bonds, 0, accum);
}

} /* namespace GarethMol */
