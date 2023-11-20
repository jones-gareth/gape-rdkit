//
// Created by Gareth Jones on 12/28/2022.
//

#define CATCH_CONFIG_MAIN

#include "catch2/catch.hpp"
#include "gape/GapeSettings.h"
#include "mol/AcceptorAtomFeature.h"
#include "mol/DonorHydrogenFeature.h"
#include "util/Reporter.h"
#include "mol/HydrogenBondingType.h"
#include "mol/Solvate.h"
#include "mol/PartialCharge.h"

using namespace Gape;
using namespace RDKit;

std::map<const Atom*, std::shared_ptr<const HydrogenBondingType>>
findDonorsOrAcceptors(HydrogenBondType bondType, RWMol& mol, const GapeSettings& settings, bool solvate = true)
{
	mol.setProp(common_properties::_Name, "Unknown molecule");
	MolOps::addHs(mol);
	MolOps::sanitizeMol(mol);
	if (solvate)
	{
		solvateMolecule(settings.getSolvationRules(), mol);
	}

	findPartialCharges(settings.getPartialCharges(), mol);
	auto donorsAndAcceptors =
		bondType == HydrogenBondType::Acceptor
			? findHydrogenBondAcceptors(settings.getHydrogenBondingTypes(), mol)
			: findHydrogenBondDonors(settings.getHydrogenBondingTypes(), mol);
	return donorsAndAcceptors;
}

bool chargedAcceptor(const Atom* const atom)
{
	return PartialCharge::getPartialCharge(atom) < AcceptorAtomFeature::maxAcceptorPartialCharge;
}

bool chargedDonor(const Atom* const atom)
{
	return PartialCharge::getPartialCharge(atom) > DonorHydrogenFeature::minDonorPartialCharge;
}

TEST_CASE("Finding hydrogen bonding types works as expected", "[hydrogenBondingType]")
{
	Reporter::setMinReportingLevel(Reporter::DEBUG);
	const GapeSettings settings;

	SECTION("Terminal Phosphate")
	{
		auto const mol = "OP(=O)(O)C"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings);
		CHECK(features.size() == 3);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 8);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			CHECK(bondingType->name == "Terminal Phosphate");
			CHECK(chargedAcceptor(atom));
		}
	}

	SECTION("Phosphinyl")
	{
		auto const mol = "OP(=O)(C)C"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings);
		CHECK(features.size() == 2);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 8);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			CHECK(bondingType->name == "Phosphinyl");
			CHECK(chargedAcceptor(atom));
		}
	}

	SECTION("Carboxylate")
	{
		auto const mol = "OC(=O)C"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings);
		CHECK(features.size() == 2);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 8);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			CHECK(bondingType->name == "Carboxylate");
			CHECK(chargedAcceptor(atom));
		}
	}

	SECTION("Asymmetric Het6 N")
	{
		auto const mol = "N1=NC=CC=C1"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings);
		CHECK(features.size() == 2);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			CHECK(bondingType->name == "Asymmetric Het6 N");
			CHECK(!chargedAcceptor(atom));
		}
	}

	SECTION("Acid O")
	{
		auto const mol = "O=C(O)C"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings, false);
		CHECK(features.size() == 2);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 8);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			if (atom->getTotalNumHs(true) == 1)
			{
				CHECK(bondingType->name == "Hydroxyl");
				CHECK(!chargedAcceptor(atom));
			}
			else
			{
				CHECK(bondingType->name == "Acid O");
				CHECK(!chargedAcceptor(atom));
			}
		}
	}

	SECTION("Symmetric Het5 N")
	{
		auto const mol = "N1=COC=C1"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings);
		CHECK(features.size() == 2);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			if (atom->getAtomicNum() == 7)
			{
				CHECK(bondingType->name == "Symmetric Het5 N");
				CHECK(!chargedAcceptor(atom));
			}
			else
			{
				CHECK(atom->getAtomicNum() == 8);
				CHECK(bondingType->name == "Het5 O Aromatic");
				CHECK(!chargedAcceptor(atom));
			}
		}
	}

	SECTION("Sulfoxide")
	{
		auto const mol = "O=S(C)C"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 8);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			CHECK(bondingType->name == "Sulfoxide");
			CHECK(!chargedAcceptor(atom));
		}
	}

	SECTION("Primary Amine")
	{
		auto const mol = "NCC"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings, false);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			CHECK(bondingType->name == "Primary Amine");
			CHECK(!chargedAcceptor(atom));
		}
	}

	SECTION("Asymmetric Het5 N")
	{
		auto const mol = "N1=NC=CO1"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings);
		CHECK(features.size() == 3);
		for (const auto& [atom, bondingType] : features)
		{
			if (atom->getAtomicNum() == 7)
			{
				CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
				CHECK(bondingType->name == "Asymmetric Het5 N");
				CHECK(!chargedAcceptor(atom));
			}
			else
			{
				CHECK(atom->getAtomicNum() == 8);
				CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
				CHECK(bondingType->name == "Het5 O Aromatic");
				CHECK(!chargedAcceptor(atom));
			}
		}
	}

	SECTION("Thiocarbonyl")
	{
		auto const mol = "S=C(C)C"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 16);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			CHECK(bondingType->name == "Thiocarbonyl");
			CHECK(!chargedAcceptor(atom));
		}
	}

	SECTION("Hydroxyl")
	{
		auto const mol = "OCC"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 8);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			CHECK(bondingType->name == "Hydroxyl");
			CHECK(!chargedAcceptor(atom));
		}
	}

	SECTION("Symmetric Het6 N")
	{
		auto const mol = "n1ccccc1"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			CHECK(bondingType->name == "Symmetric Het6 N");
			CHECK(!chargedAcceptor(atom));
		}
	}

	SECTION("Terminal Sulfate")
	{
		auto const mol = "O=S(=O)(O)C"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings);
		CHECK(features.size() == 3);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 8);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			CHECK(bondingType->name == "Terminal Sulfate");
			CHECK(chargedAcceptor(atom));
		}
	}

	SECTION("Tertiary Amine")
	{
		auto const mol = "N(C)(C)C"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings, false);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			CHECK(bondingType->name == "Tertiary Amine");
			CHECK(!chargedAcceptor(atom));
		}
	}

	SECTION("Ester O")
	{
		auto const mol = "O=C(C)OC"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings);
		CHECK(features.size() == 2);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 8);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			if (atom->getTotalDegree() == 1)
			{
				CHECK(bondingType->name == "Ester O");
				CHECK(!chargedAcceptor(atom));
			}
			else
			{
				CHECK(bondingType->name == "Ether");
				CHECK(!chargedAcceptor(atom));
			}
		}
	}

	SECTION("Carbamate O")
	{
		auto const mol = "O=C(OC)NC"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings);
		CHECK(features.size() == 3);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			if (atom->getAtomicNum() == 8 && atom->getTotalDegree() == 1)
			{
				CHECK(bondingType->name == "Carbamate O");
				CHECK(!chargedAcceptor(atom));
			}
			else if (atom->getAtomicNum() == 8 && atom->getTotalDegree() == 2)
			{
				CHECK(bondingType->name == "Ether");
				CHECK(!chargedAcceptor(atom));
			}
			else
			{
				CHECK(atom->getAtomicNum() == 7);
				CHECK(bondingType->name == "Secondary Amine");
				CHECK(!chargedAcceptor(atom));
			}
		}
	}

	SECTION("Nitrile")
	{
		auto const mol = "N#CC"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			CHECK(bondingType->name == "Nitrile");
			CHECK(!chargedAcceptor(atom));
		}
	}

	SECTION("Imine")
	{
		auto const mol = "N(C)=CC"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			CHECK(bondingType->name == "Imine");
			CHECK(!chargedAcceptor(atom));
		}
	}

	SECTION("Ketone")
	{
		auto const mol = "O=CC"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 8);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			CHECK(bondingType->name == "Ketone");
			CHECK(!chargedAcceptor(atom));
		}
	}

	SECTION("Secondary Amine")
	{
		auto const mol = "N(C)CC"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings, false);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			CHECK(bondingType->name == "Secondary Amine");
			CHECK(!chargedAcceptor(atom));
		}
	}

	SECTION("Phenol OH")
	{
		auto const mol = "Oc1ccccc1"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 8);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			CHECK(bondingType->name == "Phenol OH");
			CHECK(!chargedAcceptor(atom));
		}
	}

	SECTION("Ether")
	{
		auto const mol = "COC"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 8);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			CHECK(bondingType->name == "Ether");
			CHECK(!chargedAcceptor(atom));
		}
	}

	SECTION("Het5 O")
	{
		auto const mol = "O1C=CCC(C)(C)1"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 8);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			CHECK(bondingType->name == "Het5 O");
			CHECK(!chargedAcceptor(atom));
		}
	}

	SECTION("Het5 O Aromatic")
	{
		auto const mol = "O1C=CC=C1"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 8);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			CHECK(bondingType->name == "Het5 O Aromatic");
			CHECK(!chargedAcceptor(atom));
		}
	}

	SECTION("Nitro O")
	{
		auto const mol = "O=[N+]([O-])C"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Acceptor, *mol, settings);
		CHECK(features.size() == 2);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 8);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
			CHECK(bondingType->name == "Nitro O");
			CHECK(chargedAcceptor(atom));
		}
	}

	SECTION("Het5 NH")
	{
		auto const mol = "N1C=CC=C1"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Donor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Donor);
			CHECK(bondingType->name == "Het5 NH");
			CHECK(!chargedDonor(atom));
		}
	}

	SECTION("sp2 N+H (1) and (2)")
	{
		auto const mol = "N(C)=C(C)NC"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Donor, *mol, settings);
		CHECK(features.size() == 2);
		bool found1 = false, found2 = false;
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Donor);
			if (bondingType->name == "sp2 N+H (1)")
			{
				found1 = true;
			}
			else if (bondingType->name == "sp2 N+H (2)")
			{
				found2 = true;
			}
			CHECK(chargedDonor(atom));
		}
		CHECK(found1);
		CHECK(found2);
	}

	SECTION("sp2 N+H (3)")
	{
		auto const mol = "[nH+]1ccccc1"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Donor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Donor);
			CHECK(bondingType->name == "sp2 N+H (3)");
			CHECK(chargedDonor(atom));
		}
	}

	SECTION("Acid OH")
	{
		auto const mol = "OC=O"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Donor, *mol, settings, false);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 8);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Donor);
			CHECK(bondingType->name == "Acid OH");
			CHECK(!chargedDonor(atom));
		}
	}

	SECTION("sp2 N+H2 (1)")
	{
		auto const mol = "[NH2+]=CC"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Donor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Donor);
			CHECK(bondingType->name == "sp2 N+H2 (1)");
			CHECK(chargedDonor(atom));
		}
	}


	SECTION("sp2 N+H2 (2) ortho amino pyridine")
	{
		auto const mol = "n1c(N)cncc1"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Donor, *mol, settings);
		CHECK(features.size() == 2);
		bool found1=false, found2=false;
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Donor);
			if (bondingType->name == "sp2 N+H (3)")
			{
				found1 = true;
			}
			else if (bondingType->name == "sp2 N+H2 (2)")
			{
				found2 = true;
			}
			CHECK(chargedDonor(atom));
		}
		CHECK(found1);
		CHECK(found2);
	}

	SECTION("sp2 N+H2 (3) para amino pyridine")
	{
		auto const mol = "n1cnc(N)cc1"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Donor, *mol, settings);
		CHECK(features.size() == 2);
		bool found1=false, found2=false;
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Donor);
			if (bondingType->name == "sp2 N+H (3)")
			{
				found1 = true;
			}
			else if (bondingType->name == "sp2 N+H2 (3)")
			{
				found2 = true;
			}
			CHECK(chargedDonor(atom));
		}
		CHECK(found1);
		CHECK(found2);
	}

	SECTION("sp2 N+H2 (4) and (5)")
	{
		auto const mol = "N=C(C)N"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Donor, *mol, settings);
		CHECK(features.size() == 2);
		bool found1=false, found2=false;
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Donor);
			if (bondingType->name == "sp2 N+H2 (4)")
			{
				found1 = true;
			}
			else if (bondingType->name == "sp2 N+H2 (5)")
			{
				found2 = true;
			}
			CHECK(chargedDonor(atom));
		}
		CHECK(found1);
		CHECK(found2);
	}

	SECTION("sp3 N+H3")
	{
		auto const mol = "NC"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Donor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Donor);
			CHECK(bondingType->name == "sp3 N+H3");
			CHECK(chargedDonor(atom));
		}
	}

	SECTION("sp3 N+H2")
	{
		auto const mol = "N(C)C"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Donor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Donor);
			CHECK(bondingType->name == "sp3 N+H2");
			CHECK(chargedDonor(atom));
		}
	}

	SECTION("sp3 N+H")
	{
		auto const mol = "N(C)(C)C"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Donor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Donor);
			CHECK(bondingType->name == "sp3 N+H");
			CHECK(chargedDonor(atom));
		}
	}

	SECTION("OH")
	{
		auto const mol = "CO"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Donor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 8);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Donor);
			CHECK(bondingType->name == "OH");
			CHECK(!chargedDonor(atom));
		}
	}

	SECTION("Primary Amide NH2")
	{
		auto const mol = "NC=O"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Donor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Donor);
			CHECK(bondingType->name == "Primary Amide NH2");
			CHECK(!chargedDonor(atom));
		}
	}

	SECTION("Secondary Amide NH")
	{
		auto const mol = "N(C)C=O"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Donor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Donor);
			CHECK(bondingType->name == "Secondary Amide NH");
			CHECK(!chargedDonor(atom));
		}
	}

	SECTION("Phenyl NH")
	{
		auto const mol = "N(C)c1ccccc1"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Donor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Donor);
			CHECK(bondingType->name == "Phenyl NH");
			CHECK(!chargedDonor(atom));
		}
	}

	SECTION("Phenyl NH2")
	{
		auto const mol = "Nc1ccccc1"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Donor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Donor);
			CHECK(bondingType->name == "Phenyl NH2");
			CHECK(!chargedDonor(atom));
		}
	}

	SECTION("Imine NH")
	{
		auto const mol = "N=C"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Donor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Donor);
			CHECK(bondingType->name == "Imine NH");
			CHECK(!chargedDonor(atom));
		}
	}

	SECTION("Phenyl OH")
	{
		auto const mol = "Oc1ccccc1"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Donor, *mol, settings);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 8);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Donor);
			CHECK(bondingType->name == "Phenyl OH");
			CHECK(!chargedDonor(atom));
		}
	}

	SECTION("Primary Amine NH2")
	{
		auto const mol = "NC"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Donor, *mol, settings, false);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Donor);
			CHECK(bondingType->name == "Primary Amine NH2");
			CHECK(!chargedDonor(atom));
		}
	}

	SECTION("Secondary Amine NH")
	{
		auto const mol = "N(C)C"_smiles;
		const auto features = findDonorsOrAcceptors(HydrogenBondType::Donor, *mol, settings, false);
		CHECK(features.size() == 1);
		for (const auto& [atom, bondingType] : features)
		{
			CHECK(atom->getAtomicNum() == 7);
			CHECK(bondingType->hydrogenBondType == HydrogenBondType::Donor);
			CHECK(bondingType->name == "Secondary Amine NH");
			CHECK(!chargedDonor(atom));
		}
	}
}
