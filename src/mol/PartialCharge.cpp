#include "PartialCharge.h"
#include "util/Reporter.h"
#include <GraphMol/Substruct/SubstructMatch.h>

namespace Gape
{
	PartialCharge::PartialCharge(std::string n, const int fC, const double pC, std::string sma, const bool m) :
		name(std::move(n)), formalCharge(fC), partialCharge(pC), smarts(std::move(sma)), multiple(m)
	{
		try
		{
			query = SmartsToMol(smarts);
		}
		catch (...)
		{
			REPORT(Reporter::WARN) << "Failed to parse partial charge " << name << " pattern " << smarts;
		}
		if (query == nullptr)
		{
			REPORT(Reporter::WARN) << "Failed to parse partial charge " << name << " pattern " << smarts;
		}
	}

	void findPartialCharges(const PartialChargesList& partialChargesList, ROMol& mol)
	{
        SubstructMatchParameters params;
        params.uniquify = false;
		for (const auto& partialCharge : partialChargesList) {

			if (partialCharge->query != nullptr) {
				for (const auto& match : SubstructMatch(mol, *partialCharge->query, params)) {

				}
			}
		}
	}
} // namespace Gape
