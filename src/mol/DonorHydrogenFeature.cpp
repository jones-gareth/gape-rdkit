#include "DonorHydrogenFeature.h"

namespace Gape
{
	thread_local double DonorHydrogenFeature::hBondLen;
	thread_local double DonorHydrogenFeature::chargeFactor;
	thread_local double DonorHydrogenFeature::maxDonorAngle;
	thread_local double DonorHydrogenFeature::minDonorAngle;
	thread_local bool DonorHydrogenFeature::scoreDonorAtoms;

}
