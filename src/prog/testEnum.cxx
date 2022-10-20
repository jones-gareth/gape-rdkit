#include <iostream>
#include "RocsSdfParser.h"
#include "EnumIter.h"

enum class Color {
	Red, Green, Blue, First = Red, Last = Blue
};

int main() {
	using namespace GarethUtil;
	using namespace Difgape;

	for (auto e : EnumIter<Color>()) {
		std::cout << ((int) e) << std::endl;
	}
	for (Difgape::ScoreName e : EnumIter<Difgape::ScoreName>()) {
		std::cout << ((int) e) << " : "
				<< RocsSdfParser::scoreNameToString(e) << std::endl;
	}
}
