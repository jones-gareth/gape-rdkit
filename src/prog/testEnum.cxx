#include <iostream>
#include "util/EnumIter.h"

enum class Color {
	Red, Green, Blue, First = Red, Last = Blue
};

int main() {
	using namespace Gape;

	for (auto e : EnumIter<Color>()) {
		std::cout << ((int) e) << std::endl;
	}
}
