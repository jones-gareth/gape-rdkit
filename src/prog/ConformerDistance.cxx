/**
 * Program to investigate different distance metrics for interconformer distances
 */

#include <iostream>
#include "../mol/MultiConformerGyrationShapeEvaluation.h"
#include "../mol/MultiConformerRmsEvalution.h"

using namespace std;
using namespace GarethMol;
using namespace GarethUtil;

int main(int argc, char* argv[]) {

    assert(argc == 2);
    string file(argv[1]);

    auto molecule = MulticonformerMolecule::readAndInitializeMolecule(file);
    MultiConformerGyrationShapeEvalution shapeEval(*molecule);
    MultiConformerRmsEvalution rmsEval(*molecule);
    rmsEval.evaluateRms(false);

    ofstream out("conformer_distance.csv");

    out << "NO1" << "\t" << "NO2" << "\t"
                        << "RMS" << "\t" << "SHAPE" << endl;
    for (auto i = 0ul; i < molecule->nConformers(); i++) {
        for (auto j = i + 1; j < molecule->nConformers(); j++) {
            auto rmsDistance = rmsEval.getRmsValues().get(i, j);
            auto shapeDistance = shapeEval.getDistanceMatrix().get(i, j);
            out << to_string(i + 1) << "\t" << to_string(j + 1) << "\t"
                    << rmsDistance << "\t" << shapeDistance << endl;
        }
    }


    out.close();
    return EXIT_SUCCESS;
}
