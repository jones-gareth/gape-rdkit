# GAPE C++ Codebase with RDKit Molecular Handling

WIP.  Intended to be a C++ version of GAPE (https://github.com/Gareth Jones-gareth/gape)
with custom molecular handling replaced by RDKit.

## Future work

* Use RDKit for structure handling
* Replace TAFF with RDkit MMFF94s
* Add GAPE GA representation (fitness function, chromosome, crossover, mutation etc)

## Includes

* Steady state with no duplicates Genetic Algorithm (similar to that used in GOLD and GASP)
* C++ implementation of DIFGAPE (not validated)
* Routines for reading MOL2 and SDF file formats
* Smiles and Smarts parsers
* Substructure search
* Implementation of Tripos Associates Forcefield

## Requires

* Boost
* GTest
* Eigen

## Build

Edit CMakeLists.txt to select Boost, GTest and Eigen paths (or you can invoke cmake
with the appropriate defines). If you download GTest (googletest) make sure you build it after downloading (`cmake .; make` in the toplevel directory). Then:

Linking to RDKit cmake may select `${RDBASE}cmake-build-release/lib` as the library directory, even when
`${RDBASE}/lib` is specified in CMakeLists.txt.  Workaround is to remove cmake build directories from
RDKit source.

```sh
mkdir build
cd build
cmake ..
```

## TODO

* Island model
* Corner flipping
* User defined features
* 
