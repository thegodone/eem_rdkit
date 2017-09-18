# eem_rdkit
Implement EEM partial charges computation using RDKit toolbox

This is a simplified version that use RDKit to compute EEM partial charges similar to this:
http://www.fi.muni.cz/~xracek/neemp/ under the GNU GPLv3 licence.

it compiles on mac os sierra:


% for Lapack version
clang++ -g -Wall -I$RDBASEINCLUDE -I$BOOSTINCLUDE eem_lapack.cpp -o eem_lapack -L$RDBASELIB -L$BOOSTLIB -lRDGeneral -lGraphMol -lSmilesParse -lFileParsers -lSubstructMatch -lForceField -lForceFieldHelpers  -framework Accelerate


% for Eigen version
clang++ -g -Wall -I$RDBASEINCLUDE -I$BOOSTINCLUDE eem_eigen3.cpp -o eem_eigen3 -L$RDBASELIB -L$BOOSTLIB -lRDGeneral -lGraphMol -lSmilesParse -lFileParsers -lSubstructMatch -lForceField -lForceFieldHelpers  -I/usr/include/eigen3






current version of Eigen has memory leak :
valgrind --leak-check=full --show-leak-kinds=all -v ./eem_eigen3
