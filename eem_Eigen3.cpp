//
//  Copyright (C) 2017 Guillaume GODIN
//  conversion of neemp code to use rdkit
//
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <RDGeneral/RDLog.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <boost/foreach.hpp>
#include <Eigen/Dense>


using namespace Eigen;

/*
clang++ -g -Wall -I$RDBASEINCLUDE -I$BOOSTINCLUDE eemtest2.cpp -o eemtest2 -L$RDBASELIB -L$BOOSTLIB -lRDGeneral -lGraphMol -lSmilesParse -lFileParsers -lSubstructMatch -lForceField -lForceFieldHelpers  -I/usr/include/eigen3
clang++ -g -Wall -I$RDBASEINCLUDE -I$BOOSTINCLUDE eemtest.cpp -o eemtest -L$RDBASELIB -L$BOOSTLIB -lRDGeneral -lGraphMol -lSmilesParse -lFileParsers -lSubstructMatch -lForceField -lForceFieldHelpers  -framework Accelerate
*/

/* Symbols for chemical elements */
//static const char * const elems[] = {"??", "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr"};

/* Electronegativities of the chemical elements */
const float electronegativies[] = {0.0, 2.2, 0, 0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98, 0, 0.93, 1.31, 1.61, 1.9, 2.19, 2.58, 3.16, 0, 0.82, 1, 1.36, 1.54, 1.63, 1.66, 1.55, 1.83, 1.88, 1.91, 1.9, 1.65, 1.81, 2.01, 2.18, 2.55, 2.96, 3, 0.82, 0.95, 1.22, 1.33, 1.6, 2.16, 1.9, 2.2, 2.28, 2.2, 1.93, 1.69, 1.78, 1.96, 2.05, 2.1, 2.66, 2.6, 0.79, 0.89, 1.1, 1.12, 1.13, 1.14, 1.13, 1.17, 1.2, 1.2, 1.2, 1.22, 1.23, 1.24, 1.25, 1.1, 1.27, 1.3, 1.5, 2.36, 1.9, 2.2, 2.2, 2.28, 2.54, 2, 2.04, 2.33, 2.02, 2, 2.2, 0, 0.7, 0.9, 1.1, 1.3, 1.5, 1.38, 1.36, 1.28, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3};


// hardcoding of the parameters for the moment!
// need to find a way to read xml and create same type of arrays
const float kappa = 0.1960;

const float A1[] = {0.0,2.3594,0.0,0.0,0.0,0.0,2.4541,2.5908,2.7130,0.0,0.0,0.0,0.0,0.0,0.0,0.0,2.3833};
const float A2[] = {0.0,0.0,0.0,0.0,0.0,0.0,2.4726,2.5409,2.6766,0.0,0.0,0.0,0.0,0.0,0.0,0.0,2.4956};
const float A3[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

const float B1[] = {0.0,0.5962,0.0,0.0,0.0,0.0,0.2591,0.3316,0.5028,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.4564};
const float B2[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.2268,0.2319,0.4992,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.1493};
const float B3[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};


//extern void dspsv_(char *uplo, int *n, int *nrhs, double *ap, int *ipiv, double *b, int *ldb, int *info);

// function to retreive the atomtype value based on the highest bond type of an atom
// in the publication they don't have access to "Aromatic type"
double getAtomtype(const RDKit::ROMol *mol, const RDKit::Atom *atom) {
      double t=1.0;
      double a;
      RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
      boost::tie(nbrIdx, endNbrs) = mol->getAtomNeighbors(atom);
      while (nbrIdx != endNbrs) {
        const RDKit::Bond *bond = mol->getBondBetweenAtoms(*nbrIdx, atom->getIdx());
        a = bond->getBondTypeAsDouble();
        if (a == 1.5) {
             t = 2.0;
            }
        if (a == 2.0) {
             t = 2.0;
            }
        if (a == 3.0) {
             t = 3.0;
            }
        ++nbrIdx;
  }
  return t;
}

/* Calculate average electronegativity of a molecule (harmonic mean) */
double m_calculate_avg_electronegativity(const RDKit::ROMol *mol) {
    double hsum = 0.0;
    for(int i = 0; i < mol->getNumAtoms(); i++)
        hsum += 1.0 / electronegativies[mol->getAtomWithIdx(i)->getAtomicNum()];
    return mol->getNumAtoms() / hsum;
}

/* Calculate sum and average charge of atoms in the molecule */
double m_calculate_charge_stats(const RDKit::ROMol *mol) {
    double sum = 0.0;
    for(int i = 0; i < mol->getNumAtoms(); i++)
        sum += mol->getAtomWithIdx(i)->getFormalCharge();
    return sum;
}



/* Print matrix in packed storage format */
void print_matrix_full(const double * const A, long int n) {

    assert(A != NULL);

#define IDX(x, y) (x * n + y)
    for(long int i = 0; i < n; i++) {
        for(long int j = 0; j < n; j++)
            printf("%6.4f ", A[IDX(i, j)]);

        printf("\n");
    }
#undef IDX
}



void fill_EEM_matrix_full(double * const A, RDKit::ROMol *mol, double * dist3D) {
    const long int n = mol->getNumAtoms();
    /* Following #define works only for i <= j */
    // upper matrix definition
    #define IDX(x, y) (x * (n+1) + y)
    /* Fill the full n * n block */
    for(long int i = 0; i < n; i++) {
        int idx=mol->getAtomWithIdx(i)->getAtomicNum();
        double t = getAtomtype(mol,mol->getAtomWithIdx(i));
        double v;
        if (t== 1.0) {
                v=  B1[idx];
        }
        if (t== 2.0) {
                v = B2[idx];
        }
        if (t== 3.0) {
                v = B3[idx];
        }
        //std::cout << "v: " << v <<  "\n";

        A[IDX(i, i)] = v;
        for(long int j = i + 1; j < n; j++) {
                A[IDX(i, j)] = kappa / dist3D[i*n+j];
                A[IDX(j, i)] = A[IDX(i, j)];

            }
        }
    /* Fill last column & row */
    for(long int i = 0; i < n; i++) {
        A[IDX(n, i)] = 1.0f; // column
        A[IDX(i, n)] = -1.0f; // row

    }


    /* Set the bottom right element to zero */
    A[(n+1)*(n+1)+n] = 0.0f;
    #undef IDX
}


/* Calculate charges for a particular kappa_data structure */
void calculate_charges(RDKit::ROMol *mol, double * dist3D) {
        const int n = mol->getNumAtoms();
        void *tmp1 = NULL;
        void *tmp2 = NULL;
        posix_memalign(&tmp1, 64, ((n + 1) * (n + 1)) * sizeof(double));
        posix_memalign(&tmp2, 64, (n + 1) * sizeof(double));
        double *Af = (double *) tmp1;
        double *b = (double *) tmp2;

        fill_EEM_matrix_full(Af, mol, dist3D);
        std::cout << "CreateData:\nA:\n";

        //print_matrix_full(Af,  n+1);

        /* Fill vector b i.e. -A */
        for(int j = 0; j < n; j++) {
            double t = getAtomtype(mol,mol->getAtomWithIdx(j));
            int idx=mol->getAtomWithIdx(j)->getAtomicNum();
            if (t== 1.0) {
                b[j] = -A1[idx];
            }
            if (t== 2.0) {
                b[j] = -A2[idx];
            }
            if (t== 3.0) {
                b[j] = -A3[idx];
            }
        }

        b[n] = m_calculate_charge_stats(mol);; // sum of charges
        /*
        std::cout << "\nb:\n";

        for(int j = 0; j <= n; j++) {
                std::cout << b[j] <<  "\n";
        }
        */
        Map<MatrixXd> A(Af, n+1, n+1);
        Map<VectorXd> B(b, n+1);
        VectorXd Res(n+1);
        Res =  A.jacobiSvd(ComputeThinU | ComputeThinV).solve(B);
        //std::cout << "\nResult via Jacobian SVD:\n" << Res << std::endl;

        double* vec = Res.data();
        for(int j = 0; j <= n; j++) {
                std::cout << vec[j] <<  ", ";
        }
          std::cout << "\n";


        //std::cout << "\nResult via QR decomposition:\n"
        //<< A.colPivHouseholderQr().solve(B) << std::endl;
        //std::cout << "\nResult via normal equations:\n"
        //<< (A.transpose() * A).ldlt().solve(A.transpose() * B) << std::endl;

/*
        for(int j = 0; j <= n; j++) {
                charges[j] = (double) Res[j];
        }

    for(int j = 0; j <= n; j++) {
        std::cout << charges[j] <<  "\n";
    }

    */
   


        //A.resize(0,0);

        //B.resize(0);
        free(Af);
        free(b);


}


void eemtest() {
  std::cout << "=>start test rdf\n";

  std::string pathName = "";//getenv("RDBASE");
  std::string sdfName =
      pathName + "examples/set00.sdf";

  RDKit::SDMolSupplier reader(sdfName, true, false);

  std::string line;
  int nDone = 0;
  while (!reader.atEnd()) {
    RDKit::ROMol *m = reader.next();
    TEST_ASSERT(m);
    std::string nm;
    int confId=-1;
    // 3D distance matrix
    double* dist3D = RDKit::MolOps::get3DDistanceMat(*m, confId, false, true);
    // compute the charges
    calculate_charges(m, dist3D);

    delete m;
    ++nDone;
  }

  BOOST_LOG(rdErrorLog) << "test on : " << nDone << " molecules done"
                        << std::endl;
}

int main(int argc, char *argv[]) {
  RDLog::InitLogs();
  eemtest();
}
