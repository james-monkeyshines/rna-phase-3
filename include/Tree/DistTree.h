#ifndef DISTTREE_H
#define DISTTREE_H

#include "Tree/DistTree.h"

#include "Util/array2D.h"
#include "Util/array3D.h"
#include "Util/ParametersSet.h"

class SequenceTable;
class Model;

#include <vector>
#include <iostream>

using namespace std;

/** ****************************************************************************
 * The main class responsible for the computation of ML distance
 * In spite of its name, it is not really a tree.
 **************************************************************************** */

class DistTree{
public:

    typedef enum{
        LOWER = 0x01,
        UPPER = 0x02,
        SQUARE = 0x03                
    } MatrixFormat;
    
    /** ************************************************************************
     * DistTree
     * @semantics  constructor (empty)
     ************************************************************************ */
    DistTree();
    DistTree( ParametersSet& treeParameters );
    
    /** ************************************************************************
     * ~DistTree
     * @semantics  destructor
     ************************************************************************ */
    ~DistTree();

    /** ************************************************************************
     * initialisation
     * @input          the sequences and the model to use with this tree
     * @preconditions  no data or model already loaded
     * @semantics      store the pointers
     ************************************************************************ */
    void initialisation( SequenceTable * ptable, Model * pmodel );

    /** ************************************************************************
     * computePairwiseDist
     * @preconditions  data and model loaded
     * @semantics      launch the computation, optimize ML distance fill
     *                 dist vector. Newton-Raphson method
     *                 t_n+1 = t_n + L'(t)/L''(t)
     *                 L being the likelihood P(S1->S2|t,SubstutionModel)
     ************************************************************************ */
    void computePairwiseDist();
    
    /** ************************************************************************
     * printResults
     * @input          the stream where the results are to be sent
     * @frmt           the format for matrix output
     * @preconditions  computePairwiseDist called and successful
     * @semantics      print the matrix
     ************************************************************************ */
    void printResults( ostream & outputStream, MatrixFormat frmt );
    
protected:
    SequenceTable * ptable;
    Model * pmodel;


    /** ***********************************************************************
     * lastCalc
     * in lastCalc we store the site where the computation has been done for
     * the first and the last time, with a given identical state for all
     * descendants size(lastCalc[i]) == t_model->getNumberStates()
     * There is a vector for each model (size(lastCalc)==numberModelsUsed)
     ************************************************************************ */
    vector< vector<int> > lastCalc;
    
    /** ***********************************************************************
     * dist
     * an array to store the ML distances, the matrix is symmetric but space
     * is not an issue
     ************************************************************************ */
    array2D< double > dist;
    
    /** ***********************************************************************
     * loadLeafData
     * @input      first, if true load sequence i in partialLikelihood1,
     *             else load in partialLikelihood2
     * @semantics  fill a partialLikelihood array for the leaf i.
     *             each symbol is matched to one or more states
     *             thanks to the model equivalency table
     ************************************************************************ */
    void loadLeafData(unsigned int i, bool first);
    
    /** ***********************************************************************
     * store the equivalency tables
     ************************************************************************ */
    array2D<double>* partialLikelihood1;
    array2D<double>* partialLikelihood2;


    /** ***********************************************************************
     * fillLikArray
     * @input       double branchLength 
     * @semantics   once two sequences are loaded, fill lik/diffLik/diff2Lik
     *              for the given branchLength
     ************************************************************************ */
    void fillLikArray( double branchLength );
        
    
    /** ***********************************************************************
     * lik/diffLik/diff2Lik
     * store P(X_i->Y_i|t,M), P' and P'' for each site of a pair of sequence
     *********************************************************************** */
    vector<double>* lik;
    vector<double>* diffLik;
    vector<double>* diff2Lik;
    
    /** ***********************************************************************
     * lookup/diffLookup
     * memory space to store p(X->Y|t,r) and dp(X->Y|t,r)/dt
     * for each model declared, dim(lookup[modelID]) =
     *             numberMixtureCategories * nb_states * nb_states 
     ************************************************************************ */
    array3D<double> * lookup;
    array3D<double> * diffLookup;
    array3D<double> * diff2Lookup;            

};
#endif //DISTTREE_H
