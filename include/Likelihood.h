#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include "Phase.h"

#include <vector>

using namespace std;


class SequenceTable;
class LikelihoodTree;
class Model;

class Likelihood: public Phase{

public:
      
    void initObject();
  
    /** ***********************************************************************
     * getStateFrequencies
     * @input    speciesMask, a set of species
     * @input    symbolCategory, the model we are interested in
     * @input    rateCategoriesMask, a set of category
     * @return   double, BPP of the categories given in rateCategoriesMask
     *           (the prior is 1/k for each gamma category but the posterior
     *           can be different)
     * @output   stateFreq, composition observed for the species in speciesMask
     *           at the given rate categories. The contribution of each site
     *           to these estimates is weighted by the BPP of being in the
     *           categories
     ************************************************************************ */
    double getStateFrequencies( const vector< bool > & speciesMask,
                                unsigned int symbolCategory,
                                const vector< bool > & rateCategoriesMask,
                                vector<double>& stateFreq);
    
    /** ***********************************************************************
     * getAncestralSequence
     * @output probReconstruction, BPP for the frequencies at a site (for
     *         each node)
     * @output ancestral sequence at each node (marginal reconstruction)
     ************************************************************************ */
    void getAncestralSequences( vector< vector< vector<double> > > & probReconstruction,
                                vector< pair<InferenceNode*, string> > & ancestralSequences );
    
    /** ***********************************************************************
     * getRateCategory
     * @output probReconstruction, BPP for the rates at a site
     * @output vector< unsigned >, the MAP category chosen for each site
     ************************************************************************ */
    void getRateCategory( vector< vector<double> > & probReconstruction, vector<unsigned int> & rate);
    
    /** ***********************************************************************
     * getMaxProb
     * @input      probMatrix, a set of probability for each site in the 
     *             processed sequences
     * @output     vector< unsigned >, the MAP value chosen from probMatrix 
     *             for each site
     * @semantics  used by getAncestralSequence and getRateCategory since
     *             the two functions are doing similar things
     ************************************************************************ */
    void getMaxProb( const array2D<double> & probMatrix, vector<unsigned int> & maxProb);
    
    /** ***********************************************************************
     * ~Likelihood
     * @semantics  the destructor
     ************************************************************************ */
    virtual ~Likelihood();

        
    int run( int argc, char** argv);
    
protected:
        
    SequenceTable* table;
    Model* model;    
    LikelihoodTree* tree;
        
    /**
     * ancestralStates[n].second[symbolCategory](site,state) is the
     * probability p(state|topology, branch lengths, substitution model parameters)
     * at internal node ancestralState[n].first for the given nucleotide position
     * (symbolCategory,site) site is in the range [0 .. nb non-invariant sites
     * + invariant.size(), where invariant.size() is the number of invariant
     * patterns (when they are lumped, cf. the sequence provided by SequenceTable
     */
    list< pair<InferenceNode*, vector< array2D<double> > > > ancestralStates;


    /* rateCat[symbolCategory](site,rate) is the
     * probability p(rate|topology, branch lengths, substitution model parameters)
     * for the given nucleotide position (symbolCategory,site)
     * site is in the range [0 .. nb non-invariant sites + invariant.size()
     * invariant.size() is the number of invariant size when they are lumped */
    vector< array2D<double> > rateCat;
    
    /* siteSpecificLik[symbolCategory](site) is the
     * probability p(data|topology, branch lengths, substitution model parameters)
     * for the given nucleotide position (symbolCategory,site)
     * site is in the range [0 .. nb_non-invariant_sites + nb_invariant_sites]
     * this is the initial length for each category.
     * invariant sites are put at the end and site-specific loglik are repeated
     * for them.
     */
    vector< vector<double> > siteSpecificLik;
    
    unsigned int initialSequenceLength;
    
    /* a set of parameters for each symbol category */
    unsigned int numberCategories;
    vector< unsigned int > length;
    vector< unsigned int > nonInvLength;
    vector< unsigned int > invLength;
    vector< unsigned int > numberStates;
    vector< unsigned int > numberRates;
   
};


#endif //LIKELIHOOD_H
