#include "Tree/DistTree.h"

#include "Sequence/SequenceTable.h"

#include "Models/Model.h"

#include "Sequence/SequenceTable.h"
#include "Util/randombox.h"
#include "PatternDesign/Singleton.h"


#include <assert.h>
#include <iomanip>

#define ITER_MAX 50
#define DIST_CONV .00001

DistTree::DistTree(){}

DistTree::DistTree( ParametersSet& ){
    pmodel = NULL;
    ptable = NULL;
    lookup = NULL;
}

DistTree::~DistTree(){}

void DistTree::initialisation( SequenceTable * ptable, Model * pmodel ){
    assert ( this->ptable == NULL );
    assert ( this->pmodel == NULL );
    this->ptable = ptable;
    this->pmodel = pmodel;
    
    if (ptable->getNumberCategories() != pmodel->getNumberSymbolCategory()){
        cerr << "ERROR: " << ptable->getNumberCategories()
             << " categories of symbols have been found in your sequence table"
             << " whereas the substitution model defined can deal with "
             << pmodel->getNumberSymbolCategory() << " categories" << endl;
        exit(EXIT_FAILURE);
    }
    
    //create lookup table
    assert( lookup == NULL );
    unsigned int numberCategories = ptable->getNumberCategories();
    lookup = new array3D<double> [numberCategories];
    diffLookup = new array3D<double> [numberCategories];
    diff2Lookup = new array3D<double> [numberCategories];
    partialLikelihood1 = new array2D<double> [numberCategories];
    partialLikelihood2 = new array2D<double> [numberCategories];
    lik = new vector<double>[numberCategories];
    diffLik = new vector<double>[numberCategories];
    diff2Lik = new vector<double>[numberCategories];
    for( unsigned int cat = 0; cat < numberCategories; ++cat ){
        unsigned int numberMixtCategories = pmodel->getNumberRatesCategories(cat);
        unsigned int numberStates = pmodel->getNumberStates(cat);
        unsigned int numberSites = ptable->getSequencesLength(cat)+ptable->getInvariantsLength(cat);
        lookup[cat].resize( numberMixtCategories, numberStates, numberStates );
        diffLookup[cat].resize( numberMixtCategories, numberStates, numberStates );
        diff2Lookup[cat].resize( numberMixtCategories, numberStates, numberStates );
        lik[cat].resize( numberSites );
        diffLik[cat].resize( numberSites );
        diff2Lik[cat].resize( numberSites );
        partialLikelihood1[cat].resize( numberSites, numberStates);
        partialLikelihood2[cat].resize( numberSites, numberStates);
    }
}

void DistTree::computePairwiseDist(){
    
    Singleton < randombox > & randBox = Singleton < randombox >::instance();
    
    unsigned int numberCategories = ptable->getNumberCategories();
    
    dist.resize(ptable->getNumberSpecies(),ptable->getNumberSpecies());
    for (unsigned int i = 1; i < ptable->getNumberSpecies(); ++i){
        loadLeafData( i, true );
        for (unsigned int j = 0; j < i; ++j){
            double branch = randBox.ran()*1.0;
            loadLeafData( j, false );
            bool stopCondition = false;
            /**
             * for each pair, we maximise log( P(S1->S2|t,M) ) wrt. t
             * log( P(S1->S2|t,M) ) = sum(log P(X1_i->X2_i|t,M)) for all sites i
             * d log( P(S1->S2|t,M) ) / dt = sum( P'(X1_i->X2_i|t,M) / P(X1_i->X2_i|t,M) )
             * d2 log( P(S1->S2|t,M) ) / dt2 = sum( [P''(X1_i->X2_i|t,M)*P(X1_i->X2_i|t,M) - P'(X1_i->X2_i|t,M)^2 ]/ P(X1_i->X2_i|t,M)^2 )
             */
            for (unsigned int iteration = 0; (iteration<ITER_MAX) && !stopCondition ;++iteration){
                fillLikArray( branch );
		        double diffLogLikelihood = 0.0;
                double diff2LogLikelihood = 0.0;
                for (unsigned int cat = 0; cat < numberCategories; ++cat){
                    unsigned int numberMixtCategories = pmodel->getNumberRatesCategories(cat);
                    unsigned int varSeqLength = ptable->getSequencesLength(cat);
                    unsigned int totalSeqLength = varSeqLength + ptable->getInvariantsLength(cat);
                    for (unsigned int mixtCat = 0; mixtCat < numberMixtCategories; ++mixtCat){
                         for (unsigned int site = 0;site < totalSeqLength; ++site){
                            double diffLikelihoodSite = diffLik[cat][site] / lik[cat][site];
                            assert(lik[cat][site]!=0.0);
                            assert(!isnan(diffLikelihoodSite));
                            assert(!isinf(diffLikelihoodSite));
                            double diff2LikelihoodSite = (diff2Lik[cat][site]/lik[cat][site])-(diffLikelihoodSite*diffLikelihoodSite);
                            assert(!isnan(diff2LikelihoodSite));
                            assert(!isinf(diff2LikelihoodSite));
                            if (site<varSeqLength){
                                diffLogLikelihood += diffLikelihoodSite;
                                diff2LogLikelihood += diff2LikelihoodSite;
                            }
                            else{
                                unsigned int numberInvariants = ptable->getInvariantBases( cat )[site-varSeqLength].second;
                                diffLogLikelihood += numberInvariants * diffLikelihoodSite;
                                diff2LogLikelihood += numberInvariants * diff2LikelihoodSite;
                            }
                        }
                    }
                }
                double newBranch;
                if(diff2LogLikelihood<0){
                    newBranch=branch-diffLogLikelihood/diff2LogLikelihood;
                    if (newBranch<0.0){
                        newBranch=branch/2.0;
                        cout << 'd' << flush;
                    }
                    else{
                        cout << '.' << flush;
                    }
                }
                else{
                    cout << 'd' << flush;
                    newBranch=branch/2.0;
                }
                stopCondition = fabs(newBranch-branch)<DIST_CONV;
                branch = newBranch;
            }
            cout << endl;
            cout << "MLdist(" << ptable->species[j] << "<->" << ptable->species[i] << ") = " << branch << endl;
	        dist(i,j) = branch;
	        dist(j,i) = branch;
        }
    }
}


void DistTree::fillLikArray( double branchLength ){

    unsigned int numberCategories = ptable->getNumberCategories();
    vector<double> mixtProb;
    for (unsigned int cat = 0; cat < numberCategories; ++cat){
        unsigned int numberMixtCategories = pmodel->getNumberRatesCategories(cat);
        mixtProb.clear();
        mixtProb.reserve(numberMixtCategories);
        for (unsigned int mixtCat = 0; mixtCat < numberMixtCategories; ++mixtCat){
            mixtProb.push_back( pmodel->getRateCategoryProbability( mixtCat, cat ) );
        }
        unsigned int numberStates = pmodel->getNumberStates(cat);
        for (unsigned int mixtCat = 0; mixtCat < numberMixtCategories; ++mixtCat){
            for (unsigned int state1 = 0; state1 < numberStates; ++state1){
                for (unsigned int state2 = 0; state2 <= state1; ++state2){
                    lookup[cat](mixtCat,state1,state2) =
                       pmodel->probability( state1, state2, branchLength, mixtCat, cat );
                    diffLookup[cat](mixtCat,state1,state2) =
                       pmodel->diffProbability( state1, state2, branchLength, mixtCat, cat );
                    diff2Lookup[cat](mixtCat,state1,state2) =
                       pmodel->secondDiffProbability( state1, state2, branchLength, mixtCat, cat );
                    if (state2!=state1){
                        lookup[cat](mixtCat,state2,state1) =
                           pmodel->probability( state2, state1, branchLength, mixtCat, cat );
                        diffLookup[cat](mixtCat,state2,state1) =
                           pmodel->diffProbability( state2, state1, branchLength, mixtCat, cat );
                        diff2Lookup[cat](mixtCat,state2,state1) =
                          pmodel->secondDiffProbability( state2, state1, branchLength, mixtCat, cat );
                    }
                }
            }
        }
        
        for (unsigned int site = 0;site < ptable->getSequencesLength(cat)+ptable->getInvariantsLength(cat); ++site){
            lik[cat][site] = 0.0;
            diffLik[cat][site] = 0.0;
            diff2Lik[cat][site] = 0.0;
            for (unsigned int state1=0; state1 < numberStates; ++state1){
                if (partialLikelihood1[cat](site,state1)){
                    for (unsigned int state2=0; state2 < numberStates; ++state2){
                        if (partialLikelihood2[cat](site,state2)){
                            double fact = partialLikelihood1[cat](site,state1) *
                                               partialLikelihood2[cat](site,state2);
                            for (unsigned int mixtCat = 0; mixtCat < numberMixtCategories; ++mixtCat){
                                double fact2 = mixtProb[mixtCat] * fact * pmodel->getFrequency(state1,mixtCat,cat);
                                lik[cat][site] += fact2 * lookup[cat](mixtCat,state1,state2);
                                diffLik[cat][site] += fact2 * diffLookup[cat](mixtCat,state1,state2);
                                diff2Lik[cat][site] += fact2 * diff2Lookup[cat](mixtCat,state1,state2);
                            }
                        }
                    }
                }
            }
            assert(!isnan(lik[cat][site]));
            assert(!isnan(diffLik[cat][site]));
            assert(!isnan(diff2Lik[cat][site]));
            assert(!isinf(lik[cat][site]));
            assert(!isinf(diffLik[cat][site]));
            assert(!isinf(diff2Lik[cat][site]));
/*
            if (!lik[cat][site]){
                cout << cat << '/' << site << '/' << branchLength << endl;
                for (unsigned int state1 = 0; state1 < numberStates; ++state1){
                    for (unsigned int state2 = 0; state2 <= state1; ++state2){
                        cout << lookup[cat](numberMixtCategories-1,state1,state2) << "   ";
                    }
                    cout << endl;
                }
                for (unsigned int state=0; state < numberStates; ++state){
                    cout << (partialLikelihood1[cat](site,state)) << "  ";
                }
                cout << endl;
                for (unsigned int state=0; state < numberStates; ++state){
                    cout << (partialLikelihood2[cat](site,state)) << "  ";
                }
                cout << endl;
                for (unsigned int mixtCat = 0; mixtCat < numberMixtCategories; ++mixtCat){
                    cout << pmodel->getRateCategoryProbability( mixtCat, cat ) << endl;
                    cout << pmodel->getFrequency(0,mixtCat,cat) << "  ";
                    cout << pmodel->getFrequency(1,mixtCat,cat) << "  ";
                    cout << pmodel->getFrequency(2,mixtCat,cat) << "  ";
                    cout << pmodel->getFrequency(3,mixtCat,cat) << "  " << endl;
                }
            }
*/
            assert(lik[cat][site]);
        }
    }
}

void DistTree::loadLeafData(unsigned int i, bool first){

    array2D<double>* partialLikelihood;
    partialLikelihood = first ? partialLikelihood1 : partialLikelihood2;
    
    int symbol;
    bool error = false;
    bool stateFound;
    
    for ( unsigned int cat = 0; cat < ptable->getNumberCategories(); ++cat ){
        const array2D<string> & sequences = ptable->getSequences( cat );
        const vector< pair< string, unsigned int > > & invariants =
            ptable->getInvariantBases( cat );
        unsigned int size = sequences.numberColumns() + invariants.size();
        const array2D<double> & table = pmodel->getEquivalencyTable(cat);
        for ( unsigned int site = 0; site < size; ++site ) {
            if (site<sequences.numberColumns()){        
                symbol = pmodel->getSymbolNumber( sequences(i,site), cat );
            }
            else{
                symbol = pmodel->getSymbolNumber( invariants[site-sequences.numberColumns()].first, cat );            
            }
            if ( symbol == (int)pmodel->getNumberSymbols( cat ) ){
                if (site<sequences.numberColumns()){
                    cerr << "unrecognized symbol in taxon " << ptable->species[i] << ": " << sequences(i,site) << endl;
                }
                else{
                    cerr << "unrecognized symbol in all sequences: " << invariants[site-sequences.numberColumns()].first << endl;
                }
                error=true;
            }
            else{
                stateFound = false;
                for ( unsigned int state = 0; state < pmodel->getNumberStates(cat); ++state ){
                    if (table( symbol, state )!=0.0){
                        partialLikelihood[cat]( site, state ) =
                            table( symbol, state );
                        stateFound = true;
                    }
                    else{
                        partialLikelihood[cat]( site, state ) = 0.0;
                    }
                }
                //a symbol can be recognized but not assigned to any state in a model (eg, codon-stop).
                //we deal with this specific case here
                if (!stateFound){
                    if (site<sequences.numberColumns()){
                        cerr << "forbidden symbol in taxon " << ptable->species[i] << ": " << pmodel->getSymbol(symbol,cat) << endl;
                    }
                    else{
                        cerr << "forbidden symbol in all sequences: " << pmodel->getSymbol(symbol,cat) << endl;
                    }
                    error=true;
                }
            }
            if(error){
                cerr << "position: ";
                vector < vector< unsigned int > > pos;
                if (site<ptable->getSequencesLength(cat)){
                    pos = ptable->retrieveInitialNucleotides( cat, site );
                }
                else{
                    pos = ptable->retrieveInitialNucleotides( cat, (int)sequences.numberColumns() - (int)site - 1 );
                }
                for ( vector < vector< unsigned int > >::iterator iter = pos.begin();
                      iter != pos.end(); ++iter ){
                    for ( vector< unsigned int >::iterator iter2 = (*iter).begin();
                          iter2 != (*iter).end(); ++iter2 ){
                         if (iter2 != (*iter).begin() ){
                             cerr << ", ";
                         }
                         cerr << *iter2 + 1;
                    }
                    cerr << "; ";
                }
                cerr << endl;
                exit(EXIT_FAILURE);
            }
        } // end site
    } // end model
}

void DistTree::printResults( ostream & outputStream, MatrixFormat frmt ){
    unsigned int prec = 5;
    unsigned int spacing = prec+3;
    
    outputStream.setf(ios::fixed);
    outputStream << setprecision(prec);
    for (unsigned int i = 0; i < ptable->getNumberSpecies(); ++i){
        outputStream << setw(17) << ptable->species[i] << "  ";
        for (unsigned int j = 0; j < ptable->getNumberSpecies(); ++j){
            if (i==j){
                if(frmt==UPPER){
                    outputStream << setw(spacing) << " ";
                }
                if(frmt==SQUARE){
                    outputStream << setw(spacing) << 0.0 << " ";
                }
                //do nothing if frmt==LOWER
            }
            else{
                if (i>j){
                    //i>j and (SQUARE or LOWER) matrix => print
                    if (frmt&0x01){
                        outputStream << setw(spacing) << dist(i,j) << " ";
                    }
                    else{
                        outputStream << setw(spacing) << " " << " ";
                    }
                }
                else{
                    //i<j and (SQUARE or UPPER) matrix => print
                    if ( frmt&0x02 ){
                        outputStream << setw(spacing) << dist(i,j) << " ";
                    }
                    //do nothing if frmt==LOWER
                }
            }
        }
        outputStream << endl;
    }
    outputStream << endl;
}

