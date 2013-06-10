#include "Tree/InferenceNode.h"

#include "Tree/InferenceTree.h"
#include "Tree/TreeMap.h"

#include "Models/Model.h"

#include "PatternDesign/Singleton.h"
#include "Util/randombox.h"

#include <assert.h>

#include "Tree/Cluster.h"
#include "Tree/ClustersTree.h"
#include "Tree/ClustersTreeNode.h"

InferenceNode::InferenceNode():BasicNode(){
    lookup.resize(NUMBER_LOOKUP_TABLE);
};

InferenceNode::InferenceNode( InferenceTree* ptree ):
BasicNode(){
    //initialise the inference node like the basic node
    cnumber = -1;
    this->ptree = ptree;
    //launch the initialisation if the table was loaded
    if ( ptree->getTable() ){
        assert( ptree->getModel() );
        initialisation( false );
    }    
    lookup.resize(NUMBER_LOOKUP_TABLE);
}

InferenceNode::InferenceNode( const string& label, InferenceTree* ptree ):
BasicNode( label ) {
    //initialise the inference node like the basic node
    cnumber = -1;
    this->ptree = ptree;
    //launch the initialisation if the table was loaded
    if ( ptree->getTable() ){
        assert( ptree->getModel() );
        initialisation( true );
    }
    lookup.resize(NUMBER_LOOKUP_TABLE);
}

InferenceNode::InferenceNode( const TreeMap& treeMap, InferenceTree* ptree ):
BasicNode(){
    //if this node is an internal node
    if ( treeMap.getNumberChildren() ){
        cnumber = -1;
        this->ptree = ptree;
        //create the children
        for ( unsigned int i = 0; i < treeMap.getNumberChildren() ; ++i ) {
            addChild( new InferenceNode(treeMap.getChildMap(i), ptree),
                      treeMap.getDistance(i) );
        }
        updateNodesCount();
        //launch the initialisation if the table/model was loaded
        if ( ptree->getTable() ){
            assert( ptree->getModel() );
            initialisation( false );
        }
    }
    //if this node is a leaf
    else{
        //the basic node constructor used was for an internal node
        //correct this
        nbInternalNodes = 0;
        nbLeaves = 1;
        cnumber = -1;
        this->label = treeMap.getLabel();
        this->ptree = ptree;
        //launch the initialisation if the table/model was loaded
        if ( ptree->getTable() ){
            assert( ptree->getModel() );
            initialisation( true );
        }
    }
    lookup.resize(NUMBER_LOOKUP_TABLE);
}

InferenceNode::~InferenceNode(){
}  


void InferenceNode::initialisation( bool leafNode ) {

    identical.resize( ptree->getTable()->getNumberCategories() );
    lastCalc.resize( ptree->getTable()->getNumberCategories() );
    partialLikelihood.resize( ptree->getTable()->getNumberCategories() );
    partialLikelihoodSave.resize( ptree->getTable()->getNumberCategories() );
    

    pnodeModel = ptree->getModel();

    if ( leafNode ) {
        int sequenceId = ptree->getTable()->findSequence( label );
        if ( sequenceId == -1 ) {
            cerr << "Unable to find sequence for \"" << label << "\"" << endl;
            exit(EXIT_FAILURE);
        }
        cnumber = sequenceId + 1;
    }
    else {
        partialLikelihoodWork.resize( ptree->getTable()->getNumberCategories() );
        cnumber = -1;
    }
    
    for( unsigned int cat = 0; cat < ptree->getTable()->getNumberCategories();
         ++cat ){
        const array2D<string> & sequences =
                ptree->getTable()->getSequences( cat );
        const vector< pair< string, unsigned int > > & invariants =
            ptree->getTable()->getInvariantBases( cat );
        int size = sequences.numberColumns() + invariants.size();
        identical[cat].resize(size);
        lastCalc[cat].resize(pnodeModel->getNumberSymbols(cat));
    
        if ( leafNode ) {
            //partialLikelihood for nodes are constant (no work space)
            partialLikelihoodSave[cat] = new array3D< double >;
            partialLikelihoodSave[cat]->resize(size,
                                       pnodeModel->getNumberStates(cat),
                                       1 );
            partialLikelihood[cat] = partialLikelihoodSave[cat];
        }
        //if !leafNode
        else {
            partialLikelihoodSave[cat] = new array3D< double >;
            partialLikelihoodWork[cat] = new array3D< double >;
            partialLikelihoodSave[cat]->resize(size,
                            pnodeModel->getNumberStates(cat),
                            pnodeModel->getNumberRatesCategories(cat) );
            partialLikelihoodWork[cat]->resize(size,
                            pnodeModel->getNumberStates(cat),
                            pnodeModel->getNumberRatesCategories(cat) );
            //at the beginning, computation are invalid
            partialLikelihood[cat] = partialLikelihoodWork[cat];
        }
    }
    
    if (leafNode){
        loafLeaf();
    }
    else{
        //if children are already plugged compute the identical
        if (getNumberChildren()){
            computeIdentical();
        }
    }
}


unsigned int InferenceNode::lookupResize(){
    unsigned int multi = 0;
    for( list<BasicNode*>::iterator childrenIter = getChildrenList().begin();
                 childrenIter != getChildrenList().end(); ++childrenIter ){
        unsigned int childMulti = 0;
        if (!(*childrenIter)->isLeaf()){
            childMulti = ((InferenceNode*)(*childrenIter))->lookupResize();
        }
        multi = MAX(multi, childMulti);
    }
    lookup.resize(getNumberChildren());
    multi = MAX(multi, getNumberChildren());
    return multi;
}

InferenceNode& InferenceNode::quickCopy( const InferenceNode& src,
                                         bool withPartialLikelihood){
    //call BasicNode::operator=
    BasicNode::operator=(src);
    //both node must be initialised with same model/table or not be initialised
    assert( src.ptree->getTable() == ptree->getTable() );
    assert( src.ptree->getModel() == ptree->getModel() );
    //copy identical and partialLikelihood if necessary.
    if ( ( ptree->getTable() != NULL ) || ( ptree->getModel() != NULL ) ) {
        identical = src.identical;
        pnodeModel = src.pnodeModel;
        for( unsigned int cat = 0;
             cat < ptree->getTable()->getNumberCategories(); ++cat ){
            if ( withPartialLikelihood ){
                //copy the saved computation whatever the status
                (*partialLikelihoodSave[cat]) = *(src.partialLikelihoodSave[cat]);
                //let partialLikelihood point where it should
                if (src.partialLikelihood[cat] == src.partialLikelihoodSave[cat]){
                    partialLikelihood[cat] == partialLikelihoodSave[cat];
                }
                else{
                    partialLikelihood[cat] == partialLikelihoodWork[cat];
                }
            }
            else{
                partialLikelihood[cat] == partialLikelihoodWork[cat];
            }
        }
    }   //end if ( ptree->ptable || ptree->pmodel )
    return (*this);
}


void InferenceNode::recursiveCopy( InferenceNode * srcNode,
                                   InferenceNode * parent ) {
    if(!parent){
        pdistance = -1.0;
    }
    //the node should not have any children yet
    assert(getNumberChildren() == 0);
    
    quickCopy(*srcNode, true);
    this->parent = parent;
    
    for ( list<BasicNode*>::const_iterator iter = srcNode->getChildrenList().begin();
          iter != srcNode->getChildrenList().end(); ++iter ){
        InferenceNode* newNode;
        if (!(*iter)->isLeaf()){
            //create the new internal node and repeat the recursive process
            newNode = new InferenceNode(ptree);
            newNode->recursiveCopy( (InferenceNode*)(*iter), this );
        }
        else{
            //if it is a leaf just create it, the step contruct+init will do everything needed
            newNode = new InferenceNode( (*iter)->getLabel(), ptree );
        }
        //addChild deals with pdistance and the nbInternalNodes an nbLeaves counter
        addChild( newNode, (*iter)->getParentDistance() );
    }
}



InferenceNode* InferenceNode::getChild(unsigned int childId){
    assert( childId < getNumberChildren() );
    list< BasicNode* >::iterator childIter = getChildrenList().begin();
    while (childId>0){
        ++childIter;
        --childId;
    }
    return (InferenceNode*)(*childIter);
}

string InferenceNode::toStringNumbered( bool topologyOnly ) const{
    char d[30];
    //if leaf node
    if ( isLeaf() ) {
        sprintf(d,"%d",cnumber);
        return d;
    }
    
    //else
    string tree= "(";
    for ( list< BasicNode* >::const_iterator iter = getChildrenList().begin();
          iter != getChildrenList().end(); ++iter ){
        //comma separator
        if ( iter != getChildrenList().begin() ){
            tree += ",";
        }
        tree += ((InferenceNode*)(*iter))->toStringNumbered( topologyOnly );
        if( !topologyOnly ){
            sprintf( d, ":%.8f", (*iter)->getParentDistance() );
            tree += d;
        }
    }
    tree += ")";
    return tree;
}



void InferenceNode::likelihood() {
    
    unsigned int numberChildren = getNumberChildren();
    
    //no likelihood for a leaf node
    assert( numberChildren );

    if ( (lookup.size()<numberChildren) || (ptree->lookup.size()<numberChildren) ){
        cerr << "Sorry, maximum " << lookup.size() << " children allowed" << endl;
        exit(EXIT_FAILURE);
    }
    
    //likelihood for a parent with only 2 leaf children
    if ( (numberChildren == 2) &&
         (*getChildrenList().begin())->isLeaf() &&
         (*(--getChildrenList().end()))->isLeaf() ){
        optimizedLikelihoodLeafParent();
        return;
    }
        
    
    SequenceTable* pnodeTable = ptree->getTable();
#ifdef TESTCALC
    double save = (*(partialLikelihood[0]))(0,0,0);
#endif  
    
    // Initialize lookup table
    for ( unsigned int cat = 0; cat < pnodeTable->getNumberCategories(); ++cat){
#ifndef TESTCALC
        if (partialLikelihood[cat]==partialLikelihoodWork[cat]){
#endif  
            int numberStates = pnodeModel->getNumberStates(cat);
            int numberRateCategories =
                    pnodeModel->getNumberRatesCategories(cat);
            for ( int initialState = 0; initialState < numberStates; ++initialState ) {
                for ( int finalState = 0; finalState < numberStates; ++finalState ) {
                   for ( int rate = 0; rate < numberRateCategories; ++rate ) {
                        list<BasicNode*>::const_iterator iter = getChildrenList().begin();
                        for ( unsigned int childId = 0; childId<numberChildren;
                              ++childId){
                            ptree->lookup[childId][cat](rate, initialState, finalState ) =
                                ((InferenceNode*)(*iter))->getModel()->
                                    probability( initialState, finalState,
                                        (*iter)->getParentDistance(), rate, cat );
                            ++iter;
                        }
                    }
                }
            }
#ifndef TESTCALC
        }
#endif  
        //initialise lastCalc
        fill( lastCalc[cat].begin(), lastCalc[cat].end(), -1);
    }

    for ( unsigned int cat = 0; cat < pnodeTable->getNumberCategories(); ++cat){
#ifndef TESTCALC
        if (partialLikelihood[cat]==partialLikelihoodWork[cat]){
#endif  
            int sequenceLength = pnodeTable->getSequencesLength(cat);
            int invariantsLength = pnodeTable->getInvariantsLength(cat);
            int numberStates = pnodeModel->getNumberStates(cat);
            int numberRates = pnodeModel->getNumberRatesCategories(cat);            
            for ( unsigned int childId = 0; childId<numberChildren; ++childId){
                lookup[childId] = &(ptree->lookup[childId][cat]);
            }
            for ( int site = 0; site < sequenceLength+invariantsLength; ++site ) {
                 if ( ( identical[cat][site] == -1 ) ||
                     ( lastCalc[cat][identical[cat][site]] == -1 ) ){
                    for ( int rate = 0; rate < numberRates; ++rate ) {
                        for ( int initialState = 0; initialState < numberStates;
                              ++initialState ) {
                            (*partialLikelihood[cat])(site, initialState, rate) = 1.0;
                            unsigned int childId = 0;
                            for ( list<BasicNode*>::const_iterator iter
                                        = getChildrenList().begin();
                                  iter != getChildrenList().end(); ++iter ){
                                //if a children is a leaf rates are not in the
                                //partial Likelyhood array and endRate = 0.
                                // otherwise dim partialLikelyhood[site] =
                                // (numberStates, numberRates)
                                // and "endRate = startRate"
                                int endRate = (*iter)->isLeaf() ? 0 : rate;
                                double temp = 0.0;
                                for ( int finalState = 0; finalState < numberStates;
                                      ++finalState ) {
                                    temp += (*lookup[childId])
                                                ( rate, initialState, finalState ) *
                                (*((InferenceNode*)(*iter))->partialLikelihood[cat])
                                                (site, finalState, endRate);
                                } // for each final state
                                (*partialLikelihood[cat])(site, initialState, rate) *= temp;
                                ++childId;
                            } // for each child
                        } // for each initial state
                    } // for each rate
                    // if identical[modelId][site] save this site ID for
                    // future computation
                    if ( identical[cat][site] != -1 ) {
                        lastCalc[cat][identical[cat][site]] = site;
                    }
                }
                // if identical[site]!=-1 && already calculated
                else {
                    memcpy( &((*partialLikelihood[cat])(site,0,0)),
     &((*partialLikelihood[cat])(lastCalc[cat][identical[cat][site]],0,0)),
                      sizeof( double ) * numberStates * numberRates );
#ifdef DEBUG4
                    for ( int i = 0; i< numberStates; ++i ){
                        for ( int j = 0; j< numberRates; ++j ){
                            assert((*partialLikelihood[cat])(site, i, j) ==
    (*partialLikelihood[cat])(lastCalc[cat][identical[cat][site]], i, j) );
                        }
                    }
#endif
                } //end if "already calculated"
            } // for each site
#ifndef TESTCALC
        } //if category was invalidated
#endif
    } // for each category
    
#ifdef TESTCALC
    if (partialLikelihood[0]==partialLikelihoodSave[0])
        assert(save == (*(partialLikelihood[0]))(0,0,0));
#endif
      
}




void InferenceNode::computeIdentical() {
    unsigned int cat;

    SequenceTable* pnodeTable = ptree->getTable();

    if ( !isLeaf() ) {
#ifdef TESTIDENT
        for ( cat=0; cat<pnodeTable->getNumberCategories(); ++cat){
            unsigned int length =
                        pnodeTable->getSequencesLength(cat);
            for ( unsigned int i = 0; i < length; ++i ) {
                identical[cat][i] = -1;
            }
        }
        return;
#endif
#ifndef TESTIDENT    
        for ( cat=0; cat < pnodeTable->getNumberCategories(); ++cat){
            unsigned int length = pnodeTable->getSequencesLength(cat);
            list<BasicNode*>::const_iterator iter = getChildrenList().begin();
            //initialise identical with identical in the first child
            identical[cat] = ((InferenceNode*)(*iter))->identical[cat];
            
            
            ++iter;
            list<BasicNode*>::const_iterator iterSeek;
            //for "global" non-invariant sites, check whether the symbol
            //was invariant for all the leaves descendant from this node
            for ( unsigned int i = 0; i < length; ++i ) {
                iterSeek =  iter;
                while ( (identical[cat][i] != -1) &&
                        (iterSeek != getChildrenList().end()) ){
                    if ( identical[cat][i] !=
                         ((InferenceNode*)(*iterSeek))->identical[cat][i] ){
                        identical[cat][i] = -1;
                    }
                    ++iterSeek;
                }
            }
        }
#endif
    }
    else {
        // identical for a leaf is invariant, error if trying to calculate
        // should have been initialised when the sequence is loaded
        // (initialisation method)
        //however some node might be temporarily leaf because they
        //have no children yet
        assert( label=="" );
    }
}


#ifdef DEBUG1
void InferenceNode::checkIdentical() {

    SequenceTable* pnodeTable = ptree->getTable();

    for (unsigned int cat=0; cat < pnodeTable->getNumberCategories(); ++cat){
        const array2D<string> & sequences =
                pnodeTable->getSequences( cat );
        const vector< pair< string, unsigned int > > & invariants =
            pnodeTable->getInvariantBases( cat );
        unsigned int size = sequences.numberColumns() + invariants.size();
        if ( !isLeaf() ) {
            for ( unsigned int i = 0; i < size; ++i ) {
                if (identical[cat][i] != -1){
                    for ( list<BasicNode*>::const_iterator iter = getChildrenList().begin();
                          iter!=getChildrenList().end(); ++iter ) {
                        assert(((InferenceNode*)(*iter))->identical[cat][i]
                                == identical[cat][i]);
                    }
                }
                else{
                    list<BasicNode*>::const_iterator iter = getChildrenList().begin();
                    int check = ((InferenceNode*)(*iter))->identical[cat][i];
                    while (check != -1){
                        ++iter;
                        assert(iter!=getChildrenList().end());
                        if ( ((InferenceNode*)(*iter))->identical[cat][i]
                              != check ){
                            check = -1;
                        }
                    }
                }
            }
        }
        else {
            int sequenceId = ptree->getTable()->findSequence( label );
            for ( unsigned int site = 0; site < sequences.numberColumns();
                  ++site ) {
                assert(identical[cat][site] ==
                    pnodeModel->getSymbolNumber( sequences(sequenceId,site), cat ) );
            }
            for ( unsigned int inv = 0; inv< invariants.size(); ++inv ) {
                assert( identical[cat][inv+sequences.numberColumns()] ==
                    pnodeModel->getSymbolNumber( invariants[inv].first, cat ) );
            }
        }
    }
}
#endif


void InferenceNode::updateIdenticalRecursively(InferenceNode* stoppingPoint){
    InferenceNode* tempNode = this;
    while( tempNode != stoppingPoint ){
        assert( tempNode );
        tempNode->computeIdentical();
        tempNode = (InferenceNode*)tempNode->getParent();
    }
}

void InferenceNode::invalidateRecursively(int cat, InferenceNode* stoppingPoint){
    InferenceNode* tempNode = this;
    while( tempNode != stoppingPoint ){
#ifdef DEBUG1
        assert( tempNode );
#endif
        //invalidate all
        if (cat == -1){
            tempNode->partialLikelihood = tempNode->partialLikelihoodWork;
        }
        //invalidate just one category
        else{
            tempNode->partialLikelihood[cat] = tempNode->partialLikelihoodWork[cat];
        }
        tempNode = (InferenceNode*)tempNode->getParent();
    }
}

void InferenceNode::retrieveRecursively(int cat, InferenceNode* stoppingPoint){
    InferenceNode* tempNode = this;
    while( tempNode != stoppingPoint ){
#ifdef DEBUG1
        assert( tempNode );
#endif
        //retrieve all
        if (cat == -1){
            tempNode->partialLikelihood = tempNode->partialLikelihoodSave;
        }
        //retrieve just one category
        else{
            tempNode->partialLikelihood[cat] = tempNode->partialLikelihoodSave[cat];
        }
        tempNode = (InferenceNode*)tempNode->getParent();
    }
}

void InferenceNode::saveRecursively(int cat, InferenceNode* stoppingPoint){
    InferenceNode* tempNode = this;
    while( tempNode != stoppingPoint ){
#ifdef DEBUG1
        assert( tempNode );
#endif
        //save all
        if (cat == -1){
            for ( unsigned int cat = 0;
                  cat < ptree->getTable()->getNumberCategories(); ++cat ){
                //if work is used, swap the array with save
                if (tempNode->partialLikelihoodSave[cat] != tempNode->partialLikelihood[cat]){
                    tempNode->partialLikelihoodWork[cat] = tempNode->partialLikelihoodSave[cat];
                    tempNode->partialLikelihoodSave[cat] = tempNode->partialLikelihood[cat];
                }
            }
        }
        //retrieve just one category
        else{
            if (tempNode->partialLikelihoodSave[cat] != tempNode->partialLikelihood[cat]){
                tempNode->partialLikelihoodWork[cat] = tempNode->partialLikelihoodSave[cat];
                tempNode->partialLikelihoodSave[cat] = tempNode->partialLikelihood[cat];
            }
        }
        tempNode = (InferenceNode*)tempNode->getParent();
    }
}


void InferenceNode::loafLeaf() {

    int symbol;
    
    bool error=false;

    //identical is used for standard cases 1 symbol -> 1 state
    for( unsigned int cat = 0; cat < ptree->getTable()->getNumberCategories();
         ++cat ){
        const array2D<string> & sequences =
                ptree->getTable()->getSequences( cat );
        const vector< pair< string, unsigned int > > & invariants =
            ptree->getTable()->getInvariantBases( cat );
        unsigned int size = sequences.numberColumns() + invariants.size();
        const array2D<double> & table =
                pnodeModel->getEquivalencyTable(cat);
        for ( unsigned int site = 0; site < size; ++site ) {
            if (site<sequences.numberColumns()){        
                symbol = pnodeModel->getSymbolNumber( sequences(cnumber-1,site), cat );
            }
            else{
                symbol = pnodeModel->getSymbolNumber( invariants[site-sequences.numberColumns()].first, cat );
            }
            if ( symbol == (int)pnodeModel->getNumberSymbols( cat ) ){
                if (site<sequences.numberColumns()){
                    cerr << "unrecognized symbol in taxon " << label << ": " << sequences(cnumber-1,site) << endl;
                }
                else{
                    cerr << "unrecognized symbol in all sequences: " << invariants[site-sequences.numberColumns()].first << endl;
                }
                error=true;
            }
            else{
                identical[cat][site] = -2;
                for ( unsigned int state = 0; state < pnodeModel->getNumberStates(cat); ++state ) {
                    if (table( symbol, state )!=0.0){
                        if ( (identical[cat][site] == -2) && (table( symbol, state ) == 1.0) ){
                            identical[cat][site] = state;
                        }
                        else{
                            identical[cat][site] = -1;
                        }
                        (*partialLikelihoodSave[cat])( site, state, 0 ) = table( symbol, state );
                    }
                    else{
                        (*partialLikelihoodSave[cat])( site, state, 0 ) = 0.0;
                    }
                }
                if (identical[cat][site] == -2){
                    //a symbol can be recognized but not assigned to any state in a model.
                    //we deal with this specific case here
                    if (site<sequences.numberColumns()){
                        cerr << "forbidden symbol in taxon " << label << ": " << pnodeModel->getSymbol(symbol,cat) << endl;
                    }
                    else{
                        cerr << "forbidden symbol in all sequences: " << pnodeModel->getSymbol(symbol,cat) << endl;
                    }
                    error=true;
                }
            }
            //print position(s) in case of error
            if(error){
                cerr << "position: ";
                vector < vector< unsigned int > > pos;
                if (site<sequences.numberColumns()){
                    pos = ptree->getTable()->retrieveInitialNucleotides( cat, site );
                }
                else{
                    pos = ptree->getTable()->retrieveInitialNucleotides( cat, (int)sequences.numberColumns() - (int)site - 1);
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
        }
    }
}


void InferenceNode::optimizedLikelihoodLeafParent() {

    int left, right;
    double temp1, temp2;
    
    SequenceTable* pnodeTable = ptree->getTable();
    
    unsigned int numberCategories = pnodeTable->getNumberCategories();

    int* fastLookup = NULL;

    array3D < double > * look1, * look2;
    const array3D < double > * partialLikelihood1, * partialLikelihood2;

#ifdef TESTCALC
    double save = (*(partialLikelihood[0]))(0,0,0);
#endif  
    
    //initialise local variable for quick access to the children
    list<BasicNode*>::const_iterator iter = getChildrenList().begin();
    const InferenceNode* child1 = (InferenceNode*)*iter;
    double distance1 = child1->getParentDistance();
    const InferenceNode* child2 = (InferenceNode*)*(++iter);
    double distance2 = child2->getParentDistance();
        
    // Initialize lookup and fastLookup table
    for ( unsigned int cat = 0; cat < numberCategories; ++cat ) {
#ifndef TESTCALC
        if (partialLikelihood[cat]==partialLikelihoodWork[cat]){
#endif
            unsigned int numberRates = pnodeModel->getNumberRatesCategories( cat );
            unsigned int numberStates = pnodeModel->getNumberStates( cat );
            unsigned int numberStatesSqr = numberStates*numberStates;
            look1 = &(ptree->lookup[0][cat]);
            look2 = &(ptree->lookup[1][cat]);
            partialLikelihood1 = child1->partialLikelihood[cat];
            partialLikelihood2 = child2->partialLikelihood[cat];
            for ( unsigned int rate = 0; rate < numberRates; ++rate ) {
                for ( unsigned int initialState = 0; initialState < numberStates; ++initialState ) {
                    for ( unsigned int finalState = 0; finalState < numberStates; ++finalState ) {
                        (*look1)( rate, initialState, finalState ) =
                            child1->pnodeModel->probability( initialState, finalState,
                                                     distance1, rate, cat );
                        (*look2)( rate, initialState, finalState ) =
                            child2->pnodeModel->probability( initialState, finalState,
                                                     distance2, rate, cat );
                    }
                }
            }
            //fastLookuptable
            fastLookup = new int[numberStatesSqr];
            if (!fastLookup){
                cerr << "memory allocation error...sorry" << endl;
            }
            int * end = fastLookup + numberStatesSqr;
            for ( int *i = fastLookup; i != end; ++i ) {
                *i = -1;
            }
            int sequencesLength = pnodeTable->getSequencesLength(cat);
            int invariantsLength = pnodeTable->getInvariantsLength(cat);
                
            for ( int site = 0; site < sequencesLength + invariantsLength; ++site ){
                left = child1->identical[cat][site];
                right = child2->identical[cat][site];
                if ( (left==-1) || (right==-1) || fastLookup[(left*numberStates)+right] == -1 ) {
                    for ( unsigned int rate = 0; rate < numberRates; ++rate ) {
                        for ( unsigned int initialState = 0; initialState < numberStates;
                              ++initialState ) {
                            temp1 = 0.0;
                            temp2 = 0.0;
                            for ( unsigned int finalState = 0; finalState < numberStates;
                                  ++finalState ) {
                                temp1 += (( * look1 )(rate, initialState, finalState) *
                                     (*partialLikelihood1)(site, finalState, 0) );
                                temp2 += (( * look2 )(rate, initialState, finalState)*
                                     (*partialLikelihood2)(site, finalState, 0) );
                            }
                            (*partialLikelihood[cat])( site, initialState, rate ) =
                               temp1 * temp2;
                        } // for each initial state
                    } // for each rate
                    if ( (left!=-1) && (right!=-1) ){
                        fastLookup[(left*numberStates)+right] = site;
                    }
                }
                else {
                    memcpy(&((*partialLikelihood[cat])(site,0,0)),
                       &((*partialLikelihood[cat])(fastLookup[(left*numberStates)+right],0,0)),
                           sizeof(double)*numberRates*numberStates);
//#ifdef DEBUG4
                    for ( unsigned int i = 0; i < numberStates; ++i ){
                        for ( unsigned int j = 0; j < numberRates; ++j ){
                            assert((*partialLikelihood[cat])(site,i,j) ==
                        (*partialLikelihood[cat])(fastLookup[(left*numberStates)+right],i,j));
                        }
                    }
//#endif
                }   //if in fastlookup
            } // for each site
            delete [] fastLookup;
            fastLookup = NULL;
        } // for each model
#ifndef TESTCALC
    } // if work spspace
#endif
#ifdef TESTCALC
    if (partialLikelihood[0]==partialLikelihoodSave[0])
        assert(save == (*(partialLikelihood[0]))(0,0,0));
#endif
}


Cluster* InferenceNode::assignClustersRec( ClustersTree* clustersTree ){
    assert(ptree->getTable());
    assert(clustersTree);

    Cluster* cluster = new Cluster(ptree->getTable()->getNumberSpecies());
    for (list<BasicNode*>::const_iterator iter = getChildrenList().begin();
         iter != getChildrenList().end(); ++iter){
        if ((*iter)->isLeaf()){
            cluster->addSpecies((*iter)->cnumber-1);
        }
        else{
            Cluster* clusterChild = ((InferenceNode*)(*iter))->assignClustersRec(clustersTree);
            cluster->join(*clusterChild);
            delete clusterChild;
        }
    }
    for ( ClustersTree::const_iterator iter = clustersTree->begin();
          iter != clustersTree->end(); ++iter ){
        assert((*iter).second.first);
        if ( (*iter).second.first->matches(*cluster) ){
            (*iter).second.first->setAncestor(this);
        }
        if ( !(*iter).second.first->compatible(*cluster) ){
            cerr << "The following subtree is not compatible with the list of clusters provided" << endl;
            cerr << toString() << endl;
            cerr << "Cluster " << (*iter).second.first->getName() << '('
                 << (*iter).second.first->getKey() << ')' << endl;
            exit(EXIT_FAILURE);
        }
    }
    return cluster;
}

Cluster* InferenceNode::getCluster() const{
    assert(ptree->getTable());
    Cluster* ret = new Cluster(ptree->getTable()->getNumberSpecies());
    getClusterRec( ret );
    return ret;
}

void InferenceNode::getClusterRec( Cluster* cluster ) const{
    for (list<BasicNode*>::const_iterator iter = getChildrenList().begin();
         iter != getChildrenList().end(); ++iter){
        if ((*iter)->isLeaf()){
            cluster->addSpecies((*iter)->cnumber-1);
        }
        else{
            ((InferenceNode*)(*iter))->getClusterRec( cluster );
        }
    }
}

double InferenceNode::findClusterBranchesRec( vector< pair<InferenceNode*,double> > & branchVector,
                                              ClustersTreeNode* clustersTreeNode ){
    double sum = getParentDistance();
    branchVector.push_back( pair<InferenceNode*,double>(this,sum) );
    //termination condition
    if ( isLeaf() ){
        return sum;
    }
    //termination condition 2
    if (clustersTreeNode){
        for ( list<ClustersTreeNode*>::const_iterator iter = clustersTreeNode->begin();
              iter != clustersTreeNode->end(); ++iter ){
            if ( (*iter)->cluster->getAncestor() == this ){
                return sum;
            }
            /* NB: in the specific case where a cluster as no ancestor (it can happens
             * during ML search algorithms where clusters are build/destroyed on the fly),
             * cluster->getAncestor() == NULL means the cluster has not been assigned
             * anywhere yet. We assume that in this case its descendant are not assigned
             * as well and consequently we do not have to worry about checking for the ancestor
             * of this children cluster (they will not be assigned as well) */
        }
    }
    for (list<BasicNode*>::const_iterator iter = getChildrenList().begin();
         iter != getChildrenList().end(); ++iter){
        sum += ((InferenceNode*)(*iter))->findClusterBranchesRec( branchVector, clustersTreeNode );
    }
    return sum;
}

double InferenceNode::findClusterBranches( vector< pair<InferenceNode*,double> > & branchVector,
                                              ClustersTreeNode* clustersTreeNode ){
    if (clustersTreeNode){
        assert( clustersTreeNode->cluster->getAncestor() == this );
    }
    double sum = 0.0;
    for (list<BasicNode*>::const_iterator iter = getChildrenList().begin();
         iter != getChildrenList().end(); ++iter){
        sum += ((InferenceNode*)(*iter))->findClusterBranchesRec( branchVector, clustersTreeNode );
    }
    return sum;
}

double InferenceNode::findIncludedBranches( vector< pair<InferenceNode*,double> > & branchVector,
                                          ClustersTreeNode* clustersTreeNode ){
    //assert (ptree->findInsertion(this)==clustersTreeNode);
    return findClusterBranchesRec(branchVector, clustersTreeNode);
}

