#include "Tree/HeterogeneousMolClkMCMCTree.h"

#include "Util/ParametersSet.h"
#include "Util/FileParser.h"

#include <iostream>

//register the tree to the tree factory with the name Heterogeneous rooted MCMC tree
HeterogeneousMolClkMCMCTree HeterogeneousMolClkMCMCTree::prototype("Heterogeneous MCMC tree with molecular clock");

HeterogeneousMolClkMCMCTree::HeterogeneousMolClkMCMCTree( const string & registrationName ):
MolClkMCMCTree( registrationName ){
}

HeterogeneousMolClkMCMCTree::HeterogeneousMolClkMCMCTree( ParametersSet& parameters ):
MolClkMCMCTree( parameters ){
}

HeterogeneousMolClkMCMCTree::~HeterogeneousMolClkMCMCTree(){
    if (swapPerturbator){
        delete swapPerturbator;
    }
}

MCMCTree* HeterogeneousMolClkMCMCTree::clone(ParametersSet& parameters) const{
    return new HeterogeneousMolClkMCMCTree(parameters);
}

double HeterogeneousMolClkMCMCTree::globalTopologyChange(){

    double lnHastingRatio;

    if ( randBox.ran() < 0.67 ) {
        lnHastingRatio = MolClkMCMCTree::globalTopologyChange();
    }
    //branch model swap
    else{
        //add one to the counter of swap proposals
        swapPerturbator->addPerturbation();
        modelSave.clear();
        vector< BasicNode * >::iterator iter = nodeRefVector.begin();
        ++iter;
        while ( iter != nodeRefVector.end() ){
            //for each node, prob swap = step of swapPerturbator
            if ( randBox.ran() < swapPerturbator->getStep() ){
                //save the old model
               Model* oldModel = ((InferenceNode*)(*iter))->getModel();
               unsigned int oldModelId = modelMap[oldModel];
               modelSave.push_back( pair<InferenceNode*, Model*>
                    ( (InferenceNode*)(*iter), oldModel ) );
               //choose a new model (different from the old one)
               unsigned int newModelId = (int)(randBox.ran()*(numberModels-1));
               if ( newModelId >= oldModelId ){
                   ++newModelId;
               }
               ((InferenceNode*)(*iter))->setModel(
                       pheterogeneous->getModel( newModelId ) );
               ((InferenceNode*)(*iter)->getParent())->invalidateRecursively(-1);
            }
            ++iter;
        }
        lastGlobalChange = SWAP_PROP;
        lnHastingRatio = 0.0;
    }
    return lnHastingRatio;
}

bool HeterogeneousMolClkMCMCTree::validateTopologyChange( bool validation ){
     if (lastGlobalChange != SWAP_PROP){
        return MolClkMCMCTree::validateTopologyChange( validation );
    }
    else{
        //warn the perturbator about the outcome (for adaptation of the step)
        swapPerturbator->validatePerturbation(validation);
        if (validation){
            for ( unsigned int i = 0; i < modelSave.size(); ++i ){
                ((InferenceNode*)((modelSave[i].first)->getParent()))->saveRecursively(-1);
            }
        }
        else{
            for ( unsigned int i = 0; i < modelSave.size(); ++i ){
                (modelSave[i].first)->setModel( modelSave[i].second );
                ((InferenceNode*)((modelSave[i].first)->getParent()))->retrieveRecursively(-1);
            }
        }
        lastGlobalChange = INVALID_GLOBAL;
    }
    return true;
}


void HeterogeneousMolClkMCMCTree::initialiseMCMC(ParametersSet& parameters ){
    MolClkMCMCTree::initialiseMCMC( parameters );
    swapPerturbator = new PerturbatorBase("Swap", parameters,
                                            "Swap", .05, .01, .9,
                                            "initial probability");
}

void HeterogeneousMolClkMCMCTree::initSampling(ParametersSet& parameters, bool overwrite){
    MCMCTreeBasic::initSampling(parameters, overwrite);
    string outputFileBaseName( parameters.stringParameter( "Output file" ) );

    FileParser::confirmOpenFile( branchModelsFile,
             outputFileBaseName + ".bm", overwrite );
}

void HeterogeneousMolClkMCMCTree::sample(){
    MCMCTreeBasic::sample();
    vector < BasicNode * >::const_iterator iter = nodeRefVector.begin();
    ++iter;
    //branch models sampling
    while ( iter != nodeRefVector.end() ){
        branchModelsFile << setw(4) << modelMap[((InferenceNode*)(*iter))->getModel()] << ' ';
        ++iter;
    }
    branchModelsFile << endl;
}

void HeterogeneousMolClkMCMCTree::printPerturbationParameters(ostream& outputStream){
    MolClkMCMCTree::printPerturbationParameters( outputStream );
    swapPerturbator->printPerturbationParameters( outputStream );
}


void HeterogeneousMolClkMCMCTree::stopBurn(){
    MolClkMCMCTree::stopBurn();
    swapPerturbator->stopBurn();
}

void HeterogeneousMolClkMCMCTree::loadDataAndModel
        ( SequenceTable * psequenceTable, Model * pmodel ){
    MolClkMCMCTree::loadDataAndModel( psequenceTable, pmodel );
    pheterogeneous = dynamic_cast<Heterogeneous*>(pmodel);
    if ( !pheterogeneous ){
        cerr << "You must use a heterogeneous model with the"
             << " heterogeneous tree" << endl;
        exit(EXIT_FAILURE);
    }
    numberModels = pheterogeneous->getNumberModels();
    for ( unsigned int i = 0; i < numberModels; ++i ){
        modelMap[pheterogeneous->getModel(i)] = i;
    }
    //if the tree is already constructed
    if (nodeRefVector.size()){
        initialiseNodeModel();
    }
}

void HeterogeneousMolClkMCMCTree::update( UpdateMessage* message ){
    int cat=-1;
    int model=-1;
    assert(message->hasType(UpdateMessage::MODEL_TYPE));
    if (message->hasFlag(UpdateMessage::MODEL_PRIOR_FLAG)){
        //cout << "message prior from the model" << endl;
        return;
    }
    if (message->hasFlag(UpdateMessage::SYMBOL_CATEGORY)){
        cat = (int)message->modelMsg.symbolCategory;
    }
#ifdef DEBUG1
    cout << "Heterogeneous MolClk tree receive update from symbol " << cat << endl;
#endif
    if (message->hasFlag(UpdateMessage::MODEL_FLAG)){
        model = (int)message->modelMsg.model;
    }
#ifdef DEBUG1
    cout << "Heterogeneous MolClk tree receive update from model " << model << endl;
#endif
    if ( model == (int)pheterogeneous->getNumberModels() ){
        return;
    }
    if ( model == -1 ){
        invalidateNodes( cat );
        return;
    }
    //else
    vector<BasicNode*>::iterator iter = nodeRefVector.begin();
    ++iter;
    while ( iter != nodeRefVector.end() ){
        if ( modelMap[((InferenceNode*)(*iter))->getModel()] == model ){
            ((InferenceNode*)((*iter)->getParent()))->invalidateRecursively( cat );
        }
        ++iter;
    }
}

void HeterogeneousMolClkMCMCTree::retrieveNodesAux(){
    MCMCTreeBasic::retrieveNodesAux();
}

void HeterogeneousMolClkMCMCTree::saveNodesAux(){
    MCMCTreeBasic::saveNodesAux();
}


void HeterogeneousMolClkMCMCTree::getAllParameters( vector<double>& params ) const{
    //get the branch lengths
    MolClkMCMCTree::getAllParameters( params );
    vector< double > modelId;
    vector<BasicNode*>::const_iterator iterNode = nodeRefVector.begin();
    ++iterNode;
    //get the model ID associated with the branch
    while ( iterNode != nodeRefVector.end() ){
        map< Model*, int>::const_iterator search =
                modelMap.find( ((InferenceNode*)(*iterNode))->getModel() );
        params.push_back( (double)((*search).second) );
        ++iterNode;
    }
    assert( params.size() == getNumberTreeParameters());
}

void HeterogeneousMolClkMCMCTree::setAllParameters( const vector<double>& params){
    assert( params.size() == getNumberTreeParameters());
    vector<double>::const_iterator iter = params.begin();
    vector<double> p;
    unsigned int nbp;
    //first set basic parameters
    nbp = MolClkMCMCTree::getNumberPerturbationParameters();
    p.resize(nbp);
    for ( unsigned int i = 0; i < nbp; ++i ){
        p[i] = *iter;
        ++iter;
    }
    MolClkMCMCTree::setAllParameters( p );
    //then set the model for each node
    vector<BasicNode*>::iterator iterNode = nodeRefVector.begin();
    ++iterNode;
    while ( iterNode != nodeRefVector.end() ){
        ((InferenceNode*)(*iterNode))->setModel(
                pheterogeneous->getModel( (int)(*iter) ) );
        ++iterNode;
        ++iter;
    }
}

unsigned int HeterogeneousMolClkMCMCTree::getNumberTreeParameters() const{
    unsigned int numberParameters = MolClkMCMCTree::getNumberTreeParameters();
    //add the number of branches to store the model associated with each of
    //them
    numberParameters += getNumberBranches();
    return numberParameters;
}

void HeterogeneousMolClkMCMCTree::getAllPerturbationParameters( vector<double>& params ) const{
    MolClkMCMCTree::getAllPerturbationParameters( params );
    vector<double> p;
    swapPerturbator->getAllPerturbationParameters( p );
    params.insert( params.end(), p.begin(), p.end() );
    assert( params.size() == getNumberPerturbationParameters() );
}

void HeterogeneousMolClkMCMCTree::setAllPerturbationParameters( const vector<double>& params){
    assert( params.size() == getNumberPerturbationParameters() );
    vector<double>::const_iterator iter = params.begin();
    vector<double> p;
    unsigned int nbp = MolClkMCMCTree::getNumberPerturbationParameters();
    p.resize(nbp);
    for ( unsigned int i = 0; i < nbp; ++i ){
        p[i] = *iter;
        ++iter;
    }
    MolClkMCMCTree::setAllPerturbationParameters( p );
    nbp = swapPerturbator->getNumberPerturbationParameters();
    p.resize(nbp);
    for ( unsigned int i = 0; i < nbp; ++i ){
        p[i] = *iter;
        ++iter;
    }
    swapPerturbator->setAllPerturbationParameters ( p );
}

unsigned int HeterogeneousMolClkMCMCTree::getNumberPerturbationParameters() const{
    unsigned int numberPerturbationParameters = MolClkMCMCTree::getNumberPerturbationParameters();
    numberPerturbationParameters += swapPerturbator->getNumberPerturbationParameters();
    return numberPerturbationParameters;
}

void HeterogeneousMolClkMCMCTree::getAllPriorParameters( vector<double>& params ) const{
    MolClkMCMCTree::getAllPriorParameters(params);
}

void HeterogeneousMolClkMCMCTree::setAllPriorParameters( const vector<double>& params){
    MolClkMCMCTree::setAllPriorParameters(params);
}

unsigned int HeterogeneousMolClkMCMCTree::getNumberPriorParameters() const{
    return MolClkMCMCTree::getNumberPriorParameters();
}

double HeterogeneousMolClkMCMCTree::getLnPrior() const{
    return MolClkMCMCTree::getLnPrior();
}

void HeterogeneousMolClkMCMCTree::constructFromString( string stringTree ){
    MolClkMCMCTree::constructFromString( stringTree );
    //if the model has been loaded
    if(pheterogeneous){
        initialiseNodeModel();
    }
}

void HeterogeneousMolClkMCMCTree::constructRandomly( double maxBranchLength ){
    MolClkMCMCTree::constructRandomly( maxBranchLength );
    //if the model has been loaded
    if(pheterogeneous){
        initialiseNodeModel();
    }
}

void HeterogeneousMolClkMCMCTree::initialiseNodeModel(){
    vector< BasicNode * >::iterator iter = nodeRefVector.begin();
    while ( iter != nodeRefVector.end() ){
        if(*iter!=root){
            ((InferenceNode*)(*iter))->setModel( pheterogeneous->getModel(
                    (int)(randBox.ran()*(double)numberModels) ) );
        }
        ++iter;
    }
}
