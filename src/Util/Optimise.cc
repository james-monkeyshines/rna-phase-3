#include "Util/Optimise.h"

#include <assert.h>

#include <iostream>
#include <iomanip>
#include <numeric>

#include "Tree/OptimizerTree.h"

#include "Models/Model.h"

#define ITMAX 1000
#define ZEPS 1.0e-10
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SQR(a) (a*a)
#define EPS 1.0e-8
#define TOLX 1.0e-8
#define ALF 1.0e-8
#define STPMAX 1.0
#define SIGN(n) (n>=0.0 ? 1.0 : -1.0)


// Find the vector of partial gradients for the branch lengths (sqrt)
// branchGradientVector must have the right size
// branchesLengthSqrt must be the actual value
void Optimise::diffLoglikelihoodSqrt( OptimizerTree & tree ,
const vector < double > & branchesLengthSqrt,
vector < double > & branchGradientVector) {

    tree.diffloglikelihood(branchGradientVector) ;
    for ( unsigned int i = 0 ; i < branchGradientVector.size() ; ++i ) {
        branchGradientVector[i] *= (2.0 * branchesLengthSqrt[i]) ;
#ifdef DEBUG1
        assert( !isnan( branchGradientVector[i] ) ) ;
        assert( !isinf( branchGradientVector[i] ) ) ;
#endif
    }
}

// Find the vector of partial gradients of the model parameters
// gradientVector must have the right size
// model parameters must be the actual parameters of the model
void Optimise::modelPartialGradients( OptimizerTree & tree, Model & model,
vector < double > & modelParameters, vector<unsigned int> &optParameters, double step,
vector<double> &gradientVector ) {

    double newLikelihood;
    double newPenalty;
    
    //load the actual likelihood of the model and associated penalty
    double oldLikelihood = tree.loglikelihood();
    double oldPenalty = model.getLnPenalty();
    // for each parameter of the model that should be optimised
	for ( unsigned int j = 0; j < optParameters.size(); j++) {
    	unsigned int i = optParameters.at(j);
		
        // modify the given parameter of the model by step
        double tmp = modelParameters[i];
        modelParameters[i] += step ;
        model.setAllParameters( modelParameters );
        model.validChange();
        //and compute the new likelihood
        //tree.invalidateNodes(-1) ;
        newLikelihood = tree.loglikelihood();
        newPenalty = model.getLnPenalty();
        // put back the parameter to its old value
        modelParameters[i] = tmp;
        gradientVector[i] = ( newLikelihood - oldLikelihood + newPenalty - oldPenalty) / step ;
        assert( !isnan( gradientVector[i] ) );
    }

    // restore original parameters
    model.setAllParameters( modelParameters );
    model.validChange();
}

void Optimise::modelPenaltyPartialGradients( Model & model,
vector < double > & modelPenalty, double step,
vector<double> &gradientVector ) {
    
    double newPenalty;
    double oldPenalty = model.getLnPenalty();
    assert(modelPenalty.size());
    
    vector<double> exactGradientVector;
    model.diffLnPenalty(exactGradientVector);   
    assert(exactGradientVector.size()<=modelPenalty.size());
    
    copy(exactGradientVector.begin(),exactGradientVector.end(),gradientVector.begin());
    for ( unsigned int i = exactGradientVector.size() ; i < modelPenalty.size(); i++ ) {
        // modify the given parameter of the model by step
        double tmp = modelPenalty[i];
        modelPenalty[i] += step;
        model.setAllPenaltyParameters( modelPenalty );
        model.validChange();
        newPenalty = model.getLnPenalty();
        // put back the parameter to its old value
        modelPenalty[i] = tmp;
        gradientVector[i] = ( newPenalty - oldPenalty ) / step ;
        assert( !isnan( gradientVector[i] ) );
    }
    model.setAllPenaltyParameters( modelPenalty );
    model.validChange();
//assert(0);
}
  
double Optimise::logOptValue( OptimizerTree & tree, Model & model,
const vector<double>& sqrtInitialBranches, vector<double>& newBranches,
const vector<double>& initialModelPenalty, vector<double>& newModelPenalty,  
const vector<double>& initialParameters, vector<double>& newParameters,
vector<double>& direction, double magnitude ) {

    unsigned int branchSize = sqrtInitialBranches.size();
	bool zeroLength = false;
    // Determine new branch parameters
    for (  unsigned int i = 0 ; i < branchSize ; ++i ) {
       newBranches[i] = sqrtInitialBranches[i] + ( magnitude * direction[i] );
       newBranches[i] = newBranches[i] * newBranches[i];
	   if (newBranches[i] < TOLX) {
			zeroLength = true;
	   }
    }
	// Disallow zero-length branches; add a tiny amount to all branches as sort of pseudocount.
	if (zeroLength) {
    	for ( unsigned int i = 0 ; i < branchSize ; ++i ) {
			//newBranches[i] += 1e-6;
	   }
    }
	
    unsigned int modelPenaltySize = initialModelPenalty.size();
    //new penalty
    if(modelPenaltySize){
        for ( unsigned int i = 0; i < modelPenaltySize ; ++i ) {
            newModelPenalty[i] = initialModelPenalty[i] + ( magnitude * direction[branchSize+i] );
        }
    }
    // Determine new model parameters
    unsigned int modelSize = initialParameters.size();
    for ( unsigned int i = 0; i < modelSize ; ++i ) {
        newParameters[i] = initialParameters[i] + ( magnitude * direction[branchSize+modelPenaltySize+i] );
    }
    // return log-likelihood with new branches' length and new model parameters & penalty
    tree.setBranchVector( newBranches );
    model.setAllParameters( newParameters );
    model.setAllPenaltyParameters( newModelPenalty );
    model.validChange();
    return ( tree.loglikelihood() + model.getLnPenalty()) ;	
}


double Optimise::logOptValue( OptimizerTree & tree, Model & model,
const vector<double>& sqrtInitialBranches, vector<double>& newBranches,
const vector<double>& initialModelPenalty, vector<double>& newModelPenalty,  
vector<double>& direction, double magnitude ) {

    unsigned int branchSize = sqrtInitialBranches.size();
	bool zeroLength = false;
    // Determine new branch parameters
    for ( unsigned int i = 0 ; i < branchSize ; i++ ) {
       newBranches[i] = sqrtInitialBranches[i] + ( magnitude * direction[i] );
       newBranches[i] = newBranches[i] * newBranches[i];
	   if (newBranches[i] < TOLX) {
			zeroLength = true;
	   }
    }
	// Disallow zero-length branches; add a tiny amount to all branches as sort of pseudocount.
	if (zeroLength) {
    	for ( unsigned int i = 0 ; i < branchSize ; ++i ) {
			//newBranches[i] += 1e-6;
	   }
    }
	
    // return log-likelihood with new branches' length
    tree.setBranchVector(newBranches);
    //new penalty
    unsigned int modelPenaltySize = initialModelPenalty.size();
    if(modelPenaltySize){
        for ( unsigned int i = 0; i < modelPenaltySize ; ++i ) {
            newModelPenalty[i] = initialModelPenalty[i] + ( magnitude * direction[branchSize+i] );
        }
        model.setAllPenaltyParameters( newModelPenalty );
        model.validChange();
    } 
    return ( tree.loglikelihood() + model.getLnPenalty() ) ;
}



double Optimise::lineSearch( OptimizerTree & tree, Model & model,
const vector < double > & branchVector, const vector < double > & sqrtBranchVector,
const vector < double > & modelPenaltyVector, const vector < double > & modelVector,
vector < double > & direction, vector < double > & oppGradient, bool modelOpt ) {

    double temp , alam2 = 0.0 ,
    disc , tmplam = 0.0;
    double rhs1 = 0.0 , rhs2 = 0.0 , f2 = 0.0 , a = 0.0 , b = 0.0 ;

    //storage space to avoid repetitive memory allocation in subprocedure
    vector < double > branchJunk;
    branchJunk.resize(sqrtBranchVector.size());
    vector < double > penaltyJunk;
    penaltyJunk.resize(modelPenaltyVector.size());
    vector < double > modelJunk;
    modelJunk.resize(modelVector.size());
    // compute the slope
    double slope = 0.0 ;
    for ( unsigned int i = 0 ; i < direction.size() ; ++i ){
        slope -= (oppGradient[i] * direction[i]);
    }


    assert( slope < 0.0 );

    // Convergence testing
    // alamin = TOLX/max(direction[i]/point[i])
    double max = 0.0 ;
    for ( unsigned int i = 0 ; i < branchVector.size() ; ++i ) {
        temp = fabs( direction[i] ) / ( MAX( fabs( branchVector[i] ), 1.0 ) );
        if ( temp > max ){
            max = temp ;
        }
    }
    for ( unsigned int i = 0 ; i < modelPenaltyVector.size() ; ++i ) {
        temp = fabs( direction[branchVector.size()+i] ) / ( MAX( fabs( modelPenaltyVector[i] ), 1.0 ) );
        if ( temp > max ){
            max = temp ;
        }
    }
    for ( unsigned int i = 0 ; i < modelVector.size() ; ++i ) {
        temp = fabs( direction[branchVector.size()+modelPenaltyVector.size()+i] ) / ( MAX( fabs( modelVector[i] ), 1.0 ) );
        if ( temp > max ){
            max = temp ;
        }
    }
    double alamin = TOLX / max;
    double alam = 1.0 ;


    // Initial value, compute the likelihood in point
    double initEnergy = -tree.loglikelihood() - model.getLnPenalty();
    double f;
    //convergence testing inside the loop, it will exit with a return
    while (1) {
        // Try the full newton step first
        // compute the likelihood in point + alam*direction
        if ( modelOpt && modelVector.size() ){
            f = -logOptValue( tree, model, sqrtBranchVector, branchJunk, modelPenaltyVector, penaltyJunk,
                                           modelVector, modelJunk, direction, alam );
        }
        else{
            f = -logOptValue( tree, model, sqrtBranchVector, branchJunk, modelPenaltyVector, penaltyJunk, direction, alam );
        }

        // if alam becomes lower than the threshold : converged
        if ( alam < alamin ) {
            tree.setBranchVector( branchVector );
            if ( modelOpt && modelVector.size() ){
                model.setAllParameters( modelVector );
            }
            if(modelPenaltyVector.size()){
                model.setAllPenaltyParameters( modelPenaltyVector );
            }
            if (( modelOpt && modelVector.size() ) || modelPenaltyVector.size() ){
                model.validChange();
            }
            return ( 0.0 ) ;
        }
        // otherwise
        else {
            if ( f <= initEnergy + ALF * alam * slope ){
                tree.setBranchVector( branchVector );
                if ( modelOpt && modelVector.size() ){
                    model.setAllParameters( modelVector );
                }
                if(modelPenaltyVector.size()){
                    model.setAllPenaltyParameters( modelPenaltyVector );
                }
                if (( modelOpt && modelVector.size() ) || modelPenaltyVector.size() ){
                    model.validChange();
                }
                return alam ;
            }
            else {
                //initialisation for the first run through the loop
                if ( alam == 1.0 )
                    tmplam = -slope / ( 2.0 * ( f - initEnergy - slope ) ) ;
                else {
                    rhs1 =   f - initEnergy - ( alam * slope ) ;
                    rhs2 = f2 - initEnergy - ( alam2 * slope );
                    a = ( rhs1 / ( alam * alam ) - rhs2 / ( alam2 * alam2 ) ) /
                    ( alam - alam2 ) ;
                    b = ( -(alam2 * rhs1) / ( alam * alam ) + ( alam * rhs2 ) /
                    ( alam2 * alam2 ) ) / ( alam - alam2 ) ;
                    if ( a == 0.0 )
                        tmplam = -slope / ( 2.0 * b ) ;
                    else {
                        disc = b * b - 3.0 * a * slope ;
                        if ( disc < 0.0 )
                            tmplam = 0.5 * alam ;
                        else {
                            if ( b <= 0.0 )
                                tmplam = ( -b + sqrt( disc ) ) / ( 3.0 * a );
                            else
                                tmplam = -slope / ( b + sqrt( disc ) ) ;
                        }
                    }
                    if ( (tmplam>(0.5*alam)) || (isnan(tmplam)) || (isinf(tmplam)) ){
                        tmplam = 0.5 * alam ;
                    }
                }
            }
            alam2 = alam ;
            f2 = f ;
            alam = MAX( tmplam, ( 0.1 * alam ) );
        }
    } // End of optimization loop
} // End lineSearch



// The quasi-newton (BGFS update) optimization function
vector < double > Optimise::optimiseQuasiNewton( OptimizerTree & tree ,
Model & model , bool modelOpt , bool empFreqs , double tolerance, bool display, bool displayPoint ) {

    vector< double > modelVector;
    model.getAllParameters( modelVector );
    vector< double > modelPenaltyVector;
    model.getAllPenaltyParameters( modelPenaltyVector );
    unsigned int nbModelParameters = modelVector.size() ;
    unsigned int nbModelPenaltyParameters = modelPenaltyVector.size() ;
    
    vector < double > branchVector;
    tree.getBranchVector( branchVector );
    unsigned int nbBranches = branchVector.size();
    vector < double > sqrtBranchVector;
    sqrtBranchVector.resize(nbBranches);
    for ( unsigned int i = 0 ; i < nbBranches ; ++i ){
        sqrtBranchVector[i] = sqrt( branchVector[i] ) ;
    }
    
    
    unsigned int numberParameters = nbBranches + nbModelPenaltyParameters;
    if (modelOpt){
        numberParameters += nbModelParameters;
    }

    vector<double> startingPos(numberParameters);
    vector<double> newGradientVector(numberParameters);
    vector<double> oppGradientVector(numberParameters);
    vector<double> oldGradientVector(numberParameters);
    vector<double> modelGradientVector(nbModelParameters);
    vector<double> branchGradientVector(nbBranches);
    vector<double> modelPenaltyGradientVector(nbModelPenaltyParameters);
    vector<double> direction(numberParameters);
    vector<unsigned int> optParameters;
    
    
    array2D<double> H(numberParameters, numberParameters);
    array2D<double> update1(numberParameters, numberParameters);
    array2D<double> update2(numberParameters, numberParameters);
    array2D<double> update3(numberParameters, numberParameters);
    vector<double> delta;
    delta.resize(numberParameters);
    vector<double> omega;
    omega.resize(numberParameters);
    vector<double> u;
    u.resize(numberParameters);
    vector<double> now;
    now.resize(numberParameters);
    vector<double> tempv;
    tempv.resize(numberParameters);


    double alpha, tempx, max1, max2, den;

    tree.invalidateNodes(-1);
    double oldOptValue;
    double newOptValue = tree.loglikelihood() + model.getLnPenalty();
    assert(!isnan(newOptValue));
    assert(!isinf(newOptValue));
      
    double fac, fae , sumdg , sumxi;

      
    // Initialise H to the identity matrix
    for ( unsigned int i = 0 ; i < numberParameters ; ++i ) {
        H( i, i ) = 1.0 ;
        for ( unsigned int j = ( i + 1 ) ; j < numberParameters ; ++j ){
            H( i, j ) = H( j, i ) = 0.0 ;
        }
    }

    // Find the vector of partial gradients w.r.t to SQRT(the tree branch lengths)
    diffLoglikelihoodSqrt(tree, sqrtBranchVector, branchGradientVector);
    vector<double>::iterator last = 
            copy(branchGradientVector.begin(), branchGradientVector.end(),newGradientVector.begin()) ;
    // Find the vector of partial gradients w.r.t to the model penalty parameters
    if(nbModelPenaltyParameters){
        modelPenaltyPartialGradients( model, modelPenaltyVector, 1e-6, modelPenaltyGradientVector );
        //model.diffLnPenalty(modelPenaltyGradientVector);
        last = copy(modelPenaltyGradientVector.begin(), modelPenaltyGradientVector.end(), last);
    }
    // Find the vector of partial gradients w.r.t to the model parameters
    if ( ( modelOpt ) && ( nbModelParameters > 0 ) ){
		// Work out which parameters to optimise (ie skip frequency parameters if using empirical values)
		model.getOptimisableParameters( empFreqs, optParameters );
        modelPartialGradients( tree, model, modelVector, optParameters, 1e-6, modelGradientVector );
        copy(modelGradientVector.begin(), modelGradientVector.end(),last);
    }
    
    //initialisation done, start the optimisation loop
    for (unsigned int iter = 1; iter<ITMAX; ++iter){
        oldOptValue = newOptValue ;
        oldGradientVector = newGradientVector ; // Save the gradient vector

        // Calculate the maximum step size stepMax with sum = (sum of all squared parameters)
        double sum = 0.0;
        for ( unsigned int i = 0; i < nbBranches; ++i ){
            sum += branchVector[i];
        }
        for ( unsigned int i = 0; i < nbModelPenaltyParameters; ++i ){
            sum += ( modelPenaltyVector[i] * modelPenaltyVector[i] ) ;
        }
        if ( modelOpt ) {
            for ( unsigned int i = 0; i < nbModelParameters; ++i ){
                sum += ( modelVector[i] * modelVector[i] ) ;
            }
        }
        double stepMax = STPMAX*MAX( sqrt(sum), ( double )numberParameters );
    
        // Determine search direction and scale the direction vector if necessary
        sum = 0.0;
        for ( unsigned int i = 0 ; i < numberParameters ; ++i ) {
            direction[i] = 0.0 ;
            for ( unsigned int j = 0 ; j < numberParameters ; ++j ) {
                direction[i] += ( H( i, j ) * newGradientVector[j] );
            }
            if ( isnan(direction[i]) || isinf(direction[i]) ) {
                cerr << "Error iteration " << iter << ", search direction inf or nan" << endl ;
                cerr << "direction i = " << i << endl ;
                direction[i] = 0.0;
            }
            sum += ( direction[i] * direction[i] );
        }
        sum = sqrt( sum ) ;
                  
        //scaling if necessary    
        if ( sum > stepMax ){
            double factor = stepMax / sum;
            for ( unsigned int i = 0 ; i < numberParameters ; ++i ){
                direction[i] *= factor;
            }
        }

       
        // Perform line search, new and old gradient vector are in fact
        // opposite of the real gradient
        alpha = lineSearch( tree, model, branchVector, sqrtBranchVector, modelPenaltyVector,
                            modelVector, direction, newGradientVector, modelOpt );

        // Calculate new branch length vector
        for ( unsigned int i = 0 ; i < nbBranches ; ++i ) { // Calculate new branch parameter vector
            delta[i]   = ( alpha * direction[i] ) ;
            sqrtBranchVector[i] += delta[i];
            branchVector[i] = sqrtBranchVector[i] * sqrtBranchVector[i];
        }
        if ( nbModelPenaltyParameters ) {
            for ( unsigned int i = 0 ; i < nbModelPenaltyParameters ; ++i ) { 
                delta[i + nbBranches] = alpha * direction[i + nbBranches];
                modelPenaltyVector[i] += delta[i + nbBranches] ;
            }
        }
        // Calculate new model parameter vector
        if ( modelOpt && nbModelParameters ) {
            for ( unsigned int i = 0 ; i < nbModelParameters ; ++i ) { 
                delta[i + nbBranches + nbModelPenaltyParameters] = alpha * direction[i + nbBranches + nbModelPenaltyParameters];
                modelVector[i] += delta[i + nbBranches+ nbModelPenaltyParameters] ;
            }
        }

        // Set new values for the tree
        tree.setBranchVector( branchVector ) ;
        // Set new parameters for the model
        if ( nbModelPenaltyParameters > 0 ) {
            model.setAllPenaltyParameters( modelPenaltyVector );
            model.validChange();
        }
        if ( modelOpt && ( nbModelParameters > 0 ) ) {
            model.setAllParameters( modelVector );
            model.validChange();
        }
        // Calculate new value for the optimized quantity and print to screen
        newOptValue = tree.loglikelihood();
        if ( display ){
            cout << "Loglik=" << newOptValue;
        }
        double lnPenalty = model.getLnPenalty();
        newOptValue += lnPenalty;
        if ( display && lnPenalty ){
            cout << "      logPenalty=" << lnPenalty
                 << "      total=" << newOptValue;
        }
        
        if ( display ){
            cout << endl;
        }
        else{
            if ( displayPoint ){
                cout << "." << flush;
            }
        }
       
        if (nbModelPenaltyParameters){
            for ( unsigned int i = 0 ; i < nbModelPenaltyParameters ; ++i ) { 
                cout << modelPenaltyVector[i] << ' ';
            }
            cout << endl;
        }
        
        
        
        // Find the vector of partial gradients w.r.t to SQRT(the tree branch lengths)
        diffLoglikelihoodSqrt(tree, sqrtBranchVector, branchGradientVector);
        last = copy( branchGradientVector.begin(), branchGradientVector.end(),
                         newGradientVector.begin() ) ;
        // Find the vector of partial gradients w.r.t to the model penalty parameters
        if(nbModelPenaltyParameters){
            modelPenaltyPartialGradients( model, modelPenaltyVector, 1e-6, modelPenaltyGradientVector );
            //model.diffLnPenalty(modelPenaltyGradientVector);
            last = copy(modelPenaltyGradientVector.begin(), modelPenaltyGradientVector.end(), last);
        }
        // Find the vector of partial gradients w.r.t to the model parameters
        if ( ( modelOpt ) && ( nbModelParameters > 0 ) ){
            modelPartialGradients( tree, model, modelVector, optParameters, 1e-6, modelGradientVector );
            copy(modelGradientVector.begin(), modelGradientVector.end(),last);
        }
            
        // Since tgvector is actually the negative of the gradient
        for ( unsigned int i = 0 ; i < numberParameters ; ++i ) {
            omega[i] = oldGradientVector[i] - newGradientVector[i];
        }

        for ( unsigned int i = 0 ; i < numberParameters ; ++i ) {
            tempv[i] = 0.0 ;
            for ( unsigned int j = 0 ; j < numberParameters ; ++j ) {
                tempv[i] += ( H( i, j ) * omega[j] ) ;
            }
        }


        fac = fae = sumdg = sumxi = 0.0 ;
        for ( unsigned int i = 0; i < numberParameters; ++i ){
            fac   += ( omega[i] * delta[i] ) ;
            fae   += ( omega[i] * tempv[i] )  ;
            sumdg += ( omega[i] * omega[i] ) ;
            sumxi += ( delta[i] * delta[i] ) ;
        }


/************* Perform two types of convergence test  *******************/
        // First test for convergence in value
        max1 = 0.0;
        for ( unsigned int i = 0; i < nbBranches; ++i ) {
            tempx = fabs( delta[i] ) / MAX(fabs(sqrtBranchVector[i]), 1.0);
            if ( tempx > max1 ){
                max1 = tempx ;
            }
        }
        for ( unsigned int i = 0 ; i < nbModelPenaltyParameters ; ++i ) {
            tempx = fabs(delta[i+nbBranches]) / MAX(fabs(modelPenaltyVector[i]), 1.0);
            if ( tempx > max1 )
                max1 = tempx ;
        }
        if(modelOpt){
            for ( unsigned int i = 0; i < nbModelParameters; ++i ) {
                tempx = fabs( delta[i+nbBranches+nbModelPenaltyParameters] ) / MAX(fabs(modelVector[i]), 1.0);
                if ( tempx > max1 ){
                    max1 = tempx ;
                }
            }
        }
        // Second test for convergence in gradient
        max2 = 0.0 ;
        den = MAX( -newOptValue, 1.0 ) ;
        for ( unsigned int i = 0 ; i < nbBranches ; ++i ) {
            tempx = fabs(newGradientVector[i]) * MAX(fabs(sqrtBranchVector[i]), 1.0)/den ;
            if ( tempx > max2 )
                max2 = tempx ;
        }
        for ( unsigned int i = 0 ; i < nbModelPenaltyParameters ; ++i ) {
            tempx = fabs(newGradientVector[i+nbBranches]) * MAX(fabs(modelPenaltyVector[i]), 1.0)/den ;
            if ( tempx > max2 )
                max2 = tempx ;
        }
        if(modelOpt){
            for ( unsigned int i = 0 ; i < nbModelParameters ; ++i ) {
                tempx = fabs(newGradientVector[i+nbBranches+nbModelPenaltyParameters]) * MAX(fabs(modelVector[i]), 1.0)/den ;
                if ( tempx > max2 ){
                    max2 = tempx ;
                }
            }
        }
        if ( (max1<TOLX) || (max2<tolerance) ){
            vector<double>::iterator lastElt = 
                copy(sqrtBranchVector.begin(), sqrtBranchVector.end(),now.begin()) ;
            if( nbModelPenaltyParameters ){
                vector<double>::iterator ls = copy(modelPenaltyVector.begin(), modelPenaltyVector.end(), lastElt);
                lastElt=ls;
            }
            if ( modelOpt && nbModelParameters ) {
                copy(modelVector.begin(), modelVector.end(), lastElt);
            }
            if ( display ){
                cout << endl;
            }
            return now;
        }
// end gradient convergence test



//BGFS update
        if ( fac > sqrt( EPS * sumdg * sumxi ) ) {

            for ( unsigned int i = 0; i < numberParameters; ++i ) {
                u[i] = ( delta[i] / fac ) - ( tempv[i] / fae ) ;
            }

            for ( unsigned int i = 0; i < numberParameters; ++i ) {
                for ( unsigned int j = 0; j < numberParameters; ++j ) {
                    update3( i, j ) = u[i] * u[j] * fae ;
                }
            }

            fac = 1.0 / fac ;
            fae = 1.0 / fae ;
            for ( unsigned int i = 0; i < numberParameters; ++i ) {
                for ( unsigned int j = 0; j < numberParameters; ++j ) {
                    update1( i, j ) = delta[i] * delta[j] * fac  ;
                }
            }

            for ( unsigned int i = 0; i < numberParameters; ++i ) {
                for ( unsigned int j = 0; j < numberParameters; ++j ) {
                    update2( i, j ) = tempv[i] * tempv[j] * fae;
                }
            }
            

            // bfgs update
            for ( unsigned int i = 0; i < numberParameters; ++i ) {
                for ( unsigned int j = 0; j < numberParameters; ++j ) {
                    H( i, j ) = H( i, j ) + update1( i, j ) -
                    update2( i, j ) + update3( i, j );
                }
            }

#ifdef DEBUG1
            for ( unsigned int i = 0; i < numberParameters ; ++i) {
                for ( unsigned int j = 0; j < numberParameters ; ++j ) {
                    assert( !isnan( H( i, j ) ) );
                }
            }
#endif
        }
#ifdef DEBUG1
        for ( unsigned int i = 0; i < numberParameters ; ++i ) {
            assert( !isnan( newGradientVector[i] ) ) ;
            assert( !isinf( newGradientVector[i] ) ) ;
        }
#endif

    } //end loop
    assert(0);
    return ( now ) ;
}




