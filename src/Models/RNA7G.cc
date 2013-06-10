#include "Models/RNA7G.h"
using namespace std;

RNA7G RNA7G::prototype( "RNA7G" );
RNA7G::RNA7G( const string & registrationName ): RNA7E( registrationName ) {}
RNA7G::RNA7G( ParametersSet & parameters ): RNA7E( parameters ) {}

bool RNA7G::setParameterNumbers() {
    numberRatesRatios = 1;
    numberFrequencies = 4;
	bool bp_symmetry = true;
	return bp_symmetry;
}

// Destructor
RNA7G::~RNA7G() {
}

string RNA7G::getName( void ) const {
    string name = "RNA7G" + MatrixModel::getName();
    return ( name );
}

double RNA7G::getFrequency( unsigned int residue,
    unsigned int rateCategory, unsigned int ) const {
    switch (residue){
        case 0: case 3:
            return (*frequencies)[rateCategory][0]/2.0;
        case 1: case 4:
            return (*frequencies)[rateCategory][1]/2.0;
        case 2: case 5:
            return (*frequencies)[rateCategory][2]/2.0;
        case 6:
            return (*frequencies)[rateCategory][3];
    }
    assert("wrong residue in RNA7G:getFrequency()"==0);
    return 0.0;
}

string RNA7G::getFrequencyState( unsigned int freqState, unsigned int ) const{
    assert(freqState<4);
    if (freqState==3) return "MM";
    return getState(freqState)+'+'+getState(freqState+3);
}

Model * RNA7G::clone( ParametersSet & parameters ) const {
    return new RNA7G( parameters );
}

