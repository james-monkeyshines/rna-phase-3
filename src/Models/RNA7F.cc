#include "Models/RNA7F.h"
using namespace std;

RNA7F RNA7F::prototype( "RNA7F" );
RNA7F::RNA7F( const string & registrationName ): RNA7D( registrationName ) {}
RNA7F::RNA7F( ParametersSet & parameters ): RNA7D( parameters ) {}

bool RNA7F::setParameterNumbers() {
    numberRatesRatios = 3;
    numberFrequencies = 4;
	bool bp_symmetry = true;
	return bp_symmetry;
}

// Destructor
RNA7F::~RNA7F() {
}

string RNA7F::getName( void ) const {
    string name = "RNA7F" + MatrixModel::getName();
    return ( name );
}

double RNA7F::getFrequency( unsigned int residue,
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
    assert("wrong residue in RNA7F:getFrequency()"==0);
    return 0.0;
}

string RNA7F::getFrequencyState( unsigned int freqState, unsigned int ) const{
    assert(freqState<4);
    if (freqState==3) return "MM";
    return getState(freqState)+'+'+getState(freqState+3);
}

Model * RNA7F::clone( ParametersSet & parameters ) const {
    return new RNA7F( parameters );
}

