#include "Models/RNA7B.h"
using namespace std;

RNA7B RNA7B::prototype( "RNA7B" );
RNA7B::RNA7B( const string & registrationName ): RNA7A( registrationName ) {}
RNA7B::RNA7B( ParametersSet & parameters ): RNA7A( parameters ) {}

bool RNA7B::setParameterNumbers() {
    numberRatesRatios = 20;
    numberFrequencies = 4;
	bool bp_symmetry = true;
	return bp_symmetry;
}

// Destructor
RNA7B::~RNA7B() {}

string RNA7B::getName( void ) const {
    string name = "RNA7B" + MatrixModel::getName();
    return ( name );
}

double RNA7B::getFrequency( unsigned int residue,
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
    assert("wrong residue in RNA7B:getFrequency()"==0);
    return 0.0;
}

string RNA7B::getFrequencyState( unsigned int freqState, unsigned int ) const{
    assert(freqState<4);
    if (freqState==3) return "MM";
    return getState(freqState)+'+'+getState(freqState+3);
}

Model * RNA7B::clone( ParametersSet & parameters ) const {
    return new RNA7B( parameters );
}

