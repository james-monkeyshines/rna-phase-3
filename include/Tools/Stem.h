#ifndef STEM_H
#define STEM_H

#include <iostream>

using namespace std;

class Stem{
public:

    string name;

    //position of the 5' and 3' strand        
    pair< unsigned int, unsigned int > pos;

    //flag
    bool pseudoKnot;
    
    inline Stem( const string& name, unsigned int openingPosition ){
        this->name = name;
        pos.first = openingPosition;
        pseudoKnot = false;
    }
};


#endif
