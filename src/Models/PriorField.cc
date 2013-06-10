#include "Models/PriorField.h"

#include <algorithm>
#include <iostream>
#include <cstdio>

PriorField::PriorField( double value ){
    this->value = value;
    cst = true;
}

PriorField::PriorField( string strField ){
    cst = false;
    hyperPrior = new HyperPrior;
    //read prior name
    string::iterator iter = find(strField.begin(),strField.end(),'(');
    hyperPrior->priorName = string( strField.begin(), iter );
    if (iter==strField.end()) return;
    //get parameters/hyperprior
    string tmpString;
    ++iter;
    string::iterator endStringIter = strField.end();
    --endStringIter;
    unsigned int parenthesisCount = 1;
    bool isNotConstant = false;
    while ( iter != strField.end() ){
        if( (( *iter==',' ) || ( *iter==')' )) && (parenthesisCount==1) ){
            //if ')' check that we are at the end
            if ( *iter==')' ){
                parenthesisCount = 0;
                //if we are not at the last character then error
                if ( iter != endStringIter ){
                    cerr << "invalid prior " << strField << endl;
                    exit(EXIT_FAILURE);
                }
            }
            //if not a double, add a PriorField to the fields vector
            if(isNotConstant){
                hyperPrior->fields.push_back( new PriorField(tmpString) );
            }
            //otherwise add a ConstantValue
            else{
                char* endPtr;
                double v = strtod( tmpString.c_str(), & endPtr );
                assert(*endPtr == '\0');
                hyperPrior->fields.push_back( new PriorField(v) );
            }
            isNotConstant = false;
            tmpString.clear();
        }
        else{
           switch (*iter){
                case '(': ++parenthesisCount; break;
                case ')':
                    if ( (parenthesisCount == 0) ){
                        cerr << "invalid prior ('(' missing?): "
                             << strField << endl;
                        exit(EXIT_FAILURE);
                    }
                    --parenthesisCount;
                break;
           }
           //push the character on tmpString
           if (!isspace(*iter)){
               tmpString.push_back(*iter);
               if ( ( *iter != '.' ) && ( *iter != '-' ) && ( ( *iter<'0' ) || ( *iter>'9' ) ) ){
                   isNotConstant = true;
               }
           }
        }
        ++iter;
    }
    if(parenthesisCount){
        cerr << "invalid prior (')' missing?): " << strField << endl;
        exit(EXIT_FAILURE);
    }
}


string PriorField::toString() const{
    if (cst){
        char ret[20];
        sprintf( ret, "%.2f", value );
        return string( ret );
    }
    else{
        string ret(hyperPrior->priorName+'(');
        vector<PriorField*>::iterator iter = hyperPrior->fields.begin();
        while(iter != hyperPrior->fields.end()){
            ret += (*iter)->toString();
            ++iter;
            if (iter != hyperPrior->fields.end()) ret += ',';
        }
        ret += ')';
        return ret;
    }
}



PriorField::~PriorField(){
    if ( (!cst) && (hyperPrior) ){
        for( vector<PriorField*>::iterator iter = hyperPrior->fields.begin();
             iter != hyperPrior->fields.end(); ++iter ){
            delete *iter;
        }
        delete hyperPrior;
    }
}
