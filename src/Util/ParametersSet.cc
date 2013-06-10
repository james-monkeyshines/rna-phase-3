#include "Util/ParametersSet.h"


using namespace std;

#include <algorithm>
#include <functional>
#include <iostream>

ParametersSet::ParametersSet( const string & setName ) :
name( setName ) {
}

bool ParametersSet::findParameter( const string & parameter ) const{
    return (this->find(parameter) != end());
}

bool ParametersSet::findCategory( const string& categoryName ) const{
    for ( vector < pair<ParametersSet, bool> >::const_iterator i = categories.begin();
              i != categories.end(); ++i ) {
        if ( (*i).first.name == categoryName ) {
            return true;
        }
    }
    return false;
}


int ParametersSet::intParameter( const string & parameter ){
    char * endPtr;
    int result;
    iterator i = find( parameter );
    if ( i == end() ) {
        cerr << "missing parameter : \"" << parameter << "\" in the set "
        << name << endl;
        exit(EXIT_FAILURE);
    }
    result = strtol( ( * i ).second.first.c_str(), & endPtr, 0 );
    if ( * endPtr != '\0' ) {
        cerr << "The parameter : \"" << parameter << "\" in the set "
        << name << " is not an integer" << endl;
        exit(EXIT_FAILURE);
    }
    (*i).second.second = true;
    return result;
}

double ParametersSet::doubleParameter( const string & parameter ){
    char * endPtr;
    double result;
    iterator i = find( parameter );
    if ( i == end() ) {
        cerr << "missing parameter : \"" << parameter << "\" in the set "
        << name << endl;
        exit(EXIT_FAILURE);
    }
    result = strtod( (*i).second.first.c_str(), & endPtr );
    if ( * endPtr != '\0' ) {
        cerr << "The parameter : \"" << parameter << "\" in the set "
        << name << " is not an integer" << endl;
        exit(EXIT_FAILURE);
    }
    (*i).second.second = true;
    return result;
}

const string& ParametersSet::stringParameter( const string & parameter ){
    iterator i = find( parameter );
    if ( i == end() ) {
        cerr << "missing parameter : \"" << parameter << "\" in the set "
        << name << endl;
        exit(EXIT_FAILURE);
    }
    (*i).second.second = true;
    return (*i).second.first;
}

bool ParametersSet::boolParameter( const string & parameter ){

    bool result;

    iterator i = find( parameter );
    if ( i == end() ) {
        cerr << "missing parameter : \"" << parameter << "\" in the set "
        << name << endl;
        exit(EXIT_FAILURE);
    }
    string stringResult( ( * i ).second.first );
    transform( stringResult.begin(), stringResult.end(), stringResult.begin(), (int(*)(int))tolower );

    if ( ( stringResult == "yes" )  || ( stringResult == "y" ) ||
    ( stringResult == "true" ) || ( stringResult == "t" ) ) {
        result = true;
    }
    else {
        if ( ( stringResult == "no" ) || ( stringResult == "n" ) ||
        ( stringResult == "false" ) || ( stringResult == "f" ) ) {
            result = false;
        }
        else {
            cerr << "The parameter : \"" << parameter << "\" in the set "
            << name << " can not be interpreted as a boolean,"
            << " try yes/no values" << endl;
            exit(EXIT_FAILURE);
        }
    }
    (*i).second.second = true;
    return result;
}

ParametersSet & ParametersSet::operator () ( const string & categoryName ){
    for ( vector < pair<ParametersSet, bool> >::iterator i =
    categories.begin(); i != categories.end(); ++i ) {
        if ( ( * i ).first.name == categoryName ) {
            (*i).second = true;
            return ( (*i).first );
        }
    }
    categories.push_back(
            pair< ParametersSet, bool >(ParametersSet(categoryName), false) );
    return ( (categories.back()).first );
}


string & ParametersSet::operator [] ( const string & fieldName ){
    map< string, pair<string,bool> >::operator[](fieldName) = pair< string, bool >("",false);
    return map< string, pair<string,bool> >::operator[](fieldName).first;
}

const string & ParametersSet::getName() const {
    return name;
}

void ParametersSet::setName( const string& newName ) {
    name = newName;
}

void ParametersSet::saveToFile( ofstream & fileOutput, int indentation ) const{
    string indent;
    indent.resize( indentation, ' ' );
    for (const_iterator i = begin(); i != end(); ++i ){
        fileOutput << indent << (*i).first << " = " << (*i).second.first << endl;
    }
    for ( vector< pair<ParametersSet, bool> >::const_iterator i = categories.begin();
          i != categories.end(); ++i ){
        fileOutput << indent << "{" << (*i).first.getName() << "}" << endl;
        (*i).first.saveToFile( fileOutput, indentation+2 );
        fileOutput << indent << "{\\" << (*i).first.getName() << "}" << endl << endl;
    }
}

bool ParametersSet::checkAllUsed() const{
    bool res = true;
    for (const_iterator i = begin(); i != end(); ++i ){
        if ( !((*i).second.second) ){
            res = false;
            cerr << "Field " << (*i).first << " not used in " << name << '.' << endl;
        }
    }
    for ( vector < pair<ParametersSet, bool> >::const_iterator i = categories.begin();
            i != categories.end(); ++i ) {
        if ( !((*i).second) ){
            res = false;
            cerr << "Category " << (*i).first.name << " not used in "
                 << name << '.' << endl;
        }
        res = res && (*i).first.checkAllUsed();
    }
    return res;
}


void ParametersSet::touch(){
    for (iterator i = begin(); i != end(); ++i ){
        (*i).second.second = true;
    }
    for ( vector < pair<ParametersSet, bool> >::iterator i = categories.begin();
            i != categories.end(); ++i ) {
        (*i).second = true;
        (*i).first.touch();
    }
}

bool ParametersSet::touch( const string & parameter ){
    iterator i = find( parameter );
    if ( i != end() ) {
        (*i).second.second = true;
        return true;
    }
    //not found
    return false;
}

bool ParametersSet::touchCategory( const string & categoryName ){
    //find the category
    for ( vector < pair<ParametersSet, bool> >::iterator i =
    categories.begin(); i != categories.end(); ++i ) {
        //once the category is found, "touch" recursively
        if ( ( * i ).first.name == categoryName ) {
            (*i).second = true;
            ( * i ).first.touch();
            return true;
        }
    }
    //not found
    return false;
}
