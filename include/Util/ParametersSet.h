#ifndef PARAMETERSSET_H
#define PARAMETERSSET_H

#include <vector>
#include <string>
#include <map>
#include <fstream>

/*a class to parse files and used in many constructors with unknown parameters*/

using namespace std;

class ParametersSet : public map < string, pair < string, bool > > {
private:
    string name;
    vector < pair< ParametersSet, bool > > categories;
    
    //forbid the call to map::find
    inline iterator find( const string& parameter){
        return map < string, pair < string, bool > >::find( parameter );
    }
    inline const_iterator find( const string& parameter) const{
        return map < string, pair < string, bool > >::find( parameter );
    }
    
public:
    /** ********************************************************************
     * ParametersSet
     * @input      The name of the set to be printed on the screen in case of
     *             error
     * @semantics  Constructor of the set of parameter
     ************************************************************************ */
    ParametersSet( const string & setName );

    bool findParameter( const string& parameter ) const;
    
    bool findCategory( const string& categoryName ) const;
    
    bool boolParameter( const string & parameter );

    const string& stringParameter( const string & parameter );

    int intParameter( const string & parameter );
    
    double doubleParameter( const string & parameter );

    ParametersSet & operator() ( const string & categoryName );

    string & operator[] ( const string & fieldName );
            
    const string & getName() const;
    
    void setName( const string& newName );

    void saveToFile( ofstream & fileOutput, int indentation = 0 ) const;
    
    void touch();
    
    bool touch( const string & parameter );
    
    bool touchCategory( const string & categoryName );
           
    bool checkAllUsed() const;

};

#endif //PARAMETERSSET_H




