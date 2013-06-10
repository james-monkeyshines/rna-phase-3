#ifndef FACTORY_H
#define FACTORY_H


#include <string>
#include <iostream>
#include <map>
#include <assert.h>

class Model;
class ParametersSet;

using namespace std;

template <class T>
class Factory {

private:
    typedef typename std::map< string, const T* > registryMap;
    typedef typename registryMap::iterator mapIterator;
    registryMap registry;

protected:

    /***************************************************************************
     * Factory
     * @semantics      a protected constructor to avoid instantiation without
     *                 the Singleton
     * @preconditions  T must define a virtual clone, this method should be
     *                 virtual and redefined in
     ************************************************************************ */
    Factory(){
        registry.clear();
    };

public:
    /***************************************************************************
     * subscribe
     * @input          the name of the new prototype
     * @semantics      register a new object prototype of the class T' in
     *                 Factory<T> ( T' inherits from T ).
     *                 T' should redefine the virtual method clone() of T to
     *                 return an object T'
     * @preconditions
     ************************************************************************ */
    void subscribe( const T * prototype, string name  ){
        assert( registry.find( name ) == registry.end() );
        registry[name] = prototype;
    }

    T * create( string name, ParametersSet & parameters ){
        mapIterator i = registry.find( name );
        if ( i == registry.end() ) {
            cerr << endl;
            cerr << "Attempt to create an unknown object: " << name << endl;
            cerr << "Please check your spelling (case-sensitive) against the list of registered objects:" << endl;
            for ( i = registry.begin(); i != registry.end(); ++i ){
                cerr << i->first << endl;
            }
            cerr << endl;
            return NULL;
        }
        else {
            return dynamic_cast<T*>(( ( * i ).second )->clone( parameters ));
        }
    }
};

#endif //FACTORY_H




