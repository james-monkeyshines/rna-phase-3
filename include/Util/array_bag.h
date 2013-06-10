#ifndef ARRAY_BAG_H
#define ARRAY_BAG_H


#include <vector>
#include <fstream>

using namespace std;

class array_bag {
public :
    array_bag() ; // Constructor
    ~array_bag() ; // Destructor

    void add( vector < int > & ) ;
    void print( ostream & )  ;
    int  find( vector < int > & )    ;
    int  number_of_unique_arrays() ;

protected :
    vector < int > * arrays ;
    int bag_size ;
    int * number_of_arrays ;
    int unique_arrays ;
    void expand() ;
    void sort() ;

};
#endif //ARRAY_BAG_H




