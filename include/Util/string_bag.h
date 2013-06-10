#ifndef STRING_BAG_H
#define STRING_BAG_H

#include <string>
#include <fstream>

using namespace std;

class string_bag {
public :
    string_bag(); // Constructor
    ~string_bag(); // Destructor
    // Methods
    void add( string & ) ;
    void print( ostream & ) ;
    int  find( string & ) ;
    int  number_of_entries() ;
protected :
    string * strings ;
    int    * number_of_strings ;
    int    unique_strings ;
    int bag_size ;
    void expand() ;
    void sort() ;
};

#endif //STRING_BAG_H



