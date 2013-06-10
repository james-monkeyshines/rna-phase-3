#include "Util/array_bag.h"

#include <stdlib.h>

// Constructor
array_bag::array_bag() {
    unique_arrays = 0 ;
    bag_size = 100 ;
    arrays = new vector < int > [bag_size] ;
    number_of_arrays = new int[bag_size] ;
}

array_bag::~array_bag() {
    if ( arrays != NULL )
        delete[] arrays ;
    if ( number_of_arrays != NULL )
        delete[] number_of_arrays ;
}

// array_bag methods
void array_bag::expand() {
    vector < int > * newarrays = new vector < int > [2 * bag_size] ;
    int    * new_number_of_arrays = new int[2 * bag_size] ;
    int i ;

    for ( i = 0 ; i < bag_size ; i++ ) {
        newarrays[i] = arrays[i] ;
        new_number_of_arrays[i] = number_of_arrays[i] ;
    }
    delete[] arrays ;
    delete[] number_of_arrays ;
    arrays = newarrays ;
    number_of_arrays = new_number_of_arrays ;
    bag_size *= 2 ;
    return ;
}


void array_bag::add( vector < int > & thevector ) {
    int i = 0 ;
    bool found = false ;
    while ( ( !found ) && i < unique_arrays )
        found = ( thevector == arrays[i++] );
    if ( found )
        number_of_arrays[i - 1] ++ ;
    else {
        if ( ( unique_arrays + 1 ) > bag_size )
            expand() ;
        arrays[unique_arrays] = thevector ;
        number_of_arrays[unique_arrays] = 1 ;
        unique_arrays++ ;
    }
    return ;
}

/*
void array_bag::add(vector <int> thearray) {
int i=0 ;
bool found = false ;
while ((!found) && i<unique_arrays)
found = (thevector == arrays[i++]);
if (found)
number_of_arrays[i-1]++ ;
else {
if ((unique_arrays+1) > bag_size)
expand() ;
arrays[unique_arrays] = thevector ;
number_of_arrays[unique_arrays] = 1 ;
unique_arrays++ ;
}
return ;
}
*/

void array_bag::print( ostream & out ) {
    sort() ;
    // cout << "Unique arrays = " << unique_arrays << endl ;
    for ( int i = 0 ; i < unique_arrays ; ++i ) {
        for ( unsigned int j = 0 ; j < arrays[i].size() ; ++j ) {
            if ( arrays[i] [j] > 0 )
                out << "*" ;
            else
                out << "-" ;
        }
        out << "\t" << number_of_arrays[i] << endl ;
    }
}

void array_bag::sort() {
    int i = 0 , j ;
    int swapn ;
    vector < int > swapvector;
    //  bool swap ;
    // bubble sort
    while ( i < unique_arrays ) {
        for ( j = i + 1 ; j < unique_arrays ; j++ ) {
            if ( number_of_arrays[i] < number_of_arrays[j] ) {
                swapn = number_of_arrays[i] ;
                swapvector = arrays[i] ;
                number_of_arrays[i] = number_of_arrays[j] ;
                arrays[i] = arrays[j] ;
                number_of_arrays[j] = swapn ;
                arrays[j] = swapvector;
            }
        }
        i++ ;
    }
    return ;
}

int array_bag::find( vector < int > & search_vector ) {
    int i ;
    for ( i = 0 ; i < unique_arrays ; i++ ) {
        if ( arrays[i] == search_vector )
            return i ;
    }
    return -1 ;
}

int array_bag::number_of_unique_arrays() {
    return unique_arrays ;
}

