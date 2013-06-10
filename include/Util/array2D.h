#ifndef ARRAY2D_H
#define ARRAY2D_H

#ifdef ARRAY_CHECK
#define _USE_STRING_INLINES
#include <assert.h>
#include <fstream>
using namespace std;
#endif


// the C string header needed memcpy
#include <string.h>

#include <stdlib.h>

// Generic 2-dimensional array
template < class Q >
class array2D {
public:
    // Constructors
    array2D() {
        nrows  = ncolumns = 0 ;
        data   = NULL ;
        access = NULL ;
    };

    array2D(unsigned int numbersRows, unsigned int numbersColumns  ) {
        Q * p ;

        nrows     = numbersRows ;
        ncolumns  = numbersColumns ;
        data      = new Q[nrows * ncolumns];
        access    = new Q * [nrows] ;

        p = data ;
        for ( unsigned int i = 0 ; i < nrows ; i++ ) {
            access[i] = p ;
            p += ncolumns ;
        }
        p = NULL ;
    };

    array2D( const array2D< Q > & src ) {
        Q * p ;

        nrows    = src.nrows ;
        ncolumns = src.ncolumns ;
        data = new Q[nrows * ncolumns];
        access    = new Q * [nrows] ;

        p = data ;
        for ( unsigned int i = 0 ; i < nrows ; i++ ) {
            access[i] = p ;
            p += ncolumns ;
        }

        memcpy( data , src.data , sizeof( Q ) * nrows * ncolumns ) ;
    }

    // Destructor
    ~array2D() {
        //      cout << "Deleting 2d array" << endl ;
        if ( data != NULL )
            delete[] data;
        
        if ( access != NULL )
            delete[] access ;


        //      cout << "Deleted 2d array" << endl ;
    }

    // Members
    void resize( unsigned int nrow , unsigned int col ) {
        if ( ( nrow == nrows ) && ( ncolumns == col ) )
            return ;
        if ( access != NULL )
            delete[] access ;
        if ( data != NULL )
            delete[] data ;

        nrows    = nrow ;
        ncolumns = col ;
        data     = new Q[nrows * ncolumns] ;
        access = new Q * [nrows] ;
        for ( unsigned int i = 0 ; i < nrows ; i++ )
            access[i] = & ( data[i * ncolumns] ) ;
        return ;
    };

    unsigned int numberRows() const {
        return nrows;
    };

    unsigned int numberColumns() const {
        return ncolumns;
    };

    inline Q operator() ( unsigned int i , unsigned int j ) const {
#ifdef ARRAY_CHECK
        assert( i < nrows ) ;
        assert( j < ncolumns ) ;
#endif
        return ( access[i] [j] ) ;
    }


    inline Q & operator() ( unsigned int i , unsigned int j ) {
#ifdef ARRAY_CHECK
        assert( i < nrows ) ;
        assert( j < ncolumns ) ;
#endif
        return ( access[i] [j] ) ;
    }


    inline array2D< Q > & operator = ( const array2D< Q > & src ) {

        if ( this == & src ) // No self assignment
            return ( * this ) ;

        if ( ( nrows != src.nrows ) || ( ncolumns != src.ncolumns ) ) {
            if ( data != NULL ) {
                delete[] data ;
                delete[] access ;
            }
            nrows    = src.nrows ;
            ncolumns = src.ncolumns ;
            data = new Q[nrows * ncolumns] ;
            access = new Q * [nrows] ;
            for ( unsigned int i = 0 ; i < nrows ; i++ )
                access[i] = & ( data[i * ncolumns] ) ;
        }


        memcpy( data , src.data , sizeof( Q ) * nrows * ncolumns ) ;
        return ( * this ) ;
    }

    inline bool operator == ( const array2D< Q > & src ) {
        if ( this == & src )
            return ( true );
        if ( (nrows!=src.nrows) || (ncolumns!=src.ncolumns) ) {
            return false;
        }
        return (memcmp( data , src.data , sizeof( Q ) * nrows * ncolumns )==0);
    }

    void reclaim_space() {
        if ( data != NULL ) {
            delete[] access ;
            delete[] data;
        }
        nrows = ncolumns = 0 ;
    }

protected:
    unsigned int nrows;
    unsigned int ncolumns;
    Q * data ;
    Q * * access ;
};


#endif




