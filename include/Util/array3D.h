#ifndef ARRAY3D_H
#define ARRAY3D_H

#ifdef ARRAY_CHECK
#include <stdlib.h>
#include <assert.h>
#include <fstream>
#endif

#include <string.h>

template < class Z >
class array3D {
public:
    // Constructors
    array3D() {
        nslices = nrows = ncolumns = 0 ;
        block_slice = 0 ;
        data = NULL;
        block_index = NULL ;
    }

    array3D( int slice, int row, int columns ) {
        nslices  = slice ;
        nrows    = row ;
        ncolumns = columns ;
        data = new Z[nslices * nrows * ncolumns] ;
        block_slice = nrows * ncolumns ;
        block_index = new Z * [nslices] ;
        for ( int i = 0 ; i < nslices ; i++ )
            block_index[i] = & ( data[i * ( block_slice )] ) ;
    }

    // Destructor
    ~array3D() {
        if ( block_index != NULL )
            delete[] block_index ;
        if ( data != NULL )
            delete[] data;
    }

    array3D( const array3D < Z > & src ) {
        nslices  = src.nslices ;
        nrows    = src.nrows ;
        ncolumns = src.ncolumns ;
        data     = new Z[nslices * nrows * ncolumns] ;
        block_index = new Z * [nslices] ;
        memcpy( data , src.data , sizeof( Z ) * nslices * nrows * ncolumns ) ;
        block_slice = src.block_slice ;
        for ( int i = 0 ; i < nslices ; i++ )
            block_index[i] = & ( data[i * ( block_slice )] ) ;
    }

    inline int numberSlices() const {
        return nslices;
    }

    inline int numberRows() const {
        return nrows;
    }

    inline int numberColumns() const {
        return ncolumns;
    }

    inline Z operator() ( int s , int r , int c ) const {
#ifdef ARRAY_CHECK
        assert( s < nslices ) ;
        assert( r < nrows ) ;
        assert( c < ncolumns ) ;
#endif
        return ( block_index[s] [( r * ncolumns ) + c] ) ;
    }

    inline Z & operator() ( int s , int r , int c ) {
#ifdef ARRAY_CHECK
        assert( s < nslices ) ;
        assert( r < nrows ) ;
        assert( c < ncolumns ) ;
#endif
        return ( block_index[s] [( r * ncolumns ) + c] ) ;
    }

    array3D < Z > & operator = ( const array3D < Z > & src ) {
        int i ;
        if ( this == & src ) // No self assignment
            return ( * this ) ;

        if ( !( ( nslices == src.nslices ) && ( nrows == src.nrows ) &&
        ( ncolumns == src.ncolumns ) ) ) {
            if ( block_index != NULL )
                delete[] block_index ;
            if ( data != NULL )
                delete[] data ;
            nslices  = src.nslices ;
            nrows    = src.nrows ;
            ncolumns = src.ncolumns ;
            data = new Z[nslices * nrows * ncolumns] ;
            block_index = new Z * [nslices] ;
            block_slice = src.block_slice ;
            for ( i = 0 ; i < nslices ; i++ )
                block_index[i] = & ( data[i * ( block_slice )] ) ;
        }

        memcpy( data , src.data , sizeof( Z ) * nslices * block_slice ) ;
        return ( * this ) ;
    }

    void resize( int nslice , int nrow , int col ) {
        if ( ( nslice == nslices ) && ( nrow == nrows ) && ( ncolumns == col ) )
            return ;

        if ( block_index != NULL )
            delete[] block_index ;
        if ( data != NULL ) 	
            delete[] data ;	

        nslices  = nslice ;
        nrows    = nrow ;
        ncolumns = col ;
        data     = new Z[nslices * nrows * ncolumns] ;
        block_slice = nrows * ncolumns ;
        block_index = new Z * [nslices] ;
        for ( int i = 0 ; i < nslices ; i++ )
            block_index[i] = & ( data[i * ( block_slice )] ) ;
        return ;
    };

    bool operator == ( const array3D < Z > & src ) const{
        if (nslices == src.numberSlices()){
            if (nrows == src.numberRows()){
                if (ncolumns == src.numberColumns()){
                    return memcmp( this, &src, sizeof( Z ) * nslices * block_slice);
                }
            }
        }
        return false;
    }

    void reclaimSpace() {
        if ( block_index != NULL )
            delete[] block_index ;
        if ( data != NULL )
            delete[] data;
        nrows = ncolumns = nslices = 0 ;
    }

protected:
    int nrows;
    int ncolumns;
    int nslices;
    int block_slice ;
    Z * data ;
    Z * * block_index ;
};

#endif




