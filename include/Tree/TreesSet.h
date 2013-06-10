#ifndef TREESSET_H
#define TREESSET_H

#include <string>
#include <vector>
using namespace std;

template <class T>
class Array;

class TreeModelValue;

class TreesSet {
public :
    TreesSet() ;
    TreesSet( int numberEntries ) ;
    ~TreesSet() ;

    bool insert( const string & stringTree,
                 const vector < double > & branchLengths,
                 const vector < double > & modelParameters,
                 double maxLikelihood ) ;
    TreeModelValue * get( int index ) ;
    int size() ;
    void operator = ( const TreesSet & ) ;

protected:
    int numberEntries ;
    TreeModelValue * * entries ;
    int smallestIndex;
    int sofar;
};

#endif //TREESSET_H




