#ifndef TREEMAP_H
#define TREEMAP_H

/**
 * The TreeMap object is an intermediate representation between the string
 * representation of a tree and the final tree.
 * The TreeMap contains other TreeMap according to the recursive definition
 * of a tree.
 */
 
#include <string>
#include <iostream>
#include <algorithm>
#include <vector>
#include <assert.h>

using namespace std;

class TreeMap {
public:
    /** ************************************************************************
     * TreeMap
     * @input      The string representation of a tree
     * @semantics  Constructor of the TreeMap
     ************************************************************************ */
    TreeMap( const string& stree );

protected:
    /** ************************************************************************
     * TreeMap
     * @input      The iterator begin() end() of the string representation
     *             of the tree
     * @semantics  Protected constructor of the TreeMap, this constructor is
     *             used to build an internal TreeMap, the two string::iterator
     *             are used to delimit the portion of the string to be used.
     *             This function is recursive.
     ************************************************************************ */
    TreeMap( string::const_iterator start, string::const_iterator finish );
    
    /** ************************************************************************
     * TreeMapConstruct
     * @input      The iterator begin() end() of the string representation
     *             of the tree
     * @semantics  Constructor primitive of the TreeMap, this function does
     *             the work
     ************************************************************************ */
    void treeMapConstruct( string::const_iterator start,
                           string::const_iterator finish );

public:
    /** ************************************************************************
     * printMap
     * @input      outputStream, where to output the object in a readable
     *             format
     * @semantics  for debug purpose, should reproduce the string used for
     *             construction
     ************************************************************************ */        
    void printMap( ostream& outputStream ) const;

    /** ************************************************************************
     * getLabel
     * @semantic      return the label
     * @precondition  is a leaf
     ************************************************************************ */
    inline const string& getLabel() const{
        return label;   
    }

    /** ************************************************************************
     * getNumberChildren
     * @semantic      return the number of children of the TreeMap
     ************************************************************************ */
    inline unsigned int getNumberChildren() const{
        return childMap.size();
    }
          
    /** ************************************************************************
     * getChildMap
     * @input          childId
     * @semantic       return the child childId
     * @precondition   0 <= childId < getNumberChildren()
     ************************************************************************ */
    inline const TreeMap& getChildMap( unsigned int childId ) const{
        assert( childId < childMap.size() );
        return childMap[childId].first;
    }
    
     /** ************************************************************************
     * getChildMap
     * @input          childId
     * @semantic       return the child childId
     * @precondition   0 <= childId < getNumberChildren()
     ************************************************************************ */
    inline TreeMap& getChildMap( unsigned int childId ){
        assert( childId < childMap.size() );
        return childMap[childId].first;
    }
   
    /** ************************************************************************
     * getDistance
     * @input          childId
     * @semantic       return the distance to the child childId
     * @precondition   0 <= childId < getNumberChildren()
     ************************************************************************ */
    inline double getDistance( unsigned int childId ) const{
        assert( childId < childMap.size() );
        return childMap[childId].second;
    }
    
    /** ************************************************************************
     * getParentDistance
     * @semantic       return the distance to the parent
     * @precondition   0 <= childId < getNumberChildren()
     ************************************************************************ */
    inline double getParentDistance() const{
        return parentDistance;
    }
    
    /** ************************************************************************
     * setParentDistance
     * @semantic       set the distance to the parent
     ************************************************************************ */
    inline void setParentDistance( double distance ){
        parentDistance = distance;
    }

protected:    
    double parentDistance;
    /** if this TreeMap is an internal node it has a set of child */
    vector< pair<TreeMap,double> > childMap;    
    /** if this TreeMap is a terminal node (leaf), it must have a label */
    string label;
};

#endif




