#ifndef CLUSTERSSET_H
#define CLUSTERSSET_H

#include <iostream>
#include <fstream>
#include <map>

#include "Tree/Cluster.h"

class SequenceTable;

/**
 * a ClustersSet is a set of clusters
 * clusters are not repeated
 */


class ClustersSet : public multimap<unsigned int, pair<Cluster*, unsigned int> > {

protected:
    unsigned int totalNumberSpecies;

    /** ************************************************************************
     * insert
     * @return     false if the cluster was already in the set
     * @semantics  add the given cluster in the set, if already present
     *             return false
     ************************************************************************ */
    bool insert( Cluster* cluster);

public:
    /** ************************************************************************
     * ClustersSet
     * @semantics  default constructor
     ************************************************************************ */
    ClustersSet( unsigned int totalNumberSpecies );

    /** ************************************************************************
     * ClustersSet
     * @semantics  destructor, free the memory occupied by the clusters
     ************************************************************************ */
     ~ClustersSet();

    /** ************************************************************************
     * add
     * @return     false if the cluster was already in the set
     * @semantics  add a cluster in the set, if already present increase the
     *             counter, used to build a consensus
     ************************************************************************ */
    bool add(const Cluster* cluster);

    /** ************************************************************************
     * size
     * @return     the number of clusters in the set (counter values included)
     ************************************************************************ */
    unsigned int size();

    /** ************************************************************************
     * getCluster
     * @input      a cluster name
     * @return     the cluster
     ************************************************************************ */
    Cluster* getCluster( const string & name ) const;

    /** ************************************************************************
     * listing
     * @return     a multimap with all the clusters ranked in increasing order
     *             of times (ie number of times they have been added)
     * @semantics  useful to output all the clades found with the corresponding
     *             probability
     ************************************************************************ */
    multimap< unsigned int, Cluster* > listing() const;

    /** ************************************************************************
     * compatible
     * @return     true if compatible
     * @semantics  check whether the cluster provided is not in disagrement with
     *             clusters in the set
     ************************************************************************ */
    bool compatible( const Cluster& cluster ) const;
};

#endif //CLUSTERSSET_H

