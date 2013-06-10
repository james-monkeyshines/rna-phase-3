#include "Tree/ClustersSet.h"

#include <iostream>
#include <assert.h>
#include <algorithm>

#include "Sequence/SequenceTable.h"

using namespace std;

ClustersSet::ClustersSet(unsigned int totalNumberSpecies){
    this->totalNumberSpecies = totalNumberSpecies;
}



ClustersSet::~ClustersSet(){
    for (multimap<unsigned int, pair<Cluster*, unsigned int> >::iterator i = begin(); i != end(); ++i){
        delete ( ((*i).second).first );
    }
}

bool ClustersSet::insert(Cluster* cluster){
    if ( ( cluster->getName() != "" ) &&
          getCluster(cluster->getName()) ){
          cerr << "Error: duplicated cluster name: " << cluster->getName() << endl;
          exit(EXIT_FAILURE);
    }
    //look for a similar cluster in the map
    pair< iterator, iterator > iters = equal_range( cluster->getKey() );
    for (iterator i = iters.first; i != iters.second; ++i){
        //if clusters are equal
        if ( (((*i).second).first)->matches( *cluster ) ){
            return false;
        }
    }
    multimap<unsigned int, pair<Cluster*, unsigned int> >::insert( value_type( cluster->getKey(), mapped_type(cluster, 1) ) );
    return true;
}


bool ClustersSet::add(const Cluster* cluster){
    if ( ( cluster->getName() != "" ) &&
          getCluster(cluster->getName()) ){
          cerr << "Error: duplicated cluster name: " << cluster->getName() << endl;
          exit(EXIT_FAILURE);
    }
    //look for a similar cluster in the map
    pair< iterator, iterator > iters = equal_range( cluster->getKey() );
    for (iterator i = iters.first; i != iters.second; ++i){
        //if clusters are equal
        if ( (((*i).second).first)->matches( *cluster ) ){
            ++( ((*i).second).second);
            (((*i).second).first)->addParam( *cluster );
            return false;
        }
    }
    multimap<unsigned int, pair<Cluster*, unsigned int> >::insert( value_type( cluster->getKey(), mapped_type(cluster->dup(), 1) ) );
    return true;
}

unsigned int ClustersSet::size(){
    unsigned int ret = 0;
    for (iterator i = begin(); i != end(); ++i){
        ret += ((*i).second).second;
    }
    return ret;
}

Cluster* ClustersSet::getCluster( const string& clusterName ) const{

    string temp;
    string name = clusterName;
    transform( name.begin(), name.end(), name.begin(), (int(*)(int)) tolower );

    const_iterator iter = begin();
    while (iter != end()){
        temp = ((iter->second).first)->getName();
        transform( temp.begin(), temp.end(), temp.begin(), (int(*)(int)) tolower );
        if ( temp == name ) return (iter->second).first;
        ++iter;
    }
    return ( NULL );

}


multimap< unsigned int, Cluster* > ClustersSet::listing() const{
    multimap< unsigned int, Cluster* > content;
    for (const_iterator i = begin(); i != end(); ++i){
        multimap< unsigned int, Cluster* >::value_type elt(
                ((*i).second).second, ((*i).second).first );
        content.insert( elt );
    }
    return content;
}

bool ClustersSet::compatible( const Cluster& cluster ) const{
    for ( ClustersSet::const_iterator iter = begin(); iter != end(); ++iter ){
        if ( !(*iter).second.first->compatible(cluster) ){
            return false;
        }
    }
    return true;
}
