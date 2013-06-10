#include "Tree/Cluster.h"
#include "Tree/BasicNode.h"

#include <algorithm>
#include <assert.h>

Cluster::Cluster( unsigned int totalNumberSpecies ){
    this->totalNumberSpecies = totalNumberSpecies;
    clusterMask.resize( 1 +
            (totalNumberSpecies - 1) / (8*sizeof(unsigned int)) );
    fill(clusterMask.begin(), clusterMask.end(), 0);
    numberSpecies = 0;
    key = 0;
    commonAncestor = NULL;
}

Cluster::Cluster( unsigned int totalNumberSpecies, const string & name ){
    this->totalNumberSpecies = totalNumberSpecies;
    clusterMask.resize( 1 +
            (totalNumberSpecies - 1) / (8*sizeof(unsigned int)) );
    fill(clusterMask.begin(), clusterMask.end(), 0);
    numberSpecies = 0;
    key = 0;
    this->name = name;
    commonAncestor = NULL;
}

Cluster::Cluster(const Cluster & cluster1, const Cluster & cluster2){
    vector<unsigned int> mask = cluster1.getMask();
    vector<unsigned int>::iterator iter = mask.begin();
    vector<unsigned int>::const_iterator iterOther = cluster2.getMask().begin();
    while(iter!=mask.end()){
        *iter = *iter | *iterOther;
        ++iter;
        ++iterOther;
    }
    // this final call will initialise the mask, the key and the number
    // of species
    setMask( mask );
}

void Cluster::getSpecies( vector<unsigned int> & species ) const{
    species.clear();
    for (unsigned int i = 0; i < clusterMask.size(); ++i ){
        unsigned int x = clusterMask[i];
        unsigned int id = 8*sizeof(unsigned int)*i;
        while (x){
            if ( x & 0x01){
                species.push_back(id);
            }
            x>>=1;
            ++id;
        }
    }
}

void Cluster::join(const Cluster & other){
    vector<unsigned int> mask = clusterMask;
    vector<unsigned int>::iterator iter = mask.begin();
    vector<unsigned int>::const_iterator iterOther = other.getMask().begin();
    while(iter!=mask.end()){
        *iter = *iter | *iterOther;
        ++iter;
        ++iterOther;
    }
    // this final call will initialise the mask, the key and the number
    // of species
    setMask( mask );
}

bool Cluster::contains( unsigned int speciesId ){
    unsigned int x = clusterMask[speciesId / (8*sizeof(unsigned int))];
    unsigned int mask = 1 << speciesId % (8*sizeof(unsigned int));
    if ( x & mask ){
        return true;
    }
    else{
        return false;
    }
}

bool Cluster::contains( const Cluster & other ){
    if (getNumberSpecies() < other.getNumberSpecies()){
        return false;
    }
    vector<unsigned int>::const_iterator iter = clusterMask.begin();
    vector<unsigned int>::const_iterator iterOther = other.getMask().begin();
    while (iter!=clusterMask.end()){
        /**
         * ~ = bitwise negation, this is the mask for the species NOT
         *                       in the cluster
         * if ( ~(*iter) & (*iterOther) != 0 ) one species in other is not in
         * this cluster.
         */
        if ( ~(*iter) & (*iterOther) ){
            return false;
        }
        ++iter;
        ++iterOther;
    }
    return true;
}

bool Cluster::intersect( const Cluster & other ){
    vector<unsigned int>::const_iterator iter = clusterMask.begin();
    vector<unsigned int>::const_iterator iterOther = other.getMask().begin();
    while (iter!=clusterMask.end()){
        /**
         * if ( (*iter) & (*iterOther) != 0 ) at least one species in other is
         * in this cluster.
         */
        if ( (*iter) & (*iterOther) ){
            return true;
        }
        ++iter;
        ++iterOther;
    }
    return false;
}

#include <iostream>
bool Cluster::compatible( const Cluster & other ){
    //if other has a species not in this AND this has a specie not in other
    //then clusters are not compatible; unless completly disjoint
    if (!intersect(other)){
        return true;
    }
    bool otherDiff = false;
    bool thisDiff = false;

    vector<unsigned int>::const_iterator iter = clusterMask.begin();
    vector<unsigned int>::const_iterator iterOther = other.getMask().begin();
    while (iter!=clusterMask.end()){
       /**
         * if ( ~(*iter) & (*iterOther) != 0 ) one species in other is not in
         * this cluster.
        */
        if ( ~(*iter) & (*iterOther) ){
            if (thisDiff) return false;
            otherDiff = true;
        }
        if ( (*iter) & ~(*iterOther) ){
            if (otherDiff) return false;
            thisDiff = true;
        }
        ++iter;
        ++iterOther;
    }
    return true;
}

void Cluster::setMask( const vector<unsigned int> & newMask ){
    clusterMask = newMask;
    //count the number of species
    numberSpecies = 0;
    key = 0;
    vector<unsigned int>::const_iterator iter = clusterMask.begin();
    while (iter!=clusterMask.end()){
        unsigned int x = *iter;
        key ^= *iter;
        while (x){
            if ( x & 0x01 ){
                ++numberSpecies;
            }
            x>>=1;
        }
        ++iter;
    }
}

bool Cluster::addSpecies( unsigned int speciesId ){
    unsigned int element = speciesId / (8*sizeof(unsigned int));
    unsigned int pos = speciesId % (8*sizeof(unsigned int));
    unsigned int partialMask = 0x01 << pos;
    bool already = clusterMask[element] & partialMask;
    if (already){
        return false;
    }
    //key = key XOR partialMask
    key ^= partialMask;
    //add the species to the mask
    clusterMask[element] |= partialMask;
    ++numberSpecies;
    return true;
}


void Cluster::invert(){
    key = 0;
    numberSpecies = totalNumberSpecies - numberSpecies;
    //for the last element of clusterMask, we have to account for
    //supplementary 1
    //declare lastMask to know how the last element of clusterMask is used...
    unsigned int lastMask = getLastMask();
    for ( vector<unsigned int>::iterator iter = clusterMask.begin();
          iter != clusterMask.end(); ++iter ){
        *iter= ~(*iter);
        if (&(*iter) == &(clusterMask.back())){
            *iter &= lastMask;
        }
        key ^= *iter;
    }
}

void Cluster::invert( const Cluster& outgroup ){
    assert( outgroup.totalNumberSpecies == totalNumberSpecies );
    assert( !intersect( outgroup ) );
    key = 0;
    numberSpecies = totalNumberSpecies - outgroup.numberSpecies - numberSpecies;
    //for the last element of clusterMask, we have to account for
    //supplementary 1
    //declare lastMask to know how the last element of clusterMask is used...
    unsigned int lastMask = getLastMask();
    vector<unsigned int>::iterator iter = clusterMask.begin();
    vector<unsigned int>::const_iterator iterOther = outgroup.getMask().begin();
    while (iter != clusterMask.end()){
        *iter= ~(*iter) & ~(*iterOther);
        if (&(*iter) == &(clusterMask.back())){
            *iter &= lastMask;
        }
        key ^= *iter;
        ++iter;
        ++iterOther;
    }
}

void Cluster::remove( const Cluster& other ){
    vector<unsigned int> mask = clusterMask;
    vector<unsigned int>::iterator iter = mask.begin();
    vector<unsigned int>::const_iterator iterOther = other.getMask().begin();
    while(iter!=mask.end()){
        *iter = *iter & ~(*iterOther);
        ++iter;
        ++iterOther;
    }
    unsigned int lastMask = getLastMask();
    mask.back() &= lastMask;
    // this final call will initialise the mask, the key and the number
    // of species
    setMask( mask );
}


unsigned int Cluster::getLastMask(){
    unsigned int lastMask = 0;
    for ( unsigned int i = 0;
          i < 1+((totalNumberSpecies-1)%(8*sizeof(unsigned int))); ++i){
        lastMask<<=1;
        lastMask |= 0x01;
    }
    return lastMask;
}

void Cluster::clear(){
    fill(clusterMask.begin(), clusterMask.end(), 0);
    numberSpecies = 0;
    key=0;
}

bool Cluster::removeSpecies( unsigned int speciesId ){
    unsigned int element = speciesId / (8*sizeof(unsigned int));
    unsigned int pos = speciesId % (8*sizeof(unsigned int));
    unsigned int partialMask = 0x01 << pos;
    bool exists = clusterMask[element] & partialMask;
    if (exists){
        //key = key XOR partialMask
        key ^= partialMask;
        clusterMask[element] = clusterMask[element] & ~partialMask;
        --numberSpecies;
        return true;
    }
    return false;
}


unsigned int Cluster::intersect( const Cluster & other, Cluster& intersectCluster ){
    vector<unsigned int>::const_iterator iter = clusterMask.begin();
    vector<unsigned int>::const_iterator iterOther = other.getMask().begin();
    vector<unsigned int> intersectMask;
    while (iter!=clusterMask.end()){
        /**
         * if ( (*iter) & (*iterOther) != 0 ) at least one species in other is
         * in this cluster.
         */
        intersectMask.push_back( (*iter) & (*iterOther) );
        ++iter;
        ++iterOther;
    }
    intersectCluster.setMask( intersectMask );
    return ( intersectCluster.getNumberSpecies() );
}
