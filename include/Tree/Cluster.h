#ifndef CLUSTER_H
#define CLUSTER_H

#include <string>
#include <vector>

class BasicNode;

using namespace std;

/**
 * a Cluster is a set of species
 * Clusters are implemented with bitmask for efficiency
 */
class Cluster {
    string name;
    unsigned int totalNumberSpecies;
    vector<unsigned int> clusterMask;

    unsigned int numberSpecies;

    //the key is a XOR with the element of clusterMask, it is therefore not
    //unique but collision should be rare (2^32 different clades possible)
    //it is made for a preliminary test to save time
    unsigned int key;

    BasicNode * commonAncestor;

    double support;

public:

    /** ************************************************************************
     * Cluster
     * @semantics  default constructor
     ************************************************************************ */
    Cluster( unsigned int totalNumberSpecies );

    /** ************************************************************************
     * Cluster
     * @semantics  default constructor with a name
     ************************************************************************ */
    Cluster( unsigned int totalNumberSpecies, const string & name );

    /** ************************************************************************
     * Cluster
     * @semantics  union constructor
     ************************************************************************ */
    Cluster(const Cluster & cluster1, const Cluster & cluster2);

    /** ************************************************************************
     * ~Cluster
     * @semantics  virtual destructor
     ************************************************************************ */
    virtual ~Cluster(){}

     /** ************************************************************************
     * join
     * @semantics  union
     ************************************************************************ */
    void join(const Cluster & other);

    /** ************************************************************************
     * matches
     * @return  true if the clusters have the same set of species
     ************************************************************************ */
    inline bool matches( const Cluster & other ){
        if (getNumberSpecies() != other.getNumberSpecies()){
            return false;
        }
        return (other.getMask() == clusterMask );
    }

    /** ************************************************************************
     * getKey
     * @return  a (non-unique) key for the cluster
     ************************************************************************ */
    inline unsigned int getKey() const{
         return key;
    }

    /** ************************************************************************
     * getName
     * @return  the name given to the lustr at construction
     ************************************************************************ */
    inline const string& getName() const{
         return name;
    }

    /** ************************************************************************
     * contains
     * @return  true if this cluster contains the given species
     ************************************************************************ */
    bool contains( unsigned int speciesId );

    /** ************************************************************************
     * contains
     * @return  true if this cluster contains all the species in the other one
     ************************************************************************ */
    bool contains( const Cluster & other );

    /** ************************************************************************
     * remove
     * @semantics   remove the intersection with other from the current cluster 
     ************************************************************************ */
    void remove( const Cluster & other );

    /** ************************************************************************
     * intersect
     * @return  true if this cluster shares at least one species with
     *          the other one
     ************************************************************************ */
    bool intersect( const Cluster & other );

    /** ************************************************************************
     * intersect
     * @return  true if this cluster shares at least one species with
     *          the other one, return the cluster of shared species
     ************************************************************************ */
    unsigned int intersect( const Cluster & other, Cluster& intersectCluster );


    /** ************************************************************************
     * compatible
     * @return  true if clusters are not in disagreement
     ************************************************************************ */
    bool compatible( const Cluster & other );

    /** ************************************************************************
     * getNumberSpecies
     * @return  the number of speices in the cluster
     ************************************************************************ */
    inline unsigned int getNumberSpecies() const{
        return numberSpecies;
    }

     /** ************************************************************************
     * getSpecies
     * @return  a vector with the species ID
     ************************************************************************ */
    void getSpecies( vector<unsigned int> & species ) const;


    /** ************************************************************************
     * addSpecies
     * @return     true if the specie was not in the cluster
     * @semantics  add a species in the cluster
     ************************************************************************ */
    bool addSpecies( unsigned int speciesId );


    /** ************************************************************************
     * removeSpecies
     * @return  true if the specie was in the cluster
     * @semantics  remove a species from the cluster
     ************************************************************************ */
    bool removeSpecies( unsigned int speciesId );


     /** ************************************************************************
     * addParam
     * @input      another cluster
     * @semantics  descendant of Cluster will carry parameter with them
     *             (branch lengths, model parameters, ...)
     *             during a consensus we have to 'merge' the parameters of
     *             similar cluster, merge to this cluster the other one
     ************************************************************************ */
    virtual void addParam( const Cluster& other ){
        if (&other){}; //to avoid warning
    }

     /** ************************************************************************
     * dup
     * @semantics  duplique the current cluster to store it in a cluster set
     ************************************************************************ */
    virtual Cluster* dup() const{
        return new Cluster(*this);
    }

     /** ************************************************************************
     * invert
     * @semantics  invert the current cluster
     ************************************************************************ */
     void invert();

     /** ************************************************************************
     * invert
     * @input      an invariant cluster
     * @semantics  invert the current cluster relatively to outgroup species
     ************************************************************************* */
     void invert( const Cluster& other );

     /** ************************************************************************
     * clear
     * @semantics  empty the current cluster
     ************************************************************************ */
     void clear();

     /** ************************************************************************
     * getAncestor
     * @return  a pointer to the common ancestor
     ************************************************************************ */
    inline BasicNode* getAncestor() const{
        return commonAncestor;
    }

     /** ************************************************************************
     * setAncestor
     * @semantics  change/add a common ancestor to the cluster
     ************************************************************************ */
    inline void setAncestor( BasicNode* commonAncestor ){
        this->commonAncestor = commonAncestor;
    }


    /** ************************************************************************
     * getSupport
     * @return  the support
     ************************************************************************ */
    inline double getSupport() const{
        return support;
    }
    /** ************************************************************************
     * setSupport
     * @semantics   set the support variable
     ************************************************************************ */
    inline void setSupport(double newSupport){
        support = newSupport;
    }

    inline const vector<unsigned int> & getMask() const{
        return clusterMask;
    }

protected:
    /** ************************************************************************
     * setMask
     * @semantics   change the mask, key and number od species
     ************************************************************************ */
    void setMask( const vector<unsigned int> & newMask );

    /** ************************************************************************
     * getLastMask
     * @return      a mask with the species used in the last element of
     *              the mask
     ************************************************************************ */
    unsigned int getLastMask();


};

#endif //CLUSTER_H

