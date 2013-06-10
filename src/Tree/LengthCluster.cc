#include "Tree/LengthCluster.h"


LengthCluster::LengthCluster( unsigned int totalNumberSpecies ) :
Cluster( totalNumberSpecies ){
    length = 0.0;
    number = 0;
}

LengthCluster::LengthCluster( unsigned int totalNumberSpecies,
                     double initialLength ):Cluster( totalNumberSpecies ){
    length = initialLength;
    number = 1;
}

void LengthCluster::addParam( const Cluster& other ){
   // cout << "actual length " << length << endl;
   // cout << "actual number " << number << endl;
    length += ((LengthCluster*)(&other))->getLength();
    number += ((LengthCluster*)(&other))->getNumber();
  //  cout << "after length " << length << endl;
 //   cout << "after number " << number << endl;
}    

