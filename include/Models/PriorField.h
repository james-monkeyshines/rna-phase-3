#ifndef PRIORFIELD_H
#define PRIORFIELD_H

#include <string>
#include <vector>
#include <assert.h>

using namespace std;

class PriorField{

    bool cst;
    
    typedef struct _HyperPrior{
        string priorName;
        vector<PriorField*> fields;
    } HyperPrior;
    
    union {
        double value;
        HyperPrior* hyperPrior;
    };
    
public:
    PriorField( double value );

    PriorField( string field );

    ~PriorField();
    
    inline const string& name() const{
        assert(!cst);
        return hyperPrior->priorName;
    }
      
    string toString() const;
    
    inline bool isConstant( unsigned int index ) const{
        assert(!cst);
        assert(index<hyperPrior->fields.size());
        return hyperPrior->fields[index]->cst;
    }
    
    inline double getValue( unsigned int index ) const{
        assert(index<hyperPrior->fields.size());
        assert(hyperPrior->fields[index]->cst);
        return hyperPrior->fields[index]->value;
    }
    
    inline const PriorField& getHyperPriorField( unsigned int index ) const{
        assert(!cst);
        assert(index<hyperPrior->fields.size());
        return *(hyperPrior->fields[index]);
    }
    
    inline unsigned int getNumberParameters() const{
        return hyperPrior->fields.size();
    }    
};

#endif
