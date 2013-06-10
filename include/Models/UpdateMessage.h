#ifndef UPDATEMESSAGE_H
#define UPDATEMESSAGE_H

#include <assert.h>

class UpdateMessage{


public:

    unsigned int message;

    typedef struct{
        unsigned int symbolCategory;
        unsigned int model;
    } ModelMsg;

//    typedef struct{
//    } ParamMsg;

    ModelMsg modelMsg;

    typedef enum{
        SYMBOL_CATEGORY = 0x01,
        MODEL_FLAG = 0x02,
        MODEL_PRIOR_FLAG = 0x40
    } ModelFlags;

    typedef enum{
        FREQ = 0x0001,
        RATE = 0x0002,
        GAMMA = 0x0004,
        INVARIANT = 0x0008,
        AVERAGE_RATE = 0x0010,
        BRANCH = 0x0100,
        HYPER_PARAM = 0x1000,
        PRIOR_FLAG = 0x2000
    } ParamFlags;

    typedef enum{
        MODEL_TYPE = 0x8000,
        PARAM_TYPE = 0x4000,
    } TypeMask;

    UpdateMessage() {
        clear();
    }

    ~UpdateMessage(){
    }

    inline void clear(){
        message = 0;
    }

    inline void setFlag( ParamFlags flag ){
        assert(this->message&PARAM_TYPE);
        this->message |= flag;
    }
    inline void setFlag( ModelFlags flag ){
        assert(this->message&MODEL_TYPE);
        this->message |= flag;
    }

    inline void unsetFlag( ParamFlags flag ){
        assert(this->message&PARAM_TYPE);
        this->message &= ~flag;
    }
    inline void unsetFlag( ModelFlags flag ){
        assert(this->message&MODEL_TYPE);
        this->message &= ~flag;
    }

    inline bool hasFlag( ParamFlags flag ){
        assert(this->message&PARAM_TYPE);
        return this->message & flag;
    }
    inline bool hasFlag( ModelFlags flag ){
        assert(this->message&MODEL_TYPE);
        return this->message & flag;
    }

    inline bool hasType( TypeMask type ){
        return (this->message&type);
    }

    inline void setType( TypeMask type ){
        assert(message==0);
        this->message |= type;
    }
};

#endif
