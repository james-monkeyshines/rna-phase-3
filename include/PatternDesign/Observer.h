#ifndef OBSERVER_H
#define OBSERVER_H

#include <vector>

/**
 * pattern design Observer
 * usage:
 *        define class MyObserver : Observer< MySubject >
 *        define class MySubject : Subject< MySubject >
 *        MyObserver can inherit from many Observer<T> to observe
 *        multiple classes.
 *        MySubject : Subject<Message1>, Subject<Message2>
 *        can send different messages as well
 */
template< class T >
class Observer{
public:
    Observer() {}
    virtual ~Observer() {}
    virtual void update(T* subject) = 0;    
};

template< class T >
class Subject : public T {
public:
    Subject() {}
    virtual ~Subject() {}
    void attach( Observer< T > & observer ){
        obs.push_back( &observer );
    }
    void detach( Observer< T > & observer ){
        obs.remove( &observer );
    }
    void notify(){
        for( typename std::vector< Observer<T>* >::iterator iter = obs.begin();
             iter != obs.end(); ++iter ){
            (*iter)->update(static_cast<T*> (this) );
        }
    }
private:
    std::vector< Observer<T>*  > obs;
};

#endif
