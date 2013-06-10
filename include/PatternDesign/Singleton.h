#ifndef SINGLETON_H
#define SINGLETON_H
/* pattern design Singleton */
/* usage : Singleton<yourClass>& yourinstance =                              */
/*                   Singleton<yourClass>::instance()                        */
template <class T>
class Singleton : public T
{
private:
    static Singleton* t;

/* private constructors to prevent misuse */
private:
    Singleton(){
    } 
    
    Singleton(const Singleton&){
    }
    
public:
    static Singleton& instance();
};

template<class T> Singleton<T>* Singleton<T>::t=0;
template<class T> Singleton<T>& Singleton<T>::instance(){
    if (!t) t=new Singleton();
    return *t;
}

#endif //SINGLETON_H
