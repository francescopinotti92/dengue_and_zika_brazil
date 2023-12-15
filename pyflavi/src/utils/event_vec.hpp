//
//  event_vec.hpp
//  CrossFlavivirus
//
//  Created by MacBook Pro on 05/01/2021.
//

#ifndef event_vec_h
#define event_vec_h

#include <vector>

template <typename T>
class event_vec {
public:
    event_vec( int capacity );
    void push(T el);
    void clear();
    typename std::vector<T>::iterator begin();
    typename std::vector<T>::iterator end();
    int get_size();
    
private:
    int size;
    int capacity;
    std::vector<T> v;
};

#endif /* event_vec_h */
