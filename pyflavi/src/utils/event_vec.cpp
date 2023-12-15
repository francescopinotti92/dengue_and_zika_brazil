//
//  event_vec.cpp
//  CrossFlavivirus
//
//  Created by user on 31/07/2021.
//

#include "event_vec.hpp"


template <typename T>
event_vec<T>::event_vec(int capacity_): capacity(capacity_) {
    v.resize(0);
    v.resize(capacity);
    size = 0;
}

template <typename T>
void event_vec<T>::push(T el) {
    v[size] = el;
    size++;
}

template <typename T>
void event_vec<T>::clear() {
    size = 0;
}

template <typename T>
typename std::vector<T>::iterator event_vec<T>::begin() {
    return v.begin();
}

template <typename T>
typename std::vector<T>::iterator event_vec<T>::end() {
    return v.begin() + size;
}

template <typename T>
int event_vec<T>::get_size() {
    return size;
}


template class event_vec<int>;
template class event_vec<std::pair<int, int>>;
