//
// Created by Mingqiu Wang on 2/6/16.
//


#include "../include/Utility.h"

//template<typename T>
//bool find(std::vector<T> &a, T target) {
//    for (const auto &item : a)
//        if (item==target) return true;
//    else return false;
//}



bool find(std::vector<int> &a, int target) {
    for (const auto &item : a)
        if (item==target) return true;
    return false;
}

