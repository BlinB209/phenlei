//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      HyList.h
//! @brief     Efficient / Hybrid Structure Container.
//! @author    Bell.

#ifndef HYLIST_H
#define HYLIST_H

#include <cstdlib>
#include <iostream>
using namespace std;

//! Define the sub segment list of the HyList.
template < typename T >
class Sub_List
{
public: 
    int      NP;    //! Average No. point of sublist.
    int      num;
    T        *node;
    Sub_List *next;
    Sub_List(int nNode);
    ~Sub_List();
};


//! Define the HyList.
template < typename T > 
class HyList
{
public:
    int AverNnode;
    int num;
    Sub_List<T> *head;
    Sub_List<T> *sublist;

    //! SetAverNnode :initialize AverNnode with the input value in.
    void SetAverNnode(int in) { AverNnode = in; }

    //! push_back :add new item.
    void push_back(const T node);    //! Don't deal with overlap values.

    //! insert    :insert an item.
    void insert(const T node);       //! Deal with overlap value. 

    //! insert    :invoking exit(0).
    void insert(const int iNode, const T node);    //! Deal with overlap value. 

    //! erase     :invoking exit(0).
    void erase(const int iNode);    //! Delete the iNode th node.

    //! GetData   :get data from a Sub_List object's node property.
    T GetData(const int iNode);

    //! size      :return  num's value.
    int size();

    //! operator []: invoking GetData function.
    T operator [] (int i);
    HyList();
    ~HyList();
};

#include "HyList.hxx"
#endif