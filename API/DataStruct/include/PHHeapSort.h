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
//! @file      PHHeapSort.h
//! @brief     It defines kernel data structure of PHengLEI, named 'heapSortInplace'.
//! @author    WanYunBo.

#pragma once
using namespace std;

namespace PHSPACE
{
struct sortIdxValue
{
    RDouble val;
    int idx;
};

#define COMPARE(a, b) ((a.val) < (b.val))
    //k/2  2*k  2*k+1
#define PARENT(k) ((k) >> 1)
#define LEFT__(k) ((k) << 1)
#define RIGHT_(k) (((k) << 1) | 0x1)

template< typename T >
class heapSortInplace
{
    //heap_max root >= child for k 2*k 2*k+1

private:
    T *_data;
    int _N;

public:
    heapSortInplace(T *data, int N) : _data(data - 1), _N(N) {}

    //private:
    void swim(int k) 
    {
        while ((k > 1) && COMPARE(_data[k], _data[PARENT(k)]))
        {
            SWAP(_data[PARENT(k)], _data[k]);
            k = PARENT(k);
        }
    }

    void sink(int k, int N)
    {
        while (LEFT__(k) <= N)
        {
            int largest_child_idx = LEFT__(k);

            if (largest_child_idx < N && COMPARE(_data[largest_child_idx + 1], _data[largest_child_idx])) ++largest_child_idx;

            if (!COMPARE(_data[largest_child_idx], _data[k])) break;

            SWAP(_data[k], _data[largest_child_idx]);
            k = largest_child_idx;
        }
    }

public:
    void sort()
    {
        int N = _N;    //get_size();
        for (int i = N / 2; i > 0; -- i)
        {
            sink(i, N);
        }

        while (N > 1)
        {
            SWAP(_data[1], _data[N]);
            --N;
            sink(1, N);
        }
    }

    int get_size() const { return _N; }

    void add(const T &data) { _data[++_N] = data; }
    void set(int i, const T &data) { _data[i] = data; }
    const T & get(int i) { return _data[i]; }
};

}
