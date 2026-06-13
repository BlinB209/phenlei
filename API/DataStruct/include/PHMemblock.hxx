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
//#include "FYMemblock.h"

namespace PHSPACE
{
// Null memory block for each (template) instantiation of MemoryBlockReference
template<typename P_type>
NullMemoryBlock<P_type> MemoryBlockReference<P_type>::nullBlock_;

template<typename P_type>
void MemoryBlock<P_type>::deallocate()
{
    if (!NumericTypeTraits<T_type>::hasTrivialCtor)
    {
        for (std::size_t i = 0; i < length_; ++ i)
        {
            data_[i].~T_type();
        }
        delete [] reinterpret_cast<char*>(dataBlockAddress_);
    }
    else
    {
        delete [] dataBlockAddress_;
    }
}

template<typename P_type>
inline void MemoryBlock<P_type>::allocate(size_t length)
{
    size_t numBytes = length * sizeof(T_type);

    if (numBytes < 1024)
    {
        dataBlockAddress_ = new T_type[length];
        data_ = dataBlockAddress_;
    }
    else
    {
        // We're allocating a large array.  For performance reasons,
        // it's advantageous to force the array to start on a
        // cache line boundary.  We do this by allocating a little
        // more memory than necessary, then shifting the pointer
        // to the next cache line boundary.

        // Patches by Petter Urkedal to support types with nontrivial
        // constructors.

        const int cacheBlockSize = 128;    // Will work for 32, 16 also

        dataBlockAddress_ = reinterpret_cast<T_type*>(new char[numBytes + cacheBlockSize - 1]);
        // Shift to the next cache line boundary

        //ptrdiff_t offset = ptrdiff_t(dataBlockAddress_) % cacheBlockSize;
        ptrdiff_t offset = reinterpret_cast<ptrdiff_t>(dataBlockAddress_) % cacheBlockSize;
        ptrdiff_t shift  = (offset == 0) ? 0 : (cacheBlockSize - offset);
        data_ = reinterpret_cast<T_type *>(reinterpret_cast<char *>(dataBlockAddress_) + shift);

        // Use placement new to construct types with nontrival ctors
        for (std::size_t i = 0; i < length; ++ i)
        {
            new(&data_[i]) T_type;
        }
    }
}

}

