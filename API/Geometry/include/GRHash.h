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
//! @file      GRhash.h
//! @brief     Hash table.
//! @author    Dr. Wang, Bell.

//
// James's hash table (RDouble hashing)
// class T cannot be int
// class T needs to have the following functions:
// hash_key();
// operator ==
// operator =
//
// Note: index returned starts from 1 !!!!
// hash size must be of power of 2

// The follows is noted by Bell.
// Hash function, h(k,i) = [h1(k) + i*h2(k)] mod m.
// ii=(hash1+i_prob*hash2)%hash_size --> Double Hashing.
// iProb: is the i, probing time.
// hash_size: m, always to be 2^r
// hash1: h1(k), hash function1, always to be ODD.
// hash2: h2(k), hash function2, always to be ODD.
// alpha: load factor = n / m < 1, where, n and m are the number of data and slot, respectively.
// hash_table: hash table, the position of the element, but the real element is stored in hash_deque.
// hash_deque: stores the real element. Because of deque STL is used, which is not a random store data-structure,
//             hash_table is used to store the position, a one-to-one map is build between 
//             hash_table and the position in hash_deque.
#include "PHMpi.h"
#include "Geo_Point.h"
#include "Geo_Face.h"
#include "TK_Exit.h"
#include "PHHeader.h"
#include "Geo_LineTopo_Unstruct.h"
using namespace std;
using namespace PHSPACE;

namespace GRIDER_SPACE
{

#define HASH_MAXI_LOAD 0.8
#define HASH_INIT_SIZE 128

#define NIL 0 //! Defined as reference <Introduction to Algorithms>.

//typedef unsigned int CAulong;
typedef            int CAulong;

CAulong gethashkey(const Point3D &p);
CAulong gethashkey(const Point3D *p);

CAulong gethashkey(const Geo_LineTopo_Unstruct &line);
CAulong gethashkey(const Geo_LineTopo_Unstruct *line);

CAulong gethashkey(const Geo_Face &face);
CAulong gethashkey(const Geo_Face *face);

bool equal_hash(const Point3D &p1, const Point3D &p2);
bool equal_hash(const Point3D *p1, const Point3D *p2);

bool equal_hash(const Geo_Face &face1, const Geo_Face &face2);
bool equal_hash(const Geo_Face *face1, const Geo_Face *face2);

bool equal_hash(const Geo_LineTopo_Unstruct &line1, const Geo_LineTopo_Unstruct &line2);
bool equal_hash(const Geo_LineTopo_Unstruct *line1, const Geo_LineTopo_Unstruct *line2);

template <class T>
class Hash {
public:
    Hash(CAulong n = HASH_INIT_SIZE)
    {
        data_size  = 0;
        hash_table = 0;
        hash_size  = n/2;    // because it is doubled in resize
        hash_deque = new deque<T>; resize();
    }
    ~Hash() { delete [] hash_table; delete hash_deque; }

    int insert(const T &p);
    int insert(const T &p, bool &isExist);
    int search(const T &p);
    int erase (const T &p);

    T & operator[](const T &p) { return ((*hash_deque)[insert(p)-1]); }
    T & operator[](CAulong i) const { return ((*hash_deque)[i]); }
    T & operator[](CAulong i)  { return ((*hash_deque)[i]); }

    //CAulong size() const { return ((*hash_deque).size()); }    //! It is not suitable for deleting case.
    CAulong size() const { return data_size; }

    RDouble  load() const { return ((*hash_deque).size()/RDouble(hash_size)); }

    void    clear() 
    {
        delete [] hash_table; hash_table = NULL; hash_size = HASH_INIT_SIZE;
        (*hash_deque).clear();
        delete hash_deque;
        hash_deque = new deque<T>;
        resize();
    }
    CAulong probes() const { return i_prob; }
    void info();

    Hash(const Hash<T> &H) {                    //! Assignment constructor.
        register CAulong i; CAulong sz = H.size();
        for (i = 0; i < sz; i ++)
            (*this).insert((*(H.hash_deque))[i]);
    }

    //     Hash<T>& operator=(const Hash<T>& H) {      //! Assign operator.
    //         clear();
    //         register int i; int sz=H.size();
    //         for (i=0; i<sz; i++)
    //             (*this).insert((*(H.hash_deque))[i]);
    //         return (*this);
    //     }    
    Hash<T> &operator = (const Hash<T> &H) {      //! Assign operator.
        clear();
        hash_size = H.hash_size;
        delete [] hash_table;
        hash_table = new CAulong[hash_size];
        register CAulong i; CAulong sz = H.size();

        for (i = 0; i < hash_size; i ++) hash_table[i] = H.hash_table[i];

        for (i = 0; i < sz; i ++)
            hash_deque->push_back(H[i]);
        return (*this);
    }
private:
    CAulong   data_size;
    CAulong   hash_size;
    CAulong  *hash_table;    //! Store the position of elements, value start from '1'.
    deque<T> *hash_deque;
    //    sstring     name;

    CAulong     key, hash1, hash2, o_size, i_prob;

    void resize();

    //Hash(const Hash<T>& H) {}    // disable copy constructor
    //void operator=(const Hash<T>& H) {}    // disable assignment operator
};

template <class T>
void Hash<T>::resize()
{
    register CAulong ii, jj, kk;

    o_size = static_cast<CAulong>((*hash_deque).size());
    bool *isDeleted = new bool [o_size];
    for (int i = 0; i < o_size; ++ i) isDeleted[i] = false;

    //! Mark out which elements have been deleted.
    if (hash_table)
    {
        for (int i = 0; i < hash_size; ++ i)
        {
            if (hash_table[i] < NIL)
            {
                int position = - hash_table[i];
                isDeleted[position - 1] = true;
            }
        }
        delete [] hash_table;
        hash_table = 0;
    }

    int numberOfElement = 0;
    for (int i = 0; i < o_size; ++ i)
    {
        if (!isDeleted[i]) ++ numberOfElement;
    }

    //! Resize the hash size, multiply by two.
    hash_table = new CAulong[hash_size*=2];
    for (ii = 0; ii < hash_size; ii ++) hash_table[ii] = 0;

    //************************************************
    //! Erase the deleted element from hash_deque.
    //************************************************
    //! Firstly, re-order the non-deleted elements.
    int iStart = 0, iEnd = o_size - 1;
    while (iStart <= iEnd)
    {
        if (isDeleted[iStart])
        {
            while (isDeleted[iEnd] && iStart <= iEnd)
            {
                -- iEnd;
            }

            if (!isDeleted[iEnd] && iStart < iEnd)
            {
                (*hash_deque)[iStart] = (*hash_deque)[iEnd];

                ++ iStart;
                -- iEnd;
            }
        }
        else
        {
            ++ iStart;
        }
    }
    
    if (numberOfElement != iStart)
    {
        ostringstream oss;
        oss << "Error: numberOfElement != iStart " << numberOfElement << "!=" << iStart << endl;
        TK_Exit::ExceptionExit(oss);
    }

    //! Secondly, Erase the left elements.
    int deleteDataNumber = o_size - numberOfElement;
    int count = 0;
    while (count < deleteDataNumber)
    {
        ++ count;
        hash_deque->pop_back();
    }
    ASSERT(hash_deque->size() == numberOfElement);

    //! Re-build the image: hash_table --> hash_deque.
    o_size = static_cast<CAulong>(hash_deque->size());
    for (jj = 0; jj < o_size; jj ++)
    {
        key = gethashkey((*hash_deque)[jj]);
        hash1 = key % hash_size;
        hash2 = key % (hash_size-1) | 1;
        kk = 0;
        while (hash_table[ii=((hash1+kk++*hash2)%hash_size)] != NIL);

        hash_table[ii] = jj+1;
    }

    delete [] isDeleted;
}

/*
This function is used when some data could be erased !!!
*/
template <class T>
int Hash<T>::insert(const T& p, bool &isExist)
{
    int index = search(p);
    if (!index)
    {
        isExist = false;
        index = insert(p);
    }
    else
    {
        isExist = true;
    }

    return index;
}

/*
This function is used when NO data is erased !!!!!
insert(): 
Insert the Type p into hash table,according the hash key.
E.g. when Type == Point3D, key is the computed from (x, y, z),
The point index is the same, when the coordinates is same for different points.
*/
template <class T>
int Hash<T>::insert(const T &p)
{
    key = gethashkey(p);
    hash1 = key % hash_size;
    hash2 = key % (hash_size-1) | 1;
    register CAulong ii;

    // ii=(hash1+i_prob*hash2)%hash_size --> Double Hashing.
    for (i_prob = 0; i_prob < hash_size; i_prob ++)
    {
        if (hash_table[ii=((hash1+i_prob*hash2)%hash_size)] == NIL)
        {
            // Insert
            (*hash_deque).push_back(p);    // empty slot, insert
            hash_table[ii] = o_size = static_cast<CAulong>((*hash_deque).size());
            if (o_size > hash_size*HASH_MAXI_LOAD) resize();
            ++ data_size;
            return (o_size);
        }
        else if (hash_table[ii] < NIL)
        {
            //****************************************************************************
            //!!This is only used when insert(p, isExist) is used!!!!!
            //!!The deleted position could be re-insert, when search() is used!!!!!
            //****************************************************************************

            //! When the element at this position has been deleted,
            //! it is permission to re-insert other one into this position.
            int position = - hash_table[ii];
            hash_table[ii] = position;    //! hash_table is re-marked to be usable.

            (*hash_deque)[position-1] = p;
            ++ data_size;
            return position;
        }
        else if (equal_hash((*hash_deque)[hash_table[ii]-1], p))
        {
            return (hash_table[ii]);    //! Already exist.
        }
    }

    return (0);    //! Insert failed.
}

template <class T>
int Hash<T>::erase(const T &p)
{
    key = gethashkey(p);
    hash1 = key % hash_size;
    hash2 = key % (hash_size-1) | 1;
    register CAulong ii;

    // ii=(hash1+i_prob*hash2)%hash_size --> Double Hashing.
    for (i_prob = 0; i_prob < hash_size; i_prob ++)
    {
        if (hash_table[ii=((hash1+i_prob*hash2)%hash_size)] == NIL)
        {
            //! Data not here, do not need to erase.
            return 0;
        } 
        else if (hash_table[ii] < NIL)
        {
            //! Deleted position, continue.
            //! Because element may be exist at the position after 'ii'.
            continue;
        }
        else if (equal_hash((*hash_deque)[hash_table[ii]-1], p))
        {
            /*
            T end = (*hash_deque).back();
            (*hash_deque)[hash_table[ii]-1] = end;

            key=gethashkey(end);
            hash1=key%hash_size;
            hash2=key%(hash_size-1)|1;
            register CAulong jj;
            CAulong j_prob;
            for (j_prob=0; j_prob<hash_size; j_prob++) 
            {
                if (hash_table[jj=((hash1+j_prob*hash2)%hash_size)]==0)
                {

                } 
                //else if ((*hash_deque)[hash_table[jj]-1]==end) 
                else if (equal_hash((*hash_deque)[hash_table[jj]-1], end)) 
                {
                    hash_table[jj] = hash_table[ii];
                    break;
                }
            }
            */

            hash_table[ii] *= -1;    // This position is marked as 'Deleted'
            -- data_size;

            /*
            (*hash_deque).pop_back();    // Delete the real element from deque container.
            */

            //! Data not here, do not need to delete.
            return hash_table[ii];
        }
    }

    return 0;    // delete failed
}

template <class T>
int Hash<T>::search(const T &p)
{
    key = gethashkey(p);
    hash1 = key % hash_size;
    hash2 = key % (hash_size-1) | 1;
    register CAulong ii;

    // ii=(hash1+i_prob*hash2)%hash_size --> Double Hashing.
    for (i_prob = 0; i_prob < hash_size; i_prob ++)
    {
        if (hash_table[ii=((hash1+i_prob*hash2)%hash_size)] == NIL)
        {
            return (0);    //! Does not exist.
        }
        else if (hash_table[ii] < NIL)
        {
            //! Deleted position, continue.
            //! Because element may be exist at the position after 'ii'.
            continue;
        }
        else if (equal_hash((*hash_deque)[hash_table[ii]-1], p))
        {
            return (hash_table[ii]);    //! Found.
        }
    }

    return (0);    //! Search failed.
}

template <class T>
void Hash<T>::info()
{
    if (size() == 0) return;

    T data;
    CAulong maxp = 0, avp = 0, p;

    for (CAulong i = 0; i < size(); i ++)
    {
        data = (*this)[i];
        search(data);
        p = probes();
        avp += p;
        if (maxp < p) maxp = p;
    }
}


typedef Hash<Point3D *> PointHash;
typedef Hash<Geo_Face *> FaceHash;
}

