/****************************************************\
        High Efficient / Hybrid Structure Container
                   By BELL, 2010.08.01
                   Modified 2012.08.09
 Note:  对于模板，在包含时，应该同时包括.h & .cpp文件
\****************************************************/

#include "HyList.h"
template < typename T > Sub_List<T>::Sub_List(int nNode)
{
    num  = 0;
    NP   = nNode;
    node = new T [NP];
    next = NULL;
}

template < typename T > Sub_List<T>::~Sub_List()
{
    delete [] node;
    num = 0;
}

//****************************************************************//

template < typename T > HyList<T>::HyList()
{
    AverNnode  = 8;
    num        = 0;
    head       = new Sub_List<T>(0);
    head->next = NULL;
}

template < typename T > HyList<T>::~HyList()
{
    Sub_List< T > *p = head->next;
    while (p)
    {
        Sub_List<T> *temp;
        if (p->next)
        {
            temp = p->next;
            delete p;
            p = temp;
        }
        else
        {
            delete p;
            break;
        }
    }
    delete head;
    num = 0;
}

template < typename T > void HyList<T>::push_back(const T node)
{
    int nNode = AverNnode, i;
    Sub_List<T> *p = head->next;

    if (!p)
    {
        Sub_List<T> *temp = new Sub_List<T>(nNode);
        p = temp;
        head->next = p;
    }
    else
    {
        while (p->next)
        {
            p = p->next;
        }

        if (p->num == p->NP)
        {
            Sub_List<T> *temp = new Sub_List<T>(nNode);
            p->next = temp;
            p = p->next;
        }
    }

    i = p->num;
    p->node[i] = node;
    p->num ++;
    num ++;
}

// Before insert a node to tail of the list, chech if it is exist.
template < typename T > void HyList<T>::insert(const T node)
{
    int nNode = AverNnode, i;
    Sub_List<T> *p = head->next;

    if (!p)
    {
        Sub_List<T> *temp = new Sub_List<T>(nNode);
        p = temp;
        head->next = p;
    }
    else
    {
        while (p->next)
        {
            p = p->next;
        }

        if (p->num == p->NP)
        {
            Sub_List<T> *temp = new Sub_List<T>(nNode);
            p->next = temp;
            p = p->next;
        }
    }

    bool flag = true;
    i = 0;
    while (flag && i < num)
    {
        T data = GetData(i++);
        if (data == node)
            flag = false;
    }
    if (flag)
    {
        i = p->num;
        p->node[i] = node;
        p->num ++;
        num ++;
    }
}

// Insert "node" to the position witch is before "iNode".
template < typename T > void HyList<T>::insert(const int iNode, const T node)
{
    exit(0);
}

template < typename T > void HyList<T>::erase(const int iNode)
{
    if (iNode < 0 || iNode >= num)
    {
        exit(0);
    }
}

template < typename T > T HyList<T>::GetData(const int iNode)
{
    T data;
    if (iNode < 0)
    {
        exit(0);
    }
    else
    {
        int num = 0, i;
        Sub_List<T> *p = head->next;
        num += p->NP;
        while (num <= iNode)
        {
            p = p->next;
            num += p->NP;
        }
        i = iNode - (num - p->NP);
        data = p->node[i];
    }

    return data;
}
template < typename T > int HyList<T>::size()
{
    return num;
}

template < typename T> T HyList<T>::operator[] (int i)
{
    return GetData(i);
}