/* Copyright (C)
 * 2019 - Hu Ren, rh890127a@163.com
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */
/**
 * @file utilityContainer.hpp
 * @brief, basic container warrper and useful array manipulations
 * @author Hu Ren, rh890127a@163.com
 * @version v0.1
 * @date 2019-08-16
 */

#ifndef UTILITY_UTILITYCONTAINER_HPP
#define UTILITY_UTILITYCONTAINER_HPP

#if defined(_WIN32)
#define NOMINMAX
#endif

#include <typeinfo>

#include "utilityBasicFunction.h"
#include "utilityType.h"
#include "utilityUsingCpp.hpp"

namespace UTILITY {
/**
 * wrap some containers to UTILITY
 */
#define Word string
#define Array vector
#define List set
#define Table map
#define MultiList multiset
#define MultiTable multimap
#define HASHMAP unordered_map
#define HASHLIST unordered_set

/**
 * conversion from string to Word
 */
#define toWord to_string
#define w2f stof
#define w2d stod
#define w2ld stold
#define w2i stoi
#define w2l stol
#define w2ll stoll

/**
 * @brief count the reference pointers
 */
class RefCounted {
 private:
  int count_;  ///< the count of reference pointers
 public:
  /**
   * @brief constructor
   */
  RefCounted() : count_(1){};
  /**
   * @brief get the count
   */
  int GetRefCount() const { return count_; }
  /**
   * @brief add a reference pointer
   */
  void incRefCount() { count_++; }
  /**
   * @brief substract a reference pointer
   * @return the count of pointers
   */
  int DecRefCount() { return --count_; }
};

/**
 * @brief basic struct of AoS data struct
 * @tparam T label or scalar
 */
template <typename T>
class BasicElement {
 public:
  label num;
  T *data;
  BasicElement() {
    data = NULL;
    num = 0;
  };
  ~BasicElement() {
    num = 0;
    data = NULL;
  }
  label size() { return num; }
  T &operator[](const int i) { return data[i]; }
};

template <class T>
bool compareArray(const Array<T> &a, const Array<T> &b);

/**
 * @brief AoS data struct
 * @tparam T label or scalar
 */
template <typename T>
// class ArrayArray : public RefCounted {
class ArrayArray {
 public:
  label num;                  ///< size of structs
  label *startIdx;            ///< start index of structs
  T *data;                    ///< structs
  RefCounted *refCount;       /// count of reference pointers
  BasicElement<T> *basicEle;  /// basic elements
  /**
   * @brief constructor
   */
  ArrayArray() {
    refCount = new RefCounted();
    startIdx = NULL;
    basicEle = NULL;
    data = NULL;
    num = 0;
    // par_std_out("constructor: %d\n", refCount->GetRefCount());
  };
  /**
   * @brief copy constructor
   */
  ArrayArray(const ArrayArray<T> &arr) {
    this->num = arr.num;
    this->startIdx = arr.startIdx;
    this->data = arr.data;
    this->refCount = arr.refCount;
    this->basicEle = arr.basicEle;
    refCount->incRefCount();
    // par_std_out("copy constructor: %d\n", refCount->GetRefCount());
  };
  /**
   * @brief operator overload: =
   */
  ArrayArray<T> &operator=(const ArrayArray &arr) {
    // ArrayArray<T> tmp;
    this->num = arr.num;
    this->startIdx = arr.startIdx;
    this->basicEle = arr.basicEle;
    this->data = arr.data;
    this->refCount = arr.refCount;
    this->refCount->incRefCount();
    return *this;
  };
  /**
   * @brief deconstructor
   */
  ~ArrayArray() {
    // std::cout<<this->GetRefCount()<<std::endl;
    // par_std_out("deconstructor: %d\n", refCount->GetRefCount());
    if (refCount->DecRefCount() <= 0) {
      // if(startIdx) delete[] startIdx;
      // if(data) delete[] data;
      DELETE_POINTERS(data);
      DELETE_POINTERS(startIdx);
      DELETE_POINTER(refCount);
      DELETE_POINTERS(basicEle);
    } else {
      startIdx = NULL;
      data = NULL;
      refCount = NULL;
      basicEle = NULL;
    }
  };
  /**
   * @brief get the size of structs
   * @return the size of structs
   */
  label size() const { return num; };
  /**
   * @brief print the class to screen
   */
  void display() const {
    for (int i = 0; i < num; ++i) {
      std::cout << "Item: " << i << " (";
      for (int j = startIdx[i]; j < startIdx[i + 1]; ++j) {
        std::cout << data[j];
        if (j < startIdx[i + 1] - 1) std::cout << ",";
      }
      std::cout << ")" << std::endl;
    }
  }
  /**
   * @brief overload function of []
   */
  BasicElement<T> &operator[](const int i) {
    if (!basicEle) {
      basicEle = new BasicElement<T>[num];
      for (int i = 0; i < num; ++i) {
        basicEle[i].num = startIdx[i + 1] - startIdx[i];
        basicEle[i].data = &data[startIdx[i]];
      }
    }
    return basicEle[i];
  }
};

/**
 * @brief divide the data into two parts,
 *        the values of left part is smaller than pivot
 *        the values of right part is larger than pivot
 * @param[in] arr the unpartitioned data
 * @param[in] l the lower bound
 * @param[in] r the upper bound
 * @return the index of partition boundary
 */
template <class T>
int partition(Array<Array<T> > &arr, int l, int r) {
  int k = l, pivot = arr[r][0];
  for (int i = l; i < r; ++i) {
    if (arr[i][0] <= pivot) {
      arr[i].swap(arr[k]);
      k++;
    }
  }
  arr[k].swap(arr[r]);
  return k;
}

/**
 * @brief quick sort algorithm for Array<Array<T>>
 * @param[in,out] arr unsorted array
 * @param[in] l the lower bound
 * @param[in] r the upper bound
 */
template <class T>
void quicksortArray(Array<Array<T> > &arr, int l, int r) {
  if (l < r) {
    int pivot = partition(arr, l, r);
    quicksortArray(arr, l, pivot - 1);
    quicksortArray(arr, pivot + 1, r);
  }
}

/**
 * @brief eliminate the duplicate elements
 *        and pick the unique ones and divide them into two parts
 * @param[in,out] arr original array
 * @return label[0]: size of the unique one
 *         label[1]: size of the duplicate one
 */
template <class T>
void filterArray(Array<Array<T> > &arr, label *faceNum) {
  int num = arr.size();
  quicksortArray(arr, 0, num - 1);
  // for(int i=0;i<tmp;i++) printf("%d\n", arr[0][i]);
  sort(arr[0].begin(), arr[0].end());
  Array<Array<T> > bndArr, innArr;
  for (int i = 0; i < num; ++i) {
    sort(arr[i].begin(), arr[i].end());
  }
  int end = 0;
  bool *isInner = new bool[num];
  for (int i = 0; i < num; ++i) {
    isInner[i] = false;
  }
  while (end < num) {
    // printf("%dth elements in %d\n", end, num);
    // printf("%d, %d\n", end, num);
    if (isInner[end]) {
      end++;
      continue;
    }
    // for (int i = end+1; i < num; ++i)
    int i = end + 1;
    while (i < num && arr[i][0] == arr[end][0]) {
      if (compareArray(arr[i], arr[end])) {
        isInner[i] = true;
        isInner[end] = true;
        innArr.push_back(arr[end]);
        break;
      }
      i++;
    }
    if (!isInner[end]) bndArr.push_back(arr[end]);
    end++;
  }
  // printf("old Array Num: %d, new Array Num: %d\n", arr.size(),
  // newArr.size());
  arr.clear();
  arr.insert(arr.end(), bndArr.begin(), bndArr.end());
  arr.insert(arr.end(), innArr.begin(), innArr.end());
  printf("arr: %ld, bndArr: %ld, innArr: %ld\n", arr.size(), bndArr.size(),
         innArr.size());
  faceNum[0] = bndArr.size();
  faceNum[1] = innArr.size();
  delete[] isInner;
};

/**
 * @brief compare the two arrays if they are equal
 * @param[in] a array A
 * @param[in] b array B
 * @return equal or not
 */
template <class T>
bool compareArray(const Array<T> &a, const Array<T> &b) {
  int num_a = a.size();
  int num_b = b.size();
  num_a = std::min(num_a, num_b);
  // if(num_a!=num_b) return false;
  for (int i = 0; i < num_a; ++i) {
    // printf("%d, %d\n", a[i],b[i]);
    if (a[i] != b[i]) return false;
  }
  return true;
};

/**
 * @brief transform the Array<Array<>> to ArrayArray
 * @param[in] arr array with Array<Array<>> format
 * @param[out] res array with ArrayArray format
 */
template <class T>
void transformArray(const Array<Array<T> > &arr, ArrayArray<T> &res) {
  int cellNum = arr.size();
  res.num = cellNum;
  //防止内存泄露
  DELETE_POINTERS(res.data);
  DELETE_POINTERS(res.startIdx);
  DELETE_POINTERS(res.basicEle);
  res.startIdx = new label[cellNum + 1];
  res.basicEle = new BasicElement<T>[cellNum];
  res.startIdx[0] = 0;
  for (int i = 0; i < cellNum; ++i) {
    res.startIdx[i + 1] = res.startIdx[i] + arr[i].size();
  }
  // printf("cellNum: %d, nodeNum: %d\n", cellNum, res.startIdx[cellNum]);
  res.data = new T[res.startIdx[cellNum]];
  for (int i = 0; i < cellNum; ++i) {
    int k = 0;
    for (int j = res.startIdx[i]; j < res.startIdx[i + 1]; ++j) {
      res.data[j] = arr[i][k];
      k++;
    }
    res.basicEle[i].num = res.startIdx[i + 1] - res.startIdx[i];
    res.basicEle[i].data = &res.data[res.startIdx[i]];
  }
};

/**
 * @brief transform the ArrayArray to Array<Array<>>
 * @param[in] arr array with ArrayArray format
 * @param[out] res array with Array<Array<>> format
 */
template <class T>
void transformArray(const ArrayArray<T> &arr, Array<Array<T> > &res) {
  int cellNum = arr.size();
  res.resize(cellNum);
  // printf("%d\n", cellNum);
  for (int i = 0; i < cellNum; ++i) {
    // printf("%d, %d\n", arr.startIdx[i], arr.startIdx[i+1]);
    for (int j = arr.startIdx[i]; j < arr.startIdx[i + 1]; ++j) {
      res[i].push_back(arr.data[j]);
    }
  }
};

/**
 * @brief find if the value exists in the array
 * @param[in] arr array array
 * @param[in] value array
 * @return the index of entry
 */
template <class T>
label findArray(Array<Array<T> > &arr, Array<T> &value) {
  int num = arr.size();
  // printf("%d\n", num);
  // printf("*****************************************************\n");
  for (int i = 0; i < num; ++i) {
    // for (int j = 0; j < value.size(); ++j)
    // {
    // printf("(%d, %d)", arr[i][j], value[j]);
    // }
    // printf("\n");
    if (compareArray(arr[i], value)) {
      return i;
    };
  }
  // printf("*****************************************************\n");
  return -1;
}

/**
 * @brief find if the value exists in the sorted array
 * @param[in] arr sorted array array
 * @param[in] value array
 * @param[in] l the lower bound
 * @param[in] r the upper bound
 * @return the index of entry
 */
template <class T>
label findSortedArray(Array<Array<T> > &arr, Array<T> &value, label l,
                      label r) {
  // label num = std::min(arr.size(),value.size());
  label num = arr.size();
  bool isFinded = false;
  label m;
  while (l < r) {
    m = (l + r) / 2;
    if (l == m) {
      if (arr[r][0] == value[0]) {
        isFinded = true;
        m = r;
      }
      break;
    }
    if (arr[m][0] < value[0])
      l = m;
    else if (arr[m][0] > value[0])
      r = m;
    else {
      isFinded = true;
      break;
    }
  }
  if (!isFinded) return -1;
  int i = m;
  while (i < num && arr[i][0] == value[0]) {
    if (compareArray(arr[i], value)) return i;
    i++;
  }
  i = m;
  while (i >= 0 && arr[i][0] == value[0]) {
    if (compareArray(arr[i], value)) return i;
    i--;
  }
  // // printf("%d\n", num);
  // // printf("*****************************************************\n");
  // for (int i = 0; i < num; ++i)
  // {
  //  // printf("%d, %d\n", i, num);
  //  if(arr[i][0]!=value[0]) continue;
  //  // if(compareArray(arr[i], value)) {return i;};
  // }
  // // printf("*****************************************************\n");
  return -1;
}

/**
 * @brief sort, Sort an single array that contains data can be compared
 * @tparam D data type that can be compared with "<"
 * @param[in,out] data the data array
 * @return
 */
template <class D>
int sort(Array<D> &data) {
  MultiList<D> table;

  size_t size = data.size();
  for (size_t i = 0; i < size; i++) table.insert(data[i]);

  size_t count = 0;
  for (typename MultiList<D>::iterator iter = table.begin();
       iter != table.end(); iter++) {
    data[count] = *iter;
    count++;
  }

  if (count != size) {
    cerr << __FILE__ << " +" << __LINE__ << ":" << endl
         << __FUNCTION__ << ":" << endl;
    cerr << "Error: Sorted data is not consistent with original size!" << endl;
    exit(-1);
  }
}

/**
 * @brief sort, Sort an single array that contains data can be compared with
 * corresponding weights.
 * @tparam D data type
 * @tparam W weight type that can be compared with "<"
 * @param[in,out] data the data array
 * @param[in] weights the weight array
 * @return
 */
template <class D, class W>
int sort(Array<D> &data, const Array<W> &weights = Array<W>(0)) {
  MultiTable<W, D> table;

  size_t size = data.size();
  for (size_t i = 0; i < size; i++)
    table.insert(pair<W, D>(weights[i], data[i]));

  size_t count = 0;
  for (typename MultiTable<W, D>::iterator iter = table.begin();
       iter != table.end(); iter++) {
    data[count] = iter->second;
    count++;
  }

  if (count != size) {
    cerr << __FILE__ << " +" << __LINE__ << ":" << endl
         << __FUNCTION__ << ":" << endl;
    cerr << "Error: Sorted data is not consistent with original size!" << endl;
    exit(-1);
  }
}

/**
 * @brief compare the two value
 * @tparam W
 * @return 1: a<b; -1: a>b; 0:a=b;
 */
template <class W>
int ascend(const void *a, const void *b) {
  cout << "comparing " << *(W *)a << " and " << *(W *)b << endl;
  if (*(W *)a < *(W *)b)
    return -1;
  else if (*(W *)b < *(W *)a)
    return 1;
  else
    return 0;
}

/**
 * @brief compare the two value
 * @tparam W
 * @return 1: a>b; -1: a<b; 0:a=b;
 */
template <class W>
int descend(const void *a, const void *b) {
  cout << "comparing " << *(W *)a << " and " << *(W *)b << endl;
  if (*(W *)a > *(W *)b)
    return -1;
  else if (*(W *)b < *(W *)a)
    return 1;
  else
    return 0;
}
typedef int ComparePtr(const void *, const void *);
/**
 * @brief biSearch, search the array with binary searching
 * @tparam W the data type that can be compared with "<"
 * @param[in] data data array sorted in ascending order
 * @param[in] value the target value to search for
 * @param[in] comp the rule of bi-search
 * @return the position target value located; return data.size() if not found.
 */
template <class W>
size_t biSearch(const Array<W> &data, W value, ComparePtr comp = ascend<W>) {
  // typename Array<W>::iterator iter;
  // iter = std::lower_bound(data.begin(), data.end(), value);
  // return (size_t) (iter - data.begin() );
  size_t num = data.size();
  size_t size = sizeof(W);
  W *posi = (W *)bsearch(&value, &(data[0]), num, size, comp);
  // cout<<"find value: "<<*posi<<endl;
  return (W *)posi - (W *)&(data[0]);
}

}  // namespace UTILITY

#endif  // UTILITY_CONTAINER_HPP
