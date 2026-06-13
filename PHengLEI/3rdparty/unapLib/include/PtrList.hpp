/// \file PtrList.hpp
/// \author Hanfeng GU
/// \version 1.0
/// \date 2019-01-30
/// \brief 对象指针列表，支持增删

#ifndef PTRLIST_HPP
#define PTRLIST_HPP

#include <set>

#include "unapMPI.hpp"

namespace UNAP {
/// \brief list of object pointers, fast add/delete supported
template <typename T>
class PtrList {
 private:
  /// \brief 指针列表大小
  label size_;

  /// \brief 指针列表首地址
  T **ptrs_;

  Array<int> reUseIndex_;

 public:
  /// \brief null constructor
  PtrList();

  /// \brief construct with size specified
  /// \param[in] label 指针列表大小
  PtrList(const label);

  /// \brief copy constructor
  PtrList(const PtrList<T> &);

  /// \brief destructor
  ~PtrList() {
    if (ptrs_) {
      forAll(i, size_) {
        if (ptrs_[i] && reUseIndex_[i] == 0) {
          delete ptrs_[i];
          ptrs_[i] = NULL;
        }
      }

      delete[] ptrs_;
      ptrs_ = NULL;
    }
  }

  /// \brief return the number of elements in the PtrList
  inline label size() const;

  /// \brief set size of PtrList
  /// \param[in] i 指针列表大小
  inline void SET_size(const label i);

  /// \brief return true if the PtrList is empty (ie, size() is zero)
  inline bool isEmpty() const;

  /// \brief deep copy
  PtrList &operator=(const PtrList<T> &);

  /// \brief return element const reference
  inline T &operator[](const label i) const {
    if (!ptrs_[i]) {
      COUT << "\nError: null pointer detected!" << ENDL;
      ERROR_EXIT;
    } else {
      return *ptrs_[i];
    }
  }
  /// \brief make the ptrs_[i] point to the object
  void setLevel(const label i, T &obj, bool reUse = false);

  /// \brief delete a level
  void removeLevel(const label i);
};

template <class T>
PtrList<T>::PtrList() : size_(50), ptrs_(NULL) {
  ptrs_ = new T *[size_];
  reUseIndex_.resize(size_);
  forAll(i, size_) {
    ptrs_[i] = NULL;
    reUseIndex_[i] = 0;
  }
}

template <class T>
PtrList<T>::PtrList(const label size) : size_(size), ptrs_(NULL) {
  ptrs_ = new T *[size_];
  reUseIndex_.resize(size_);
  forAll(i, size_) {
    ptrs_[i] = NULL;
    reUseIndex_[i] = 0;
  }
}

/*! 深拷贝，新分配了空间给每个指针 */
template <class T>
PtrList<T>::PtrList(const PtrList<T> &oldObj) : size_(-1), ptrs_(NULL) {
  size_ = oldObj.size();
  ptrs_ = new T *[size_];
  reUseIndex_.resize(size_);
  forAll(i, size_) {
    ptrs_[i] = new T();
    T &newObj = *ptrs_[i];
    newObj = oldObj[i];
    reUseIndex_[i] = 0;
  }
}

template <class T>
PtrList<T> &PtrList<T>::operator=(const PtrList<T> &oldObj) {
  if (this == &oldObj) {
    return *this;
  } else {
    //- delete existed
    if (this->ptrs_) {
      forAll(i, this->size_) { DELETE_POINTER(this->ptrs_[i]) }
      delete[] this->ptrs_;
      this->ptrs_ = NULL;
    }

    //- copy
    this->size_ = oldObj.size();
    reUseIndex_.resize(size_);
    this->ptrs_ = new T *[this->size_];
    forAll(i, this->size_) {
      this->ptrs_[i] = new T(oldObj[i]);
      reUseIndex_[i] = 0;
    }
    return *this;
  }
}

template <class T>
inline label PtrList<T>::size() const {
  return size_;
}

template <class T>
inline void PtrList<T>::SET_size(const label newSize) {
#ifdef DEBUG
  if (newSize < 0) {
    POUT << "Error in PtrList SET_size: "
         << "bad new set size " << newSize << ENDL;

    ERROR_EXIT;
  }
#endif
  label oldSize = size();

  if (newSize == 0) {
    this->~PtrList();
  } else if (oldSize == newSize) {
    //- nothing to do
  } else {
    T **ptrsOld = ptrs_;
    ptrs_ = new T *[newSize];
    reUseIndex_.resize(newSize);
    forAll(i, newSize) {
      if (i < oldSize) {
        ptrs_[i] = ptrsOld[i];
      } else {
        ptrs_[i] = NULL;
        reUseIndex_[i] = 0;
      }
    }
    size_ = newSize;
    delete[] ptrsOld;
  }
}

template <class T>
inline bool PtrList<T>::isEmpty() const {
  if (size_ == 0) {
    return true;
  } else {
    return false;
  }
}

/*! 删除了该level对应的空间 a*/
template <class T>
void PtrList<T>::removeLevel(const label leveli) {
  label oldSize = size_;
  if (leveli < size_) {
    delete ptrs_[leveli];
    for (label i = leveli; i < size_ - 1; i++) {
      ptrs_[i] = ptrs_[i + 1];
    }
    ptrs_[size_ - 1] = NULL;
    size_--;
  } else {
    POUT << "Error: removeLevel failed!" << ENDL;
    POUT << "Number of levels is " << oldSize
         << ", while the level tended to delete is " << leveli << ENDL;
    ERROR_EXIT;
  }
  reUseIndex_.erase(reUseIndex_.begin() + leveli);
}

/*! 浅拷贝，未分配空间 */
template <typename T>
void PtrList<T>::setLevel(const label leveli, T &obj, bool reUse) {
  if (leveli < size_) {
    ptrs_[leveli] = &obj;
    if (reUse) reUseIndex_[leveli] = 1;
  } else {
    POUT << "Error: setLevel failed!" << ENDL;
    POUT << "Number of levels is " << size_
         << ", while the level tended to set is " << leveli << ENDL;
    ERROR_EXIT;
  }
}

}  // namespace UNAP

#endif  //- PTRLIST_HPP
