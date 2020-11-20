/*************************************************************************
 *  Copyright (c) 2016.
 *  All rights reserved.
 *  This file is part of the CLAIRE library.
 *
 *  CLAIRE is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  CLAIRE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with CLAIRE. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/
#pragma once

#ifndef _MEMORYUTILS_HPP_
#define _MEMORYUTILS_HPP_

#include "TypeDef.hpp"

namespace reg {

/*! throw error (PETSc interface) defined elsewhere */
extern PetscErrorCode ThrowErrorMsg(std::bad_alloc&, int, const char*);
/*
template<class T> inline PetscErrorCode AllocateArray(T*& ptr, size_t N) {
  try {
    ptr = new T[N];
  } catch (std::bad_alloc& err) {
    return reg::ThrowError(err);
  }
  return 0;
}
template<class A, class T, class ... Args> inline PetscErrorCode Allocate(T*& ptr, Args ... args) {
  try {
    ptr = new A(args...);
  } catch (std::bad_alloc& err) {
    return reg::ThrowError(err);
  }
  return 0;
}
template<class T, class ... Args> inline PetscErrorCode Allocate(T*& ptr, Args ... args) {
  return Allocate<T>(ptr, args...);
}
*/
template<class T> inline PetscErrorCode AllocateArrayOnce(T*& ptr, size_t N) {
  if (ptr == nullptr) {
    try {
      ptr = new T[N];
    } catch (std::bad_alloc& err) {
      return ThrowError(err);
    }
  }
  return 0;
}
template<class A, class T, class ... Args> inline PetscErrorCode AllocateOnce(T*& ptr, Args ... args) {
  if (ptr == nullptr) {
    try {
      ptr = new A(args...);
    } catch (std::bad_alloc& err) {
      return ThrowError(err);
    }
  }
  return 0;
}
template<class T, class ... Args> inline PetscErrorCode AllocateOnce(T*& ptr, Args ... args) {
  return AllocateOnce<T, T, Args...>(ptr, args...);
}
template<class T> inline PetscErrorCode AllocateMemoryOnce(T*& ptr, size_t size) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  if (ptr == nullptr) {
#ifdef REG_HAS_CUDA
    ierr = cudaMalloc(reinterpret_cast<void**>(&ptr), size); CHKERRCUDA(ierr);
#else
    ptr = reinterpret_cast<T*>(claire_alloc(size));
#endif
  }
  PetscFunctionReturn(ierr);
}

template<class T> inline PetscErrorCode Free(T*& ptr) {
  if (ptr != nullptr)
    delete ptr;
  ptr = nullptr;
  return 0;
}
template<class T> inline PetscErrorCode FreeArray(T*& ptr) {
  if (ptr != nullptr)
    delete[] ptr;
  ptr = nullptr;
  return 0;
}
template<class T> inline PetscErrorCode FreeMemory(T*& ptr) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  if (ptr != nullptr) {
#ifdef REG_HAS_CUDA
    ierr = cudaFree(ptr); CHKERRCUDA(ierr);
#else
    claire_free(ptr);
#endif
  }
  ptr = nullptr;
  PetscFunctionReturn(ierr);
}

template<class T>
class ManagedMemory {
  public:
    typedef ManagedMemory Self;

    ManagedMemory();
    ManagedMemory(size_t);
    virtual ~ManagedMemory();

    PetscErrorCode Resize(size_t);

    explicit operator T*();
    T& operator [](size_t);

    const T* ReadHost();
    T* WriteHost();
    T* ReadWriteHost();

    const T* ReadDevice();
    T* WriteDevice();
    T* ReadWriteDevice();

    PetscErrorCode CopyHostToDevice();
    PetscErrorCode CopyDeviceToHost();

    PetscErrorCode AllocateHost();
    PetscErrorCode AllocateDevice();

  protected:
    bool m_HostValid;
    bool m_DeviceValid;

    T *m_HostPtr;
    T *m_DevicePtr;

    T* m_CurrentPtr;

    size_t m_N;

  private:
    PetscErrorCode Initialize();
    PetscErrorCode ClearMemory();
};

template<class T> ManagedMemory<T>::ManagedMemory() {
  this->Initialize();
}
template<class T> ManagedMemory<T>::ManagedMemory(size_t N) {
  this->Initialize();
  this->Resize(N);
}
template<class T> ManagedMemory<T>::~ManagedMemory() {
  this->ClearMemory();
}

template<class T> PetscErrorCode ManagedMemory<T>::Resize(size_t N) {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  if (N != this->m_N) {
    ierr = this->ClearMemory(); CHKERRQ(ierr);
    this->m_N = N;
    this->m_HostValid = false;
    this->m_DeviceValid = false;
  }

  PetscFunctionReturn(ierr);
}

template<class T> const T* ManagedMemory<T>::ReadHost() {
  this->AllocateHost();
  this->CopyDeviceToHost();

  return this->m_HostPtr;
}
template<class T> T* ManagedMemory<T>::WriteHost() {
  this->AllocateHost();
  this->m_HostValid = true;
  this->m_DeviceValid = false;

  this->m_CurrentPtr = this->m_HostPtr;

  return this->m_HostPtr;
}
template<class T> T* ManagedMemory<T>::ReadWriteHost() {
  this->AllocateHost();
  this->CopyDeviceToHost();
  this->m_HostValid = true;
  this->m_DeviceValid = false;

  return this->m_HostPtr;
}

template<class T> const T* ManagedMemory<T>::ReadDevice() {
  this->AllocateDevice();
  this->CopyHostToDevice();

  return this->m_DevicePtr;
}
template<class T> T* ManagedMemory<T>::WriteDevice() {
  this->AllocateDevice();
  this->m_HostValid = false;
  this->m_DeviceValid = true;

  this->m_CurrentPtr = this->m_DevicePtr;

  return this->m_DevicePtr;
}
template<class T> T* ManagedMemory<T>::ReadWriteDevice() {
  this->AllocateDevice();
  this->CopyHostToDevice();
  this->m_HostValid = false;
  this->m_DeviceValid = true;

  return this->m_DevicePtr;
}

template<class T> PetscErrorCode ManagedMemory<T>::CopyHostToDevice() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

#ifdef REG_HAS_CUDA
  if (this->m_DevicePtr) {
    if (!this->m_DeviceValid && this->m_HostPtr) {
      ierr = cudaMemcpy(static_cast<void*>(this->m_DevicePtr),
        static_cast<const void*>(this->m_HostPtr),
        this->m_N*sizeof(T),
        cudaMemcpyHostToDevice); CHKERRCUDA(ierr);
      this->m_DeviceValid = true;
    }
  } else {
    ierr = this->AllocateDevice(); CHKERRQ(ierr);
    this->m_DeviceValid = true;
    if (this->m_HostPtr) {
      ierr = cudaMemcpy(static_cast<void*>(this->m_DevicePtr),
        static_cast<const void*>(this->m_HostPtr),
        this->m_N*sizeof(T),
        cudaMemcpyHostToDevice); CHKERRCUDA(ierr);
    }
  }
#endif
  this->m_CurrentPtr = this->m_DevicePtr;

  PetscFunctionReturn(ierr);
}
template<class T> PetscErrorCode ManagedMemory<T>::CopyDeviceToHost() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

#ifdef REG_HAS_CUDA
  if (this->m_HostPtr) {
    if (!this->m_HostValid && this->m_DevicePtr) {
      ierr = cudaMemcpy(static_cast<void*>(this->m_HostPtr),
        static_cast<const void*>(this->m_DevicePtr),
        this->m_N*sizeof(T),
        cudaMemcpyDeviceToHost); CHKERRCUDA(ierr);
      this->m_HostValid = true;
    }
  } else {
    ierr = this->AllocateHost(); CHKERRQ(ierr);
    this->m_HostValid = true;
    if (this->m_DevicePtr) {
      ierr = cudaMemcpy(static_cast<void*>(this->m_HostPtr),
        static_cast<const void*>(this->m_DevicePtr),
        this->m_N*sizeof(T),
        cudaMemcpyDeviceToHost); CHKERRCUDA(ierr);
    }
  }
#endif
  this->m_CurrentPtr = this->m_HostPtr;

  PetscFunctionReturn(ierr);
}

template<class T> PetscErrorCode ManagedMemory<T>::AllocateHost() {
  PetscFunctionBegin;

  if (this->m_HostPtr == nullptr)
    this->m_HostPtr = reinterpret_cast<T*>(claire_alloc(this->m_N*sizeof(T)));

  this->m_CurrentPtr = this->m_HostPtr;

  PetscFunctionReturn(0);
}
template<class T> PetscErrorCode ManagedMemory<T>::AllocateDevice() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  if (this->m_DevicePtr == nullptr) {
#ifdef REG_HAS_CUDA
    ierr = cudaMalloc(reinterpret_cast<void**>(&this->m_DevicePtr), this->m_N*sizeof(T)); CHKERRCUDA(ierr);
#else
    ierr = this->AllocateHost(); CHKERRQ(ierr);
    this->m_DevicePtr = this->m_HostPtr;
#endif
  }

  this->m_CurrentPtr = this->m_DevicePtr;

  PetscFunctionReturn(ierr);
}

template<class T> PetscErrorCode ManagedMemory<T>::Initialize() {
  PetscFunctionBegin;

  this->m_HostPtr = nullptr;
  this->m_DevicePtr = nullptr;
  this->m_CurrentPtr = nullptr;
  this->m_N = 0;
  this->m_HostValid = false;
  this->m_DeviceValid = false;

  PetscFunctionReturn(0);
}
template<class T> PetscErrorCode ManagedMemory<T>::ClearMemory() {
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;

  if (this->m_HostPtr)
    claire_free(this->m_HostPtr);
#ifdef REG_HAS_CUDA
  if (this->m_DevicePtr) {
    ierr = cudaFree(this->m_DevicePtr); CHKERRCUDA(ierr);
  }
#endif
  this->m_HostPtr = nullptr;
  this->m_DevicePtr = nullptr;
  this->m_CurrentPtr = nullptr;
  this->m_N = 0;
  this->m_HostValid = false;
  this->m_DeviceValid = false;

  PetscFunctionReturn(ierr);
}

template<class T> ManagedMemory<T>::operator T* () {
  return this->m_CurrentPtr;
}

template<class T> T& ManagedMemory<T>::operator [] (size_t idx) {
  return this->m_CurrentPtr[idx];
}

} // namespace reg

#endif
