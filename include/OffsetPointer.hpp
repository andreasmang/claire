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


#ifndef _OFFSETPOINTER_HPP_
#define _OFFSETPOINTER_HPP_

namespace reg {

template<typename T>
class OffsetPointer {
 public:
  OffsetPointer() : m_RawPointer(nullptr), m_Offset(0) {};
  virtual ~OffsetPointer() {}
  
  T*& raw() { return m_RawPointer; }
  T* ptr() { return m_RawPointer + m_Offset; }
  
  T& operator [] (IntType idx) { return m_RawPointer[m_Offset + idx]; }
  
  operator T*() const { return m_RawPointer + m_Offset; }
  operator T**() { return &m_RawPointer; }
  
  void SetOffest(size_t offset) { m_Offset = offset; }
 protected:
  T* m_RawPointer;
  size_t m_Offset;
};

} // namespace reg

#endif
