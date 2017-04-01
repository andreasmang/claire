
// Written by Dhairya Malhotra
// modified by Amir Gholami
#include <iostream>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <stdlib.h>
#include <cstdlib>
#include <cstdio>
#include <stdint.h>
#include <vector>
#include <stack>
#include <map>
#include <set>
#include <omp.h>
typedef     long  Integer;  // bounded numbers < 32k
typedef  int64_t     Long;  // problem size

namespace pvfmm{
  class MemoryManager;
#define MEM_ALIGN 64
#define GLOBAL_MEM_BUFF 0LL //in MB
#define PVFMM_WARN(msg) \
do { \
  std::cerr<<"\n\033[1;31mWarning:\033[0m "<<msg<<'\n'; \
}while(0)

#define PVFMM_ERROR(msg) \
do { \
  std::cerr<<"\n\033[1;31mError:\033[0m "<<msg<<'\n'; \
  abort(); \
}while(0)

#define PVFMM_ASSERT_MSG(cond, msg) \
do { \
  if (!(cond)) PVFMM_ERROR(msg); \
}while(0)

#define PVFMM_ASSERT(cond) \
do { \
  if (!(cond)) { \
    fprintf (stderr, "\n%s:%d: %s: Assertion `%s' failed.\n",__FILE__,__LINE__,__PRETTY_FUNCTION__,#cond); \
    abort(); \
  } \
}while(0)

#include <iostream>

#ifdef PVFMM_MEMDEBUG
  template <class ValueType> class ConstIterator{

    template<typename T> friend class ConstIterator;

    template<typename T> friend class Iterator;

    public:
      typedef std::random_access_iterator_tag iterator_category;
      typedef const ValueType& reference;
      typedef Long difference_type;
      typedef ValueType value_type;
      typedef const ValueType* pointer;

    protected:
      char* base;
      difference_type len, offset;
      Long alloc_ctr;
      void* mem_head;
      static const Long ValueSize = sizeof(ValueType);

    public:

      ConstIterator(void* base_=NULL) {
        base=(char*)base_;
        len=0;
        offset=0;
        alloc_ctr=0;
        mem_head=NULL;
      }

      template <size_t LENGTH> ConstIterator(ValueType (&base_)[LENGTH]) { // DEPRECATED
        PVFMM_ASSERT(false);
      }

      ConstIterator(const ValueType* base_, difference_type len_, bool dynamic_alloc=false);

      template <class AnotherType> explicit ConstIterator(const ConstIterator<AnotherType>& I) {
        this->base=I.base;
        this->len=I.len;
        this->offset=I.offset;
        PVFMM_ASSERT_MSG((uintptr_t)(this->base+this->offset) % alignof(ValueType) == 0, "invalid alignment during pointer type conversion.");
        this->alloc_ctr=I.alloc_ctr;
        this->mem_head=I.mem_head;
      }

      //value_type* like operators
      reference operator*() const;

      const value_type * operator->() const;

      reference operator[](difference_type off) const;

      //Increment / Decrement
      ConstIterator& operator++(){
        offset+=(Long)sizeof(ValueType);
        return *this;
      }

      ConstIterator operator++(int){
        ConstIterator<ValueType> tmp(*this);
        ++*this;
        return tmp;
      }

      ConstIterator& operator--(){
        offset-=(Long)sizeof(ValueType);
        return *this;
      }

      ConstIterator operator--(int){
        ConstIterator<ValueType> tmp(*this);
        --*this;
        return tmp;
      }

      //Arithmetic
      ConstIterator& operator+=(difference_type i){
        offset+=i*(Long)sizeof(ValueType);
        return *this;
      }

      ConstIterator operator+(difference_type i) const{
        ConstIterator<ValueType> tmp(*this);
        tmp.offset+=i*(Long)sizeof(ValueType);
        return tmp;
      }

      friend ConstIterator operator+(difference_type i, const ConstIterator& right){
        return (right+i);
      }

      ConstIterator& operator-=(difference_type i){
        offset-=i*(Long)sizeof(ValueType);
        return *this;
      }

      ConstIterator operator-(difference_type i) const{
        ConstIterator<ValueType> tmp(*this);
        tmp.offset-=i*(Long)sizeof(ValueType);
        return tmp;
      }

      difference_type operator-(const ConstIterator& I) const{
        if(base!=I.base) PVFMM_WARN("comparing two unrelated memory addresses.");
        Long diff=((ValueType*)(base+offset))-((ValueType*)(I.base+I.offset));
        PVFMM_ASSERT_MSG(I.base+I.offset+diff*(Long)sizeof(ValueType)==base+offset, "invalid memory address alignment.");
        return diff;
      }

      //Comparison operators
      bool operator== (const ConstIterator& I) const{
        return (base+offset==I.base+I.offset);
      }

      bool operator!= (const ConstIterator& I) const{
        return !(*this==I);
      }

      bool operator< (const ConstIterator& I) const{
        if(base!=I.base) PVFMM_WARN("comparing two unrelated memory addresses.");
        return (base+offset)<(I.base+I.offset);
      }

      bool operator<= (const ConstIterator& I) const{
        if(base!=I.base) PVFMM_WARN("comparing two unrelated memory addresses.");
        return (base+offset)<=(I.base+I.offset);
      }

      bool operator> (const ConstIterator& I) const{
        if(base!=I.base) PVFMM_WARN("comparing two unrelated memory addresses.");
        return (base+offset)>(I.base+I.offset);
      }

      bool operator>= (const ConstIterator& I) const{
        if(base!=I.base) PVFMM_WARN("comparing two unrelated memory addresses.");
        return (base+offset)>=(I.base+I.offset);
      }

      friend std::ostream &operator<<(std::ostream &out, const ConstIterator &I) {
        out<<"("<<(long long)I.base<<"+"<<I.offset<<":"<<I.len<<")";
        return out;
      }
  };

  template <class ValueType> class Iterator : public ConstIterator<ValueType>{

    public:
      typedef std::random_access_iterator_tag iterator_category;
      typedef ValueType& reference;
      typedef Long difference_type;
      typedef ValueType value_type;
      typedef ValueType* pointer;

    public:

      Iterator(void* base_=NULL) : ConstIterator<ValueType>(base_){}

      template <size_t LENGTH> Iterator(ValueType (&base_)[LENGTH]) : ConstIterator<ValueType>(base_){}

      Iterator(ValueType* base_, difference_type len_, bool dynamic_alloc=false) : ConstIterator<ValueType>(base_, len_, dynamic_alloc) {}

      template <class AnotherType> explicit Iterator(const ConstIterator<AnotherType>& I) : ConstIterator<ValueType>(I) {}

      //value_type* like operators
      reference operator*() const;

      value_type* operator->() const;

      reference operator[](difference_type off) const;

      //Increment / Decrement
      Iterator& operator++(){
        this->offset+=(Long)sizeof(ValueType);
        return *this;
      }

      Iterator operator++(int){
        Iterator<ValueType> tmp(*this);
        ++*this;
        return tmp;
      }

      Iterator& operator--(){
        this->offset-=(Long)sizeof(ValueType);
        return *this;
      }

      Iterator operator--(int){
        Iterator<ValueType> tmp(*this);
        --*this;
        return tmp;
      }

      // Arithmetic
      Iterator& operator+=(difference_type i){
        this->offset+=i*(Long)sizeof(ValueType);
        return *this;
      }

      Iterator operator+(difference_type i) const{
        Iterator<ValueType> tmp(*this);
        tmp.offset+=i*(Long)sizeof(ValueType);
        return tmp;
      }

      friend Iterator operator+(difference_type i, const Iterator& right){
        return (right+i);
      }

      Iterator& operator-=(difference_type i){
        this->offset-=i*(Long)sizeof(ValueType);
        return *this;
      }

      Iterator operator-(difference_type i) const{
        Iterator<ValueType> tmp(*this);
        tmp.offset-=i*(Long)sizeof(ValueType);
        return tmp;
      }

      difference_type operator-(const ConstIterator<ValueType>& I) const{
        return static_cast<const ConstIterator<ValueType>&>(*this) - I;
      }
  };

  template <class ValueType, Long DIM> class StaticArray : public Iterator<ValueType>{

    public:
      StaticArray();

      ~StaticArray();

    private:
      StaticArray(const StaticArray&);
      StaticArray& operator=(const StaticArray&);

      Iterator<ValueType> arr;
  };

#define      PVFMM_PTR2ITR(type, ptr, len) pvfmm::     Iterator<type>((      type*)ptr,len)
#define PVFMM_PTR2CONSTITR(type, ptr, len) pvfmm::ConstIterator<type>((const type*)ptr,len)
#else
  template <typename ValueType> using Iterator = ValueType*;
  template <typename ValueType> using ConstIterator = const ValueType*;
  template <typename ValueType, Integer DIM> using StaticArray = ValueType[DIM];

#define      PVFMM_PTR2ITR(type, ptr, len) (      type*)ptr
#define PVFMM_PTR2CONSTITR(type, ptr, len) (const type*)ptr
#endif

  /**
   * \brief Identify each type uniquely.
   */
  template <class T> class TypeTraits{

    public:

      static uintptr_t ID();

      static bool IsPOD();
  };
  /**
   * \brief MemoryManager class declaration.
   */
  class MemoryManager{

    public:

      static const char init_mem_val=42;

      /**
       * \brief Header data for each memory block.
       */
      struct MemHead{
        Long n_indx;
        Long n_elem;
        Long type_size;
        Long alloc_ctr;
        uintptr_t type_id;
        unsigned char check_sum;
      };

      /**
       * \brief Constructor for MemoryManager.
       */
      MemoryManager(Long N);

      /**
       * \brief Constructor for MemoryManager.
       */
      ~MemoryManager();

      static MemHead& GetMemHead(char* p);

      static void CheckMemHead(const MemHead& p);

      Iterator<char> malloc(const Long n_elem=1, const Long type_size=sizeof(char), const uintptr_t type_id=TypeTraits<char>::ID()) const;

      void free(Iterator<char> p) const;

      void print() const;

      static void test();

      // Check all free memory equals init_mem_val
      void Check() const;

      // A global MemoryManager object. This is the default for aligned_new and aligned_free
      static MemoryManager& glbMemMgr(){
        static MemoryManager m(GLOBAL_MEM_BUFF*1024LL*1024LL);
        return m;
      }

    private:

      // Private constructor
      MemoryManager();

      // Private copy constructor
      MemoryManager(const MemoryManager& m);

      /**
       * \brief Node structure for a doubly linked list, representing free and
       * occupied memory blocks. Blocks are split, merged or state is changed
       * between free and occupied in O(1) time given the pointer to the MemNode.
       */
      struct MemNode{
        bool free;
        Long size;
        char* mem_ptr;
        Long prev, next;
        std::multimap<Long, Long>::iterator it;
      };

      /**
       * \brief Return index of one of the available MemNodes from node_stack or
       * create new MemNode by resizing node_buff.
       */
      Long new_node() const;

      /**
       * \brief Add node index for now available MemNode to node_stack.
       */
      void delete_node(Long indx) const;

      char* buff; // pointer to memory buffer.
      Long buff_size; // total buffer size in bytes.
      Long n_dummy_indx; // index of first (dummy) MemNode in link list.

      mutable std::vector<MemNode> node_buff; // storage for MemNode objects, this can only grow.
      mutable std::stack<Long> node_stack; // stack of available free MemNodes from node_buff.
      mutable std::multimap<Long, Long> free_map; // pair (MemNode.size, MemNode_id) for all free MemNodes.
      mutable omp_lock_t omp_lock; // openmp lock to prevent concurrent changes.
      mutable std::set<void*> system_malloc; // track pointers allocated using system malloc.
  };

  inline uintptr_t align_ptr(uintptr_t ptr){
    const uintptr_t ALIGN_MINUS_ONE=MEM_ALIGN-1;
    const uintptr_t NOT_ALIGN_MINUS_ONE=~ALIGN_MINUS_ONE;
    return ((ptr+ALIGN_MINUS_ONE) & NOT_ALIGN_MINUS_ONE);
  }
  /**
   * \brief Aligned allocation as an alternative to new. Uses placement new to
   * construct objects.
   */
  template <class ValueType> Iterator<ValueType> aligned_new(Long n_elem=1, const MemoryManager* mem_mgr=&MemoryManager::glbMemMgr());
  /**
   * \brief Aligned de-allocation as an alternative to delete. Calls the object
   * destructors. Not sure which destructor is called for virtual classes, this
   * is why we also match the TypeTraits<T>::ID()
   */
  template <class ValueType> void aligned_delete(Iterator<ValueType> A, const MemoryManager* mem_mgr=&MemoryManager::glbMemMgr());
  /**
   * \brief Wrapper to memcpy. Also checks if source and destination pointers are
   * the same.
   */
  template <class ValueType> Iterator<ValueType> memcopy(Iterator<ValueType> destination, ConstIterator<ValueType> source, Long num);

  template <class ValueType> Iterator<ValueType> memset(Iterator<ValueType> ptr, int value, Long num);
}

#include <compact_mem_mgr.txx>
// end header
