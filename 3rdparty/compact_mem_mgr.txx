
// Written by Dhairya Malhotra
// Modified by Amir Gholami
namespace pvfmm{
#ifdef PVFMM_MEMDEBUG
  template <class ValueType> inline ConstIterator<ValueType>::ConstIterator(const ValueType* base_, difference_type len_, bool dynamic_alloc) {
    this->base=(char*)base_;
    this->len=len_*(Long)sizeof(ValueType);
    this->offset=0;
    PVFMM_ASSERT_MSG((uintptr_t)(this->base+this->offset) % alignof(ValueType) == 0, "invalid alignment during pointer type conversion.");
    if (dynamic_alloc) {
      MemoryManager::MemHead& mh=*&MemoryManager::GetMemHead((char*)this->base);
      MemoryManager::CheckMemHead(mh);
      alloc_ctr = mh.alloc_ctr;
      mem_head=&mh;
    }
    else mem_head=NULL;
  }

  template <class ValueType> inline static void IteratorAssertChecks(Long offset, Long len, void* mem_head, Long alloc_ctr){
    PVFMM_ASSERT_MSG(offset>=0 && offset+(Long)sizeof(ValueType)<=len, "access to pointer [B"<<(offset<0?"":"+")<<offset<<",B"<<(offset+(Long)sizeof(ValueType)<0?"":"+")<<offset+(Long)sizeof(ValueType)<<") is outside of the range [B,B+"<<len<<").");
    if(mem_head){
      MemoryManager::MemHead& mh=*(MemoryManager::MemHead*)(mem_head);
      PVFMM_ASSERT_MSG(mh.alloc_ctr==alloc_ctr , "invalid memory address or corrupted memory.");
    }
  }

  template <class ValueType> inline typename ConstIterator<ValueType>::reference ConstIterator<ValueType>::operator*() const{
    IteratorAssertChecks<ValueType>(this->offset, this->len, this->mem_head, this->alloc_ctr);
    return *(ValueType*)(base+offset);
  }

  template <class ValueType> inline const typename ConstIterator<ValueType>::value_type* ConstIterator<ValueType>::operator->() const{
    IteratorAssertChecks<ValueType>(this->offset, this->len, this->mem_head, this->alloc_ctr);
    return (ValueType*)(base+offset);
  }

  template <class ValueType> inline typename ConstIterator<ValueType>::reference ConstIterator<ValueType>::operator[](difference_type j) const{
    IteratorAssertChecks<ValueType>(this->offset+j*(Long)sizeof(ValueType), this->len, this->mem_head, this->alloc_ctr);
    return *(ValueType*)(base+offset+j*(Long)sizeof(ValueType));
  }


  template <class ValueType> inline typename Iterator<ValueType>::reference Iterator<ValueType>::operator*() const{
    IteratorAssertChecks<ValueType>(this->offset, this->len, this->mem_head, this->alloc_ctr);
    return *(ValueType*)(this->base+this->offset);
  }

  template <class ValueType> inline typename Iterator<ValueType>::value_type* Iterator<ValueType>::operator->() const{
    IteratorAssertChecks<ValueType>(this->offset, this->len, this->mem_head, this->alloc_ctr);
    return (ValueType*)(this->base+this->offset);
  }

  template <class ValueType> inline typename Iterator<ValueType>::reference Iterator<ValueType>::operator[](difference_type j) const{
    IteratorAssertChecks<ValueType>(this->offset+j*(Long)sizeof(ValueType), this->len, this->mem_head, this->alloc_ctr);
    return *(ValueType*)(this->base+this->offset+j*(Long)sizeof(ValueType));
  }


  template <class ValueType, Long DIM> inline StaticArray<ValueType,DIM>::StaticArray() {
    arr=aligned_new<ValueType>(DIM);
    Iterator<ValueType>::operator=(arr);
  }

  template <class ValueType, Long DIM> inline StaticArray<ValueType,DIM>::~StaticArray(){
    aligned_delete<ValueType>(arr);
  }

  template <class ValueType, Long DIM> inline StaticArray<ValueType,DIM>::StaticArray(const StaticArray& I) {
    arr=aligned_new<ValueType>(DIM);
    Iterator<ValueType>::operator=(arr);
    for(Long i=0;i<DIM;i++) (*this)[i]=I[i];
  }

  template <class ValueType, Long DIM> inline StaticArray<ValueType,DIM>& StaticArray<ValueType,DIM>::operator=(const StaticArray& I) {
    for(Long i=0;i<DIM;i++) (*this)[i]=I[i];
    return *this;
  }

#endif


  template <class T> inline uintptr_t TypeTraits<T>::ID(){
    return (uintptr_t)&ID;
  }

  template <class T> inline bool TypeTraits<T>::IsPOD(){
    return false;
  }

#define PVFMMDefinePOD(type) template<> bool inline TypeTraits<type>::IsPOD(){return true;};
  PVFMMDefinePOD(char);
  PVFMMDefinePOD(float);
  PVFMMDefinePOD(double);
  PVFMMDefinePOD(int);
  PVFMMDefinePOD(long long);
  PVFMMDefinePOD(unsigned long);
  PVFMMDefinePOD(char*);
  PVFMMDefinePOD(float*);
  PVFMMDefinePOD(double*);
#undef PVFMMDefinePOD


  inline MemoryManager::MemoryManager(Long N){
    buff_size=N;
    { // Allocate buff
      assert(MEM_ALIGN <= 0x8000);
      Long alignment=MEM_ALIGN-1;
      char* base_ptr=(char*)::malloc(N+2+alignment);
      PVFMM_ASSERT_MSG(base_ptr, "memory allocation failed.");
      buff=(char*)((uintptr_t)(base_ptr+2+alignment) & ~(uintptr_t)alignment);
      ((uint16_t*)buff)[-1] = (uint16_t)(buff-base_ptr);
    }
    { // Initialize to init_mem_val
#ifdef PVFMM_MEMDEBUG
#pragma omp parallel for
      for(Long i=0;i<buff_size;i++){
        buff[i]=init_mem_val;
      }
#endif
    }
    n_dummy_indx=new_node();
    Long n_indx=new_node();
    MemNode& n_dummy=node_buff[n_dummy_indx-1];
    MemNode& n=node_buff[n_indx-1];

    n_dummy.size=0;
    n_dummy.free=false;
    n_dummy.prev=0;
    n_dummy.next=n_indx;
    n_dummy.mem_ptr=&buff[0];
    assert(n_indx);

    n.size=N;
    n.free=true;
    n.prev=n_dummy_indx;
    n.next=0;
    n.mem_ptr=&buff[0];
    n.it=free_map.insert(std::make_pair(N,n_indx));

    omp_init_lock(&omp_lock);
  }

  inline MemoryManager::~MemoryManager(){
    Check();
    MemNode* n_dummy=&node_buff[n_dummy_indx-1];
    MemNode* n=&node_buff[n_dummy->next-1];
    if (!n->free || n->size!=buff_size || node_stack.size()!=node_buff.size()-2 || !system_malloc.empty()) {
      PVFMM_WARN("memory leak detected.");
    }
    omp_destroy_lock(&omp_lock);

    { // free buff
      assert(buff);
      ::free(buff-((uint16_t*)buff)[-1]);
    }
  }

  inline MemoryManager::MemHead& MemoryManager::GetMemHead(char* I){
    PVFMM_ASSERT_MSG(I!=NULL, "NULL pointer exception.");
    static uintptr_t alignment=MEM_ALIGN-1;
    static uintptr_t header_size=(uintptr_t)(sizeof(MemHead)+alignment) & ~(uintptr_t)alignment;
    return *(MemHead*)(((char*)I)-header_size);
  }

  inline void MemoryManager::CheckMemHead(const MemHead& mem_head){ // Verify header check_sum
#ifdef PVFMM_MEMDEBUG
    Long check_sum=0;
    const unsigned char* base_=(const unsigned char*)&mem_head;
    for(Integer i=0;i<sizeof(MemHead);i++){
      check_sum+=base_[i];
    }
    check_sum-=mem_head.check_sum;
    check_sum=check_sum & ((1UL << (8*sizeof(mem_head.check_sum)))-1);
    PVFMM_ASSERT_MSG(check_sum==mem_head.check_sum, "invalid memory address or corrupted memory.");
#endif
  }

  inline Iterator<char> MemoryManager::malloc(const Long n_elem, const Long type_size, const uintptr_t type_id) const{
    if(!n_elem) return NULL;
    static uintptr_t alignment=MEM_ALIGN-1;
    static uintptr_t header_size=(uintptr_t)(sizeof(MemHead)+alignment) & ~(uintptr_t)alignment;

    Long size=n_elem*type_size+header_size;
    size=(uintptr_t)(size+alignment) & ~(uintptr_t)alignment;
    char* base=NULL;

    omp_set_lock(&omp_lock);
    static Long alloc_ctr = 0; alloc_ctr++;
    Long head_alloc_ctr = alloc_ctr;
    std::multimap<Long, Long>::iterator it=free_map.lower_bound(size);
    Long n_indx=(it!=free_map.end()?it->second:0);
    if(n_indx){ // Allocate from buff
      Long n_free_indx=(it->first>size?new_node():0);
      MemNode& n=node_buff[n_indx-1];
      assert(n.size==it->first);
      assert(n.it==it);
      assert(n.free);

      if(n_free_indx){ // Create a node for the remaining free part.
        MemNode& n_free=node_buff[n_free_indx-1];
        n_free=n;
        n_free.size-=size;
        n_free.mem_ptr=(char*)n_free.mem_ptr+size;
        { // Insert n_free to the link list
          n_free.prev=n_indx;
          if(n_free.next){
            Long n_next_indx=n_free.next;
            MemNode& n_next=node_buff[n_next_indx-1];
            n_next.prev=n_free_indx;
          }
          n.next=n_free_indx;
        }
        assert(n_free.free); // Insert n_free to free map
        n_free.it=free_map.insert(std::make_pair(n_free.size,n_free_indx));
        n.size=size; // Update n
      }

      n.free=false;
      free_map.erase(it);
      base = n.mem_ptr;
    }
    omp_unset_lock(&omp_lock);
    if(!base){ // Use system malloc
      Long end_padding=8; // to check for out-of-bound writes
      char* p = (char*)::malloc(size+2+alignment+end_padding);
      PVFMM_ASSERT_MSG(p, "memory allocation failed.");
      #ifdef PVFMM_MEMDEBUG
      { // system_malloc.insert(p)
        omp_set_lock(&omp_lock);
        system_malloc.insert(p);
        omp_unset_lock(&omp_lock);
      }
      { // set p[*] to init_mem_val
        #pragma omp parallel for
        for(Long i=0;i<size+2+alignment+end_padding;i++) p[i]=init_mem_val;
      }
      #endif
      { // base <-- align(p)
        base = (char*)((uintptr_t)(p+2+alignment) & ~(uintptr_t)alignment);
        ((uint16_t*)base)[-1] = (uint16_t)(base-p);
      }
    }

    { // Check out-of-bounds write
#ifdef PVFMM_MEMDEBUG
      if(n_indx){
#pragma omp parallel for
        for(Long i=0;i<size;i++) PVFMM_ASSERT_MSG(base[i]==init_mem_val, "memory corruption detected.");
      }
#endif
    }

    MemHead& mem_head=*(MemHead*)base;
    { // Set mem_head
      #ifdef PVFMM_MEMDEBUG
      for(Integer i=0;i<sizeof(MemHead);i++) base[i]=init_mem_val;
      #endif
      mem_head.n_indx=n_indx;
      mem_head.n_elem=n_elem;
      mem_head.type_size=type_size;
      mem_head.alloc_ctr=head_alloc_ctr;
      mem_head.type_id=type_id;
    }
    { // Set header check_sum
#ifdef PVFMM_MEMDEBUG
      Long check_sum=0;
      unsigned char* base_=(unsigned char*)base;
      mem_head.check_sum=0;
      for(Integer i=0;i<sizeof(MemHead);i++) check_sum+=base_[i];
      check_sum=check_sum & ((1UL << (8*sizeof(mem_head.check_sum)))-1);
      mem_head.check_sum=check_sum;
#endif
    }
    #ifdef PVFMM_MEMDEBUG
    return Iterator<char>(base+header_size, n_elem*type_size,true);
    #else
    return base+header_size;
    #endif
  }

  inline void MemoryManager::free(Iterator<char> p) const{
    if(p==NULL) return;
    static uintptr_t alignment=MEM_ALIGN-1;
    static uintptr_t header_size=(uintptr_t)(sizeof(MemHead)+alignment) & ~(uintptr_t)alignment;

    MemHead& mem_head=GetMemHead(&p[0]);
    Long n_indx=mem_head.n_indx;
    Long n_elem=mem_head.n_elem;
    Long type_size=mem_head.type_size;
    char* base=(char*)&mem_head;

    { // Verify header check_sum; set array to init_mem_val
      #ifdef PVFMM_MEMDEBUG
      CheckMemHead(mem_head);
      Long size=mem_head.n_elem*mem_head.type_size;
      #pragma omp parallel for
      for(Long i=0;i<size;i++) p[i]=init_mem_val;
      for(Integer i=0;i<sizeof(MemHead);i++) base[i]=init_mem_val;
      #endif
    }

    if(n_indx==0){ // Use system free
      assert(base<&buff[0] || base>=&buff[buff_size]);
      char* p_;
      { // p_ <-- unalign(base)
        p_=(char*)((uintptr_t)base-((uint16_t*)base)[-1]);
      }
      #ifdef PVFMM_MEMDEBUG
      { // Check out-of-bounds write
        base[-1] = init_mem_val;
        base[-2] = init_mem_val;

        Long size=n_elem*type_size+header_size;
        size=(uintptr_t)(size+alignment) & ~(uintptr_t)alignment;
        Long end_padding=8; // to check for out-of-bound writes
        #pragma omp parallel for
        for(Long i=0;i<size+2+alignment+end_padding;i++) {
          PVFMM_ASSERT_MSG(p_[i]==init_mem_val, "memory corruption detected.");
        }
      }
      { // system_malloc.erase(p_)
        omp_set_lock(&omp_lock);
        PVFMM_ASSERT_MSG(system_malloc.erase(p_)==1, "double free or corruption.");
        omp_unset_lock(&omp_lock);
      }
      #endif
      return ::free(p_);
    }else{
      assert(n_indx<=node_buff.size());
      omp_set_lock(&omp_lock);
      MemNode& n=node_buff[n_indx-1];
      assert(!n.free && n.size>0 && n.mem_ptr==base);
      if(n.prev!=0 && node_buff[n.prev-1].free){
        Long n_prev_indx=n.prev;
        MemNode& n_prev=node_buff[n_prev_indx-1];
        n.size+=n_prev.size;
        n.mem_ptr=n_prev.mem_ptr;
        n.prev=n_prev.prev;
        free_map.erase(n_prev.it);
        delete_node(n_prev_indx);

        if(n.prev){
          node_buff[n.prev-1].next=n_indx;
        }
      }
      if(n.next!=0 && node_buff[n.next-1].free){
        Long n_next_indx=n.next;
        MemNode& n_next=node_buff[n_next_indx-1];
        n.size+=n_next.size;
        n.next=n_next.next;
        free_map.erase(n_next.it);
        delete_node(n_next_indx);

        if(n.next){
          node_buff[n.next-1].prev=n_indx;
        }
      }
      n.free=true; // Insert n to free_map
      n.it=free_map.insert(std::make_pair(n.size,n_indx));
      omp_unset_lock(&omp_lock);
    }
  }

  inline void MemoryManager::print() const{
    if(!buff_size) return;
    omp_set_lock(&omp_lock);

    Long size=0;
    Long largest_size=0;
    MemNode* n=&node_buff[n_dummy_indx-1];
    std::cout<<"\n|";
    while(n->next){
      n=&node_buff[n->next-1];
      if(n->free){
        std::cout<<' ';
        largest_size=std::max(largest_size,n->size);
      }
      else{
        std::cout<<'#';
        size+=n->size;
      }
    }
    std::cout<<"|  allocated="<<round(size*1000.0/buff_size)/10<<"%";
    std::cout<<"  largest_free="<<round(largest_size*1000.0/buff_size)/10<<"%\n";

    omp_unset_lock(&omp_lock);
  }

  inline void MemoryManager::test(){
    Long M=2000000000;
    { // With memory manager
      Long N=M*sizeof(double)*1.1;
      double tt;
      Iterator<double> tmp;

      std::cout<<"With memory manager: ";
      MemoryManager memgr(N);

      for(Integer j=0;j<3;j++){
        tmp=(Iterator<double>)memgr.malloc(M*sizeof(double));
        PVFMM_ASSERT(tmp!=NULL);
        tt=omp_get_wtime();
#pragma omp parallel for
        for(Long i=0;i<M;i+=64) tmp[i]=(double)i;
        tt=omp_get_wtime()-tt;
        std::cout<<tt<<' ';
        memgr.free((Iterator<char>)tmp);
      }
      std::cout<<'\n';
    }
    { // Without memory manager
      double tt;
      double* tmp;

      std::cout<<"Without memory manager: ";
      for(Integer j=0;j<3;j++){
        tmp=(double*)::malloc(M*sizeof(double));
        PVFMM_ASSERT(tmp!=NULL);
        tt=omp_get_wtime();
#pragma omp parallel for
        for(Long i=0;i<M;i+=64) tmp[i]=(double)i;
        tt=omp_get_wtime()-tt;
        std::cout<<tt<<' ';
        ::free(tmp);
      }
      std::cout<<'\n';
    }
  }

  inline void MemoryManager::Check() const{
#ifdef PVFMM_MEMDEBUG
    //print();
    omp_set_lock(&omp_lock);
    MemNode* curr_node=&node_buff[n_dummy_indx-1];
    while(curr_node->next){
      curr_node=&node_buff[curr_node->next-1];
      if(curr_node->free){
        char* base=curr_node->mem_ptr;
#pragma omp parallel for
        for(Long i=0;i<curr_node->size;i++){
          PVFMM_ASSERT_MSG(base[i]==init_mem_val, "memory corruption detected.");
        }
      }
    }
    omp_unset_lock(&omp_lock);
#endif
  }

  inline Long MemoryManager::new_node() const{
    if(node_stack.empty()){
      node_buff.resize(node_buff.size()+1);
      node_stack.push(node_buff.size());
    }

    Long indx=node_stack.top();
    node_stack.pop();
    assert(indx);
    return indx;
  }

  inline void MemoryManager::delete_node(Long indx) const{
    assert(indx);
    assert(indx<=node_buff.size());
    MemNode& n=node_buff[indx-1];
    n.free=false;
    n.size=0;
    n.prev=0;
    n.next=0;
    n.mem_ptr=NULL;
    node_stack.push(indx);
  }


  template <class ValueType> inline Iterator<ValueType> aligned_new(Long n_elem, const MemoryManager* mem_mgr){
    if(!n_elem) return NULL;

    static MemoryManager def_mem_mgr(0);
    if(!mem_mgr) mem_mgr=&def_mem_mgr;
    Iterator<ValueType> A=(Iterator<ValueType>)mem_mgr->malloc(n_elem, sizeof(ValueType));
    PVFMM_ASSERT_MSG(A!=NULL, "memory allocation failed.");

    if(!TypeTraits<ValueType>::IsPOD()){ // Call constructors
      //printf("%s\n", __PRETTY_FUNCTION__);
#pragma omp parallel for
      for(Long i=0;i<n_elem;i++){
        ValueType* Ai=new(&A[i]) ValueType();
        assert(Ai==(&A[i]));
      }
    }else{
#ifdef PVFMM_MEMDEBUG
      static Long random_init_val=1;
      Iterator<char> A_=(Iterator<char>)A;
#pragma omp parallel for
      for(Long i=0;i<n_elem*sizeof(ValueType);i++){
        A_[i]=random_init_val+i;
      }
      random_init_val+=n_elem*sizeof(ValueType);
#endif
    }

    return A;
  }

  template <class ValueType> inline void aligned_delete(Iterator<ValueType> A, const MemoryManager* mem_mgr){
    if (A==NULL) return;

    if(!TypeTraits<ValueType>::IsPOD()){ // Call destructors
      //printf("%s\n", __PRETTY_FUNCTION__);
      MemoryManager::MemHead& mem_head=MemoryManager::GetMemHead((char*)&A[0]);
#ifdef PVFMM_MEMDEBUG
      MemoryManager::CheckMemHead(mem_head);
      //PVFMM_ASSERT_MSG(mem_head.n_elem==1 || mem_head.type_id==TypeTraits<ValueType>::ID(), "pointer to aligned_delete has different type than what was used in aligned_new.");
#endif
      Long n_elem=mem_head.n_elem;
      for(Long i=0;i<n_elem;i++){
        A[i].~ValueType();
      }
    }else{
#ifdef PVFMM_MEMDEBUG
      MemoryManager::MemHead& mem_head=MemoryManager::GetMemHead((char*)&A[0]);
      MemoryManager::CheckMemHead(mem_head);
      //PVFMM_ASSERT_MSG(mem_head.type_id==TypeTraits<ValueType>::ID(), "pointer to aligned_delete has different type than what was used in aligned_new.");
      Long size=mem_head.n_elem*mem_head.type_size;
      Iterator<char> A_=(Iterator<char>)A;
#pragma omp parallel for
      for(Long i=0;i<size;i++){
        A_[i]=0;
      }
#endif
    }

    static MemoryManager def_mem_mgr(0);
    if(!mem_mgr) mem_mgr=&def_mem_mgr;
    mem_mgr->free((Iterator<char>)A);
  }

  template <class ValueType> inline Iterator<ValueType> memcopy(Iterator<ValueType> destination, ConstIterator<ValueType> source, Long num){
    if(destination!=source && num){
#ifdef PVFMM_MEMDEBUG
      destination[num-1];
      source[num-1];
#endif
      if(TypeTraits<ValueType>::IsPOD()){
        memcpy ( &destination[0], &source[0], num*sizeof(ValueType) );
      }else{
        for(Long i=0;i<num;i++) destination[i] = source[i];
      }
    }
    return destination;
  }

  template <class ValueType> inline Iterator<ValueType> memset( Iterator<ValueType> ptr, int value, Long num){
    if(num){
#ifdef PVFMM_MEMDEBUG
      ptr[0];
      ptr[num-1];
#endif
      ::memset ( &ptr[0], value, num*sizeof(ValueType) );
    }
    return ptr;
  }
}
