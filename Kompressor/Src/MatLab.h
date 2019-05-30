//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
// Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
//

#include <vector>
#if !defined(_WIN32)
#include <sys/socket.h> // Needed for the socket functions
#include <netdb.h>      // Needed for the socket functions
#endif
#include <cstdlib>      // for exit 

#include <string>
#include <fstream>

#ifndef MatLab_H
#define MatLab_H

#if defined(_WIN32)

#else
#include <stddef.h>
#endif 
 
class MatlabData{
public:
 enum {Double, Char, Int8, Int16, Int32, Uint8, Uint16, Uint32};
 std::string typeName;
 size_t ndims;
 std::vector<ptrdiff_t> dsizes;
 unsigned fullsize;
 void SetType(const std::string& t){ 
	 typeName = t;
 };
 void SetNDims(const unsigned & n){
 	ndims = n;  
        dsizes.resize(n);
        std::fill (dsizes.begin(),dsizes.end(),0);
 }
#if !defined(_WIN32)
  void recv_base(int new_sd){
  typeName = "ABCD";
  recv(new_sd,&typeName[0],4,0);
  unsigned int n= 0;
  recv(new_sd,(char *)&n,sizeof(unsigned int),0);
  if (n > 3){exit(1);}
  this->SetNDims(n);
  recv(new_sd,&dsizes[0],sizeof(int)*n,0);
 };
  void send_base(int new_sd){
    send(new_sd, "double\n", 7, 0);
    send(new_sd,(char *)& ndims,sizeof(unsigned int),0);
    send(new_sd,&dsizes[0],sizeof(int)*ndims,0);
  };
#endif
};

template<class T>
class MatLabDataMatrix: public MatlabData{
  bool internal_allocation;
  std::vector<T> internalStorage;
public:
 MatLabDataMatrix(){
   internal_allocation = false;
   SetInternalAllocation(true);
 }
 T *Data;
 void SetInternalAllocation(bool ia){
   
 if(internal_allocation == ia)  return;
   if(!ia)  {
      Data = NULL;
      internalStorage.resize(0);
   };
 internal_allocation = ia;
 }
 
 
 void allocate(){
  unsigned d=1;
  for(unsigned i=0; i < ndims; ++i){
	d*=dsizes[i];
	}
  fullsize = d;
  if(internal_allocation){
    internalStorage.resize(d);
    Data = &internalStorage[0];
    std::fill (internalStorage.begin(),internalStorage.end(),T(0));
  }
 }
 void recv(int new_sd){
   recv_base(new_sd);
   this->allocate();
   ::recv(new_sd,(void *) this->Data,sizeof(T)*fullsize,MSG_WAITALL);
 }
 void send(int new_sd){
    send_base(new_sd);
    ::send(new_sd,(void *) this->Data,sizeof(T)*fullsize,0);
 }
 T& operator()(unsigned i,unsigned j){
    return Data[dsizes[0]*j+i];
   
  }
  void ReadFromFStream(std::fstream& infile){
    
    this->SetNDims(2);
    infile.read ((char*)&this->dsizes[0], sizeof(int));
    infile.read ((char*)&this->dsizes[1], sizeof(int));
    this->allocate();
    infile.read ((char*) this->Data, sizeof(T)*this->fullsize);
  }
  int LD(){ return dsizes[0]; };
  
};

#endif //MatLab_H
