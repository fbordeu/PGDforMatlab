//
// C++ Interface: inputoutput
//
// Description: Class to read the option of the command line
//
//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
// Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
//
#ifndef GEM_inputoutput_h
#define GEM_inputoutput_h

#include "sstream"
#include "iostream"
#include <vector>
#include <cstring>
#include <cassert>

namespace xpgd{
//
/// function to transform numbers to string
template<class T>
std::string tostr(const T& input){
  std::ostringstream out;
  out << input;
  return out.str();
}
//
template<class T>
std::string tostr(const std::vector<T>& input){
  std::ostringstream out;
  for(int  i=0; i < input.size(); ++i ){
    if(i != 0 )
      out << " ";
    out << input[i];  
  }
  return out.str();
}
//
template<class T>
std::istringstream &line_input ( std::istringstream &s, T &val ){
    s >> val;
    return s;
}
//
template<typename T1, typename T2>
std::istringstream &line_input ( std::istringstream &s, std::pair<T1, T2> &val ){
  
    
    char str2[2000];
    std::strcpy(str2,s.str().c_str());
  
    char * pch;
    pch = strtok (str2,"(), ");
    bool first=1;
    while (pch != NULL){
      std::istringstream istreamtemp ( pch );
      //T val;
      if(first) {
	line_input ( istreamtemp, val.first ) ;
	first =0;
      }else {
	line_input ( istreamtemp, val.second ) ;
      };
      //vec.push_back(val);
      pch = strtok (NULL, "(), ");
    
    }
    return s;
}
//
template<typename T>
void input (const char* str, T &val ){
    std::istringstream istreamtemp ( str );
    line_input ( istreamtemp, val );
}
//
template<typename T>
void input (const char* str, std::vector<T> &vec ){
  unsigned intsize = vec.size();
  
  //std::string line(str);
  char str2[2000];
  std::strcpy(str2,str);
  
  char * pch;
  pch = strtok (str2," ");
  unsigned cpt = 0;
  while (pch != NULL)
  {
   
    std::istringstream istreamtemp ( pch );
    T val;
    line_input ( istreamtemp, val );
    if(intsize){
      vec[cpt++] = val;
      if(cpt == intsize) break;
    } else {
      vec.push_back(val);
    }
    pch = strtok (NULL, " ");
  }
   //if(istreamtemp.good()) std::cout << "Error reading the string" << std::endl;
}
//
}
#endif //