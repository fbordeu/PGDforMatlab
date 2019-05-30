//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
// Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
//


#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstring>

#include "MatLab.h"

#include "cQuotientCore.h"
#define MAIN
#include "CommandLineParser.h"

template<typename T>
void wrapper(std::fstream& infile, xpgd::CommandLineParser& maindata, PGD_Options& options, std::string& outfile_name){

      ///***** Reading Data ****************************************************
      std::cout << "Reading file from disk" << std::endl;
      
      /// number of dimensions
      int ndims;
      infile.read ((char*)&ndims, sizeof(int));
      std::cout << "Number of dims " << ndims << std::endl;
      std::vector<MatLabDataMatrix<T> > FF1;
      FF1.resize(ndims);
      
      for(int i =0; i < ndims; ++i){
        std::cout << "Reading dim " << i << std::endl;
        FF1[i].ReadFromFStream(infile);
      }
      
      /// number of dimensions
      infile.read ((char*)&ndims, sizeof(int));
      std::cout << "Number of dims " << ndims << std::endl;
      std::vector<MatLabDataMatrix<T> > FF2;
      FF2.resize(ndims);
      
      for(int i =0; i < ndims; ++i){
        std::cout << "Reading dim " << i << std::endl;
        FF2[i].ReadFromFStream(infile);
      }
      
      infile.close();
      
      std::vector<MatLabDataMatrix<T> > sol;
      
      ///***** Recompact ************************************************************************
      
      if(maindata.get_option("+NIM"))
        options.improve_modes = 0;
      if(maindata.has_atribute("-MAM")){
        maindata.get_atribute("-MAM", options.max_added_modes);
      } else {
        options.max_added_modes = std::max(FF1[0].dsizes[1],FF2[0].dsizes[1]);  
      }
      
      if(maindata.has_atribute("-lastImproveModesLoop")){
        maindata.get_atribute("-lastImproveModesLoop", options.lastImproveModesLoop);
      };
          
      maindata.get_atribute("-Tol", options.res_reduc);
      
      options.improve_modes_max = std::max(10, (int) round(std::max(FF1[0].dsizes[1],FF2[0].dsizes[1])/10) );
      
      Quotient( FF1,FF2, sol, options);
      ///***** Writing data  ************************************************************************
      std::cout << "Writing data to disk" << std::endl;
      std::cout << "  to file : " <<  outfile_name.c_str()  << std::endl;
      
      std::fstream outfile;
      outfile.open (outfile_name.c_str(), std::ios::binary | std::ios::out);
      
      outfile.write ((char*)&ndims, sizeof(int)); 
      
      for(int i =0; i < ndims; ++i){
        outfile.write ((char*)&sol[i].dsizes[0], sizeof(int));
        outfile.write ((char*)&sol[i].dsizes[1], sizeof(int));
      //}
       
      //for(int i =0; i < ndims; ++i){
         std::cout << "Writing dim " << i << std::endl;
         outfile.write((char*) sol[i].Data, sizeof(T)*sol[i].fullsize);
      }
      outfile.close();

}


int main( int argc, char**argv ){
  
  PGD_Options options;
  
  std::cout << "test" << argc << std::endl;
  ///***** Parsing Command line *****
  xpgd::CommandLineParser maindata(argc,argv);   
 
  maindata.addHelp("file1", "File with the field to recompact");
  maindata.addHelp("file2", "Output file name (if not present file out_file1 is created)");
  maindata.addHelp("+NIM", "No Improve Modes (default activated)");
  maindata.addHelp("-MAM", "Max Added Modes ");
  maindata.addHelp("+float", "to treat all the data as floats");
  maindata.addHelp("-Tol", "Exit tolerance (default 1e-8) ");
  maindata.addHelp("-lastImproveModesLoop", "How many times we improve modes after the last enrichement");
  maindata.doneReading(); 
  
  
  std::string infile_name;
  if(!maindata.get_file(infile_name)) {
    std::cout << "ERROR: I need exactly a file " << std::endl; 
    exit(1);
  }
  
  size_t pos = infile_name.find_last_of('/');
  std::string outfile_name;
  if( pos > std::string::npos )
     outfile_name = "out_" + infile_name;
  else{ 
     outfile_name = infile_name.substr(0,pos+1);
     outfile_name += "out_" + infile_name.substr(pos+1,infile_name.size());
  }

  maindata.get_file(outfile_name,1);  
  
  
  ///***** opening file *****
  std::fstream infile;
  infile.open (infile_name.c_str(), std::ios::binary | std::ios::in);
  
  if (infile.is_open()) {
    
        if(maindata.get_option("+float")){
        wrapper<float>(infile,maindata, options,  outfile_name);
        } else {
        wrapper<double>(infile,maindata, options, outfile_name);
        }
        
      
  } else {
      std::cout << "ERROR: Error opening file " << infile_name << std::endl;
  }
  return 0;
}