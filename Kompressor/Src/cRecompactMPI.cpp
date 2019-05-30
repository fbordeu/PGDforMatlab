//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
// Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
//

#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <multithread_save.h>
#include <multithread_save.cpp>
 
#include "MatLab.h"
#include "cRecompactCore.h"

#define MAIN
#include "CommandLineParser.h"




template<typename T>    
void wrapper(std::string& infile_name, xpgd::CommandLineParser& maindata,MPI_DATA& mpi_data, std::string& outfile_name){
  
  
  MPI::File infile;
  infile = MPI::File::Open(MPI_COMM_WORLD, infile_name.c_str(), MPI_MODE_RDONLY,MPI::INFO_NULL);
  
  ///***** Reading Data ****************************************************
  std::cout << "Reading file from disk" << std::endl;
  
  /// number of dimensions
  int ndims;
  infile.Read_at(0,(&ndims), 1, MPI_INT);
  MPI::Offset offset=sizeof(int);
  std::cout << "Number of dims " << ndims << std::endl;
  
  /// distributed version of FF
  std::vector<MatLabDataMatrix<T> > FF;
  FF.resize(ndims);
  
  /// Global size of FF
  std::vector<std::vector<int> > FFGlobalSizes; // s1 and s2
  FFGlobalSizes.resize(ndims);
      
  std::vector<std::vector<int> > FFLocalSizes; // s1  for every rank
  FFLocalSizes.resize(ndims);
  std::vector<MPI_Datatype> file_type; // the vector with the size of each block
  file_type.resize(ndims);
  
  
  
  for (int i =0;i  < ndims; i++){
      FF[i].SetNDims(2);
      FFGlobalSizes[i].resize(2);
      FFLocalSizes[i].resize(mpi_data.Get_size());
      
      /// the partition  ********************************************
      /// we cut only the columns 
      std::vector<int> array_of_psizes;
      array_of_psizes.resize(2);
      array_of_psizes[0]= mpi_data.Get_size();
      array_of_psizes[1]= 1;

      std::cout << "array_of_psizes "  << array_of_psizes[0] << " " << array_of_psizes[1] << std::endl;
  
      int array_of_dargs[2] = {MPI_DISTRIBUTE_DFLT_DARG,  MPI_DISTRIBUTE_DFLT_DARG};
      int array_of_distribs[2] =  { MPI_DISTRIBUTE_BLOCK,  MPI_DISTRIBUTE_NONE};

      ///************************************************************
  
      int s1;
      int s2;
      std::cout << " ---------------------------------- " << std::endl;
  //for(int i =0; i < ndims; ++i){
      std::cout << "offset " << offset << std::endl;
      MPI_File_set_view(infile, offset, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
      infile.Read(&s1, 1, MPI_INT);
      offset+= (1)*sizeof(int);
      
      std::cout << "offset " << offset << std::endl;
      MPI_File_set_view(infile, offset, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
      infile.Read(&s2, 1, MPI_INT);
      offset+= (1)*sizeof(int);
    
    std::cout << " original size is " << s1 << " , " << s2 << std::endl;  
    
    FFGlobalSizes[i][0] = s1;
    FFGlobalSizes[i][1] = s2;
    
    /// partitioning the file into the different ranks
    /// generation of the partition for every matrix

    MPI_Type_create_darray( mpi_data.Get_size(),
                            mpi_data.Get_rank(),
                           2,
                           FFGlobalSizes[i].data(),
                           array_of_distribs,
                           array_of_dargs,
                           array_of_psizes.data(),
                           MPI_ORDER_FORTRAN,
                           TypeIsDouble<T>::MPI_TYPE,
                           &(file_type[i]));
    MPI_Type_commit(&(file_type[i]));
    
    // Explore the returned type 
    //MPI_Aint file_type_extent; //only to exlore the file_type
    //MPI_Type_extent(file_type[i], &file_type_extent);
    int file_type_size;  // only to exlore the file_type
    MPI_Type_size(file_type[i], &file_type_size);
    //std::cout << " file_type_size   = " << file_type_size << std::endl;
    //std::cout << " file_type_extent   = " << file_type_extent << std::endl;
    s1 = file_type_size/sizeof(T)/s2;
    std::cout <<  " Size is " << s1 << std::endl;
    
    /// allocation for the data
    FF[i].dsizes[0] = s1;
    FF[i].dsizes[1] = s2; 
    FF[i].allocate(); 
    
  //}
      
  //int offset = (1+2*ndims)*sizeof(int);
  
  //for(int i =0; i < ndims; ++i){
     
     std::cout << "Reading dim " << i << std::endl;
     
     //std::cout << "offset " <<  offset << std::endl;
     MPI_File_set_view(infile, offset, TypeIsDouble<T>::MPI_TYPE, file_type[i], "native", MPI_INFO_NULL);
     //std::cout << "reading " << FFGlobalSizes[i][1] << " elements" << std::endl;
     
     MPI::Status status;
     infile.Read_at(0,(char*) FF[i].Data, FF[i].fullsize, TypeIsDouble<T>::MPI_TYPE, status);
     offset+= FFGlobalSizes[i][0]*FFGlobalSizes[i][1]*sizeof(T);
     //int count;
     //count  = status.Get_count(TypeIsDouble<T>::MPI_TYPE);
     //std::cout << " *** read "  <<  count << " double " << std::endl;;
     
  }
  infile.Close();
      
  std::vector<MatLabDataMatrix<T> > sol;
    
  ///***** Recompact ************************************************************************
  
  ///work in progress
  //MPI_Type_vector( NN,1, INCY, MPI_DOUBLE, &everydim);
  //MPI_Type_commit(&everydim);
  

 
  PGD_Options options; 
  if(maindata.get_option("+NIM"))
        options.improve_modes = 0;
  
  if(maindata.has_atribute("-MAM")){
     maindata.get_atribute("-MAM", options.max_added_modes);
  } else {
     options.max_added_modes = FF[0].dsizes[1];  
  }
   
  if(maindata.has_atribute("-lastImproveModesLoop")){
     maindata.get_atribute("-lastImproveModesLoop", options.lastImproveModesLoop);
   };
  maindata.get_atribute("-Tol", options.res_reduc);
    
  options.improve_modes_max = std::max(10, (int) round(FF[0].dsizes[1]/10) );

  Recompact( FF, sol, options);
  ///***** Writing data  ********************************************************************
  std::cout << "Writing data to disk" << std::endl;
  std::cout << "  to file : " <<  outfile_name.c_str()  << std::endl;
   
  MPI::File outfile;
  outfile = MPI::File::Open(MPI_COMM_WORLD, outfile_name.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE ,MPI::INFO_NULL);
  
  
  outfile.Write(&ndims, 1, MPI_INT);
  offset = (1)*sizeof(int);
  for(int i =0; i < ndims ; i++){
    std::cout << "Writing dim " << i << std::endl;
    
    MPI_File_set_view(outfile, offset, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
    outfile.Write(&FFGlobalSizes[i][0], 1, MPI_INT);
    offset += 1*sizeof(int);
    
    MPI_File_set_view(outfile, offset, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
    outfile.Write(&sol[i].dsizes[1], 1, MPI_INT);
    offset += 1*sizeof(int); 

    MPI_File_set_view(outfile, offset, TypeIsDouble<T>::MPI_TYPE, file_type[i], "native", MPI_INFO_NULL);

    MPI::Status status;
    outfile.Write_at(0,(char*) sol[i].Data, sol[i].fullsize, TypeIsDouble<T>::MPI_TYPE, status);
    offset+= FFGlobalSizes[i][0]*sol[i].dsizes[1]*sizeof(T);
 }
  
  outfile.Close();
  
}

int main(int argc, char *argv[]){
  
  /// MPI init ****************************************************************
  MPI_DATA mpi_data; 
  mpi_data.Init(argc, argv);

  ///***** Parsing Command line ************************************************
  xpgd::CommandLineParser maindata(argc,argv);   

  /// Log file name **********************************************************
  std::string coutfile;
  coutfile = "cRecompactMPI";
  std::ostringstream temp;
  temp << "_" << mpi_data.Get_rank();
  temp << ".log";
  coutfile += temp.str();
  
  maindata.coutToFile(coutfile);
  
  maindata.addHelp("file1", "File with the field to recompact");
  maindata.addHelp("file2", "Output file name (if not present file out_file1 is created)");
  maindata.addHelp("+NIM", "No Improve Modes (default activated)");
  maindata.addHelp("-MAM", "Max Added Modes ");
  maindata.addHelp("+float", "to treat all the data as floats");
  maindata.addHelp("-Tol", "Exit tolerance (default 1e-8) ");
  maindata.addHelp("-lastImproveModesLoop", "How many times we improve modes after the last enrichement");
  maindata.doneReading();
  
  /// input file ****************************************************************
  std::string infile_name;
  if(!maindata.get_file(infile_name)) {
    std::cout << "ERROR: I need exactly a file " << std::endl; 
    exit(1);
  }
  /// output file name generation ****************************************************************
  size_t pos = infile_name.find_last_of('/');
  std::string outfile_name;
  if( pos > std::string::npos )
     outfile_name = "out_" + infile_name;
  else{ 
     outfile_name = infile_name.substr(0,pos+1);
     outfile_name += "out_" + infile_name.substr(pos+1,infile_name.size());
  }

  maindata.get_file(outfile_name,1);  
  std::cout << " output file " << outfile_name << std::endl;
  

  if(maindata.get_option("+float")){
    wrapper<float>(infile_name,maindata, mpi_data,  outfile_name);
  } else {
    wrapper<double>(infile_name,maindata, mpi_data, outfile_name);
  }
  
  
  mpi_data.Finalize();;
    
  return 0;
  //maindata destruction
}


