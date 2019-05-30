//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
// Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
//


#include <climits>      // Needed for memset
#include "cRecompactCore.h"
#include "MatLab.h"
#include <cstring>      // Needed for memset
#include "iostream"
#include <math.h>


#ifdef USE_MPI
  #include "mpi.h"
  //std::vector<double> mpitempvec;
  const MPI_Datatype TypeIsDouble<double>::MPI_TYPE = MPI_DOUBLE;
  
  /// this is work in progress
MPI_Datatype everydim;
MPI_Op myOp; 
// void mySum( void *in_, void *inout_, int *len, MPI_Datatype *dptr ) { 
//   int i; 
//   double *in = (double*)in_;
//   double *inout = (double*)inout_;
//   for (i=0; i< *len; ++i) { 
//     *inout =*in+*inout;
//     in++;
//     inout++; 
//   } 
// } 
template<typename T>
void mySum( void *in_, void *inout_, int *len, MPI_Datatype *dptr ) { 
  int i; 
  T *in = (T*)in_;
  T *inout = (T*)inout_;
  //printf( "len %d \n ",*len);
  
  int file_type_size;  // only to exlore the file_type
  MPI_Type_size(*dptr, &file_type_size);
  //printf( " file_type_size   = %d \n ",file_type_size );
      
  MPI_Aint file_type_extent; //only to exlore the file_type
  MPI_Type_extent(*dptr, &file_type_extent);
  //printf( " file_type_extent   = ‰d \n" ,file_type_extent);
  int step = file_type_extent/ file_type_size;
  int ll = file_type_size/sizeof(T);
  //printf( "step %d \n" ,step );
  //printf( "ll %d \n" ,ll );
  for (i=0; i< ll*(*len); ++i) { 
    //printf( "( %d , )" ,*in , *inout) ;
    *inout =*in+*inout;
    in = in + step;
    inout= inout +step; 
  } 
  //printf("\n");
} 
  
  
#endif
//to activate the random generator for c++11
//#include <random>
//
#include "ExtraFunctions.h"
//
template<typename T>
void UpdateTerm(const PGD_Options& options, unsigned& pfix_cpt,
                ptrdiff_t& nDims,
                unsigned int& nInitialModes,
                MatLabDataMatrix<T >& Ct,
                unsigned& enri_cpt,
                MatLabDataMatrix<T>& C,
                std::vector<T>& weight,
                std::vector<T>& wk,
                std::vector<T>& w,
                std::vector<ptrdiff_t>& dimSize,
                std::vector<MatLabDataMatrix<T> >&FF,
                std::vector<T*>& RS,
                std::vector<MatLabDataMatrix<T> >&sol,
                std::vector<std::vector<T> >& RS_old, bool print =  true);

template<typename T>
void UpdateTerm(const PGD_Options& options, unsigned& pfix_cpt,
                ptrdiff_t& nDims,
                unsigned int& nInitialModes,
                MatLabDataMatrix<T >& Ct,
                unsigned& enri_cpt,
                MatLabDataMatrix<T>& C,
                std::vector<T>& weight,
                std::vector<T>& wk,
                std::vector<T>& w,
                std::vector<ptrdiff_t>& dimSize,
                std::vector<MatLabDataMatrix<T> >&FF,
                std::vector<T*>& RS,
                std::vector<MatLabDataMatrix<T> >&sol,
                std::vector<std::vector<T> >& RS_old, bool print) {


    
  
    while(pfix_cpt < options.fp_max_iter) {
        if(print){
          if(options.verbose ==true) printf( "." );
          if(options.verbose ==true) std::cout.flush();
        }



        /// for each dimension
        for(ptrdiff_t dim =0; dim < nDims; ++dim ) {
            /// put ones in the row of the Ct matrix
            for(unsigned j = 0; j < nInitialModes; ++j) {
                Ct.Data[dim+j*nDims] = 1;
            }
            
            
            // puting ones (or alphas in the case of the last dimension) in the row of the C matrix
            //if( dim < nDims -1){
            // for(int i = 0; i < enri_cpt; ++i){
            //   C.data[dim+i*dim] = 1;
            // }
            //}else {
            for(unsigned k = 0; k < enri_cpt; ++k) {
                C.Data[dim+k*nDims] = weight[k];
            }
            
            //}

            /// calcul of the wk
            std::fill(wk.begin(), wk.end(),(T)1);
            for(ptrdiff_t i=0; i < nDims; ++i) {
                for(unsigned j=0; j < nInitialModes; ++j) {
                    wk[j] *= Ct.Data[j*nDims+i];
                }
            }
            
            /// calcul of the w
            std::fill(w.begin(), w.end(),(T)1);
            for(ptrdiff_t i=0; i < nDims; ++i) {
                for(unsigned k=0; k < enri_cpt; ++k) {
                    w[k] *= C.Data[k*nDims+i];
                }
            }

            /// constribution of FF 
            cblas_Xgemv(CblasColMajor,CblasNoTrans, dimSize[dim], (int)FF[dim].dsizes[1], (T)(1.), &FF[dim].Data[0],  dimSize[dim], &wk[0], 1, (T)0.0, RS[dim], 1 );
            /// constribution of the new solution
            
            cblas_Xgemv(CblasColMajor,CblasNoTrans, dimSize[dim], enri_cpt, (T)(-1.), &sol[dim].Data[0], dimSize[dim], &w[0], 1, (T)1.0, RS[dim], 1 );


            ///normalization
          #ifdef USE_MPI
            T norm_ = cblas_Xdot( dimSize[dim], RS[dim], (ptrdiff_t)1, RS[dim], (ptrdiff_t)1 );
            T norm;
            MPI_Allreduce (&norm_, &norm, 1, TypeIsDouble<T>::MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);     
            norm =  sqrt(norm);
          #else
            T norm = sqrt(cblas_Xdot( dimSize[dim], RS[dim], (ptrdiff_t)1, RS[dim], (ptrdiff_t)1 ));
          #endif /* MPI_VERSION */ 
            T factor = 1/norm;
            if( dim == nDims -1) {
                weight[enri_cpt] =  norm;
            }
            cblas_Xscal(dimSize[dim],factor,RS[dim] ,(ptrdiff_t)1);

            ///update of C Ct
            update_C(Ct,RS[dim], FF[dim],FF[dim].dsizes[1] ,dim);
            update_C(C,RS[dim], sol[dim],sol[dim].dsizes[1],dim);
            
        }/// end for every dimension
        T error = 2;
        for(ptrdiff_t i=0; i < nDims; ++i) {
          #ifdef USE_MPI
            T localError = cblas_Xdot( dimSize[i], RS[i], (ptrdiff_t)1, &RS_old[i][0], (ptrdiff_t)1 );
            MPI_Allreduce (&localError, &localError, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);     
            error *= localError;
          #else
            error *= cblas_Xdot( dimSize[i], RS[i], (ptrdiff_t)1, &RS_old[i][0], (ptrdiff_t)1 );
          #endif  
        }

        if (options.fp_tol > (2-error)) {
            if(print)
              if(options.verbose ==true) printf("CR %f \n",weight[enri_cpt]/weight[0]);
            break;
        } else {
            for(ptrdiff_t i =0; i < nDims; ++i) {
                for(int j =0; j < dimSize[i]; ++j) {
                    RS_old[i][j] = *(RS[i]+j);
                }
            }
        }

        ++pfix_cpt;
    }/// end for every iteration



}



void DummyFunciton(){

  std::vector<MatLabDataMatrix<double> > FFD;
  std::vector<MatLabDataMatrix<double> > solD;
  Recompact(FFD,solD, PGD_Options());

  std::vector<MatLabDataMatrix<float> > FFS;
  std::vector<MatLabDataMatrix<float> > solS;
  Recompact(FFS,solS, PGD_Options());
}


template<typename T>
void Recompact(std::vector<MatLabDataMatrix<T> >& FF,std::vector<MatLabDataMatrix<T> >& sol,const PGD_Options& options) {

    ptrdiff_t nDims = FF.size();

  #ifdef USE_MPI

    int array_of_sizes[2] = {(int)nDims, 1};
    int array_of_subsize[2] = {1,1} ;
    int array_of_starts[2] = {0,0};
    MPI_Type_create_subarray(2, array_of_sizes, array_of_subsize, 
        array_of_starts, MPI_ORDER_FORTRAN, TypeIsDouble<T>::MPI_TYPE, &everydim) ;
    MPI_Type_commit(&everydim);
  
    MPI_Op_create( mySum<T>, true, &myOp ); 
  
  #endif /* USE_MPI */ 

    
    if(options.verbose == true) printf( " *** cRecompact *** \n");
    options.PrintOptions();
    if(options.verbose == true) printf( "     NumberOfDims  : %d", nDims );
    if (nDims == 0){ 
      printf( "Nothing to do exiting...\n");
      return;
      
    }

    
    unsigned nInitialModes = FF[0].dsizes[1];
    if(options.verbose ==true) printf( "     nInitialModes : %d \n" ,nInitialModes);

    std::vector<ptrdiff_t> dimSize;
    dimSize.resize(nDims);

    for(ptrdiff_t i =0; i < nDims; ++i) {
        dimSize[i] = FF[i].dsizes[0];
        if(options.verbose ==true) printf( "     Size of dim ( %d ) : %d \n ", i ,dimSize[i]);
    }
    
    if (nDims == 1){ 
        printf( "Only One dimension. Solution with only one mode.\n" );
        sol.resize(nDims);
        sol[0].SetNDims(2);
        sol[0].dsizes[0] = dimSize[0];
        sol[0].dsizes[1] = 1; 
        sol[0].allocate();
        std::fill(sol[0].Data, sol[0].Data+dimSize[0],(T)0.0);
        // sum of the terms
        for(ptrdiff_t j =0; j < FF[0].dsizes[1]; ++j){
          for(ptrdiff_t i =0; i < FF[0].dsizes[0]; ++i){
            *(sol[0].Data+i) = *(sol[0].Data+i) +  *(FF[0].Data+i+j*FF[0].dsizes[0]);
          }
        }
        
        return;
    }
    
    std::vector<T> wk;
    wk.resize(nInitialModes); // vector to store the prod(Ct,1)

    std::vector<T> w;
    w.resize(options.max_added_modes); // vector to store the prod(C,1)


    std::vector<T> weight;
    weight.resize(options.max_added_modes+1); // we need and extra space to work for the weight when improve_modes is one
    std::fill(weight.begin(),weight.end(),(T)1.);

    /// solution allocation

    sol.resize(nDims);
    for(ptrdiff_t i =0; i < nDims; ++i) {
        sol[i].SetNDims(2);
        sol[i].dsizes[0] = dimSize[i];
        sol[i].dsizes[1] = options.max_added_modes;
        sol[i].allocate();
        sol[i].dsizes[1] = 0;
        sol[i].fullsize = 0;
        std::fill(sol[i].Data, sol[i].Data+options.max_added_modes*dimSize[i],(T)0.0);
    }


    /// matrix to calculate the weight
    //MatLabDataMatrix<T> kweight;
    //kweight.SetNDims(2);
    //kweight.dsizes[0] = options.max_added_modes;;
    //kweight.dsizes[1] = options.max_added_modes;
    //kweight.allocate();


    // pointer to RS
    std::vector<T*> RS;
    RS.resize(nDims);

    std::vector<std::vector<T> > RS_old;
    RS_old.resize(nDims);
    for(ptrdiff_t i =0; i < nDims; ++i) {
        RS_old[i].resize(dimSize[i]);
        std::fill(RS_old[i].begin(), RS_old[i].end(),(T)0);
    }

    MatLabDataMatrix<T> C;
    C.SetNDims(2);
    C.dsizes[0] = nDims;
    C.dsizes[1] = options.max_added_modes;
    C.allocate();
    C.dsizes[1] =0;

    MatLabDataMatrix<T> Ct;
    Ct.SetNDims(2);
    Ct.dsizes[0] = nDims;
    Ct.dsizes[1] = nInitialModes;
    Ct.allocate();

    //
    T weight0;

    //
#ifdef  _RANDOM_H
    std::default_random_engine generator;
    std::normal_distribution<T> distribution(0.0,1.0);
#endif
    //T one =1;
    // for each enrichment
    unsigned enri_cpt = 0;
    while(enri_cpt < options.max_added_modes) {
        if(options.verbose ==true) printf( " E %d ",  enri_cpt) ;

        for(ptrdiff_t i =0; i < nDims; ++i) {
            RS[i] = sol[i].Data+dimSize[i]*enri_cpt;
        }


        for(ptrdiff_t i =0; i < nDims; ++i) {
            for(int j=0; j < dimSize[i]; ++j) {
                
                #ifdef  _RANDOM_H
                *(RS[i]+j) = distribution(generator);
                #else
                *(RS[i]+j) = rand() / (T)RAND_MAX;
                #endif
                //*(RS[i]+j) = 1;
                //printf( " %f", *(RS[i]+j) );
            }
        }

        for(ptrdiff_t i = 0; i < nDims; ++i) {
            update_C(Ct,RS[i], FF[i],FF[i].dsizes[1] ,i);
            update_C(C,RS[i], sol[i],sol[i].dsizes[1],i);
        }


        /// for each point fix iteration

        unsigned pfix_cpt=0;

        UpdateTerm( options,  pfix_cpt,
                    nDims,
                    nInitialModes,
                    Ct,
                    enri_cpt,
                    C,
                    weight,
                    wk,
                    w,
                    dimSize,
                    FF,
                    RS,
                    sol,
                    RS_old);


        if(enri_cpt == 0) {
            weight0 = weight[0];
        }

        if (pfix_cpt == options.fp_max_iter ) {
            if(options.verbose == true) printf( "CNR %f \n", weight[enri_cpt]);
        }
        if(weight[enri_cpt] != weight[enri_cpt]) {
             if(options.verbose == true) printf( " NAN detected Quiting...\n");
            break;
        }
        ++enri_cpt;
        // we take into acount the new mode
        for(ptrdiff_t i =0; i < nDims; ++i) {
            sol[i].dsizes[1] = enri_cpt;
            sol[i].fullsize = sol[i].dsizes[0]*sol[i].dsizes[1];
        }
#ifdef MATLAB_MEX_FILE   
        if(utIsInterruptPending()) {
            utSetInterruptPending(false);
            mexPrintf("Ctrl-C Detected. END\n\n");
            break;
        }
#endif /* MATLAB_MEX_FILE */

        if(weight[enri_cpt-1]/weight0<options.res_reduc) {
            break;
        }

        if(options.improve_modes) {
          int from_modes_to_improves = std::max(signed (0),signed(enri_cpt)-signed(options.improve_modes_max));
          int to_modes_to_improves = enri_cpt-1;

          if(from_modes_to_improves != to_modes_to_improves ){
            if(options.verbose ==true) printf( "improvement of modes %d  to %d " , from_modes_to_improves, to_modes_to_improves) ;//<< std::endl;
            for(int jj= 0; jj <std::max(1,options.lastImproveModesLoop*(enri_cpt==options.max_added_modes)); ++jj){
            unsigned pfix_cpt_in_max =0;
            
            for(int improve = from_modes_to_improves; improve <= to_modes_to_improves ; ++improve ) {
                //we recover the pointer to the rs to be enriched
                for(ptrdiff_t i =0; i < nDims; ++i) {
                    RS[i] = &sol[i].Data[dimSize[i]*improve];
                }
                // set the weight for that rs equal to zero
                weight[improve]=0;
                for(ptrdiff_t i =0; i < nDims; ++i) {
                    update_C(Ct,RS[i], FF[i],FF[i].dsizes[1] ,i);
                }

                unsigned pfix_cpt_in=0;
                UpdateTerm( options,  pfix_cpt_in,
                            nDims,
                            nInitialModes,
                            Ct,
                            enri_cpt,
                            C,
                            weight,
                            wk,
                            w,
                            dimSize,
                            FF,
                            RS,
                            sol,
                            RS_old , false);
                if(options.verbose ==true){
                    if (pfix_cpt_in == options.fp_max_iter ) {
                         printf( "*" );
                         std::cout.flush();
                    } else {
                        printf( "." );
                        std::cout.flush();
                    }
                }
                
                pfix_cpt_in_max = std::max(pfix_cpt_in_max,pfix_cpt_in);
                
                // we put back the correct weight to the term
                weight[improve]=weight[enri_cpt];
                weight[enri_cpt] = 1.0;

            }
            if(options.verbose ==true) printf( "\n");
            
            //if(pfix_cpt_in_max <= 2 ) break;
            }
            }
        }

    }/// end for every enrichmenet

    // to put the correct size in the matrix (the allocated vector is bigger)
    for(ptrdiff_t i =0; i < nDims; ++i) {
        sol[i].fullsize = sol[i].dsizes[0]*sol[i].dsizes[1];
    }
    
    // rescale putting the weight in the first dimension
    if(options.verbose ==true) printf( " weight are : \n"  );
    for(unsigned a =0; a < enri_cpt; ++a) {
        if(options.verbose ==true) printf( "%f  ", weight[a] ) ;
        cblas_Xscal(dimSize[0],weight[a],& sol[0](0,a) ,(ptrdiff_t)1);
    }
    if(options.verbose ==true) printf("\n");

};




