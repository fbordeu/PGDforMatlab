//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
// Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
//


#include <vector>
#include <math.h>
#include "cQuotientCore.h"


//#ifdef __APPLE__
//#include   <Accelerate/Accelerate.h>
//#else
//#if defined(_WIN32)
//   #include <blas.h>
//#else
//#include <cblas.h>
//#endif
//#endif
#include "ExtraFunctions.h"

#include <math.h>

template<typename T>
void UpdateTermWithFF2(const PGD_Options& options, unsigned& pfix_cpt,
                       ptrdiff_t& nDims,
                       unsigned int& nInitialModes,
                       unsigned& enri_cpt,
                       std::vector<ptrdiff_t>& dimSize,
                       std::vector<MatLabDataMatrix<T> >&FF1,
                       std::vector<MatLabDataMatrix<T> >&FF2,
                       
                       MatLabDataMatrix<T >& Ck,
                       std::vector<T>& wk,
                       std::vector< MatLabDataMatrix< T > > D,
                       std::vector< std::vector< T > > Cu,
                       std::vector< T > wu,
                       MatLabDataMatrix<T>& Cl,
                       std::vector<T>& wl,
                       

                       std::vector<T*>& RS,
                       std::vector<std::vector<T> >& RS_old,
                       std::vector<MatLabDataMatrix<T> >&sol,
                       std::vector<T>& weight
                       , bool print =  true ) {
  
   
    std::vector<T> workingplace;
    for(int i = 0; i < nDims; ++i) {
      workingplace.resize(std::max(dimSize[i],(ptrdiff_t)workingplace.size()));
    }
    
    int NumberOfOps = FF2[0].dsizes[1];
    int Dterms = enri_cpt*NumberOfOps;

    while(pfix_cpt < options.fp_max_iter) {
        if(print) {
            if(options.verbose ==true) std::cout << "." ;
            if(options.verbose ==true) std::cout.flush();
        }
        
        
        weight[enri_cpt] = 1;
        /// for each dimension
        for(unsigned dim =0; dim < nDims; ++dim ) {
            /// put ones in the row of the Ct matrix
            for(int j = 0; j < nInitialModes; ++j) {
                Ck.Data[dim+j*nDims] = 1;//weight[enri_cpt];
            }

            /// puting ones (or alphas in the case of the last dimension) in the row of the C matrix
            for (int l = 0; l < NumberOfOps; ++l) {
              for(int k = 0; k < enri_cpt; ++k) {
                    Cu[dim][l*enri_cpt+k] = weight[k];
              }
              //std::fill(Cu[dim].begin()+l*enri_cpt,Cu[dim].begin()+(l+1)*enri_cpt+,weight[l])
              Cl.Data[dim+l*nDims] = weight[enri_cpt];
            }

            /// calcul of the wk & wl
            std::fill(wk.begin(), wk.end(),1);
            std::fill(wl.begin(), wl.end(),1);
            for(int i=0; i < nDims; ++i) {
                for(int j=0; j < nInitialModes; ++j) {
                    wk[j] *= Ck.Data[j*nDims+i];
                }
                for(int j=0; j < NumberOfOps; ++j) {
                    wl[j] *= Cl.Data[j*nDims+i];
                }
            }
            /// calcul of the wu
            std::fill(wu.begin(), wu.end(),1);
            for(int i=0; i < nDims; ++i) {
                for(int k=0; k < Dterms; ++k) {
                    wu[k] *= Cu[i][k];
                }
            }
            
            /// contribution of FF1 (FF1(dim)*wk
            cblas_Xgemv(CblasColMajor,CblasNoTrans, dimSize[dim], FF1[dim].dsizes[1], (T)1., FF1[dim].Data,  dimSize[dim], &wk[0], 1, (T)0.0, RS[dim], 1 );

            
            /// contribution of the operators
            /// contribution of the new solution  (-D(dim)*wu)
            int one =1;
            for( int i=0; i <enri_cpt ; ++i){
                for( int j=0; j <NumberOfOps ; ++j){
                  for( int k=0; k <dimSize[dim] ; ++k){
                    *(RS[dim]+k) -= wu[enri_cpt*j+i]*FF2[dim](k,j)*sol[dim](k,i);
                    
                  }
                }
            }
            //cblas_Xgemv(CblasColMajor,CblasNoTrans, dimSize[dim], Dterms, (T)-1., D[dim].Data, dimSize[dim], &wu[0], 1, (T)1.0, RS[dim], 1 );

            
            ///division by the diag of the operator
            cblas_Xgemv(CblasColMajor,CblasNoTrans, dimSize[dim], FF2[dim].dsizes[1], (T)1., FF2[dim].Data,  dimSize[dim], &wl[0], 1, (T)0.0, &workingplace[0], 1 );

            for(int i =0; i < dimSize[dim];++i){
              RS[dim][i] /= workingplace[i];
            }

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
            weight[enri_cpt] *=  norm;
            cblas_Xscal(dimSize[dim],factor,RS[dim] ,(ptrdiff_t)1);
            
            ///update of Ck,D,Cl,Cu
            update_C(Ck,RS[dim], FF1[dim],FF1[dim].dsizes[1] ,dim);
            update_D(D[dim], FF2[dim], RS[dim]); // OK
            
            update_C(Cl,RS[dim],D[dim],D[dim].dsizes[1] ,dim);
            update_CU(Cu[dim],sol[dim],D[dim]);
        }  /// end for every dimension
        T error = 2;
        for(int i=0; i < nDims; ++i) {
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
                if(options.verbose ==true) std::cout  << "CR " << weight[enri_cpt]/weight[0] << std::endl;
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


void DummyFunciton() {

    std::vector<MatLabDataMatrix<double> > FFD1;
    std::vector<MatLabDataMatrix<double> > FFD2;
    std::vector<MatLabDataMatrix<double> > solD;
    Quotient(FFD1,FFD2,solD, PGD_Options());

    std::vector<MatLabDataMatrix<float> > FFS1;
    std::vector<MatLabDataMatrix<float> > FFS2;
    std::vector<MatLabDataMatrix<float> > solS;
    Quotient(FFS1,FFS2,solS, PGD_Options());
}


template<typename T>
void Quotient(std::vector<MatLabDataMatrix<T> >& FF1,std::vector<MatLabDataMatrix<T> >& FF2,std::vector<MatLabDataMatrix<T> >& sol,const PGD_Options& options) {

    ptrdiff_t nDims = FF1.size();
    ptrdiff_t dDims = FF2.size();

    if(options.verbose == true) std::cout << " *** cQuotient ***"  << std::endl;
    options.PrintOptions();
    if(options.verbose == true) std::cout << "     NumberOfDims  : " << nDims << std::endl;
    if(options.verbose == true) std::cout << "     NumberOfDims  : " << dDims << std::endl;
    if (nDims == 0 || dDims == 0 ) {
        std::cout << "Nothing to do exiting..." << std::endl;
        return;
    }

    if (dDims !=nDims) {
        std::cout << "ERROR: N D don't have the same number of dimensions..." << std::endl;
        return;
    }


    unsigned nInitialModes = FF1[0].dsizes[1];
    if(options.verbose == true) std::cout << "     nInitialModes : " << nInitialModes  << std::endl;

    unsigned dInitialModes = FF2[0].dsizes[1];
    if(options.verbose == true) std::cout << "     dInitialModes : " << dInitialModes  << std::endl;

    std::vector<ptrdiff_t> dimSize;
    dimSize.resize(nDims);

    for(ptrdiff_t i =0; i < nDims; ++i) {
        dimSize[i] = FF1[i].dsizes[0];
        if(options.verbose ==true) std::cout << "     Size of dim (" << i<< " ) : " << dimSize[i]  << std::endl;
        if (FF1[i].dsizes[0] != FF2[i].dsizes[0] ) {
            std::cout << "ERROR: dimension " << i << " in N D don't have the same number of values." << std::endl;
            return ;
        }
    }

    std::vector<T> weight;
    std::vector<T*> RS;
    std::vector<std::vector<T> > RS_old;
    MatLabDataMatrix<T> Cl;
    MatLabDataMatrix<T> Ck;
    std::vector<MatLabDataMatrix<T> > D;
    std::vector<std::vector<T> > Cu;
    std::vector<T> wl;
    std::vector<T> wk;
    std::vector<T> wu;
    T weight0;
    
  try{
    
    
    weight.resize((options.max_added_modes+1)*dInitialModes); // we need and extra space to work for the weight when "improve modes" is one
    std::fill(weight.begin(),weight.end(),1);


    /// solution allocation
    sol.resize(nDims);
    for(ptrdiff_t i =0; i < nDims; ++i) {
        sol[i].SetNDims(2);
        sol[i].dsizes[0] = dimSize[i];
        sol[i].dsizes[1] = options.max_added_modes;
        sol[i].allocate();
        std::fill(sol[i].Data, sol[i].Data+sol[i].fullsize,0.0);
        sol[i].dsizes[1] = 0;
        sol[i].fullsize = 0;
    }

    /// RS allocation
    // pointer to RS
    
    RS.resize(nDims);

    
    RS_old.resize(nDims);
    for(ptrdiff_t i =0; i < nDims; ++i) {
        RS_old[i].resize(dimSize[i]);
        std::fill(RS_old[i].begin(), RS_old[i].end(),0);
    }

    /// to treat the left hand side contribution
    
    Cl.SetNDims(2);
    Cl.dsizes[0] = nDims;
    Cl.dsizes[1] = dInitialModes;
    Cl.allocate();
    //Cl.dsizes[1] =0;
    
    
    wl.resize(dInitialModes); // vector to store the prod(Cl,1)

    /// to treat the contribution of the nominator
    
    Ck.SetNDims(2);
    Ck.dsizes[0] = nDims;
    Ck.dsizes[1] = nInitialModes;
    Ck.allocate();

    
    wk.resize(nInitialModes); // vector to store the prod(Ck,1)

    /// to treat the contribution of the denominator
    
    D.resize(nDims);
    for(int i = 0; i < nDims; ++i) {
        D[i].SetNDims(2);
        D[i].dsizes[0] = FF2[i].dsizes[0];
        D[i].dsizes[1] = FF2[i].dsizes[1];
        D[i].allocate();
    }

    //MatLabDataMatrix<T> Cu;
    //Cu.SetNDims(2);
    //Cu.dsizes[0] = nDims;
    //Cu.dsizes[1] = options.max_added_modes*dInitialModes;
    //Cu.allocate();

    
    Cu.resize(nDims);
    for(ptrdiff_t i =0; i < nDims; ++i) {
      Cu[i].resize(options.max_added_modes*dInitialModes);
    }
    
    
    wu.resize(options.max_added_modes*dInitialModes); // vector to store the prod(Cu,1)

    
    
  } catch(std::bad_alloc&) {
    std::cout << "ERROR : Not Enough Mermory to treat this problem" << std::endl ;
    for(size_t i =0; i < sol.size(); ++i) {
      sol[i].dsizes[0] = dimSize[i];
      sol[i].dsizes[1] = 0;
      sol[i].allocate();
    }
    return;
  }
  
    //
#ifdef  _RANDOM_H
    std::default_random_engine generator;
    std::normal_distribution<T> distribution(0.0,1.0);
#endif
    //for(ptrdiff_t i =0; i < nDims; ++i) {
      //PrintM(FF1[i])
      //PrintM(FF2[i])
    //}
    
    unsigned enri_cpt = 0;
    while(enri_cpt < options.max_added_modes) {
        if(options.verbose ==true) std::cout << " E " << enri_cpt << " " ;

        /// Pointer to the solution vector
        for(ptrdiff_t i =0; i < nDims; ++i) {
            RS[i] = sol[i].Data+dimSize[i]*enri_cpt;
        }

        /// Random init
        for(ptrdiff_t i =0; i < nDims; ++i) {
            for(int j=0; j < dimSize[i]; ++j) {
#ifdef  _RANDOM_H
                *(RS[i]+j) = distribution(generator);
#else
                *(RS[i]+j) = rand() / (T)RAND_MAX;
#endif
            }
        }

        /// Update Ck for every dimension 
        /// Update D  for every dimension 
        /// Update Cl  for every dimension 
        /// Update Cu  for every dimension 
        for(int i = 1; i < nDims; ++i) {
            update_C(Ck,RS[i], FF1[i],FF1[i].dsizes[1] ,i);
            update_D(D[i], FF2[i], RS[i]); // OK
            update_C(Cl,RS[i],D[i],D[i].dsizes[1] ,i);
            update_CU(Cu[i],sol[i],D[i]);
        }

        /// for each point fix iteration
        unsigned pfix_cpt=0;

        UpdateTermWithFF2( options,  pfix_cpt,
                           nDims,
                           nInitialModes,
                           enri_cpt,
                           dimSize,
                           FF1,
                           FF2,
                           
                           Ck,
                           wk,
                           D,
                           Cu,
                           wu,
                           Cl,
                           wl,

                           RS,
                           RS_old,
                           sol,
                           weight
                           );

        if(enri_cpt == 0) {
            weight0 = weight[0];
        }

        if (pfix_cpt == options.fp_max_iter ) {
            if(options.verbose == true) std::cout << "CNR "<< weight[enri_cpt] << std::endl;
        }
        if(weight[enri_cpt] != weight[enri_cpt]) {
            if(options.verbose == true) std::cout << " NAN detected Quiting..." << std::endl;
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
            if(options.verbose ==true) std::cout << "improvement of modes " << from_modes_to_improves <<" to "<< to_modes_to_improves ;//<< std::endl;
            for(int jj= 0; jj <std::max(1,options.lastImproveModesLoop*(enri_cpt==options.max_added_modes)); ++jj){
            unsigned pfix_cpt_in_max =0;
            
            for(int improve = from_modes_to_improves; improve <= to_modes_to_improves ; ++improve ) {
                //we recover the pointer to the rs to be enriched
                for(ptrdiff_t i =0; i < nDims; ++i) {
                    RS[i] = &sol[i].Data[dimSize[i]*improve];
                }
                // set the weight for that rs equal to zero
                weight[improve]=0;
                for(int i =1; i < nDims; ++i) {
                    update_C(Ck,RS[i], FF1[i],FF1[i].dsizes[1] ,i);
                    update_D(D[i], FF2[i], RS[i]); // OK
                    update_C(Cl,RS[i],D[i],D[i].dsizes[1] ,i);
                    update_CU(Cu[i],sol[i],D[i]);
                }

                unsigned pfix_cpt_in=0;
               
                        UpdateTermWithFF2( options,  pfix_cpt_in,
                           nDims,
                           nInitialModes,
                           enri_cpt,
                           dimSize,
                           FF1,
                           FF2,
                           
                           Ck,
                           wk,
                           D,
                           Cu,
                           wu,
                           Cl,
                           wl,

                           RS,
                           RS_old,
                           sol,
                           weight
                           , false);
                
                if(options.verbose ==true){
                    if (pfix_cpt_in == options.fp_max_iter ) {
                         std::cout << "*" ;
                         std::cout.flush();
                    } else {
                        std::cout << ".";
                        std::cout.flush();
                    }
                }
                
                pfix_cpt_in_max = std::max(pfix_cpt_in_max,pfix_cpt_in);
                
                // we put back the correct weight to the term
                weight[improve]=weight[enri_cpt];
                weight[enri_cpt] = 1.0;

            }
            if(options.verbose ==true) std::cout << std::endl;
            
            //if(pfix_cpt_in_max <= 2 ) break;
            }
            }
        }
        
        
        

    }/// end for every enrichmenet

    /// to put the correct size in the matrix (the allocated vector is bigger)
    for(ptrdiff_t i =0; i < nDims; ++i) {
        sol[i].fullsize = sol[i].dsizes[0]*sol[i].dsizes[1];
        //PrintM(sol[i]);
    }
    /// rescale putting the weight in the first dimension
    if(options.verbose ==true) std::cout << " weight are : "  << std::endl;
    for(size_t a =0; a < enri_cpt; ++a) {
        if(options.verbose ==true) std::cout << weight[a] << " " ;
        cblas_Xscal(dimSize[0],weight[a],& sol[0](0,a) ,(ptrdiff_t)1);
    }
    if(options.verbose ==true) std::cout << std::endl;
}















