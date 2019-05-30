//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
// Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
//

#include <iostream> 

#ifdef MATLAB_MEX_FILE 
#include "mex.h"
#ifdef __cplusplus 
    extern "C" bool utIsInterruptPending();
    extern "C" void utSetInterruptPending(bool);
#else
    extern bool utIsInterruptPending();
    extern void utSetInterruptPending(bool);
#endif
#endif /* MATLAB_MEX_FILE */

//if in matlab we use the matlab blas library
#ifdef MATLAB_MEX_FILE 
    #include <blas.h>
    inline double cblas_ddot( ptrdiff_t& a, double*& b,const int& c, double* d,const int& e){
        ptrdiff_t cc = c;
        ptrdiff_t ee = e;
        return ddot(&a,b,&cc,d,&ee);
   }
   enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102 };
   enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113,
                         AtlasConj=114};
    
template<typename T5,typename T6,typename T7,typename T8,typename T9,typename T12>
void cblas_dgemv(const CBLAS_ORDER& a, const CBLAS_TRANSPOSE& b, const ptrdiff_t & c, const ptrdiff_t& d,const T5& e,T6 f,T7& g,T8 h,const T9& i,const double& j,double* k,const T12& l){
    char bb= 'n';
    if (b ==CblasTrans) bb = 't';
    ptrdiff_t cc = c;
    ptrdiff_t dd = d;
    double ee = e;
    ptrdiff_t ii = i;
    double jj = j;
    ptrdiff_t ll  =  l;
    dgemv( &bb,  &cc,    &dd,    &ee,    f,    &g,    h,    &ii,    &jj,    k,   &ll);
};
      
template<typename T5,typename T6,typename T7, typename T8,typename T9,typename T12>
void cblas_sgemv(const CBLAS_ORDER& a, const CBLAS_TRANSPOSE& b, const ptrdiff_t& c, const ptrdiff_t& d, const T5& e,T6 f,T7& g,T8 h,const  T9& i,const float& j,float* k,const T12& l){
    char bb= 'n';
    if (b ==CblasTrans) bb = 't';
    ptrdiff_t cc = c;
    ptrdiff_t dd = d;
    float ee = e;
    ptrdiff_t ii = i;
    float jj = j;
    ptrdiff_t ll  =  l;

    sgemv( &bb,  &cc,    &dd,    &ee,    f,    &g,    h,    &ii,    &jj,    k,   &ll);
};

 void cblas_dgemm(const CBLAS_ORDER& a, const CBLAS_TRANSPOSE& b, const CBLAS_TRANSPOSE& c, ptrdiff_t d, ptrdiff_t e, ptrdiff_t f ,double g, double* h, ptrdiff_t i,double* j , ptrdiff_t k,double l,double* m,ptrdiff_t n){
    char bb= 'n';
    if (b ==CblasTrans) bb = 't';
    char cc= 'n';
    if (c ==CblasTrans) cc = 't';

    dgemm( &bb, &cc, &d, &e, &f, &g, h, &i, j, &k, &l,m, &n);
}

template <typename pdiff>
void  cblas_dscal (pdiff& a,double& b, double* c,const pdiff& d){
    pdiff dd = d;
    dscal(&a,&b,c,&dd);
}

template <typename pdiff>
void  cblas_sscal (pdiff& a,float& b, float* c,const pdiff& d){
    pdiff dd = d;
    sscal(&a,&b,c,&dd);
}

template <typename pdiff>
inline float cblas_sdot( pdiff& a, float*& b,const pdiff& c, float* d,const pdiff& e){
    pdiff cc = c;
    pdiff ee = e;
    return sdot(&a,b,&cc,d,&ee);
}
#else
// if not matlab we use the 'system' blas

#ifdef __APPLE__
   #include   <Accelerate/Accelerate.h>
#else
#if defined(_WIN32)
   /// for now no compilation can be done in windows without matlab
#else
   /// linux
   #include <cblas.h>
   #include <vector>
#endif
#endif
#endif /* MATLAB_MEX_FILE */


template <typename ptrdiff, typename T>
T cblas_Xdot( ptrdiff&, T*&,const ptrdiff&, T*,const ptrdiff& );


template <typename ptrdiff>
double cblas_Xdot( ptrdiff& a, double*& b,const ptrdiff& c, double* d,const ptrdiff& e){
return cblas_ddot(a,b,c,d,e);
}

template <typename ptrdiff>
float cblas_Xdot( ptrdiff& a, float*& b,const ptrdiff& c, float* d,const ptrdiff& e){
return cblas_sdot(a,b,c,d,e);
}


//template<typename T>
//void cblas_Xgemv(const CBLAS_ORDER&  , const CBLAS_TRANSPOSE&  ,int&   ,int&  ,const T&       ,     T*  , int&  ,T*      , const int&   ,const T&      , T*       ,const int&);
//

template <typename ptrdiff>
void cblas_Xgemv(const CBLAS_ORDER& a, const CBLAS_TRANSPOSE& b,const int& c,const int& d,const double& e,double* f, ptrdiff& g,double* h,const int& i ,const double& j,double* k,const int& l){
            cblas_dgemv(a,b,c,d,e,f,g,h,i,j,k,l);
};

template <typename ptrdiff>
void cblas_Xgemv(const CBLAS_ORDER& a,const CBLAS_TRANSPOSE& b,const int& c,const int& d,const float& e,float* f,ptrdiff& g,float* h,const  int& i,const float& j,float* k,const int& l){
            cblas_sgemv(a,b,c,d,e,f,g,h,i,j,k,l);
};


//template<typename T1,typename T2, typename T3,typename T4> 
//void cblas_Xscal(T1&,T2&,T3,const T4&);

//template<typename T1> 
void  cblas_Xscal(ptrdiff_t& a,double& b,double* c,const ptrdiff_t& d){
        cblas_dscal(a,b,c ,d);
}

//template<typename T1> 
void  cblas_Xscal(ptrdiff_t& a,float& b, float* c,const ptrdiff_t& d){
        cblas_sscal(a,b,c,d);
}

template <typename T>
void update_C(MatLabDataMatrix <T >& C,T*  r, MatLabDataMatrix <T >& MM,ptrdiff_t& N, const ptrdiff_t& dim);

template<>
void update_C(MatLabDataMatrix <double >& C,double*  r, MatLabDataMatrix <double >& MM,ptrdiff_t& NN, const ptrdiff_t& dim) {
    int N = NN;
    int M = MM.dsizes[0];
    int INCY = C.dsizes[0];
    cblas_dgemv(CblasColMajor,CblasTrans, M, N, 1.0, &MM.Data[0], MM.dsizes[0],      r, 1,      0, &C.Data[dim], C.dsizes[0] );
    #ifdef USE_MPI
       MPI_Allreduce(MPI_IN_PLACE, &C.Data[dim], NN, everydim, myOp, MPI_COMM_WORLD);  
    #endif /* USE_MPI */ 
};

template<>
void update_C(MatLabDataMatrix <float >& C,float*  r, MatLabDataMatrix <float >& MM,ptrdiff_t& NN, const ptrdiff_t& dim) {
    int N = NN;
    int M = MM.dsizes[0];
    int INCY = C.dsizes[0];
    cblas_sgemv(CblasColMajor,CblasTrans, M, N, (float)1.0, &MM.Data[0], MM.dsizes[0],      r, 1,      0, &C.Data[dim], C.dsizes[0] );
    #ifdef USE_MPI
       MPI_Allreduce(MPI_IN_PLACE, &C.Data[dim], NN, everydim, myOp, MPI_COMM_WORLD);   
    #endif /* USE_MPI */
};

template<typename T>
void update_D(MatLabDataMatrix <T >& D, MatLabDataMatrix <T >& FF, T* r){
  int cpt = -1;
  for(int j =0; j < FF.dsizes[1]; ++j){
    for(int i =0; i < FF.dsizes[0]; ++i){  
        D.Data[++cpt] = FF(i,j)*r[i];
    }
  }
};




//template <typename T>
//void update_CWithFF2(MatLabDataMatrix <T >& C,T*  r, MatLabDataMatrix <T >& MM);
//,unsigned& N, const unsigned& dim, MatLabDataMatrix <T >& DD, std::vector<T> &workingplace);

//template<>
void update_CU(std::vector<double >& C,MatLabDataMatrix <double >&  sol, MatLabDataMatrix <double >& D) {
  
    if( sol.dsizes[1] == 0) return;
  
    int Sr,Sc,Dr,Dc;
    Sr = sol.dsizes[0];
    Sc = sol.dsizes[1];
    
    Dr = D.dsizes[0];
    Dc = D.dsizes[1];
    
    ///C(Cr,Cc) = sol(Sr,Sc)'*D(Dr,Dc);
    
    
    //printf (" Initializing data for matrix multiplication C=sol'*D for matrix \n"
    //        " sol(%ix%i) and matrix D(%ix%i)\n\n", Sr, Sc, Dr, Dc);
    
    //std::cout << "sol.size() : " << sol.fullsize << " =? " <<Sr*Sc << std::endl;  
    //std::cout << "D.size() : " << D.fullsize << " =? " <<Dr*Dc << std::endl;  
    
    
   //std::cout << "--------------------------------------" << std::endl;
   
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 
                Sc,Dc, Dr, 1.0, sol.Data, sol.LD(), D.Data , D.LD(), 0, &C[0], sol.dsizes[1]);
    
    
    //std::cout << " ------------------------- DONE ------------------------" << std::endl;
    #ifdef USE_MPI
       MPI_Allreduce(MPI_IN_PLACE, &C[0], C.size(), everydim, myOp, MPI_COMM_WORLD);  
    #endif /* USE_MPI */ 
};
//
void update_CU(std::vector<float >& C,MatLabDataMatrix <float >&  sol, MatLabDataMatrix <float >& D) {
//#warning to correct for the double implementation 
    //cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, sol.dsizes[1], D.dsizes[1], sol.dsizes[0], 1.0, sol.Data, sol.dsizes[0], D.Data , D.dsizes[0], 0, &C[0],sol.dsizes[0]);
    #ifdef USE_MPI
       MPI_Allreduce(MPI_IN_PLACE, &C[0], C.size(), everydim, myOp, MPI_COMM_WORLD);  
    #endif /* USE_MPI */ 
};

#define PrintM(X) std::cout << #X " : "  << std::endl << "  "; \
                  for(int ii =0; ii < X.dsizes[0] ; ++ii){\
                    for(int jj =0; jj < X.dsizes[1] ; ++jj){\
                      std::cout << X(ii,jj) << " ";\
                    }\
                    std::cout << std::endl << "  "; \
                  }\
                  std::cout << std::endl;
                  
#define PrintV(X) std::cout << #X " : "  << std::endl << "  "; \
                  for(int ii=0; ii < X.size(); ++ii){\
                    std::cout << X[ii] << " "; \
                  } \
                  std::cout << std::endl;
                  
// template<>
// void update_CWithFF2(MatLabDataMatrix <float >& C,float*  r, MatLabDataMatrix <float >& MM,unsigned& NN, const unsigned& dim, MatLabDataMatrix <float >& DD, std::vector<float> &workingplace) {    int N = NN;
//     ///TODO function to be coded
//     int M = MM.dsizes[0];
//     int LDA = MM.dsizes[0];
//     int INCY = C.dsizes[0];
//     cblas_sgemv(CblasColMajor,CblasTrans, M, N, 1.0, &MM.Data[0], LDA,      r, 1,      0, &C.Data[dim], C.dsizes[0] );
//     #ifdef USE_MPI
//        MPI_Allreduce(MPI_IN_PLACE, &C.Data[dim], NN, everydim, myOp, MPI_COMM_WORLD);   
//     #endif /* USE_MPI */ 
//  
// };

