//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
// Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
//


#include <vector>
#include "MatLab.h"


#ifdef USE_MPI
#include "mpi.h"
#endif

#include "PGD_Options.h"

template< class T >
struct TypeIsDouble{
    static const bool value = false;
    #ifdef USE_MPI
    static const MPI_Datatype MPI_TYPE;
    #endif
};

#ifdef USE_MPI
template< class T >
const MPI_Datatype TypeIsDouble<T>::MPI_TYPE = MPI_FLOAT;
#endif

template<>
struct TypeIsDouble< double >{
    static const bool value = true;
    #ifdef USE_MPI
    static const MPI_Datatype MPI_TYPE ;
    #endif
};

template <typename T>
void Quotient(std::vector<MatLabDataMatrix<T> >& FF1, std::vector<MatLabDataMatrix<T> >& FF2, std::vector<MatLabDataMatrix<T> >& sol,const  PGD_Options& options);
