//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
// Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
//

#ifndef __PGD_OPTIONES__
#define __PGD_OPTIONES__

#include <climits>
#include <iostream>

struct PGD_Options {
    bool  residual_minimization;
    bool  reweight_modes;
    double fp_tol;
    double res_reduc;
    unsigned max_added_modes;
    unsigned fp_max_iter;
    bool improve_modes;
    //unsigned improve_modes_first;
    unsigned improve_modes_max;
    //iter_dims = zeros(1,dim);
    //Niter_dims = zeros(1,dim);
    //heavy_dims = zeros(1,dim);
    //heavy_weight = 1;
    //BB_fct_update_frequency = 0;
    int lastImproveModesLoop;
    bool verbose;
    PGD_Options();
    void PrintOptions()const ;
};

#endif //__PGD_OPTIONES__

