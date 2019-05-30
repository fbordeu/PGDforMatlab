//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
// Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
//

#include "PGD_Options.h"

PGD_Options::PGD_Options() {
        residual_minimization = false;
        reweight_modes = true;
        fp_tol = 1e-8;
        res_reduc = 1e-8;
        max_added_modes = 10000;
        fp_max_iter= 200;
        improve_modes = true;
        //improve_modes_first = 1;
        improve_modes_max = INT_MAX;
        //iter_dims = zeros(1,dim);
        //Niter_dims = zeros(1,dim);
        //heavy_dims = zeros(1,dim);
        //heavy_weight = 1;
        //BB_fct_update_frequency = 0;
        lastImproveModesLoop = 1;
        verbose = true;
    }
//
void PGD_Options::PrintOptions()const {
if(this->verbose ==true) std::cout << "Using options : " << std::endl;
if(this->verbose ==true) std::cout << " residual_minimization : " << residual_minimization << std::endl;
if(this->verbose ==true) std::cout << " reweight_modes    : " << reweight_modes << std::endl;
if(this->verbose ==true) std::cout << " fp_tol            : " << fp_tol << std::endl;
if(this->verbose ==true) std::cout << " res_reduc         : " << res_reduc << std::endl;
if(this->verbose ==true) std::cout << " max_added_modes   : " << max_added_modes << std::endl;
if(this->verbose ==true) std::cout << " fp_max_iter       : " << fp_max_iter << std::endl;
if(this->verbose ==true) std::cout << " improve_modes     : " << improve_modes << std::endl;
if(this->verbose ==true) std::cout << " improve_modes_max : " << improve_modes_max << std::endl;
if(this->verbose ==true) std::cout << " lastImproveModesLoop : " << lastImproveModesLoop << std::endl;
if(this->verbose ==true) std::cout << " verbose           : " << verbose << std::endl;

};