//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
// Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
//

// to compile mex -lblas quotientmex.cpp


#include "mex.h"
#include "MatLab.h"
#include "iostream"
#include "math.h"
#if defined(_WIN32)
inline double round(double d){
return floor(d+0.5);
}
#endif
#include "cQuotientCore.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if(nrhs==0) {
      std::cout << "quotientmex Quotient of a Separated Solution" << std::endl;
      std::cout << "" << std::endl;
      std::cout << "a = quotientmex(FF1,FF2)" << std::endl;
      std::cout << "a = quotientmex(FF1,FF2, ops)" << std::endl;
      std::cout << "" << std::endl;
      std::cout << "FF1 : is a cell containing the terms" << std::endl;
      std::cout << "FF2 : is a cell containing the terms" << std::endl;
      std::cout << "ops : are the opntion in the folowing order" << std::endl;
      std::cout << "[bool  reweight_modes" << std::endl;
      std::cout << " double fp_tol" << std::endl;
      std::cout << " double res_reduc" << std::endl;
      std::cout << " unsigned max_added_modes" << std::endl;
      std::cout << " unsigned fp_max_iter" << std::endl;
      std::cout << " bool improve_modes" << std::endl;
      std::cout << " unsigned improve_modes_max ]" << std::endl;
      return;
  }
    
  /* Check for proper number of arguments. */
  if(nrhs>3 || nrhs<2) {
    mexErrMsgIdAndTxt( "MATLAB:quotientmex:invalidNumInputs",
            "Two input required. Three optional");
  } else if(nlhs>1) {
    mexErrMsgIdAndTxt( "MATLAB:quotientmex:maxlhs",
            "Too many output arguments.");
  } else if (nlhs==0){
      mexErrMsgIdAndTxt( "MATLAB:quotientmex:maxlhs",
            "Too few output arguments.");
  }

  /* check if cell. */
  if (! mxIsCell (prhs[0]))
        mexErrMsgIdAndTxt( "MATLAB:recompactmex",
            "Input must be a cell.");
        
  /* check if cell. */
  if (! mxIsCell (prhs[1]))
        mexErrMsgIdAndTxt( "MATLAB:recompactmex",
            "Input must be a cell.");        
        
  mwSize nn = mxGetNumberOfElements (prhs[0]);
  mwSize nd = mxGetNumberOfElements (prhs[1]);
  
  PGD_Options options;
  
  /// we fill the FF with the incoming data
  std::vector<MatLabDataMatrix<double> > FF1;
  FF1.resize(nn);
  
  std::vector<MatLabDataMatrix<double> > FF2;
  FF2.resize(nd);
  for(int i =0; i < nn; ++i){
    FF1[i].SetInternalAllocation(false);
    FF1[i].SetNDims(2);
    mxArray *mat = mxGetCell(prhs[0], i);
    FF1[i].dsizes[0] =  mxGetM(mat);
    FF1[i].dsizes[1] =  mxGetN(mat);;
    FF1[i].allocate(); 
    FF1[i].Data = mxGetPr(mat);
   

    FF2[i].SetInternalAllocation(false);
    FF2[i].SetNDims(2);
    mat = mxGetCell(prhs[1], i);
    FF2[i].dsizes[0] =  mxGetM(mat);
    FF2[i].dsizes[1] =  mxGetN(mat);;
    FF2[i].allocate(); 
    FF2[i].Data = mxGetPr(mat);
  }
 
  // we prepare the sol 

  std::vector<MatLabDataMatrix<double> > sol;
  sol.resize(nn);
  
  options.max_added_modes = FF1[0].dsizes[1];  
  options.improve_modes_max = std::max(10, (int) round(FF1[0].dsizes[1]/10) );
  
  if(nrhs==3){
      if(mxGetNumberOfElements (prhs[2]) < 9){
              mexErrMsgIdAndTxt( "MATLAB:recompactmex:invalidNumInputs",
            "Options not of the correct size");
       };
    options.reweight_modes       = *(mxGetPr(prhs[2])+0);
    options.fp_tol               = *(mxGetPr(prhs[2])+1);  
    options.res_reduc            = *(mxGetPr(prhs[2])+2);
    options.max_added_modes      = *(mxGetPr(prhs[2])+3);
    options.fp_max_iter          = *(mxGetPr(prhs[2])+4);
    options.improve_modes        = *(mxGetPr(prhs[2])+5);
    options.improve_modes_max    = *(mxGetPr(prhs[2])+6);
    options.verbose              = *(mxGetPr(prhs[2])+7);
    options.lastImproveModesLoop = *(mxGetPr(prhs[2])+8); 
  }

  Quotient( FF1,FF2, sol, options);
  
  mxArray *msol = mxCreateCellArray (1, &nn);
  for(int i =0; i < nn; ++i){
        mxArray * s= mxCreateDoubleMatrix(sol[i].dsizes[0], sol[i].dsizes[1], mxREAL);
        std::copy(sol[i].Data,sol[i].Data+sol[i].fullsize,mxGetPr(s));
        mxSetCell(msol,i, s);
  }
  plhs[0] = msol;
}
