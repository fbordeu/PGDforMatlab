//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
// Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
//

// to compile mex -lblas recompactmex.cpp


#include "mex.h"
#include "MatLab.h"
#include "iostream"
#include "math.h"
#if defined(_WIN32)
inline double round(double d){
return floor(d+0.5);
}
#endif
#include "cRecompactCore.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if(nrhs==0) {
      printf("recompactmex Recompact a Separated Solution \n");
      printf(" \n");
      printf("a = recompactmex(FF) \n");
      printf("a = recompactmex(FF, ops) \n");
      printf(" \n");
      printf("FF : is a cell containing the terms \n");
      printf("ops : are the opntion in the folowing order \n");
      printf("[bool  reweight_modes \n");
      printf(" double fp_tol \n");
      printf(" double res_reduc \n");
      printf(" unsigned max_added_modes \n");
      printf(" unsigned fp_max_iter \n");
      printf(" bool improve_modes \n");
      printf(" unsigned improve_modes_max ] \n");
      printf(" bool verbose ] \n");
      printf(" unsigned lastImproveModesLoop ] \n");
      
      return;
  }
    
  /* Check for proper number of arguments. */
  if(nrhs>2) {
    mexErrMsgIdAndTxt( "MATLAB:recompactmex:invalidNumInputs",
            "One input required. Two optional");
  } else if(nlhs>1) {
    mexErrMsgIdAndTxt( "MATLAB:recompactmex:maxlhs",
            "Too many output arguments.");
  } else if (nlhs==0){
      mexErrMsgIdAndTxt( "MATLAB:recompactmex:maxlhs",
            "Too few output arguments.");
  }

  /* check if cell. */
  if (! mxIsCell (prhs[0]))
        mexErrMsgIdAndTxt( "MATLAB:recompactmex",
            "Input must be a cell.");
        
  mwSize n = mxGetNumberOfElements (prhs[0]);
  
  PGD_Options options;
  
  /// we fill the FF with the incoming data
  std::vector<MatLabDataMatrix<double> > FF;
  FF.resize(n);
  for(int i =0; i < n; ++i){
    FF[i].SetInternalAllocation(false);
    FF[i].SetNDims(2);
    
    mxArray *mat = mxGetCell(prhs[0], i);
    
    FF[i].dsizes[0] =  mxGetM(mat);
    FF[i].dsizes[1] =  mxGetN(mat);;
    FF[i].allocate(); 
    FF[i].Data = mxGetPr(mat);
  }
 
  // we prepare the sol 

  std::vector<MatLabDataMatrix<double> > sol;
  sol.resize(n);
  
  options.max_added_modes = FF[0].dsizes[1];  
  options.improve_modes_max = std::max(10, (int) round(FF[0].dsizes[1]/10) );
  
  if(nrhs==2){
      if(mxGetNumberOfElements (prhs[1]) < 9){
              mexErrMsgIdAndTxt( "MATLAB:recompactmex:invalidNumInputs",
            "Options not of the correct size");
       };
    options.reweight_modes       = (*(mxGetPr(prhs[1])+0) != 0)? true : false;
    options.fp_tol               = *(mxGetPr(prhs[1])+1);  
    options.res_reduc            = *(mxGetPr(prhs[1])+2);
    options.max_added_modes      = (unsigned)*(mxGetPr(prhs[1])+3);
    options.fp_max_iter          = (unsigned)*(mxGetPr(prhs[1])+4);
    options.improve_modes        = (*(mxGetPr(prhs[1])+5) !=0 )? true : false;
    options.improve_modes_max    = (unsigned)*(mxGetPr(prhs[1])+6);
    options.verbose              = (*(mxGetPr(prhs[1])+7)!=0);
    options.lastImproveModesLoop = (unsigned)*(mxGetPr(prhs[1])+8); 
  }

  Recompact( FF, sol, options);
  mxArray *msol = mxCreateCellArray (1, &n);
  for(int i =0; i < n; ++i){
        mxArray * s= mxCreateDoubleMatrix(sol[i].dsizes[0], sol[i].dsizes[1], mxREAL);
        std::copy(sol[i].Data,sol[i].Data+sol[i].fullsize,mxGetPr(s));
        mxSetCell(msol,i, s);
  }
  plhs[0] = msol;
}
