function [FF] = PGD_SYM(AA,BB,GG,N_NT,Dirichlet,Max_terms,Max_fp_iter,epsilon,epsilon_tilde)
    %General purpose PGD function with residual minimization, not optimized for efficiency
    %This specific implementation relies on the fact that the residual
    %minimization algorithm can be implemented as the original PGD
    %algorithm (as in PGD.m) applied on a symetrized problem.
    %Outputs:
    %       FF :  A cell array of matrices. The d^th element is a matrix whose
    %       columns contains all the enrichment functions computed for that
    %       dimension.
    %Inputs:
    %       AA: cell array of matrices. We assume that the differential operator
    %       is written as a sum of simple multidimensional operators that can be
    %       written as the tensorial product of operators along each of the PGD
    %       dimensions.
    %       The i^th column of AA contains the discretized form (matrix) of
    %       those operator whose tensorial product yield the i^th simple
    %       multidimensional operator.
    %
    %       BB: separated form of the right hand side (structure similar to FF)
    %       GG: separated form of the a priori known terms of the PGD solution (structure similar to FF)
    %       N_NT: cell array of matrices that are used to compute the
    %       appropriate L2 norm along each PGD dimension
    %       Dirichlet: cell array of vectors indicating for each dimension the
    %       nodal values to which one has to apply homogeneous Dirichlet
    %       boundary conditions
    %       Max_terms: Maximum number of enrichments
    %       Max_fp_iter: Maximum number of fixed point iterations
    %       epsilon: termination criterion for the fixed point (Eq. (2.8))
    %       epsilon_tilde: termination criterion used for the enrichment process (Eq. (2.26))
    %
    %
    %Copyright (c) 2013, Francisco Chinesta (Ecole Centrale de Nantes), Roland Keunings (Universite catholique de Louvain), Adrien Leygue(CNRS)
    %Author: Adrien Leygue.
    %All rights reserved.
    %See License file for more details
    
    
    %Definition of some useful quantities
    %Number of separated PGD dimensions
    Ndims = size(AA,1);
    %Number of simple multidimensional operators
    Nops = size(AA,2);
    %Number of variables along each PGD dimension
    Nvars = zeros(1,Ndims);
    for d = 1:Ndims
        Nvars(d) = size(AA{d,1},1);
    end
    %mask is used identify the variables subject to Dirichlet BC
    mask =  cellfun(@(dummy) true(dummy,1),num2cell(Nvars),'uniformoutput',0);
    for d = 1:Ndims
        mask{d}(Dirichlet{d}) = false;
    end
    
    %We transform the problem to have the residual minimization form with
    %no Dirichlet BC: all known enrichments are passed to te right hand
    %side and the dofs subject to a homogeneous Dirichlet BC are removed
    %from the problem
    
    %New problem initialization
    AA_new = cell(Ndims,Nops^2);
    BB_new = cell(Ndims,1);
    BB_tmp = cell(Ndims,1);
    GG_new = cell(Ndims,1);
    N_NT_new = cell(Ndims,1);
    Dirichlet_new = cell(Ndims,1);
    
    %removal of Dirichlet dofs in BB and N_NT
    for d = 1:Ndims
        BB_tmp{d} = BB{d}(mask{d},:);
        N_NT_new{d} = N_NT{d}(mask{d},mask{d});
    end
    %GG is incorporated to BB. The operator is applied on GG
    if(size(GG{1},2)>0)
        for i=1:Nops % specific treatment of the first separated dimension to incorporate the minus sign
            BB_tmp{1} = [BB_tmp{1} -AA{1,i}(mask{1},:)*GG{1}];
        end
        for d=2:Ndims
            for i=1:Nops
                BB_tmp{d} = [BB_tmp{d} AA{d,i}(mask{d},:)*GG{d}];
            end
        end
    end
    
    %Symetrization of the problem. The tranpose of the operators is applied
    %on BB and on the operators themselves
    for d = 1:Ndims
        BB_new{d} = zeros(nnz(mask{d}),0);
        for i=1:Nops
            BB_new{d} = [BB_new{d} AA{d,i}(mask{d},mask{d})'*BB_tmp{d}];
            for j=1:Nops
                AA_new{d,(i-1)*Nops+j} = AA{d,i}(mask{d},mask{d})'*AA{d,j}(mask{d},mask{d});
            end
        end
    end
    
    %PGD solution of the modified problem
    FF_new = PGD(AA_new,BB_new,GG_new,N_NT_new,Dirichlet_new,Max_terms,Max_fp_iter,epsilon,epsilon_tilde);
    
    %The solution is expanded to add the previously removed Dirichlet dofs.
    FF = cell(Ndims,1);
    Nterms = size(FF_new{1},2);
    for d=1:Ndims
        tmp = zeros(Nvars(d),Nterms);
        tmp(mask{d},:) = FF_new{d};
        FF{d} = [GG{d} tmp;];
    end
end