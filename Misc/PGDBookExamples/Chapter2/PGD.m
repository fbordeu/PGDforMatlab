function [FF] = PGD(AA,BB,GG,N_NT,Dirichlet,Max_terms,Max_fp_iter,epsilon,epsilon_tilde)
    %General purpose PGD function, not optimized for efficiency
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
    %Number of terms in the separated representation of the RHS
    Nb = size(BB{1},2);
    %Number of a priori known terms in the solution
    N0 = size(GG{1},2);
    %Number of variables along each PGD dimension
    Nvars = zeros(1,Ndims);
    for i=1:Ndims
        Nvars(i) = size(AA{i,1},1);
    end
    
    %auxiliary function to compute the multidimensional integral of f*g where
    %both f and g have a separated representation with one term
    aux_int = @(f,g) prod(cellfun(@(a,b,c) a'*b*c,f,N_NT,g));
    
    
    %Solution initialization
    FF = cell(Ndims,1);
    %Add GG to FF
    for i=1:Ndims
        FF{i} = GG{i};
    end
    
    %Enrichment loop
    for term = (N0+1):Max_terms
        %RS INITIALIZATION + Application of the Dirichlet boundary conditions
        RS = cell(Ndims,1);
        for i=1:Ndims
            RS{i} = randn(Nvars(i),1);
            RS{i}(Dirichlet{i}) = 0;
        end
        
        %FIXED POINT ITERATIONS
        for iter=1:Max_fp_iter
            RS_old = RS;
            %UPDATE EACH ALONG EACH DIMENSION
            for i=1:Ndims
                %mask: DIMENSIONS OTHER THAN i
                mask=[ (1:i-1) (i+1:Ndims)];
                %CONSTRUCTION OF THE LHS
                K = sparse(Nvars(i),Nvars(i));
                %ACCUMULATION OF THE SIMPLE OPERATORS
                for j=1:Nops
                    %COMPUTATION OF THE PRE-FACTOR
                    factor_K = 1;
                    for tmp_dim=mask
                        factor_K = factor_K * (RS{tmp_dim}'*AA{tmp_dim,j}*RS{tmp_dim});
                    end
                    %ASSEMBLY
                    K = K + factor_K*AA{i,j};
                end
                %CONSTRUCTION OF THE RHS
                %INITIAL RHS (BB)
                factors_b = ones(1,Nb);
                for tmp_dim=mask
                    factors_b = factors_b .* (RS{tmp_dim}'*BB{tmp_dim});
                end
                b = BB{i}*factors_b.';
                
                %KNOWN MODES PASSED TO THE RHS (FF)
                if(term>1)
                    %EACH SIMPLE OPERATOR IS APPLIED TO FF
                    for j=1:Nops
                        factors_b = ones(1,term-1);
                        for tmp_dim=mask
                            factors_b = factors_b .* ( RS{tmp_dim}'*AA{tmp_dim,j}*FF{tmp_dim});
                        end
                        b = b -(AA{i,j}*FF{i})*factors_b.';
                    end
                end
                %APPLICATION OF THE DIRICHLET BOUNDARY CONDITIONS
                tmp = true(Nvars(i),1);
                tmp(Dirichlet{i}) = false;
                %SOLUTION
                RS{i}(tmp) = K(tmp,tmp)\b(tmp);
            end
            %fixed point termination test
            if (sqrt(aux_int(RS,RS)+aux_int(RS_old,RS_old)-2*aux_int(RS,RS_old))/ sqrt(aux_int(RS,RS)))<epsilon
                break;
            end
        end
        
        %Update of FF WITH RS
        for i=1:Ndims
            FF{i} = [FF{i} RS{i}];
        end
        %the first enrichment is stored for the enrichment loop termination
        %test
        if term==(N0+1)
            First_enrichment = cellfun(@(a) a(:,term),FF,'uniformoutput',false);
        end
        %enrichment loop termination test
        if (sqrt(aux_int(RS,RS))/ sqrt(aux_int(First_enrichment,First_enrichment)))<epsilon_tilde
            break;
        end
    end
end
