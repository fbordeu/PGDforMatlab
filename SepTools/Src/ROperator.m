classdef (InferiorClasses = {?RTensor}) ROperator 
    properties
        data
        is_sym = true;
    end
    
    methods
        % Constructors
        function obj = ROperator(arg1)
            obj.data = arg1;     
        end
        % Inspection functions ********************************

        % Unitary operators *************************************************
        function result = conj(arg1)
          result = cellfun(@conj, arg1.data,'UniformOutput',false);
          result = ROperator(result);
        end
        % Binary operators **************************************************
        function res = mtimes(arg1,arg2)
           ndim = size(arg1.data,1);
           if isa(arg1,'ROperator')
              % right multiplication
              ntM = size(arg1.data,2);
              ntV = size(arg2.data{1},2);
              res = cell(ndim,1);
              for d =1:ndim
                 ress = zeros(size(arg2.data{d},1),ntM*ntV);
                 for kk = 1:size(arg1.data,2)
                    ress(:,(1+(kk-1)*ntV):(ntV*kk)) = arg1.data{d,kk}*arg2.data{d};
                 end
                 res{d} = ress;
              end
              res = RTensor(res);
           else
               error('not yet sorry');
           end
        end

      function res = diag(arg1)
          ndims = size(arg1.data,1);
          res = ROperator(cell(ndims,1));
          for d =1: ndims;
              res.data{d,1} = full(cell2mat(cellfun(@diag, arg1.data(d,:),'UniformOutput',false)));
          end
      end
      function result = mldivide(arg1,arg2)
          
          ndims = size(arg1.data,1);
          nops = size(arg1.data,2);
          
          [problem, N_NT, options] = PGD_skeleton(ndims,nops);
          
          problem.AA = arg1.data;
          problem.BB = arg2.data;
          
          options.residual_minimization = ~arg1.is_sym;
          
          solution = PGD_v7b(problem,N_NT,options);
          
          result = ScaleTerms(RTensor(solution.FF),solution.alpha');
      end
      

    end

end
