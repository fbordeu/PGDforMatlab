classdef RTensor 
   properties 
     data;
   end
   
   methods  
      %% Constructors
      function obj = RTensor(arg1)
         if(exist('arg1','var'))
            if isa(arg1,'RTensor')
                obj.data = arg1.data;  
            else
                obj.data = arg1(:);     
            end
         else
            obj.data = {[]};
         end
      end
      %
      function result = ones(arg1)
          result = RTensor(cellfun(@(arg)(ones(size(arg,1),1)), arg1.data,'UniformOutput',false));
      end
      %
      function result = zeros(arg1)
          result = RTensor(cellfun(@(arg)(zeros(size(arg,1),1)), arg1.data,'UniformOutput',false));
      end
      %% Inspection functions ********************************
      function save(arg1)
          name = inputname(1);
          if numel(name) >= 0
              eval([ name ' = arg1;'])
              save(name,name);
          else
              disp('You must save a variable and not an expression')
          end
      end
      
      function disp(arg1)
        if numel(arg1) > 1
            for i = 1: numel(arg1)
                disp(arg1(i))
            end
            return
        end
        disp('  RTensor : Rank 1 Canonical Tensor ')
        %disp([ '  '  num2str(size(arg1,1))  ' Dimensions' ])
        if dsize(arg1,1) > 0
            disp([ '  with '  num2str(size(arg1.data{1},2))  ' terms' ])
        end
        flag = {'' 'NANs detected!!!'};
        for i = 1: dsize(arg1,1)
            disp(['    Size of dim ( ' num2str(i) ' ) : ' num2str(size(arg1.data{i},1)) ' ' flag{any(isnan(arg1.data{i}(:)))+1} ])
        end
        
        disp(' ')
        
      end
      %
      function plot(arg1,arg2,varargin)
          plot3D = false;
          if numel(varargin) >0 
            if isa(varargin{1},'RTensor')
                plot3D = true;
                arg3 = varargin{1};
                if numel(varargin) >1 
                    npoints = varargin{2};        
                else
                    npoints = 1000;      
                end
            else 
              npoints = varargin{1};  
            end
          else 
              npoints = 1000;
          end
                  
          sdims = cellfun(@(arg)(size(arg,1)), arg1.data )';
          ndims = numel(sdims);
          index = ceil(bsxfun(@times, rand(npoints,ndims), sdims));
          index = unique(index,'rows');
          a1 = arg1.eval(index);
          a2 = arg2.eval(index);
          
          if plot3D 
              TRI = delaunay(a1,a2);
              a3 = arg3.eval(index);
              %trisurf(TRI,a1,a2,a3)
              plot3(a1,a2,a3,'.');
              xlabel(inputname(1))
              ylabel(inputname(2))
              zlabel(inputname(3))
          else
              plot(a1,a2,'.');
              xlabel(inputname(1))
              ylabel(inputname(2))
          end
          
          
      end
      %
      function res = NumberOfTerms(arg1)
          
            res = size(arg1.data{1},2);

      end
      %
      function result = dsize(arg1,varargin)
          if numel(arg1.data)
            result = size(arg1.data,varargin{:}); 
          else
              result =0;
         end
      end
      %
      function res = GetAlphas(arg1)
        res = getalphas(arg1.data);
      end
      %
      function result = norm(arg1)
          result = sepnorm(arg1.data);
      end
      %
      function result = export(arg1)
          result = arg1.data;
      end
      
      function res= extractTerm(arg1,n)
          res = RTensor(cellfun(@(arg) (arg(:,n)), arg1.data, 'UniformOutput',false));
      end
      function res = eval(arg1,where)
          res = sepeval(arg1.data,where);
      end
      function res = PartFor(arg1,varargin)
          res = arg1;
          for i= 1:numel(varargin)
              index = varargin{i};
              if isnumeric(index)
                res.data{i} = arg1.data{i}(index,:);
              end
          end
      end
      %% no modification *******************************************************
      function res= SortTermsByAlphas(arg1)
          [~,I] = sort(GetAlphas(arg1));
          I = fliplr(I);
          res = RTensor(cellfun(@(arg) (arg(:,I)), arg1.data, 'UniformOutput',false));
      end
      %% modifcation ******************************************************
      function res = ScaleTerms(arg1, alphas)
        res = RTensor(setalphas(arg1.data,alphas));
      end
      %
      function res = Prune(arg1, tol)
          a = getalphas(arg1.data);
          m = max(a);
          I = (a/m) > tol;
          res = RTensor(cellfun(@(arg) (arg(:,I)), arg1.data, 'UniformOutput',false)   );
      end
      %% Unitary operators *************************************************
      function result = uminus(arg)
          result = arg.mtimes(-1);
      end
      %
      function result = recompact(arg1,varargin) 
          result = recompact(arg1.data,varargin{:});
          result = RTensor(result);
      end
      %
      function result = sign(arg1,varargin) % experimental
          %from "Tensor Spaces and Numerical Tensor Calculus", 
          % Wolfgang Hackbusch 
          
          % we do an extra normalization for the algorithm
          scaling = 0.5*sqrt(sepmax(export(arg1*arg1)));
          result = arg1/scaling;
          difold =0;
          for i = 1:100 ; 
              result1 = 0.5*(result)*(3-result*result); 
              %opt = PopRecompactOptions(false);
              %opt.recompact  = false;
              %PushRecompactOptions(opt);
              %  result1 = result.polyval([-0.5 0 1.5 0]);
              %PopRecompactOptions();
              %result1 = recompact(result1);
              
              % Recover the initial dif
              if  i == 1
                  n0 =  norm(result1-result);
              end
              % we do 10 iteration mininum 
              if  i > 9 
                dif = norm(result1-result); 
                % check of dif is contracting or if the diff is very small
                if  i > 10  && (dif > difold || dif/n0 < 1e-5 )
                  break;
                end
                difold = dif ;
              end
              result = result1;
              %plot(arg1,result,100000)
              %pause;
          end
      end
      % Binary operators **************************************************
      function result = plus(a,b)
            
            %overloads the sum function
            if isnumeric(b)
                result = sepsum(a.data,b);
            else
                if isnumeric(a)
                    result = sepsum(a,b.data);
                else
                    result = sepsum(a.data,b.data);
                end
            end
            result = RTensor(result);
      end
      %
      function result = times(arg1,arg2)
          result =  mtimes(arg1,arg2);
      end
      %
      function result = mtimes(arg1,arg2)
            %overloads the sum function
            if isnumeric(arg2)
                result = sepprod(arg1.data,arg2);
            else
                if isnumeric(arg1)
                   result = sepprod(arg1,arg2.data);
                else
                   result = sepprod(arg1.data,arg2.data);
                end
            end
            result = RTensor(result);
      end
      %
      function result = mrdivide(arg1,arg2)
          if isnumeric(arg2)
              result = sepfrac(arg1.data,arg2);
          else
              if isnumeric(arg1)
                result = sepfrac(arg1,arg2.data);
              else
                result = sepfrac(arg1.data,arg2.data);
              end
          end
          
          result = RTensor(result);
      end
      function result = mldivide(arg1,arg2)
          if isnumeric(arg1)
              result = arg2/arg1;
          else
              error('RTensor:mldivide', 'Cant do mldivide ')
          end
      end
      function result = minus(a,b)
          result = plus(a,-b);
      end
      function result = dot(arg1,arg2)
          result = sepdot(arg1.data,arg2.data);
      end
      % complex functions *************************************************
      function result = power(arg,opt)
          result = mpower(arg,opt);
      end
          
      function result = mpower(arg1,arg2)
          assert(isscalar(arg2),'RTensor cant do power using this exponent');
          
          assert(arg2 > 0 ,'RTensor cant do power with negative number');
          if round(arg2) == arg2 
            poly = zeros(1,arg2+1);
            poly(1) = 1;
            result = seppolyval(poly, arg1.data);
          else
            result = sepfracpow(arg1.data,arg2 );
          end
          
          result = RTensor(result);
      end
      %
      function result = polyval(arg1,arg2)
          result = RTensor(seppolyval(arg2,arg1.data));
      end
      %
      function result = exp(arg1)
          result = sepexppade(arg1.data);
          result = RTensor(result);
      end
      %
      function result = expForNegativeNumbers(arg1)
          result = sepmexppade(arg1.data);
          result = RTensor(result);
          %result = ones(arg1);
          %for i = 1: arg1.NumberOfTerms()
          %    ops = PopRecompactOptions(false); % we recover the options
          %    ops.recompact =  false; 
          %    PushRecompactOptions( ops );
          %        pade1term = RTensor(sepmexppade(arg1.extractTerm(i)));
          %    PopRecompactOptions();
          %    result = result*pade1term;
          %end
              
          
      end
      %
      function result = padeval(arg1,pade)
        result = seppadeval(pade, arg1.data);
        result = RTensor(result);
      end 
   end
end