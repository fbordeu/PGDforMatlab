function varargout = validateAA(AA,nonSquare)
%[Ndims,Ndofs_in,Ndofs_out,Nterms] = validateAA(AA,true)
%[Ndims,Ndofs,Nterms] = validateAA(AA,false)
%[Ndims,Ndofs,Nterms] = validateAA(AA)
%Checks if AA is a valid separated representation of an operator
%if the nonSquare optional argument is set to true, 
%allows for non-square matrices in AA.
%
%Return values:
% Ndims: number of dimensions
% Ndofs: number of degrees of freedom per dimension
% Ndofs_in: number of degrees of freedom per dimension for the source tensor
%(i.e. number of columns of the elements of AA)
% Ndofs_out: number of degrees of freedom per dimension for the image tensor
%(i.e. number of lines of the elements of AA)
% Nterms: number of terms in the separated representation
%
%Author: Adrien Leygue.

%do we have a non-empty cell array?
if nargin==1
    nonSquare = true;
end
validateattributes(nonSquare,{'logical'},{'scalar'},2);

validateattributes(AA,{'cell'},{'2d','nonempty'},'','AA');

Ndims = size(AA,1);
Nterms = size(AA,2);
Ndofs_in = zeros(Ndims,1);
Ndofs_out = zeros(Ndims,1);

%is the content of each cell a nonempty 2d array of floating point numbers?
for i=1:Ndims
    %accumulate the size of the matrices to check if they are the same over
    %each dimension
    tmp_in  = zeros(1,Nterms);
    tmp_out = zeros(1,Nterms);
    for j=1:Nterms
        validateattributes(AA{i,j},{'double','single'},{'nonempty','2d'},'',['AA{' num2str(i) ',' num2str(j) '}']);
        [tmp_out(j), tmp_in(j)] = size(AA{i,j});
    end
    err_msg = ['Inconsistent operator sizes for dimension ' num2str(i) ];
    %check number of columns
    assert(all(tmp_in==tmp_in(1)),err_msg);
    %check number of lines
    assert(all(tmp_out==tmp_out(1)),err_msg);
    Ndofs_in(i)  = tmp_in(1);
    Ndofs_out(i) = tmp_out(1);    
end
%check if all matrices are square if nonSquare is false.
%then build the return values.
if ~nonSquare
    assert(all(Ndofs_in==Ndofs_out),'Non square operators found in AA');
    Ndofs = Ndofs_in;
    varargout = {Ndims,Ndofs,Nterms};
else
    varargout = {Ndims,Ndofs_in,Ndofs_out,Nterms};
end
end