%% Examples of the SYM package (tested with Matlab R2012b only)
%% Basics
% The SYM package allows you to manipulate simple symbolic expressions. It
% can currently handle the following entities:
%%
% * Zero : SYM_ZERO
% * One : SYM_ONE
% * Real numbers (double precision) : SYM_REAL
% * Named constants : SYM_CONSTANT
% * Coordinates : SYM_COORD
% * Functions / Fields of one or several variables : SYM_FIELD
% * Products of those objects : SYM_PROD
% * Sums of those objects : SYM_SUM
% * Arrays of those objects.
%
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Adrien Leygue (Adrien.Leygue@ec-nantes.fr)
%

% A first example
clear all
E = SYM_REAL(2.718);
C = SYM_CONSTANT('c');
RESULT = E*C + SYM_ONE - SYM_REAL(4)
%%
% SYM_ZERO has the appropriate properties in additions and multiplications.
RESULT = SYM_ZERO*C + SYM_CONSTANT('nu')
%%
% SYM_ONE has the appropriate properties in multiplications.
RESULT = C*SYM_ONE
%%
% Products of sums are expanded
RESULT = (C+SYM_ONE)*(SYM_CONSTANT('nu')+E)
%%
% Doubles in symbolic expressions are automatically converted to SYM_REAL.
RESULT = 3*(C+1)
%%
% In a sum, terms are automatically collected.
RESULT = C + 3*C
%%
% Operands are reordered automatically. It's not a bug, it's a feature!
RESULT = SYM_CONSTANT('nu')+C
%%
% Arrays work as expected in Matlab
A = SYM_ONE(3)
B = C*A
RESULT = B + SYM_REAL(rand(3))
V = RESULT*SYM_ONE(3,1)
W = V .* [2;C;SYM_ZERO]

%% Manipulating functions
% As this code is oriented towards Finite Element and weak formulation, on
% can define three types of functions:
%%
% * known functions
% * unknown functions (solution of some PDE)
% * test functions, identified by the sign ° following the function name
% When the function is created, one has to explicitely state the coordinate
% set.
clear all
x = SYM_COORD('x');
y = SYM_COORD('y');
t = SYM_COORD('t');

F_known = SYM_FIELD('F',[x y])
U_unknown = SYM_FIELD('U',[t x y],'unknown')
V_unknown = SYM_FIELD('V',[t x y],'unknown')
V_test = SYM_FIELD('V',[t x y],'test')

%%
% On can compute derivatives, gradient, divergence of the functions or of
% more complex expressions
diff(F_known,x)
dU = gradient(U_unknown,[x y])
Laplacian = divergence(dU,[x y])
RESULT = gradient((3*U_unknown+V_unknown)*F_known,[x y])
RESULT(1)
%%
% Finally, it can be convenient to define a test function as the variation of
% an unknown function
U_unknown
U_unknown.test_fct
test_fct(U_unknown)
test_fct(F_known)
test_fct(F_known*U_unknown)
RESULT = test_fct(U_unknown*diff(V_unknown,x))
%%
% Expressions can be converted to LaTex
RESULT.display_tex

%% Advanced topics
% The subs function allow one to make advanced substitutions of objects or
% of entire sub-expression. In that last case the sub-expression has to appear explicitely.
clear all
C = SYM_CONSTANT('c');
N = SYM_CONSTANT('nu');
a = 3*C + 4*C*C
RESULT_1 = subs(a,C,SYM_REAL(5))
RESULT_2 = subs(a,C*3,N)
%%
% Finally given an expression representing the weak form of a linear PDE,
% it is possible to identify all the terms that have to be assembled respectively in the
% left hand side and right hand side of the problem.
clear all
x = SYM_COORD('x');
y = SYM_COORD('y');

U = SYM_FIELD('U',[x y],'unknown');
F = SYM_FIELD('F',[x y]);
dU = gradient(U,[x y]);
K = [SYM_CONSTANT('Kxx') SYM_CONSTANT('Kxy') ; SYM_CONSTANT('Kxy') SYM_CONSTANT('Kyy')]

WEAK = dU.test_fct'*K*dU + F*U.test_fct
[LHS,RHS] = extract_LHS_RHS(WEAK);
disp('LHS terms')
LHS(:)
disp('RHS terms')
RHS(:)
%%
% When the problem will be solved using the PGD, it is possible to
% specify a partition of the coordinates.
clear all
x = SYM_COORD('x');
y = SYM_COORD('y');
p = SYM_COORD('p');

U = SYM_FIELD('U',[p x y],'unknown');
F = SYM_FIELD('F',[p x y]);
dU = gradient(U,[x y]);
K = [SYM_CONSTANT('Kxx') SYM_CONSTANT('Kxy') ; SYM_CONSTANT('Kxy') SYM_CONSTANT('Kyy')]

WEAK = dU.test_fct'*K*dU + F*U.test_fct
disp('Variable separation: (x,p) and y')
[LHS,RHS] = WEAK.extract_separated_LHS_RHS({[p x],y})






