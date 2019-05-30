%This script calls the functions Newton and Picard implementing two different 
%PGD-based non-linear solvers for the problem 6.1 in order to generate
%Figure 6.1
%
%Copyright (c) 2013, Francisco Chinesta (Ecole Centrale de Nantes), Roland Keunings (Universite catholique de Louvain), Adrien Leygue(CNRS)
%Author: Adrien Leygue.
%All rights reserved.
%See License file for more details


clear all
close all
%Mesh definition
Nx = 151;
Nt = 361;
k = 0.1;
x = linspace(0,1,Nx)';
t = linspace(0,1,Nt)';
%grid size
hx = x(2)-x(1);
dt = t(2)-t(1);

%FE Mass matrix for x
%elemental mass matrix
Mex = reshape([hx.*(1.0./3.0),hx.*(1.0./6.0),hx.*(1.0./6.0),hx.*(1.0./3.0)],[2,2]);
Mx = sparse(Nx,Nx);
%FE assembly of the FE mass matrix
for el = 1:Nx-1
    Mx([el el+1],[el el+1]) = Mx([el el+1],[el el+1]) + Mex;
end
%"Mass matrix" along t: trapezoidal rule weights on the diagonal
wt = dt*[0.5; ones(Nt-2,1);0.5];
Mt = spdiags(wt,0,Nt,Nt);

%reference solution: Matlab built-in ODE solver
%Computation of the FD differentiation matrix
e = ones(Nx,1);
D2 = k*spdiags([e -2*e e]/(hx^2),[-1 0 1],Nx,Nx);
D2([1 Nx],:) = 0;
%definition of the ODE
udot = @(in_t,u) D2*u - u.^2+[0; ones(Nx-2,1);0];
%solution
options = odeset('reltol',1e-8,'abstol',1e-12);
[t, U_ref] = ode113(udot,t,zeros(Nx,1),options);
U_ref = U_ref';



%Parameters explored: Three cases
%Number of PGD enrichments for each iteration
Max_terms_VALUES= [1  2   4];
%For each case:Number of iterations of the non-linear solver
NL_iter_VALUES  = [15 10  5];

handles = [];
legends = {};
colors = {'b','m','r'};
figure;
hold on;
set(gca,'yscale','log','fontsize',14);
for n=1:numel(Max_terms_VALUES)
%Non-linear solver parameter
%Number of non-linear iterations
NL_iter = NL_iter_VALUES(n);
%PGD parameters
%Maximum number of PGD enrichments at each non linear iteration
Max_terms = Max_terms_VALUES(n);
%Enrichment continuation tolerance=0 to enforce Max_terms enrichments
epsilon_tilde = 0;
%Fixed point parameters
Max_fp_iter = 50;
epsilon = 1e-8;

%Solution with Picard iterations
%all_FF_Picard is a cell array where each entry contains the PGD solution of a non-linear it  
all_FF_Picard = Picard(x,t,NL_iter,Max_terms,Max_fp_iter,epsilon,epsilon_tilde,k);
%Solution with Newton iterations
%all_FF_Picard is a cell array where each entry contains the PGD solution of a non-linear iteration  
all_FF_Newton = Newton(x,t,NL_iter,Max_terms,Max_fp_iter,epsilon,epsilon_tilde,k);


%Computation of the error as a function of the number of
%containers for the solutions
FF_Picard = cell(2,1);
FF_Newton = cell(2,1);
%Container for the error levels
E_Picard = zeros(1,NL_iter);
E_Newton = zeros(1,NL_iter);

for i=1:NL_iter
    %concatenation of the i^h increment with the current solution
    FF_Picard = cellfun(@horzcat,FF_Picard,all_FF_Picard{i},'uniformoutput',false);
    FF_Newton = cellfun(@horzcat,FF_Newton,all_FF_Newton{i},'uniformoutput',false);    
    %reconstruction of the solution
    U_Picard = FF_Picard{1}*FF_Picard{2}';
    U_Newton = FF_Newton{1}*FF_Newton{2}';
    %computation of the errors using the trapezoidal rule
    E_Picard(i) = trapz(t,trapz(x,(U_ref-U_Picard).^2,1),2);
    E_Newton(i) = trapz(t,trapz(x,(U_ref-U_Newton).^2,1),2);
end
%plot for Figure 6.1
h1 =line('xdata',1:NL_iter,'ydata',E_Picard,'linewidth',1.5,'color',colors{n},'linestyle','-.');
h2 = line('xdata',1:NL_iter,'ydata',E_Newton,'linewidth',1.5,'color',colors{n},'linestyle','-');
handles = [handles h1 h2];
legends = [legends {['$u_{k,' num2str(Max_terms_VALUES(n)) '}^{\mathrm{P}}$'],['$u_{k,' num2str(Max_terms_VALUES(n)) '}^{\mathrm{N}}$']}];
end
[legend_h,object_h,plot_h,text_strings]=legend(handles,legends,'interpreter','latex');
xlabel('Non-linear iterations (k)');
ylabel('Error level');
