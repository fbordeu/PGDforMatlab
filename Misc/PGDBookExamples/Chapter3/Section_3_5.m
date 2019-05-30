%This script calls the function PGD which implements the PGD
%solution of a general problem in order to solve 2D image compression
%example of section 3.5 and generate the corresponding figure.
%
%Copyright (c) 2013, Francisco Chinesta (Ecole Centrale de Nantes), Roland Keunings (Universite catholique de Louvain), Adrien Leygue(CNRS)
%Author: Adrien Leygue.
%All rights reserved.
%See License file for more details


close all
clear all
IMG = imread('IMAGE.png');
IMG = rgb2gray(IMG);

B = double(IMG);

[Ny, Nx] = size(B);

%Container for the differential operator
AA = cell(2,1);
%Container for the image to compress
BB = cell(2,1);
%Container for the a priori known terms of the PGD solution
GG = cell(2,1);
%Container to store the nodes subject to a Dirichlet condition
Dirichlet = cell(2,1);
%Container to store the mass matrix of each dimension
N_NT = cell(2,1);


%Here there is no differential operator
AA{1,1} = speye(Nx);
AA{2,1} = speye(Ny);

%Trivial separated representation of B
BB{1} = eye(Nx,Nx);
BB{2} = B;

N_NT{1} = speye(Nx);
N_NT{2} = speye(Ny);

%Fixed point loop exit criterion
epsilon = 1e-8;
%PGD enrichment tolerance
epsilon_tilde = 1e-8;
%Max number of fixed point iterations
N_fp = 20;
%Maximum dumber of enrichments
Max_terms = 50;

%Compute the separated approximation
[FF] = PGD(AA,BB,GG,N_NT,Dirichlet,Max_terms,N_fp,epsilon,epsilon_tilde);

%Figure 3.1
figure;
subplot(2,2,1),imshow(IMG);
subplot(2,2,1),caxis([min(IMG(:)) max(IMG(:))]);
set(gca,'fontsize',14),title('Original Image','interpreter','latex','fontsize',14)

N = 1;
APPROX = FF{2}(:,1:N)*FF{1}(:,1:N)';
%convert to integer grayscale values
APPROX = uint16(APPROX);
subplot(2,2,2),imshow(APPROX);
subplot(2,2,2),caxis([min(APPROX(:)) max(APPROX(:))]);
set(gca,'fontsize',14),title('N=1','interpreter','latex','fontsize',14)

N = 15;
APPROX = FF{2}(:,1:N)*FF{1}(:,1:N)';
%convert to integer grayscale values
APPROX = uint16(APPROX);
subplot(2,2,3),imshow(APPROX);
subplot(2,2,3),caxis([min(APPROX(:)) max(APPROX(:))]);
set(gca,'fontsize',14),title('N=15','interpreter','latex','fontsize',14)

N = 50;
APPROX = FF{2}(:,1:N)*FF{1}(:,1:N)';
%convert to integer grayscale values
APPROX = uint16(APPROX);
subplot(2,2,4),imshow(APPROX);
subplot(2,2,4),caxis([min(APPROX(:)) max(APPROX(:))]);
set(gca,'fontsize',14),title('N=50','interpreter','latex','fontsize',14)
