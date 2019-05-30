close all
clear all

%% création des indices: un nom et un range
i = Eindex('i',2);
j = Eindex('j',2);
k = Eindex('k',2);
l = Eindex('l',2);

%% création de tenseurs
%ordre 0
A0 = Etensor(3)
%%ordre 1 : données (vecteur colonne) et un indice, e.g. i.
A1 = Etensor(rand(2,1),i)
B1 = Etensor(rand(2,1),i)

%% changement d'indice
C1 = B1.idx(j);

%% ordre 2 : une matrice et deux indices
A2 = Etensor(rand(2,2),i,j)
B2 = Etensor(rand(2,2),i,j)

%% changement d'indice
C2 = B2.idx(k,l)
%% permutation des indices
D2 = B2.permute(j,i)

%% l'identité
I = Etensor(eye(2),i,j)

%% quelques opérations
%A0 * B1_i
contraction(A0,B1)

%% A1_i * B1_i
contraction(A1,B1) %ici A1 et B1 dépendent déjà de i

%% A1_i * B1_j
contraction(A1,B1.idx(j)) %on change d'indice à la volée

%% autre exemple
contraction(A2.idx(i,j),B2.idx(j,k))

%% autre façon
contraction(A2.idx(i,j),B2.idx(k,l),I.idx(j,k))

%% en passant par un tenseur d'ordre 4
tmp = contraction(A2.idx(i,j),B2.idx(k,l))
contraction(tmp,I.idx(j,k))

%% si on a des objets pour lesquels l'opérateur de différentiation est
%défini...

A0 = Etensor(sym('x*y+sin(x)'))
coords = Etensor([sym('x');sym('y')],i)
A1 = Etensor([sym('x*y+sin(x)');sym('x+y')],i)

%% gradient d'un scalaire
nabla(A0,coords)

%% gradient d'un vecteur
nabla(A1,coords.idx(j))

%% divergence
nabla(A1.idx(l),coords.idx(l))

%% plus complexe
tmp = contraction(A0,A1)
tmp = nabla(tmp.idx(i),coords.idx(j))
contraction(tmp,A1.idx(i))



