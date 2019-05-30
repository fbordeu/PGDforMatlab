1)Contributors (in alphabetical order):

    Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
    Domenico Borzacchiello (Domenico.Borzacchiello@ec-nantes.fr)
    Adrien Leygue (Adrien.Leygue@ec-nantes.fr)

2) installation

****** To checkout ****** 

* for the p code version:

> svn checkout https://svn.ec-nantes.fr/ecn/PGD/trunk/MatlabTools MatlabTools

* for the developer version:

> svn checkout https://svn.ec-nantes.fr/ecn/PGDEV/trunk/MatlabTools MatlabTools

****** To update the library ****** 
(inside the MatlabTools folder)
> svn update


3) Configuration

in matlab :

>> addpath('/home/users/mater/fbordeu/PGD/MatlabTools')
>> setMatlabTools()

if you don't want to modifie your path permanently 

>> setMatlabTools('nosave') 

4) for verification and testing the integrity of the library

>> TestAll()

4) for generation of the p code version of the library 



