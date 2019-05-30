#!/bin/bash

PGDPublic_PATH='PGDPublic'

/opt/matlab2013a/bin/matlab -nodesktop -nosplash -r "restoredefaultpath();addpath(pwd);setMatlabTools('nosave');copytopublic('$PGDPublic_PATH',12345);exit();"




