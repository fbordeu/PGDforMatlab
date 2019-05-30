--------------------------------------------------------------------
1) this is the header to put in every file.

for matlab code 

%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : first  last (mail)
%

of for c++

//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
// Principal developer : first  last (mail)
//

for python code 

#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
# Principal developer : first  last (mail)
#
--------------------------------------------------------------------

2) Time bomb, we use clock

%the time bomb
if all([ (clock >= [2015 06 0 0 0 0])  (builtin('clock') >= [2015 06 0 0 0 0])] )
   error('License Expired!!!');
end

functions with time bomb:

MatlabTools/PGDSolvers/Src/V3/SPTF_FV_v3.m
MatlabTools/PGDSolvers/Src/V3/SPTF_Min_v3.m
MatlabTools/PGDSolvers/Src/V4/SPTF_FV_v4.m
MatlabTools/PGDSolvers/Src/V4/SPTF_Min_v4.m
MatlabTools/PGDSolvers/Src/V5/SPTF_Min_v5.m
MatlabTools/PGDSolvers/Src/V6/PGD_v6.m
MatlabTools/PGDSolvers/Src/V7b/PGD_v7b.m

MatlabTools/Kompressor/Src/recompact.m
MatlabTools/Kompressor/Src/quotient.m

MatlabTools/SymTools/Src/SYM_COORD.m

