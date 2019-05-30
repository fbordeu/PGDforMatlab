#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
# Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
#

type := $(shell uname)

all : cRecompactServer_$(type) cRecompactServer_debug_$(type) cRecompactFile_$(type) cRecompactMPI_$(type) cQuotientFile_$(type) matlab

matlab : recompactmex quotientmex

cRecompactFile_$(type): cRecompactFile.cpp MatLab.h  cRecompactCore.*  PGD_Options.* ExtraFunctions.h
	g++ -o cRecompactFile_$(type) -D${type} cRecompactFile.cpp cRecompactCore.cpp  CommandLineParser.cpp PGD_Options.cpp -lblas -O3 -march=native -std=c++0x -std=gnu++0x

cRecompactServer_$(type): cRecompactServer.cpp MatLab.h NetWork.h cRecompactCore.*  PGD_Options.h ExtraFunctions.h
	g++ -o cRecompactServer_$(type) cRecompactServer.cpp cRecompactCore.cpp PGD_Options.cpp -lblas  -O3  -std=c++0x -std=gnu++0x

cRecompactServer_debug_$(type): cRecompactServer.cpp MatLab.h NetWork.h cRecompactCore.*  PGD_Options.h ExtraFunctions.h
	g++ -o cRecompactServer_debug_$(type) cRecompactServer.cpp cRecompactCore.cpp PGD_Options.cpp -lblas -ggdb -g -std=c++0x -std=gnu++0x


cRecompactMPI_$(type): cRecompactMPI.* cRecompactCore.*  PGD_Options.h
	mpic++ -O3 -IMPI_includes -DUSE_MPI -o cRecompactMPI_$(type) cRecompactMPI.cpp PGD_Options.cpp CommandLineParser.cpp cRecompactCore.cpp -lblas -std=c++0x -std=gnu++0x

cQuotientFile_$(type): cQuotientFile.* cQuotientCore.* ExtraFunctions.h
	g++ -o cQuotientFile_$(type) cQuotientFile.cpp cQuotientCore.cpp PGD_Options.cpp CommandLineParser.cpp -lblas -ggdb -g -std=c++0x -std=gnu++0x

recompactmex:
	/opt/matlab2013a/bin/mex -lmwblas -lut recompactmex.cpp cRecompactCore.cpp  PGD_Options.cpp

quotientmex:
	/opt/matlab2013a/bin/mex -lmwblas -lut quotientmex.cpp cQuotientCore.cpp PGD_Options.cpp


clean:
	-rm cRecompactFile_$(type)
	-rm cRecompactServer_$(type)
	-rm cRecompactServer_debug_$(type)
	-rm cRecompactMPI_$(type)
