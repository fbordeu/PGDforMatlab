//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
// Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
//

#include "iostream"
#include <cstring>      // Needed for memset

#include <cmath>

#include <algorithm>

//#include <Eigen/Dense>
//#typedef Matrix<double, Dynamic, Dynamic> MatrixXd;

#include <signal.h>


#include "MatLab.h"
#include "NetWork.h"
#include "cRecompactCore.h"



int main( int argc, char**argv ){
    signal(SIGCHLD, SIG_IGN);
    std::cout << "*****************************" << std::endl;
    std::cout << "* cRecompact routine        *" << std::endl;
    std::cout << "* By Felipe Bordeu          *" << std::endl;
    std::cout << "* GeM Ecole Centrale Nantes *" << std::endl;
    std::cout << "* 2013                      *" << std::endl;
    std::cout << "*****************************" << std::endl;

    std::cout << "Connecting to port 5556" << std::endl;
    NetWork mynet("5556");

    pid_t ID; 
    while(true) {

        struct sockaddr_storage their_addr;
        std::cout << "Waiting for a client..." << std::endl;
        int new_sd = mynet.Accept(their_addr);
        ID = fork();
        
        if(ID == 0) {
              
            std::cout << "Waiting to recieve data..."  << std::endl;
            MatLabDataMatrix<double> flag;

            flag.recv(new_sd);
            while(flag(0,0) >0 ) {
              
                // we recive a option structure
                
                
                MatLabDataMatrix<double> ops;
                bool opt = false;
                if(flag(0,0) > 1) {
                  std::cout << "reciving options" << std::endl;
                  opt = true;          
                  ops.recv(new_sd);
                  for(int ii=0; ii< 9 ; ++ii){
                    std::cout << "op " << ii << " " << ops(0+ii,0)<< std::endl;
                  }
                }
              
                std::cout << "waintig for a job " << std::endl;
                MatLabDataMatrix<double> NumberOfDims;
                NumberOfDims.recv(new_sd);
                std::cout << "Number of Dimensions : " <<   NumberOfDims(0,0) << std::endl;

                std::vector<MatLabDataMatrix<double> >FF;
                FF.resize(NumberOfDims(0,0));
                for (unsigned i =0; i < NumberOfDims(0,0); ++i) {
                    std::cout << "reciving data " << i << "..";
                    FF[i].recv(new_sd);
                    std::cout << ".. Done" << std::endl;
                }
                std::cout << "Data in C " << std::endl;
                std::vector<MatLabDataMatrix<double> > sol;
                if(NumberOfDims(0,0)!=0){
                    PGD_Options options;

                    
                    if(opt){
                      options.reweight_modes       = ops(0+0,0);
                      options.fp_tol               = ops(0+1,0); 
                      options.res_reduc            = ops(0+2,0);
                      options.max_added_modes      = ops(0+3,0);
                      options.fp_max_iter          = ops(0+4,0);
                      options.improve_modes        = ops(0+5,0);
                      options.improve_modes_max    = ops(0+6,0);    
                      options.verbose              = ops(0+7,0);
                      options.lastImproveModesLoop = ops(0+8,0); 
                    }else {
                      options.fp_tol = 1e-8;
                      options.res_reduc = 1e-8;                   
                      options.improve_modes_max = std::max(10, (int) round(FF[0].dsizes[1]/10) );
                      options.max_added_modes = FF[0].dsizes[1];                      
                    }
                    //sleep(10);
                    Recompact(FF, sol,options);
                    }

                for (unsigned i =0; i < NumberOfDims(0,0); ++i) {
                    sol[i].send(new_sd);
                }

                // zero to quit or 1 to compress a new
                flag.recv(new_sd);
            }
            std::cout << "Clossing connection" << std::endl;
            close(new_sd);
            break;
        } else {
          std::cout << " Runing new process for client id " << ID << std::endl;
          close(new_sd);
        }
    }
};