//
// C++ Implementation:
//
// Description:
//
//
//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
// Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
//
#ifndef GEM_GEM_CommandLineParser_cpp
#define GEM_GEM_CommandLineParser_cpp

#include "CommandLineParser.h"
#include <stdlib.h>
#include <iomanip>

namespace xpgd {
//
__attribute__((noinline)) std::string CommandLineParser::get_lib_compilation_date(){ return std::string("Library compiled ")+ __DATE__ +" at " +__TIME__ ;};

void CommandLineParser::init(int _argc, char **_argv) {
    argc = _argc;
    argv = _argv;
};
//
void CommandLineParser::addHelp(const std::string& com, const std::string help){
  man.push_back(manint(com,help));
}
void CommandLineParser::doneReading(){
  
    addHelp("--CompileDate","To get the compilation date.");
    addHelp("--cerr file","Send the cerr output to the file");
    addHelp("--cout file","Send the cout output to the file");
    
    if( get_option("--CompileDate") ){
      /// please add #def MAIN before the #include "commandLineParser.h" in your main
      std::cout << get_exec_compilation_date() << std::endl;
      std::cout << get_lib_compilation_date() << std::endl;
      exit(0);
    }
    if (has_atribute("--help")  || has_atribute("-help")  || has_atribute("-h")) {
      std::size_t  si = 0;
      for(unsigned i =0; i < man.size(); ++i){
        si = std::max(man[i].command.size(),si);
      }
      for(unsigned i =0; i < man.size(); ++i){
        std::cout <<std::left << std::setw(si) << man[i].command << " : "  ;
        std::size_t  found = man[i].help.find_first_of("\n");
        std::size_t  oldfound =  0;
        while (found!=std::string::npos){
          std::cout << man[i].help.substr(oldfound,found-oldfound) << std::endl;
          oldfound = found+1;
          found= man[i].help.find_first_of("\n",found+1);
        }
        if(oldfound >1){
          std::cout << std::left << std::setw(si+3) <<" "  ;
        }
        std::cout  << man[i].help.substr(oldfound,found-oldfound) << std::endl;
      }
      exit(0);
    } 
      
    std::string file;
    if(get_atribute("--cerr", file) ){
      cerrToFile(file);
    }
    if(get_atribute("--cout", file) ){
      coutToFile(file);
    }
};
//
CommandLineParser::CommandLineParser():  cerr_strbuf(0), debug (0) , cout_strbuf(0) ,output (0)  {};
//
CommandLineParser::CommandLineParser(int _argc, char **_argv): cerr_strbuf(0), debug (0) , cout_strbuf(0) ,output (0) {
    init(_argc,_argv);
};
//
bool CommandLineParser::get_file(std::string &file, unsigned int f_number) {
    unsigned cpt=0;
    for (int i=1; i<argc; i++) {
        if (!strncmp("+",&(argv[i][0]),1)) continue;
        if (!strncmp("-",&(argv[i][0]),1)) {
            i++;
            continue;
        };
        if (cpt == f_number) {
            file = argv[i];
            return 1;
        } else {
            cpt++;
            continue;
        }
    }
    return 0;
}
//
bool CommandLineParser::get_option(const char *key){
  
    for (int i=0; i<argc; i++) {
        if (!strcmp(argv[i], key)) {
            return 1;
        }
    }
    return 0;
};
//
void CommandLineParser::get_option(const char *key, bool& value) {
  value = get_option(key);
}
//
CommandLineParser::~CommandLineParser(){
  if(this->debug){
    std::cerr.rdbuf(cerr_strbuf);
    delete this->debug;
    this->debug = 0;
  }
  if(this->output){
    std::cout.rdbuf(cout_strbuf);
    delete this->output;
    this->output = 0;
  }
  man.resize(0);
}    
//
bool CommandLineParser::cerrToFile(const std::string filename){
      debug = new std::ofstream(filename.c_str());  
      if((*debug).is_open()) {
       // Overriding std::cerr to a file.
        cerr_strbuf = std::cerr.rdbuf();
        std::cerr.rdbuf(debug->rdbuf());
	std::cout << " **** sending all the cerr outputs to file :" << filename << " **** " << std::endl;
	return 0;
      } else {
	std::cerr << " **** Error opening file " << filename << ", cerr output not modified **** " << std::endl;
        std::cerr << " **** cerr output on the screen **** " << std::endl;
	return 1;
      };
}
//
bool CommandLineParser::coutToFile(const std::string filename){
      output = new std::ofstream(filename.c_str());  
      if((*output).is_open()) {
       // Overriding std::cerr to a file.
        std::cout << " **** sending all the cout outputs to file :" << filename << " **** " << std::endl;
        cout_strbuf = std::cout.rdbuf();
        std::cout.rdbuf(output->rdbuf());
        std::cout << " **** sending all the cout outputs to file :" << filename << " **** " << std::endl;
	return 0;
      } else {
	std::cerr << " **** WARNING: Error opening file " << filename << ", cout output not modified **** " << std::endl;
        std::cerr << " **** cout output on the screen **** " << std::endl;
	return 1;
      };
}
//
}
//
#endif //GEM_CommandLineParser_cpp
