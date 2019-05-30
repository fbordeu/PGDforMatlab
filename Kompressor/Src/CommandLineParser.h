//
// C++ Interface: CommandLineParser
//
// Description: Class to read the option of the command line
//
//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
// Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
//

#ifndef COMMANDLINEPARSER_H
#define COMMANDLINEPARSER_H

#include <cstring>
#include <fstream>
#include <signal.h>

#include "xPGDInputOutput.h"

namespace xpgd {
/* \defgroup xPGD_IO  "xPGD io" */
/**
\brief class to extract the information from the command line in a very easy way. and redirect the cerr and/or cout to a file
\ingroup xpgdio
\verbatim
>$app   (+letter1)   (-letter2 variable2)  (file1) (file2)  --cerr log.cerr
\endverbatim
the option --CompileDate to print application and library  compilation date and time.

example:

\verbatim
>$./myappr --CompileDate
Exec compiled May  9 2012 at 10:17:33
Library compiled May  9 2012 at 10:17:05
\endverbatim

example 
\verbatim
$app -l 4 +R toto.xml index.html
\endverbatim
this is equivalent to :
\verbatim
$app toto.xml +R index.html -l 4
\endverbatim

the code in c++

\verbatim
int main(int argc,char **argv){

    CommandLineParser maindata(argc,argv);   
    
    int val_l;
    maindata.get_atribute("-l", val_l);
    std::cout << "-l = " << val_l << std::endl;

    bool opt_R=0;
    maindata.get_option("+R",opt_R);
    std::cout << "+R = " << opt_R << std::endl;
    
    std::string file;
    if(!maindata.get_file(file))
	std::cout << "first file : "<<file << std::endl;
    else 
	std::cout << "file not present" << std::endl;

    
    if(!maindata.get_file(file,1))
	std::cout << "second file : "<<file << std::endl;
    else 
	std::cout << "second file not present" << std::endl;
  
    
   return 0;
};
\endverbatim
*/

typedef void (*handler)(void);
void* signal(int signum, handler);


//void sigusr1(int sig){
//  std::cerr << " *** Calculation Probleme  *** (all files are correctly closed)" << std::endl;
//  std::cout << " *** Calculation Probleme  *** (all files are correctly closed)" << std::endl;
//  exit(0);
// }
//
//    signal(SIGTERM, term);      // register a SIGTERM handler // raise(SIGTERM); // will cause term() to to exit normaly
//    signal(SIGINT, term);   
//
//     struct sigaction {
//		  void (*sa_handler)(int);
//		  void (*sa_sigaction)(int, siginfo_t *, void *);
//		  sigset_t sa_mask;
//		  int sa_flags;
//		  void (*sa_restorer)(void);
//	      }
//    sigaction myaction;
//    myaction.sa_andler = sigusr1;
//    myaction.sa_flags = 0;
//    sigaction(SIGUSR1,);
// 
//
class CommandLineParser {
    /// please add #def MAIN before the #include "commandLineParser.h" in you main file
    __attribute__((noinline)) std::string get_exec_compilation_date();
    __attribute__((noinline)) std::string get_lib_compilation_date();
public:
    void init(int _argc, char **_argv);
    //
    CommandLineParser();
    //
    CommandLineParser(int _argc, char **_argv);
    //
    template<class TV>
    bool get_atribute(const char *key, TV& value, unsigned num = 0) {
        for (int i=0; i<argc; i++) {
            if (!std::strcmp(argv[i], key)) {
	        if(num){
		  --num;
		} else {
		  if(argc < (i+1))
		    return 0;
		  input(argv[i+1],value);
		  return 1;
		}
            }
        }
        return 0;
    }
    //
    bool has_atribute(const char *key) {
        for (int i=0; i<argc; i++) {
            if (!strcmp(argv[i], key)) {
                return 1;
            }
        }
        return 0;
    }
    //
    bool get_file(std::string &file, unsigned int f_number = 0);
    //
    bool get_option(const char *key);
    //
    void get_option(const char *key, bool& value);
    //
    inline void get_option(const std::string key, bool& value){get_option(key.c_str(), value);};
    //
    bool cerrToFile(const std::string filename);
    //
    bool coutToFile(const std::string filename);
    //
    ~CommandLineParser();
    void doneReading();
    void addHelp(const std::string& com, const std::string help);
private:
    int argc;
    char **argv;
    std::streambuf* cerr_strbuf;                    ///< just to save the buffer for cerr.
    std::ofstream* debug;                           ///< The buffer for the std::cerr
    std::streambuf* cout_strbuf;                    ///< just to save the buffer for cout.
    std::ofstream* output;                           ///< The buffer for the std::cout
    struct manint{
      std::string command;
      std::string help;
      manint(const std::string com,const  std::string& he): command(com), help(he){ };
      manint(): command(), help(){ };
    };
    std::vector<manint> man;
};

#ifdef MAIN
/// to get the compilation date of the program (not the library)
__attribute__((noinline)) std::string CommandLineParser::get_exec_compilation_date(){ return std::string("Exec compiled ")+ __DATE__ +" at " +__TIME__ ;};
#endif
}

#endif //COMMANDLINEPARSER_H
