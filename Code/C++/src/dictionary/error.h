#ifndef error_H
#define error_H

#include <string>
using namespace std;

 class error
 {
    const string ERROR_ = "\n\t ERROR: ";
    const string newl = "\n";
   public:
    error();
    error(string);
    error(string,string);
    error(string,string,string);
 };

#endif
