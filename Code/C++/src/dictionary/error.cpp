
 #include <iostream>
 #include <stdlib.h>
 #include "error.h"

 error::error()
 {
   exit(EXIT_FAILURE);
 }

 error::error(string fileName)
 {
   cout << ERROR_ << "unable to open file '" << fileName << "'" << newl << endl;
   exit(EXIT_FAILURE);
 }

 error::error(string fileName,string entry)
 {
   cout << ERROR_ << "entry '" << entry << "' not found in file '" << fileName << "'" << newl << endl;
   exit(EXIT_FAILURE);
 }

 error::error(string fileName, string entry, string terminator)
 {
   cout << ERROR_ << "'" << entry << "' is expected to be terminated by '" << terminator << "' in file '" << fileName << "'" << newl << endl;
   exit(EXIT_FAILURE);
 }

