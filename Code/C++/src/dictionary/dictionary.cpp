#include <iostream>
#include <sstream>
#include <fstream>
#include "dictionary.h"

 dictionary::dictionary (const vector<string>& entryList)
 :
 fileName("inputDict"),
 entryNames(entryList),
 entryValues(entryList.size())
 {
   read(dictToStr);

   for (int it=0; it < entryList.size(); it++)
        findEntry(entryValues[it],dictToStr,entryList[it]);
 }

 dictionary::dictionary (const dictionary& dict, string subDictName, const vector<string>& entryList)
 :
 entryNames(entryList),
 entryValues(entryList.size())
 {
   subDictExist_ = dict.getSubDict(dictToStr, subDictName);

   for (int it=0; it < entryList.size() && subDictExist_; it++)
        findEntry(entryValues[it],dictToStr,entryList[it]);
 }

  void dictionary::squeezeString(string& str)
 {
   int count=0;

   for (int i=0; str[i]; i++)
   {
       if (str[i] != ' ' && str[i] != '\t')
        str[count++]=str[i];
   }
   str.resize(count);
 }

  void dictionary::read(string& dictStr)
 {
   ifstream idict(fileName);
   if (idict.is_open())
   {
     string line;
     while ( getline(idict,line) )
     {
	squeezeString(line);
	dictStr += line;
     }
     idict.close();
   }
   else
      error msg(fileName);
 }

  void dictionary::findEntry(string& entry, const string& dictStr,string entryName)
 {
    size_t pos = dictStr.find(entryName); 
    if (pos == string::npos)
	error msg(fileName,entryName);

    pos += entryName.length();
    size_t terminator = dictStr.find(";",pos);
    if ( terminator == string::npos )
	error msg(fileName,entryName,";");

    entry=dictStr.substr(pos,terminator-pos);
 }

  string dictionary::get(int i) const //string entry)
 {
/*
   int it=0;
   while (entryNames[it].compare(entry) != 0) // && it < entryNames.size())
	it++;
*/
   return entryValues[i];
 }

  double dictionary::getfl(int i) const // string entry)
 {
   stringstream ss(entryValues[i]); 
   double x; ss >> x;

   return x;
 }

/*
 subdictionary::subdictionary(string& parentDictStr, string name)
 :
 dictionary(parentDictStr,name)
*/
 bool dictionary::getSubDict(string& subDictStr, string name) const
 {
   name = name + "{";
   size_t pos = dictToStr.find(name);
   if (pos != string::npos)
   {
     // error msg(fileName,name);

     pos += name.length();
     size_t terminator = dictToStr.find("}",pos);
     if (terminator == string::npos)
        error msg(fileName,name,"}");

     subDictStr=dictToStr.substr(pos,terminator-pos);

     return true;
   }

   return false;
//   dictToStr.erase(pos-name.length(),terminator-pos+name.length()+1);
//   findEntry(type,subDictStr,"name");

//   cout << "\n\t" << name << ": " << type << endl;
 }

