#ifndef UTILS_H
#define UTILS_H

#include "stl.h"
#include "data_type.h"
#include <ctype.h>
#include <sstream>

// test if a value is nan
template<typename T>
inline bool is_nan(T x) 
{
    return x!=x;
}

// find index of a string in a string vector (fisrt occurence only)
template<typename T>
inline int which(const vector<T> &data, const T query)
{
	for (int i=0; i<(int)data.size();i++){
		if (data[i]==query) 
			return i;
	}
	return -1;
}

/*inline int which(const vector<int> &data, const int query)
{
        for (int i=0; i<(int)data.size();i++){
                if (data[i]==query)
                        return i;
        }
        return -1;
}

inline int which(const vector<string> &data, const string query)
{
        for (int i=0; i<(int)data.size();i++){
                if (data[i]==query) 
                        return i;
        }
        return -1;
}*/

/* -------------------string operation--------------------- */
// convert number to string
template <typename T>
inline string num2str(T num)
{
	ostringstream ss;
	ss << num;
	return ss.str();	
}

// convert string to uppercase
inline string str2upper(const string &str)
{
	string _str = str;
	int len = (int)_str.size();
	for (int i=0; i<len; i++){
		_str[i] = toupper(_str[i]);	
	}
	return _str;
}
// convert string to lowercase
inline string str2lower(const string &str)
{
        string _str = str;
        int len = (int)_str.size();
        for (int i=0; i<len; i++){
                _str[i] = tolower(_str[i]);
        }
        return _str;
}

// test if a string is lowercase (only test the first char)
inline bool islower(const string &str)
{
	if (str.size()==0){
		cerr << "Warning: size of str is 0" << endl;
		return false;
	}
	if (str[0]>='a' && str[0]<='z'){
		return true;
	}else{
		if (str[0]>='A' && str[0]<='Z'){
			return false;
		}else{
			cerr << "Warning: str is not letter" << endl;
			return false;
		}
	}
}
// split string 
inline vector<NtSeq> split_NtSeq(const string &s, char delim, bool rm_empty=true) {
    stringstream ss(s);
    string item;
    vector<NtSeq> tokens;
    while (getline(ss, item, delim)) {
	if (rm_empty){
        	if (item!="") tokens.push_back(str2NtSeq(item));
	}else{
		tokens.push_back(str2NtSeq(item));
	}
    }
    return tokens;
}

inline vector<string> split(const string &s, char delim, bool rm_empty=true) {
    stringstream ss(s);
    string item;
    vector<string> tokens;
    while (getline(ss, item, delim)) {
	if (rm_empty){
        	if (item!="") tokens.push_back(item);
	}else{
		tokens.push_back(item);
	}
    }
    return tokens;
}

inline vector<int> split_int(const string &s, char delim) {
    stringstream ss(s);
    string item;
    vector<int> tokens;
    while (getline(ss, item, delim)) {
        if (item!="")
                tokens.push_back(atoi(item.c_str()));
    }
    return tokens;
}

// count 

inline map<string, int> strcount(const vector<string> &str_vec, string empty_key="-", bool ignore_case=true, bool is_prob=false)
{
	map<string, int> freq;
	map<string, int>::iterator it = freq.end();
	for (int i=0; i<(int)str_vec.size(); i++){
		string key = str_vec[i]=="" ? "-" : str_vec[i] ;
		if (ignore_case == true){
			for (int j=0; j<(int)key.size(); j++)
                        	key[j] = toupper(key[j]);
		}
			
		it = freq.find(key);
		if (it==freq.end()){
			freq[key] = 1;
		}else{
			it->second++;
		}
	}
	return freq;
}

#endif





