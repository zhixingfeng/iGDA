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

inline vector<int64_t> split_int64_t(const string &s, char delim) {
    stringstream ss(s);
    string item;
    vector<int64_t> tokens;
    while (getline(ss, item, delim)) {
        if (item!="")
            tokens.push_back(stoll(item));
    }
    return tokens;
}


inline vector<double> split_double(const string &s, char delim) {
    stringstream ss(s);
    string item;
    vector<double> tokens;
    while (getline(ss, item, delim)) {
        if (item!="")
            tokens.push_back(stod(item));
    }
    return tokens;
}

// split lines of a file (suffix starts from 0)
inline void split_file(string infile, string outprefix, size_t nlines)
{
    ifstream fs_infile; open_infile(fs_infile, infile);
    ofstream fs_outfile;
    int64_t cur_id = -1;
    size_t num = 0;
    while (true) {
        string buf;
        getline(fs_infile, buf);
        if (fs_infile.eof())
            break;
        int64_t id = int64_t(num / nlines);
        if (id != cur_id){
            fs_outfile.close();
            open_outfile(fs_outfile, outprefix + ".part." + to_string(id));
            cur_id = id;
        }
        if (cur_id == -1)
            throw runtime_error("split_file(): cur_id == -1");
        fs_outfile << buf << endl;
        ++num;
    }
    
    fs_infile.close();
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





