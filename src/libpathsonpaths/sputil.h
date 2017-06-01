#ifndef SPUTIL_H
#define SPUTIL_H


#include <ctime>
#include <cstdlib>

#include <iostream>
#include <sstream>
#include <string>
#include <exception>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>


using boost::lexical_cast;

template<unsigned N>
inline double pow(double x)
	{
	double res = 1.0;

	for (unsigned i=0; i<N; i++)
		res *= x;

	return res;
	}


inline void get_time_short(std::string & smess)
	/*get the current time and put only the day number and time in smess*/
	{
	smess = lexical_cast<std::string>(time(NULL));
	}

// split string into two halves at position pos and cast results to T
template<typename T>
void splitStr(const std::string & str, int pos, T & x1, T & x2)
	{
	x1 = lexical_cast<T>(str.substr(0, pos));
	x2 = lexical_cast<T>(str.substr(pos+1, str.size()-pos-1));
	}

// split string at separator sep and feed pieces into iterator iter
template<typename ITER>
void splitStr(const std::string & str, char sep, ITER & iter)
	{
	typedef typename ITER::container_type::value_type value_type;
	size_t last = 0;
	for (size_t i=0; i<=str.size(); i++)
		if (i==str.size() || str[i]==sep)
			{
			*iter = lexical_cast<value_type>(str.substr(last, i-last));
			last = i + 1;
			}
	}

// read next non-blank line
inline void skip_space(std::istream & inp, std::string & str)
	{
	while(getline(inp, str) && boost::all(str, boost::is_space()));
	}

template<typename T> 
bool get_value (std::istream & inp_file, T & value)
	{
	static std::string str;

	skip_space(inp_file, str);
	if (!inp_file) return false;

	std::istringstream tmp_s;
	tmp_s.str(str);
	tmp_s >> value;
	if (tmp_s.fail()) return false;

	return true;
	}


#endif // SPUTIL_H
