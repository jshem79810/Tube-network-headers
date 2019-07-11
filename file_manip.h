#ifndef FILE_MANIP_H
#define FILE_MANIP_H

#include <sys/stat.h>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

inline bool file_exists(const std::string &name)
{
	struct stat buffer;   
	return (stat (name.c_str(), &buffer) == 0); 
}

inline bool check_infile(const std::string &name)
{
	std::ifstream file;
	file.open(name);
	if (file.good())
	{
		return false;
	}
	std::cout << "Error, could not read input " << name << " \n";    //check file
	return true;
}

inline bool check_outfile(const std::string &name)
{
	std::ofstream file;
	file.open(name);
	if (file.good())
	{
		return false;
	}
	std::cout << "Error, could not read input " << name << " \n";    //check file
	return true;
}

inline size_t count_lines(const std::string &name)
{
	std::ifstream file;
	file.open(name);
	std::string s;
	size_t count = 0;
	while(std::getline(file, s))
	{
		++count;
	}
	file.close();
	return count;
}

inline std::vector<std::string> get_all_lines(const std::string &name)
{
	size_t Nlines = count_lines(name);
	std::vector<std::string> output;
	output.resize(Nlines);
	std::ifstream file;
	file.open(name);
	for(size_t iline=0; iline<Nlines; iline++)
	{
		 std::getline(file, output[iline]);
	}
	file.close();
	return output;
}

template <typename T>
inline T StringToNumber(const std::string &Text)
{
	std::istringstream ss(Text);
	T result;
	return ss >> result ? result : 0;
}

template <typename T>
inline std::vector<std::vector<T>> parse_csv_file(const std::string &name)
{
	std::ifstream file;
	std::string s;
	std::stringstream ss;
	std::vector<std::vector<T>> data;

	file.open(name);
	unsigned count = 0;
	while(getline(file, s))
	{
		ss.str(s);
		std::string nc;
		unsigned n = 0;
		data.resize(count+1);
		while(getline(ss, nc, ','))
		{
			data[count].push_back(StringToNumber<T>(nc));
			n++;
		}
		ss.clear();
		count++;
	}
	file.close();

	return data;
}

inline std::vector<std::string> string_split(const std::string & str, const std::string & delim)
{
    std::vector<std::string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);
        if (pos == std::string::npos) pos = str.length();
        std::string token = str.substr(prev, pos-prev);
        tokens.push_back(token);
        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());
    return tokens;
}

#endif
