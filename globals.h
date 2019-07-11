#ifndef GLOBALS_H
#define GLOBALS_H

#include <iostream>
#include <vector>

inline void abort_on_failure()
{
	exit(EXIT_FAILURE);
}

template <typename T>
inline size_t sizet_round(const T & num)
{
	size_t nfloor = size_t(num);
	if(num - ((T) nfloor) >= 0.5)
	{
		return nfloor + 1;
	}
	else
	{
		return nfloor;
	}

}

template <typename T>
std::vector<T> vectorise(const T * value_array, const int & N)
{
	std::vector<T> v;
	for(int n = 0; n < N; n++)
	{
		v.push_back(value_array[n]);
	}
	return v;
}



#endif