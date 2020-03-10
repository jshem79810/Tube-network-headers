#ifndef POSITION_H
#define POSITION_H

#if defined(_WIN32) || defined(_WIN64)
	#ifndef _USE_MATH_DEFINES
		#define _USE_MATH_DEFINES
	#endif
#endif

#include <math.h>
#include <string>
#include <sstream>

//class for storing and manipulating 3D vectors
namespace pos
{
	template<typename T=double> class Position 
	{
	public:
		T x[3];

		Position(){};
		Position(const T & xi, const T & yi, const T & zi)
		{
			x[0] = xi;
			x[1] = yi;
			x[2] = zi;
		}
		inline T magnitude() const 
		{ 
			Position<T> copy(*(this));
			return sqrt(this->dot(copy)); 
		}

		inline T distance(const Position<T> &pos) const 
		{ 		
			return sqrt((x[0] - pos.x[0]) * (x[0] - pos.x[0])
					  + (x[1] - pos.x[1]) * (x[1] - pos.x[1]) 
					  + (x[2] - pos.x[2]) * (x[2] - pos.x[2]));
		}

		inline T dot(const Position<T> &pos) const 
		{
			return (this->x[0] * pos.x[0] + this->x[1] * pos.x[1] + this->x[2] * pos.x[2]);
		}

		inline T angle(const Position<T> &pos) const 
		{
			return (acos(this->dot(pos) / (this->magnitude()*pos.magnitude())));
		}

		inline Position<T> cross(const Position<T> &pos) const 
		{
			return Position(this->x[1]*pos.x[2] - this->x[2]*pos.x[1],
				            this->x[2]*pos.x[0] - this->x[0]*pos.x[2],
							this->x[0]*pos.x[1] - this->x[1]*pos.x[0]);
		}

		inline std::string toString() const
		{
			std::stringstream ss;
			ss << "(" << this->x[0] << "," << this->x[1] << "," << this->x[2] << ")";
			return (ss.str());
		}

		inline void normalise()
		{
			T norm = magnitude();
			for(int i = 0; i<3; i++)
			{
				x[i] /= norm;
			}
		}

		Position<T> operator-(const Position<T> &pos) const;
		Position<T> operator+(const Position<T> &pos) const;
		bool operator==(const Position<T> &pos) const;
		bool operator!=(const Position<T> &pos) const;
		Position<T> operator*(const T &fac) const;
		Position<T> operator/(const T &fac) const;
	};

	template<typename T> Position<T> Position<T>::operator-(const Position<T> &pos) const 
	{
		Position<T> out;
		for(int i = 0; i<3; i++)
		{
			out.x[i] = x[i] - pos.x[i];
		}
		return out;
	}

	template<typename T> Position<T> Position<T>::operator+(const Position<T> &pos) const 
	{
		Position<T> out;
		for(int i = 0; i<3; i++)
		{
			out.x[i] = x[i] + pos.x[i];
		}
		return out;
	}

	template<typename T> Position<T> Position<T>::operator*(const T &fac) const 
	{
		Position<T> out;
		for(int i = 0; i<3; i++)
		{
			out.x[i] = x[i] * fac;
		}
		return out;
	}

	template<typename T> Position<T> Position<T>::operator/(const T &fac) const 
	{
		Position<T> out;
		for(int i = 0; i<3; i++)
		{
			out.x[i] = x[i] / fac;
		}
		return out;
	}

	template<typename T> bool Position<T>::operator==(const Position<T> &pos) const 
	{
		if(this->x[0] == pos.x[0] && this->x[1] == pos.x[1] && this->x[2] == pos.x[2]) return true;
		return false;
	}

	template<typename T> bool Position<T>::operator!=(const Position<T> &pos) const 
	{
		return (!(*(this) == pos));
	}

}

#endif
