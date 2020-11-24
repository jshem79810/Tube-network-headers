#ifndef LIST_TEMPLATE_H
#define LIST_TEMPLATE_H

#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include "file_manip.h"
#include "globals.h"

//classes for reading and storing options and parameters

namespace inlist
{
	//generic template for options and params lists
	template <class T1, class T2>
	class List
	{
	protected:
		std::unordered_map<std::string, T1> dict1;
		std::unordered_map<std::string, T2> dict2;
	public:
		//functions
		List(){};

		void read_file(const std::string &);
		//parse template function
		bool parse(const std::string &name, const std::string &nc);
		inline void add(const std::string & name, T1 entry){ this->dict1[name] = entry; }
		inline void add(const std::string & name, T2 entry){ this->dict2[name] = entry; }
		void get(const std::string & name, T1 & t) const { t = this->dict1.at(name);}
		void get(const std::string & name, T2 & t) const { t = this->dict2.at(name);}
	};
	
	template<class T1, class T2> 
	void List<T1,T2>::read_file(const std::string &filename)  //constructer with filename argument
	{
		char buff[1000];   //stores line
		unsigned buffsize = 1000;   
		std::string temp, nc;   //stores parsed info
		std::stringstream ss;
		std::cout << "Reading file\n";
		if (!check_infile(filename))
		{
			std::ifstream infile;
			infile.open(filename);
			while (infile.good())                            //read
			{
				infile.getline(buff, buffsize);             //read line from file
				ss.str(buff);
				temp.clear();
				ss >> std::skipws >> temp;                //parse line
				if (temp.size() > 0 && temp[0] != '%')
				{
					ss >> std::skipws >> nc;                     //read value into string
					std::cout << temp << ' ' << nc << std::endl;
					if (!parse(temp, nc))    //send to parse option
					{
						std::cout << "Did not recognise option/parameter " << temp << ".\n";
					}
				}
				ss.str("");
				ss.clear();				//throw away rest of line after
			}
			infile.close();
		}
	}

	template <class T1, class T2>
	inline bool List<T1, T2>::parse(const std::string &name, const std::string &nc)
	{
		//check it is a valid option
		if(dict1.find(name) != dict1.end()) 
		{
			dict1[name]->read(nc);
			return true;
		}
		if(dict2.find(name) != dict2.end()) //check it is a valid option
		{
			dict2[name]->read(nc);
			return true;
		}

		return false;
	}

	//template for generic input, can be optino or parameter
	template <typename T>
	class Input{
	protected:
		T Default_value, value;
		std::string name;
	public:
		Input(){};
		Input(const T & val, const std::string & nam)
		{
			value = val;
			Default_value = val;
			name = nam;
		}

		T get_value(){ return this->value; }
		inline std::string get_name(){ return this->name; }
		inline void set_default_value(const T & val)
		{ 
			Default_value = val;
		}
		virtual void set_to_default()
		{ 
			value = Default_value;
		}
		virtual void update_value( const T & v ){ value = v; }
		virtual void read(const std::string &){};
		virtual std::string print() const {return "Undefined";}
	};

	//option class
	template <typename T>
	class Option: public Input<T>{
	protected:
		std::vector<T> possible_values;
		std::vector<std::string> value_names;
	public:
		std::string print() const 
		{
			for (unsigned i = 0; i < ((unsigned) possible_values.size()); i++)
			{
				if (this->value == possible_values[i])
				{
					return value_names[i];
				}
			}
			return "unknown";
		}

		Option<T>():Input<T>(){};
		Option<T>(const T & val, const std::string & nam):Input<T>(val, nam){};
		Option<T>(const T & val, const std::string & nam, const T *option_list,
				  const std::string *option_name_list, const int & option_count):Input<T>(val, nam)
		{
			possible_values = vectorise(option_list, option_count);
			value_names = vectorise(option_name_list, option_count);
		}
		virtual void read(const std::string & code){};
		std::string get_value_name() const 
		{
			for(int i = 0; i < this->possible_values.size(); i++)
			{
				if(this->value == this->possible_values[i]) return this->value_names[i];
			}
			return "Error";
		}
	};

	//character option -- for multiple choice options
	template<> inline void Option<char>::read(const std::string & code)
	{
		bool error = true;
		std::cout << code << std::endl;
		std::cout << possible_values.size() << std::endl;
		for (unsigned i = 0; i < ((unsigned) possible_values.size()); ++i)
		{
			if(code[0] == possible_values[i])
			{
				this->value = possible_values[i];
				error = false;
			}
		}
		if(error)
		{
			std::cout << name << " option code " << code << " not recognised\n";
		}
	}

	//boolian option
	template<> inline void Option<bool>::read(const std::string & code)
	{
		std::cout << code << std::endl;
		switch(code[0])
		{
			case 't':
			case 'T':
			{
				this->value = true;
			} break;
		
			case 'f':
			case 'F':
			{
				this->value = false;
			} break;
		
			default:
			{
				std::cout << name << " option code " << code << " not recognised\n";
			}
		}
	}

	//parameter -- for storing values that may have physical units
	template <typename T>
	class Parameter: public Input<T>
	{
	public:
		Parameter<T>():Input<T>(){};
		Parameter(const T & val, const std::string & nam):Input<T>(val, nam){};
		virtual void read(const std::string & code){};

		virtual void set_conversion( T * ){};
		virtual bool isOK() const { return true; }
		virtual T get_phys_value() const { return this->value; }
		virtual T get_SI_value() const { return this->value; }
		virtual void set_phys_value( const T & val ){ this->value = val; }
		virtual void calc_sim_from_phys_value(){};
		virtual std::string phys_value_string() const { return ""; }
	};

	//integer parameters
	template<> inline void Parameter<int>::read(const std::string &c){ value = atoi(c.c_str());}  //this is common to all params of type double
	
	//double parameters
	template<> inline void Parameter<double>::read(const std::string &c){ value = atof(c.c_str()); }

	//template for list of two types of options
	template<typename T1, typename T2> class OptionList: public List<std::shared_ptr<Option<T1>>,
	                                                                 std::shared_ptr<Option<T2>>>
	{
	protected:
		std::unordered_map<std::string, std::vector<std::string>> filenames;
	public:

		OptionList(){};   //default constructor is blank
		template<typename T3> T3 get_option_value(const std::string & name) const
		{
			T3 o;
			if(typeid(T3) == typeid(T1) || typeid(T3) == typeid(T2))
			{
				this->get_option_value(name, o);
			}
			else
			{
				std::cerr << "Error: option type not recognised." << std::endl;
			}
			return o;
			
		}
		void get_option_value(const std::string & name, T1 & o) const { o = this->dict1.at(name)->get_value(); }
		void get_option_value(const std::string & name, T2 & o) const { o = this->dict2.at(name)->get_value(); }
		std::string get_option_name(const std::string & name) const
		{
			if(this->dict1.find(name) != this->dict1.end())  this->dict1.at(name)->get_value_name();
			if(this->dict2.find(name) != this->dict2.end())  this->dict2.at(name)->get_value_name();
			return "";
		}
		
		inline void add_filename(const std::string & code, const std::string & fname)
		{ 
			if(!this->filename_exists(code)) this->filenames[code] = std::vector<std::string>(); 
			this->filenames[code].push_back(fname);
		}
		inline size_t count_files_with_ext(const std::string & code) const { return this->filenames.at(code).size(); }
		inline std::string get_filename(const std::string & code) const 
		{ 
			if(this->filenames.at(code).size() > 1) std::cout << "Warning, more than one file with extension " << code << '\n';
			return this->filenames.at(code)[0]; 
		}
		inline std::string get_filename(const std::string & code, const size_t & n) const { return this->filenames.at(code)[n]; }
		inline bool filename_exists(const std::string & code) const { return (this->filenames.find(code) != this->filenames.end()); }
		void get_filenames_from_args(const std::vector<std::string> & extensions, const int &argc, char** argv);
	};

	template<typename T1, typename T2> void OptionList<T1,T2>::get_filenames_from_args
		              (const std::vector<std::string> & extensions, const int &argc, char** argv) 
	{
		//arguments sorted into ordered map, where keyword is their file extension
		for(int n = 1; n < argc; n++)
		{
			std::string fname = argv[n];
			size_t i = fname.rfind('.', fname.length());
			std::string ext = fname.substr(i+1, fname.length()-1);   //get file extension
			bool sorted = false;
			for (unsigned m = 0; m < extensions.size(); m++)
			{
				if(extensions[m] == ext)  //if filename extension matches
				{
					this->add_filename(ext, fname);
					sorted = true;
				}
			}
			if (!sorted)
			{
				std::cout << "Input " << fname << " not used -- extension not recognised.\n";  //error message
			}
		}
		//tested
	}

	//template for two types of parameters
	template<typename T1, typename T2> class ParameterList: public List<std::shared_ptr<Parameter<T1>>,
	                                                                    std::shared_ptr<Parameter<T2>>>
	{
	protected:
		std::unordered_map<std::string, double> conversions_phys_to_sim;
		inline void check_OK()
		{
			for (auto it = this->dict1.begin(); it!= this->dict1.end(); ++it)
			{
				if(!it->second->isOK()) it->second->set_to_default();
			}
			for (auto it =  this->dict2.begin(); it!= this->dict2.end(); ++it)
			{
				if(!it->second->isOK()) it->second->set_to_default();
			}
		}
	public:
		ParameterList(){};    //default constructor uses pre-defined default params
		template<typename T3> T3 get_param_value(const std::string & name) const
		{
			T3 p;
			if(typeid(T3) == typeid(T1) || typeid(T3) == typeid(T2))
			{
				this->get_param_value(name, p);
			}
			else
			{
				std::cerr << "Error: param type not recognised." << std::endl;
			}
			return p;
		}
		void get_param_value(const std::string & name, T1 & p) const { p = this->dict1.at(name)->get_value(); }
		void get_param_value(const std::string & name, T2 & p) const { p = this->dict2.at(name)->get_value(); }
		void set_param_phys_value(const std::string & name, const T1 & p)
		{
			if(this->dict1.find(name) != this->dict1.end())  this->dict1[name]->set_phys_value(p);
		}
		void set_param_phys_value(const std::string & name, const T2 & p)
		{
			if(this->dict2.find(name) != this->dict2.end())  this->dict2[name]->set_phys_value(p);
		}
		void set_param_to_default(const std::string & name)
		{
			if(this->dict1.find(name) != this->dict1.end())  this->dict1[name]->set_to_default();
			if(this->dict2.find(name) != this->dict2.end())  this->dict2[name]->set_to_default();
		}
		
		
		inline void set_conversion(const std::string & name, const double & val){ this->conversions_phys_to_sim[name] = val; }
		inline double get_conversion(const std::string & name) const { return (this->conversions_phys_to_sim.at(name)); }
		inline double* get_conversion_ptr(const std::string & name) const { return &(this->conversions_phys_to_sim.at(name)); }
		virtual void check_and_convert(OptionList<char, bool> *){ this->check_OK(); }
	};

}

#endif
