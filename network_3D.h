#ifndef NETWORK_3D_H
#define NETWORK_3D_H

#if defined(_WIN32) || defined(_WIN64)
	#ifndef _USE_MATH_DEFINES
		#define _USE_MATH_DEFINES
	#endif
#endif

#include "position.h"
#include "globals.h"
#include "file_manip.h"
#include <vector>
#include <math.h>
#include <map>
#include <unordered_map>
#include <iostream>
#include <iomanip>

//name of file extensions for node, edge and term node files
#define NODE_FILE_EXT "nodes"
#define BRANCH_FILE_EXT "branches"
#define TERM_NODE_FILE_EXT "termnodes"

namespace network
{
	typedef pos::Position<double> Position;

	//for storing geometry of a tube that can have inner and outer areas
	class TubeGeometry {
	private:
		//-------Member variables-------//
		double radius, length, outer_radius, i_area, o_area, i_vol, o_vol;
		inline double area_calc(const double & rad){ return M_PI * rad * rad; }
		inline double vol_calc(const double & area, const double & len){ return area * len; }
	public:
		//-------Member functions-------//
		TubeGeometry(){};
		TubeGeometry::TubeGeometry(const double & rad, const double & len)
		{
			radius = rad;
			length = len;
			outer_radius = rad;
			i_area = area_calc(radius);
			i_vol = vol_calc(i_area, length);
			o_area = area_calc(outer_radius);
			o_vol = vol_calc(o_area, length);
		}
		TubeGeometry::TubeGeometry(const double & rad, const double & len, const double & o_rad)
		{
			radius = rad;
			length = len;
			outer_radius = o_rad;
			i_area = area_calc(radius);
			i_vol = vol_calc(i_area, length);
			o_area = area_calc(outer_radius);
			o_vol = vol_calc(o_area, length);
		}

		inline double get_inner_radius() const { return radius; }
		inline double get_outer_radius() const { return outer_radius; }
		inline double get_length() const { return length; }
		inline double inner_area() const { return i_area; }
		inline double inner_volume() const { return i_vol; }
		inline double outer_area() const { return o_area; }
		inline double outer_volume() const { return o_vol; }

		inline void update_inner_radius( const double & rad )
		{
			radius = rad;
			i_area = area_calc(radius);
			i_vol = vol_calc(i_area, length);
		}
		inline void update_outer_radius( const double & rad )
		{
			outer_radius = rad;
			o_area = area_calc(outer_radius);
			o_vol = vol_calc(o_area, length);
		}
		inline void update_length( const double & len ){ this->length = len; }
	};

	//base class for a node in the network
	class Node {
	protected:
		//-------Member variables-------//
		Position pos; //, *network_pos;
		double Npts, volume;
		Position *original_pos;
		double *original_Npts;
		double seed_rad;
		inline void copy_node_vals(Node* n)
		{
			this->pos = n->get_pos();
			this->Npts = n->point_count();
			this->volume = n->get_volume();
			this->seed_rad = n->get_seed_rad();
			if(n->has_pos_changed())
			{
				delete (this->original_pos);
				original_pos = new Position(*(n->original_pos));
			}
			if(n->has_npts_changed()) 
			{
				delete (this->original_Npts);
				original_Npts = new double(*(n->original_Npts));
			}
		}
		inline void initialise()
		{
			this->volume = 0;
			this->Npts = 1;
			this->seed_rad = -1;
			this->original_pos = &(this->pos);
			this->original_Npts = &(this->Npts);
		}
	public:
		//-------Member functions-------//
		Node(){};
		Node(const double &x, const double &y, const double &z)
		{
			pos.x[0] = x;
			pos.x[1] = y;
			pos.x[2] = z;
			this->initialise();
		}
		Node(const Position & p)
		{
			this->pos = p;
			this->initialise();
		}

		inline Position get_pos() const { return (this->pos); }
		inline double get_pos(const int & p) const { return (this->pos.x[p]); }
		inline double point_count() const { return (this->Npts); }
		inline Position get_original_pos() const { return (*(this->original_pos)); }
		inline bool has_pos_changed() const { return (this->original_pos != &(this->pos)); }
		inline double get_original_npts() const { return (*(this->original_Npts)); }
		inline bool has_npts_changed() const { return (this->original_Npts != &(this->Npts)); }
		virtual double get_volume() const { return (this->volume); }
		inline bool has_seed_rad_changed() const { return (this->seed_rad > 0); }
		inline double get_seed_rad() const { return (this->seed_rad); }

		inline void set_pos(const Position & p){ this->pos = p; }
		inline void set_point_count(const double & p){ this->Npts = p; }
		inline void set_original_pos(Position *p){ this->original_pos = p; }
		inline void set_original_npts(double *c){ this->original_Npts = c; }
		inline void update_volume(const double & v){ this->volume = v; }
		inline void set_seed_rad(const double & r){ this->seed_rad = r; }
		
		virtual void copy_all_vals(Node* n){ this->copy_node_vals(n); };
	};

	//class for an edge in the network
	template<class NodeType> class Edge
	{
	protected:
		NodeType *node_in, *node_out;
		double Nbranches;
		TubeGeometry *geom;
		TubeGeometry *original_geom;
		inline void set_nodes(NodeType *ni, NodeType *no)
		{
			this->node_in = ni;
			this->node_out = no;
		}
		inline void copy_edge_vals(Edge<NodeType> *edge)
		{
			*(this->geom) = *(edge->get_geom());
			this->Nbranches = edge->branch_count();
			if(this->has_geom_changed())
			{
				delete (this->original_geom);
				original_geom = new TubeGeometry(*(edge->get_original_geom()));
			}
		}
	public:
		//-------Member functions-------//
		Edge(){};
		Edge(NodeType *ni, NodeType *no)
		{
			this->set_nodes(ni, no);
		}
		Edge(NodeType *ni, NodeType *no, const double & Nb, const double & rad)
		{
			this->set_nodes(ni, no);
			double len = ni->get_pos().distance(no->get_pos());
			this->geom = new TubeGeometry(rad, len);
			this->original_geom = this->geom;
			this->Nbranches = Nb;
		}

		inline NodeType* get_node_in() const { return (this->node_in); }
		inline NodeType* get_node_out() const { return (this->node_out); }
		inline double branch_count() const { return (this->Nbranches); }
		inline TubeGeometry* get_geom() const { return (this->geom); }
		inline TubeGeometry* get_original_geom() const { return (this->original_geom); }
		
		inline bool has_geom_changed() const { return (this->original_geom != this->geom); }
		inline void set_original_geom(TubeGeometry *g){ this->original_geom = g; }

		virtual double get_inner_volume() const { return ((this->Nbranches)*(this->geom->inner_volume())); }
		virtual double get_outer_volume() const { return ((this->Nbranches)*(this->geom->outer_volume())); }

		inline void update_length() { this->geom->update_length( this->node_in->get_pos().distance(this->node_out->get_pos())); }

		virtual void update_geometry(const TubeGeometry &tg){ *(this->geom) = tg; };
		virtual void update_geometry(const double & rad, const double & len){ *(this->geom) = TubeGeometry(rad,len); }
		

		virtual void copy_all_vals(Edge<NodeType> *edge){ this->copy_edge_vals(edge); }
	};

	//template class for network
	template<class NodeType, class EdgeType> class Network 
	{
	private:
		void initialise_from_maps(std::map<long int, NodeType*> &node, std::map<long int, EdgeType*> &edge);
		int read_network_files(std::map<long int, NodeType*> &node, std::map<long int, EdgeType*> &edge, 
								const std::string & node_fname, const std::string & edge_fname, 
								const std::string & term_fname, const double & l_scale);
	protected:
		//storage
		std::vector<EdgeType*> EdgeVec;
		std::vector<NodeType*> NodeVec;
		std::unordered_map<NodeType*, size_t> node_index_map;
		std::unordered_map<EdgeType*, size_t> edge_index_map;

		//for navigating 
		size_t term_start;
		std::vector<std::vector<size_t>> edge_in_indices,  edge_out_indices;
		std::vector<std::vector<size_t>> edges_sorted_by_horsfield, edges_sorted_by_weibel;
		std::vector<size_t> node_in_indices, node_out_indices;
		std::vector<size_t> edge_horsfield_order, edge_weibel_order;

		int copy_tree_vals(Network<NodeType,EdgeType> *);
		inline void setup()
		{
			this->update_node_edge_maps();
			if(this->reorder_network()) abort_on_failure();   //returns 1 on error
		}
	public:
		//constructors
		Network(){};
		Network(std::vector<NodeType*> &node, std::vector<EdgeType*> &edge)
		{
			this->NodeVec = node;
			this->EdgeVec = edge;
			this->setup();
		}
		Network(std::map<long int, NodeType*> &node, std::map<long int, EdgeType*> &edge);
		Network(const std::string & node_fname, const std::string & edge_fname, const std::string & term_fname, const double & l_scale);

		//void build_network_matrices();
		int reorder_network();
		void update_node_edge_maps();
		void calc_edge_orders();
		int print_files_for_input(const std::string &fhead, const double & length_scale) const;
		virtual int print_vtk(const std::string &fhead, const double & length_scale) const;
		void copy_structure(Network<NodeType,EdgeType> *);
		virtual int copy_vals(Network<NodeType,EdgeType> *t){ return (this->copy_tree_vals(t)); }

		inline size_t get_first_term_index() const { return (this->term_start); }
		inline size_t count_nodes() const { return (this->NodeVec.size()); }
		inline size_t count_term_nodes() const { return (this->NodeVec.size() - term_start); }
		inline size_t count_edges() const { return (this->EdgeVec.size()); }
		inline size_t count_weibel_orders() const { return (this->edges_sorted_by_weibel.size()); }
		inline size_t count_horsfield_orders() const { return (this->edges_sorted_by_horsfield.size()); }
		inline size_t count_edges_in(const size_t & node_no) const { return (this->edge_in_indices[node_no].size()); }
		inline size_t count_edges_out(const size_t & node_no) const { return (this->edge_out_indices[node_no].size()); }
		inline size_t count_edges_in_weibel_order(const size_t & wo) const { return (this->edges_sorted_by_weibel[wo].size()); }
		inline size_t count_edges_in_horsfield_order(const size_t & ho) const { return (this->edges_sorted_by_horsfield[ho].size()); }
	
		inline size_t get_node_index(NodeType* n) const { return (this->node_index_map.at(n)); }
		inline size_t get_edge_index(EdgeType* e) const { return (this->edge_index_map.at(e)); }
		inline NodeType* get_node(const size_t & k) const { return (this->NodeVec[k]); }
		inline EdgeType* get_edge(const size_t & j) const { return (this->EdgeVec[j]); }
		inline NodeType* get_entry_node() const { return (this->NodeVec[0]); }
		inline size_t get_edge_in_index(const size_t & node_no, const size_t & count) const { return (this->edge_in_indices[node_no][count]); }
		inline size_t get_edge_out_index(const size_t & node_no, const size_t & count) const { return (this->edge_out_indices[node_no][count]); }
		inline size_t get_edge_index_from_weibel_order(const size_t & wo, const size_t & jo) const { return (this->edges_sorted_by_weibel[wo][jo]); }
		inline size_t get_edge_index_from_horsfield_order(const size_t & ho, const size_t & jo) const { return (this->edges_sorted_by_horsfield[ho][jo]); }
		inline size_t get_node_in_index(const size_t & edge_no) const { return (this->node_in_indices[edge_no]); }
		inline size_t get_node_out_index(const size_t & edge_no) const { return (this->node_out_indices[edge_no]); }
		inline size_t get_horsfield_order(const size_t & edge_no) const { return (this->edge_horsfield_order[edge_no]); }
		inline size_t get_weibel_order(const size_t & edge_no) const { return (this->edge_weibel_order[edge_no]); }

		inline bool node_exists(NodeType* n) const { return (this->node_index_map.find(n) != this->node_index_map.end()); }
		inline bool edge_exists(EdgeType* n) const { return (this->edge_index_map.find(n) != this->edge_index_map.end()); }

		void remove_nodes(std::vector<std::size_t> & node_indices);
		void remove_edges(std::vector<std::size_t> & edge_indices);
		double get_edge_volume() const;
		double get_node_volume() const;
		double return_tot_volume() const;
		double count_branches_in(const size_t & node_no) const;
	};

	//constructor from input files
	template<class NodeType, class EdgeType> Network<NodeType,EdgeType>::Network(const std::string & node_fname, const std::string & edge_fname, 
																				 const std::string & term_fname, const double & l_scale)
	{
		std::map<long int,NodeType*> node;     //stores indexed list of nodes -- from file
		std::map<long int,EdgeType*> edge;     //stores indexed list of edges

		if(this->read_network_files(node, edge, node_fname, edge_fname, term_fname, l_scale)) abort_on_failure();
		this->initialise_from_maps(node, edge);
	}

	//read input files
	template<class NodeType, class EdgeType> 
	int Network<NodeType,EdgeType>::read_network_files(std::map<long int, NodeType*> &node, std::map<long int, EdgeType*> &edge, 
														const std::string & node_fname, const std::string & edge_fname, 
														const std::string & term_fname, const double & l_scale)
	{
		//----------get all node positions-----------//
		//check node file and parse
		if(check_infile(node_fname)) return 1;
		std::vector<std::string> node_data = get_all_lines(node_fname);
		size_t Nnodes = node_data.size();   //size of node array

		//assign node array
		for (size_t row = 0; row < Nnodes; ++row)
		{
			long int node_no;
			std::vector<std::string> colon_split = string_split(node_data[row],":");
			//different args separated by colons
			for(size_t ipart=0; ipart<colon_split.size(); ipart++)
			{
				//values separated by commas
				std::vector<std::string> part_split = string_split(colon_split[ipart],",");
				if(ipart == 0)   //first part is core info
				{
					if(part_split.size() < 4)
					{
						std::cout << "Error, not enough columns for node " << row << '\n';
						abort_on_failure();
					}
					else
					{
						node_no = StringToNumber<long int>(part_split[0]);
						//node positions given in mm
						node[node_no] = new NodeType(l_scale*StringToNumber<double>(part_split[1]),
							                         l_scale*StringToNumber<double>(part_split[2]),
													 l_scale*StringToNumber<double>(part_split[3]));   //convert positions to sim units
					}
				}
				else  //other parts are optional extras
				{
					char arg = part_split[0][part_split[0].find_first_not_of(' ')];   //first value is character determining what is to follow
					switch(arg)
					{
					case 'p':   //original node position to follow
						{
							if(part_split.size() < 4)
							{
								std::cout << "Error, not enough columns for node " << row << '\n';
								abort_on_failure();
							}
							else
							{
								node[node_no]->set_original_pos(
									   new Position(l_scale*StringToNumber<double>(part_split[1]),
							                        l_scale*StringToNumber<double>(part_split[2]),
											        l_scale*StringToNumber<double>(part_split[3])));
							}

						} break;

					default:
						{
							std::cout << "Do not recognise node file option " << arg << '\n';
						} break;
					}
				}
			}
		}

		//parse term node file
		//--------get terminal node list-----------//
		if(check_infile(term_fname)) return 1;
		std::vector<std::string> term_node_data = get_all_lines(term_fname);
		size_t Ntermnodes = term_node_data.size();   //size of node array

		//assign termnode array
		for (size_t row = 0; row < Ntermnodes; ++row)
		{
			long int ni;
			std::vector<std::string> colon_split = string_split(term_node_data[row],":");
			//different args separated by colons
			for(size_t ipart=0; ipart<colon_split.size(); ipart++)
			{
				//values separated by commas
				std::vector<std::string> part_split = string_split(colon_split[ipart],",");
				if(ipart == 0)   //first part is core info
				{
					if(part_split.size() < 1)
					{
						std::cout << "Error, not enough columns for node " << row << '\n';
						abort_on_failure();
					}
					else
					{
						ni = StringToNumber<long int>(part_split[0]);
					}
				}
				else  //other parts are optional extras
				{
					char arg = part_split[0][part_split[0].find_first_not_of(' ')];   //first value is character determining what is to follow
					switch(arg)
					{
					case 'p':   //original node position to follow
						{
							if(part_split.size() < 4)
							{
								std::cout << "Error, not enough columns for node " << row << '\n';
								abort_on_failure();
							}
							else
							{
								node[ni]->set_original_pos( 
									new Position(l_scale*StringToNumber<double>(part_split[1]),
							                     l_scale*StringToNumber<double>(part_split[2]),
											     l_scale*StringToNumber<double>(part_split[3])));
							}

						} break;

					case 'n':   //original npts to follow
						{
							if(part_split.size() < 2)
							{
								std::cout << "Error, not enough columns for node " << row << '\n';
								abort_on_failure();
							}
							else
							{
								node[ni]->set_original_npts(new double(StringToNumber<double>(part_split[1])));
							}

						} break;

					case 'r':   //seed rad to follow
						{
							if(part_split.size() < 2)
							{
								std::cout << "Error, not enough columns for node " << row << '\n';
								abort_on_failure();
							}
							else
							{
								node[ni]->set_seed_rad(l_scale*StringToNumber<double>(part_split[1]));
							}

						} break;

					default:
						{
							std::cout << "Do not recognise term node file option " << arg << '\n';
						} break;
					}
				}
			}
		}

		//parse branch file
		//--------get branch list-----------//
		if(check_infile(edge_fname)) return 1;
		std::vector<std::string> branch_data = get_all_lines(edge_fname);
		size_t Nbranches = branch_data.size();   //size of node array

		//assign termnode array
		for (size_t row = 0; row < Nbranches; ++row)
		{
			long int branch_no;
			long int ki;
			long int ko;
			NodeType *n_in;
			NodeType *n_out;
			double radius, orig_rad;
			bool orad = false;
			double Nb = 1;
			std::vector<std::string> colon_split = string_split(branch_data[row],":");
			//different args separated by colons
			for(size_t ipart=0; ipart<colon_split.size(); ipart++)
			{
				//values separated by commas
				std::vector<std::string> part_split = string_split(colon_split[ipart],",");
				if(ipart == 0)   //first part is core info
				{
					if(part_split.size() < 4)
					{
						std::cout << "Error, not enough columns for node " << row << '\n';
						abort_on_failure();
					}
					else
					{
						branch_no = StringToNumber<long int>(part_split[0]);
						ki = StringToNumber<long int>(part_split[1]);
						ko = StringToNumber<long int>(part_split[2]);
						radius = l_scale*StringToNumber<double>(part_split[3]);
						n_in = node[ki];
						n_out = node[ko];
					}
				}
				else  //other parts are optional extras
				{
					char arg = part_split[0][part_split[0].find_first_not_of(' ')];   //first value is character determining what is to follow
					switch(arg)
					{
					case 'n':   //number branches to follow
						{
							if(part_split.size() < 2)
							{
								std::cout << "Error, not enough columns for edge " << row << '\n';
								abort_on_failure();
							}
							else
							{
								Nb = StringToNumber<double>(part_split[1]);
							}

						} break;

					case 'r':   //orig rad to follow
						{
							if(part_split.size() < 2)
							{
								std::cout << "Error, not enough columns for edge " << row << '\n';
								abort_on_failure();
							}
							else
							{
								orig_rad = l_scale*StringToNumber<double>(part_split[1]);
								orad = true;
							}

						} break;

					default:
						{
							std::cout << "Do not recognise term node file option " << arg << '\n';
						} break;
					}
				}
			}
			edge[branch_no] = new EdgeType(n_in, n_out, Nb, radius);
			n_out->set_point_count(Nb);
			if(orad)
			{
				edge[branch_no]->set_original_geom(new TubeGeometry(orig_rad, edge[branch_no]->get_geom()->get_length()));
			}
		}

		return 0;
	}

	//template constructor from maps
	template<class NodeType, class EdgeType> Network<NodeType,EdgeType>::Network(std::map<long int, NodeType*> &node, std::map<long int, EdgeType*> &edge)
	{
		this->initialise_from_maps(node,edge);
	}

	//initialise using maps
	template<class NodeType, class EdgeType> void Network<NodeType,EdgeType>::initialise_from_maps(std::map<long int, NodeType*> &node, std::map<long int, EdgeType*> &edge)
	{
		for (auto it = node.begin(); it != node.end(); ++it)   //iterate over nodes
		{
			this->NodeVec.push_back(it->second);
		}
		for (auto it = edge.begin(); it != edge.end(); ++it)  //make Edge vector
		{
			this->EdgeVec.push_back(it->second);
		}
		this->setup();
	}

	//template network sort
	template<class NodeType, class EdgeType> int Network<NodeType,EdgeType>::reorder_network()   //check everything is connected correctly, and order edge and node vectors so that term edges and nodes at end, and trachea at start
	{
		if(this->EdgeVec.size() > 0)
		{
			size_t count_entrance_nodes = 0;
			size_t trachea_node_index;
			//for storing indices
			std::vector<size_t> trachea_edge_indices, term_node_indices, all_other_node_indices, term_edge_indices, all_other_edge_indices;
			term_node_indices.reserve(this->NodeVec.size()-1);
			all_other_node_indices.reserve(this->NodeVec.size()-1);
			term_edge_indices.reserve(this->EdgeVec.size()-1);
			all_other_edge_indices.reserve(this->EdgeVec.size()-1);

			for(size_t k = 0; k < this->count_nodes(); k++)   //reorder nodes
			{
				if(this->count_edges_in(k) == 0)    //tracheal node
				{
					trachea_node_index = k;
					trachea_edge_indices = this->edge_out_indices[k];
					count_entrance_nodes++;
				}
				else
				{
					if(this->count_edges_out(k) == 0)    //terminal node, goes at end
					{
						term_node_indices.push_back(k);
						term_edge_indices.insert(term_edge_indices.end(), this->edge_in_indices[k].begin(), this->edge_in_indices[k].end());
					}
					else    //all other nodes put at start
					{

						all_other_node_indices.push_back(k);
						for(size_t ni=0; ni<this->edge_in_indices[k].size(); ni++)
						{
							if(this->count_edges_in(this->get_node_in_index(this->get_edge_in_index(k,ni))) > 0)   //add to all other edge indices if not trach edge
							{
								all_other_edge_indices.push_back(this->get_edge_in_index(k,ni));
							}
						}
					}
				}
			}

			//check network is valid
			if(count_entrance_nodes == 0)
			{
				std::cout << "Error, no entrance node.\n";
				return 1;
			}
			if(count_entrance_nodes > 1)
			{
				std::cout << "Error, more than one entrance node.\n";
				return 1;
			}

			//copy old list of nodes and edges
			std::vector<NodeType*> old_node_vec = this->NodeVec;
			std::vector<EdgeType*> old_edge_vec = this->EdgeVec;
			//reset node and edge vector
			this->NodeVec.clear();
			this->NodeVec.resize(old_node_vec.size());
			this->EdgeVec.clear();
			this->EdgeVec.resize(old_edge_vec.size());

			//re-fill node vector -- trachea, then intermediate, then terminal
			this->NodeVec[0] = old_node_vec[trachea_node_index];
			size_t count = 1;
			for(size_t ni = 0; ni < all_other_node_indices.size(); ni++)
			{
				this->NodeVec[ni + count] = old_node_vec[all_other_node_indices[ni]];
			}
			count += all_other_node_indices.size();
			this->term_start = count;
			for(size_t nt = 0; nt < term_node_indices.size(); nt++)
			{
				this->NodeVec[nt + count] = old_node_vec[term_node_indices[nt]];
			}

			//re-fill edge vector -- trachea, then intermediate, then terminal
			count = 0;
			for(size_t j0 = 0; j0 < trachea_edge_indices.size(); j0++)
			{
				this->EdgeVec[j0] = old_edge_vec[trachea_edge_indices[j0]];
				count++;
			}
			for(size_t ji = 0; ji < all_other_edge_indices.size(); ji++)
			{
				this->EdgeVec[ji + count] = old_edge_vec[all_other_edge_indices[ji]];
			}
			count += all_other_edge_indices.size();
			for(size_t jt = 0; jt < term_edge_indices.size(); jt++)
			{
				this->EdgeVec[jt + count] = old_edge_vec[term_edge_indices[jt]];
			}

			this->update_node_edge_maps();
			this->calc_edge_orders();
		}

		return 0;
	}

	//template calculate weibel and horsfield ordering
	template<class NodeType, class EdgeType> void Network<NodeType,EdgeType>::calc_edge_orders()
	{
		this->edge_horsfield_order.clear();
		this->edge_horsfield_order.resize(this->count_edges());
		this->edge_weibel_order.clear();
		this->edge_weibel_order.resize(this->count_edges());

		this->edges_sorted_by_horsfield.clear();
		this->edges_sorted_by_weibel.clear();

		//stores next set of edges to loop over
		std::vector<size_t> next_edges;    
		next_edges.resize(this->count_nodes() - this->term_start);   //starts at terminal edges

		std::vector<bool> edges_done;   //tracks which edges are accounted for
		edges_done.resize(this->count_edges());
		fill(edges_done.begin(),edges_done.end(), false);

		//fill with terminal edges and add to list
		for(size_t k = this->term_start; k < this->count_nodes(); k++)
		{
			next_edges[k - this->term_start] = this->get_edge_in_index(k,0); //list of edges in
		}

		//calculate horsfield order
		size_t order = 0;
		while(next_edges.size()>0)
		{
			std::vector<size_t> new_edges;
			new_edges.reserve(((size_t) (0.5 * next_edges.size())));
			this->edges_sorted_by_horsfield.push_back(next_edges);
			for(size_t ne = 0; ne < next_edges.size(); ne++)
			{
				size_t node_in_index = this->get_node_in_index(next_edges[ne]);
				this->edge_horsfield_order[next_edges[ne]] = order;  //fill in horsfield order
				edges_done[next_edges[ne]] = true;   //mark edge as done

				while(this->count_edges_out(node_in_index) == 1 && this->count_edges_in(node_in_index) > 0 &&
					  this->get_edge(this->get_edge_in_index(node_in_index,0))->branch_count()
											   == this->get_edge(next_edges[ne])->branch_count())  //skip past edges on same branch
				{
					next_edges[ne] = this->get_edge_in_index(node_in_index,0);
					this->edge_horsfield_order[next_edges[ne]] = order;
					this->edges_sorted_by_horsfield[order].push_back(next_edges[ne]);
					edges_done[next_edges[ne]] = true;   //mark edge as done
					node_in_index = this->get_node_in_index(next_edges[ne]);
				}

				if(node_in_index != 0)   //terminate on reaching entry node
				{
					bool can_add = true;  
					for (size_t jo = 0; jo < this->count_edges_out(node_in_index); jo++)  //check if this edge can be updated next
					{
						//edge in can only be updated if all edges out are done
						if(!edges_done[this->get_edge_out_index(node_in_index, jo)]) can_add = false;
					}
					if(can_add) 
					{
						new_edges.push_back(this->get_edge_in_index(node_in_index, 0));  //add to list for next updates
					}
				}
			}
			next_edges = new_edges;
			order++;
		}

		//calculate weibel order
		next_edges.clear();
		next_edges.resize(this->count_edges_out(0));
		for(size_t jo = 0; jo < this->count_edges_out(0); jo++)
		{
			next_edges[jo] = this->get_edge_out_index(0,jo);
		}
		order = 0;
		while(next_edges.size()>0)
		{
			std::vector<size_t> new_edges;
			new_edges.reserve(2*next_edges.size());
			this->edges_sorted_by_weibel.push_back(next_edges);
			for (size_t ne = 0; ne < next_edges.size(); ne++)
			{
				size_t node_out_index = this->get_node_out_index(next_edges[ne]);
				this->edge_weibel_order[next_edges[ne]] = order;
				while(this->count_edges_out(node_out_index) == 1 &&
					  this->get_edge(this->get_edge_out_index(node_out_index,0))->branch_count()
											   == this->get_edge(next_edges[ne])->branch_count())  //skip past edges on same branch
				{
					next_edges[ne] = this->get_edge_out_index(node_out_index,0);
					this->edge_weibel_order[next_edges[ne]] = order;
					this->edges_sorted_by_weibel[order].push_back(next_edges[ne]);
					node_out_index = this->get_node_out_index(next_edges[ne]);
				}
				for (size_t jo = 0; jo < this->count_edges_out(node_out_index); jo++)
				{
					new_edges.push_back(this->get_edge_out_index(node_out_index, jo));
				}
			}
			next_edges = new_edges;
			order++;
		}
	}

	//template to fill all structure info from node and edge vectors
	template<class NodeType, class EdgeType> void Network<NodeType,EdgeType>::update_node_edge_maps()
	{
		//fill map from nodes to indices
		this->node_index_map.clear();
		for(size_t k = 0; k < this->NodeVec.size(); k++)
		{
			this->node_index_map[this->get_node(k)] = k;
		}

		this->edge_index_map.clear();    //clear edge index map

		this->node_in_indices.clear();     //clear node in index map
		this->node_in_indices.resize(this->EdgeVec.size());

		this->node_out_indices.clear();    //clear node out index map
		this->node_out_indices.resize(this->EdgeVec.size());

		this->edge_in_indices.clear();   //clear edge in index map
		this->edge_in_indices.resize(this->NodeVec.size());

		this->edge_out_indices.clear();    //clear edge out index map
		this->edge_out_indices.resize(this->NodeVec.size());

		for(size_t j = 0; j < this->EdgeVec.size(); j++)
		{
			this->edge_index_map[this->EdgeVec[j]] = j;   //map from edge to index

			this->node_in_indices[j] = this->node_index_map[this->EdgeVec[j]->get_node_in()];       //map j to node in index
			this->node_out_indices[j] = this->node_index_map[this->EdgeVec[j]->get_node_out()];    //map j to node out index

			this->edge_in_indices[this->node_out_indices[j]].push_back(j);      //map node out index to j
			this->edge_out_indices[this->node_in_indices[j]].push_back(j);      //map node in index to j
		}
	}

	//template to print files to import network
	template<class NodeType, class EdgeType> int Network<NodeType,EdgeType>::print_files_for_input(const std::string &fhead, const double & length_scale) const 
	{
		std::stringstream ss;
		ss << fhead << '.' << NODE_FILE_EXT;
		std::ofstream node_file;
		node_file.open(ss.str().c_str());
		//double lconv = 1.0 / p->conversions_phys_to_sim[LENGTH_KEY];

		for(size_t k = 0; k < this->NodeVec.size(); k++)
		{
			node_file << k << ", " << this->get_node(k)->get_pos(0) * length_scale << ", " << this->get_node(k)->get_pos(1) * length_scale << ", " << this->get_node(k)->get_pos(2) * length_scale << '\n';
		}
		node_file.close();

		std::ofstream edge_file;
		ss.clear();
		ss.str("");
		ss << fhead << '.' << BRANCH_FILE_EXT;
		edge_file.open(ss.str().c_str());
		for(size_t j = 0; j < this->EdgeVec.size(); j++)
		{
			edge_file << j << ", " << this->get_node_index(this->get_edge(j)->get_node_in()) << ", "
					  << this->get_node_index(this->get_edge(j)->get_node_out()) << ", " << this->get_edge(j)->get_geom()->get_inner_radius() * length_scale
					  << ": n, " << this->get_edge(j)->branch_count();
			if(this->get_edge(j)->has_geom_changed())
			{
				edge_file << ": r, " << this->get_edge(j)->get_original_geom()->get_inner_radius();
			}

			edge_file << '\n';
		}
		edge_file.close();

		std::ofstream tnode_file;
		ss.clear();
		ss.str("");
		ss << fhead << '.' << TERM_NODE_FILE_EXT;
		tnode_file.open(ss.str().c_str());
		for(size_t k = this->get_first_term_index(); k < this->NodeVec.size(); k++)
		{
			tnode_file << k;
			if(this->get_node(k)->has_pos_changed())   //term nodes may have been averaged, give info about original term nodes if so
			{
				tnode_file << ": p";
				for(size_t np = 0; np < 3; np++) tnode_file <<  ", " << this->get_node(k)->get_original_pos().x[np];
			}
			if(this->get_node(k)->has_npts_changed())
			{	
				tnode_file <<  ": n, " <<  this->get_node(k)->get_original_npts();
			}
			if(this->get_node(k)->has_seed_rad_changed())
			{
				tnode_file << ": r, " << this->get_node(k)->get_seed_rad();
			}
			tnode_file << '\n';
		}

		tnode_file.close();

		return 0;
	}

	template<class NodeType, class EdgeType> int Network<NodeType,EdgeType>::print_vtk(const std::string &fhead, const double & length_scale) const 
	{
		std::stringstream ss;
		ss << fhead << ".vtk";
		std::string filename = ss.str().c_str();
		if(check_outfile(filename)) return 1;   //check if file can be written to
		std::ofstream output;
		output.open(filename);
		output << std::fixed << std::setprecision(12);
		
		size_t Npts =  this->count_nodes();
		output <<  "# vtk DataFile Version 3.0\nNetwork_data\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		output <<  "POINTS " << Npts << " float\n";

		for (size_t k = 0; k < Npts; k++)
		{
			output << this->get_node(k)->get_pos(0) * length_scale << ' ' << this->get_node(k)->get_pos(1) * length_scale << ' ' << this->get_node(k)->get_pos(2) * length_scale << '\n';
		}
		//could define extra points and edges here for shapes for each acinus

		size_t Ncells = this->count_edges();

		output <<  "\n\nCELLS " << Ncells << ' ' << 3*Ncells << '\n';
		for (size_t j = 0; j < Ncells; j++)
		{
			output << 2 << ' ' << this->get_node_in_index(j) << ' ' << this->get_node_out_index(j) << '\n';  //node numbers are order they are printed to vtk
		}

		output <<  "\n\nCELL_TYPES " << Ncells << '\n';  //all cell types are edges
		for (size_t j = 0; j < Ncells; j++) output << "3\n";

		output <<  "\n\nCELL_DATA " << Ncells << "\n";
		output <<  "\nSCALARS Nbranches float\nLOOKUP_TABLE default\n";
		for (size_t j = 0; j < Ncells; j++)
		{
			output << this->get_edge(j)->branch_count() << '\n';  //node numbers are order they are printed to vtk
		}

		output <<  "\nSCALARS Inner_rad float\nLOOKUP_TABLE default\n";
		for (size_t j = 0; j < Ncells; j++)
		{
			output << this->get_edge(j)->get_geom()->get_inner_radius() * length_scale << '\n';  //node numbers are order they are printed to vtk
		}

		output.close();

		return 0;
	}

	template<class NodeType, class EdgeType> void Network<NodeType,EdgeType>::copy_structure(Network<NodeType,EdgeType> *tree_in)
	{
		this->NodeVec.clear();
		this->NodeVec.resize(tree_in->count_nodes());
		for(size_t k = 0; k < tree_in->count_nodes(); k++)
		{
			this->NodeVec[k] = new NodeType();    //create new node   
			this->get_node(k)->copy_structure(tree_in->get_node(k), false);   //copy the node structure
		}

		this->EdgeVec.clear();
		this->EdgeVec.resize(tree_in->count_edges());
		for(size_t j = 0; j < tree_in->count_edges(); j++)
		{
			size_t kin = tree_in->get_node_in_index(j);
			size_t kout = tree_in->get_node_out_index(j);
			this->EdgeVec[j] = new EdgeType(this->get_node(kin), this->get_node(kout));   //create new edge from new nodes
			this->get_edge(j)->copy_structure(tree_in->get_edge(j));
		}	
		this->node_index_map = tree_in->node_index_map;
		this->edge_index_map = tree_in->edge_index_map;
		this->term_start = tree_in->term_start;
		this->edge_in_indices = tree_in->edge_in_indices;
		this->edge_out_indices = tree_in->edge_out_indices;
		this->edges_sorted_by_horsfield = tree_in->edges_sorted_by_horsfield;
		this->edges_sorted_by_weibel = tree_in->edges_sorted_by_weibel;
		this->node_in_indices = tree_in->node_in_indices;
		this->node_out_indices = tree_in->node_out_indices;
		this->edge_horsfield_order = tree_in->edge_horsfield_order;
		this->edge_weibel_order = tree_in->edge_weibel_order;
	}

	template<class NodeType, class EdgeType> int Network<NodeType,EdgeType>::copy_tree_vals(Network<NodeType,EdgeType> *tree_in)
	{
		if(this->count_nodes() == tree_in->count_nodes() && this->count_edges() == tree_in->count_edges())
		{
			for(size_t k = 0; k < this->count_nodes(); k++) 
			{
				this->get_node(k)->copy_all_vals(tree_in->get_node(k));
			}
			for(size_t j = 0; j < this->count_edges(); j++)
			{
				this->get_edge(j)->copy_all_vals(tree_in->get_edge(j));
			}

			return 0;
		}
		else
		{
			std::cout << "Cannot copy acini with different structures.\n";
			return 1;
		}
	}

	template<class NodeType, class EdgeType> void Network<NodeType,EdgeType>::remove_nodes(std::vector<std::size_t> & node_indices)
	{
		//sort into descending order
		std::sort(node_indices.begin(), node_indices.end(), std::greater<size_t>());
		//remove nodes
		for(size_t kd = 0; kd < node_indices.size(); kd++)
		{
			this->NodeVec.erase(this->NodeVec.begin() + node_indices[kd]);
		}

	}
	
	template<class NodeType, class EdgeType> void Network<NodeType,EdgeType>::remove_edges(std::vector<std::size_t> & edge_indices)
	{
		//sort into descending order
		std::sort(edge_indices.begin(), edge_indices.end(), std::greater<size_t>());
		//remove edges
		for(size_t jd = 0; jd < edge_indices.size(); jd++)
		{
			this->EdgeVec.erase(this->EdgeVec.begin() + edge_indices[jd]);
		}
	}

	template<class NodeType, class EdgeType> double Network<NodeType,EdgeType>::get_edge_volume() const 
	{
		double vol = 0;
		for(size_t j = 0; j < this->EdgeVec.size(); j++)
		{
			vol += this->EdgeVec[j]->get_outer_volume();
		}
		return vol;
	}

	template<class NodeType, class EdgeType> double Network<NodeType,EdgeType>::get_node_volume() const 
	{
		double vol = 0;
		for(size_t k = 0; k < this->NodeVec.size(); k++)
		{
			vol += this->get_node(k)->get_volume();
		}
		return vol;
	}

	template<class NodeType, class EdgeType> double Network<NodeType,EdgeType>::return_tot_volume() const 
	{
		double vol = this->get_edge_volume();
		vol += this->get_node_volume();
		return vol;
	}

	template<class NodeType, class EdgeType> double Network<NodeType,EdgeType>::count_branches_in(const size_t & node_no) const 
	{
		double N = 0;
		for (size_t j = 0; j < this->count_edges_in(node_no); j++)
		{
			N += this->get_edge(this->get_edge_in_index(node_no, j))->branch_count();
		}
		return N; 
	}

}


//void AirwayNetwork::build_network_matrices()
//{
//	size_t Nnodes = NodeVec.size();
//	AdjacencyCond = new Eigen::SparseMatrix<int>(Nnodes, Nnodes); //declare matrices
//	DegreeCond = new Eigen::SparseMatrix<int>(Nnodes, Nnodes);
//
//	std::vector<Eigen::Triplet<int>> adj_triplet, deg_triplet;   //continue HERE
//	deg_triplet.reserve(Nnodes);
//	adj_triplet.reserve(2*EdgeVec.size());
//	for (size_t k = 0; k < Nnodes; k++) //update edge references of nodes
//	{
//		int degree = 0;
//		for (size_t ki = 0; ki < NodeVec[k]->n_edges_in(); ki++)   //add edges going into node
//		{
//			degree += NodeVec[k]->get_edge_in(ki)->branch_count();
//			long int ni_index = NodeVec[k]->get_edge_in(ki)->get_node_in_index();
//			adj_triplet.push_back(Eigen::Triplet<int>(k, ni_index, NodeVec[k]->get_edge_in(ki)->branch_count()));
//			//cout << adj_triplet[adj_triplet.size()-1].col() << ' ' << adj_triplet[adj_triplet.size()-1].row() << ' ' << adj_triplet[adj_triplet.size()-1].value() << '\n';
//		}
//		for (size_t ko = 0; ko < NodeVec[k]->n_edges_out(); ko++)  //add edges coming out of node
//		{
//			degree += NodeVec[k]->get_edge_out(ko)->branch_count();
//			long int no_index = NodeVec[k]->get_edge_out(ko)->get_node_out_index();
//			adj_triplet.push_back(Eigen::Triplet<int>(k, no_index, NodeVec[k]->get_edge_out(ko)->branch_count()));
//			//cout << adj_triplet[adj_triplet.size()-1].col() << ' ' << adj_triplet[adj_triplet.size()-1].row() << ' ' << adj_triplet[adj_triplet.size()-1].value() << '\n';
//		}
//		deg_triplet.push_back(Eigen::Triplet<int>(k,k,degree));
//		//cout << deg_triplet[deg_triplet.size()-1].col() << ' ' << deg_triplet[deg_triplet.size()-1].row() << ' ' << deg_triplet[deg_triplet.size()-1].value() << '\n';
//	}
//
//	AdjacencyCond->setFromTriplets(adj_triplet.begin(), adj_triplet.end());
//	DegreeCond->setFromTriplets(deg_triplet.begin(), deg_triplet.end());
//}

#endif
