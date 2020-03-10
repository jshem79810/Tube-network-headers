# Tube-network-headers
Header files for reading, writing and storing tube networks

Can be used stand-alone, and forms the basis of other projects in this account.

There are 2 key header files:
  - network_3d.h: contains all objects in the network:: namespace.
  - list_template.h contains all objects in the inlist:: namespace.
Some of the features of network_3d.h depend on inclusion of list_template.h. 

In network_3d.h there are three key objects defined:
  1. network::Node
  2. network::Edge
  3. network::Network
  
A guide to each is given below.

1. Node object.

  The node object stores the 3D position of a node in the network (see position.h for a definition of the position object).     
  
    There are 4 options for constructing a node:
    - Blank constructor: Should be avoided if possible, none of the object variables are initialised.
    - Node(const double & x, const double & y, const double & z): position initialised as (x,y,z).
    - Node(const Position & p): position initialised as p.
    - The default copy constructor.
    
  There is a floating point node property -- Npts -- which is used in some codes, and is set to 1 unless specified otherwise, in which case it is used to codify that the node represents Npts in the same spatial location.

    Functions to access node properties are
    - const Position & get_pos(): returns position for access and copying
    - double get_pos(const int & i): returns position in dimension 0 <= i <= 2 (0=x, 1=y, 2=z)
    - double get_point_count(): returns value of Npts
 
    Functions to update node properties are
    - void set_pos(const Position & p): sets node position to p
    - void set_pos(const double & x, const double & y, const double & z): sets node position to (x,y,z)
    - void set_point_count(const double & n): sets Npts to n
 
2. Edge object.

  The edge object stores smart pointers (shared_ptr) to two node object, the entry (node_in) and exit (node_out), and the geometric properties of the tube(s) connecting them (in a network::Tube_Geometry object). Note the edge object is a template that takes and argument of the form <NodeType>, which can be network::Node or any class derived from it.
  
    Contructors available:
    - Blank: no properties initialised.
    - Egde(std::shared_ptr<NodeType> ni, std::shared_ptr<NodeType> no): specifies and shares ownership of the node in and node out. Shared pointer to node must have been created elsewhere, but even if it goes out of scope before Edge, this will still work. Other properties are initialised to default parameters. Length of tube is calculated from node positions.
    - Egde(std::shared_ptr<NodeType> ni, std::shared_ptr<NodeType> no, const double & rad): same as above but radius of tube is initialised as rad.
    - Egde(std::shared_ptr<NodeType> ni, std::shared_ptr<NodeType> no, const double & Nb, const double & rad): Same as above but now edge represents Nb tubes in parallel (stored in Nbranches -- default otherwise is 1).
    - Egde(std::shared_ptr<NodeType> ni, std::shared_ptr<NodeType> no, const double & Nb, const double & rad, const double & orad): Same as above but tube geometry has "outer radius" differs from inner radius.
    
    Access functions are:
    - 


