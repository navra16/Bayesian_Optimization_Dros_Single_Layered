//#include <curand.h>
#include <set>
#include <list>
#include <vector>
#include <memory>
#include "System.h"
#include "SystemBuilder.h"
#include "SystemStructures.h"
# define M_PI 3.14159265358979323846  /* pi */

// Constructor of SystemBuilder that sets the time step and solve time.
SystemBuilder::SystemBuilder(double timestep_, int solve_time_) {
	dt = timestep_/1000.0;
	solve_time = solve_time_;
};

// Function to add a node to the system with the specified coordinates (x, y, z).
void SystemBuilder::addNode(double x, double y, double z) {
	hostSetInfoVecs.nodeLocX.push_back(x);
	hostSetInfoVecs.nodeLocY.push_back(y);
	hostSetInfoVecs.nodeLocZ.push_back(z);

	hostSetInfoVecs.isNodeFixed.push_back(false);
}

// Function to add non-dimensional node data (x1 - x9) for the node. This is likely related to a specific simulation case.
void SystemBuilder::addNndata(double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8, double x9) {
	hostSetInfoVecs.nndata1.push_back(x1);
	hostSetInfoVecs.nndata2.push_back(x2);
	hostSetInfoVecs.nndata3.push_back(x3);
	hostSetInfoVecs.nndata4.push_back(x4);
	hostSetInfoVecs.nndata5.push_back(x5);
	hostSetInfoVecs.nndata6.push_back(x6);
	hostSetInfoVecs.nndata7.push_back(x7);
	hostSetInfoVecs.nndata8.push_back(x8);
	hostSetInfoVecs.nndata9.push_back(x9);
}

// Function to add an edge between two nodes with the specified IDs (idL and idR) and calculate its initial length.
void SystemBuilder::addEdge(int idL, int idR) {
	hostSetInfoVecs.edges2Nodes_1.push_back(idL);
	hostSetInfoVecs.edges2Nodes_2.push_back(idR);

	double xL = hostSetInfoVecs.nodeLocX[idL];
	double yL = hostSetInfoVecs.nodeLocY[idL];
	double zL = hostSetInfoVecs.nodeLocZ[idL];
	double xR = hostSetInfoVecs.nodeLocX[idR];
	double yR = hostSetInfoVecs.nodeLocY[idR];
	double zR = hostSetInfoVecs.nodeLocZ[idR];
	double dist = std::sqrt((xL - xR) * (xL - xR) + (yL - yR) * (yL - yR) + (zL - zR) * (zL - zR));
	hostSetInfoVecs.edge_initial_length.push_back(dist);
  hostSetInfoVecs.edge_rest_length.push_back(dist);
}

// Function to add an edge between two nodes with the specified IDs (idL and idR) and specified initial length (edge_initial_length).
void SystemBuilder::addEdge(int idL, int idR, double edge_initial_length) {
	hostSetInfoVecs.edges2Nodes_1.push_back(idL);
	hostSetInfoVecs.edges2Nodes_2.push_back(idR);
	hostSetInfoVecs.edge_initial_length.push_back(edge_initial_length);
   
}

// Function to add an element (triangle) to the system using three node IDs (idA, idB, and idC).
void SystemBuilder::addElement(int idA, int idB, int idC) {
	hostSetInfoVecs.triangles2Nodes_1.push_back(idA);
	hostSetInfoVecs.triangles2Nodes_2.push_back(idB);
	hostSetInfoVecs.triangles2Nodes_3.push_back(idC);
}

// Function to add an element (triangle) to the system using three edge IDs (idA, idB, and idC).
void SystemBuilder::addElement2Edge(int idA, int idB, int idC) {
	hostSetInfoVecs.triangles2Edges_1.push_back(idA);
	hostSetInfoVecs.triangles2Edges_2.push_back(idB);
	hostSetInfoVecs.triangles2Edges_3.push_back(idC);
}

// Function to add an edge to element mapping by providing element ID (idA) and two edge IDs (idB and idC).
void SystemBuilder::addEdge2Elem(int idA, int idB) {
	hostSetInfoVecs.edges2Triangles_1.push_back(idA);
	hostSetInfoVecs.edges2Triangles_2.push_back(idB);
}

// Function to fix a node with the specified ID (id).
void SystemBuilder::fixNodes(int id) {
	hostSetInfoVecs.isNodeFixed[id] = true;
}

// Function to add a capsid node to the system with the specified coordinates (x, y, z).
void SystemBuilder::addCapsidNode(double x, double y, double z) {
	hostSetInfoVecs.capsidNodeLocX.push_back(x);
	hostSetInfoVecs.capsidNodeLocY.push_back(y);
	hostSetInfoVecs.capsidNodeLocZ.push_back(z);
}

// Function to create the system and return a shared pointer to it.
std::shared_ptr<System> SystemBuilder::createSystem() {
	// Create a shared pointer to the System object.
	std::shared_ptr<System> host_ptr_System = std::make_shared<System>();

	// Set individual parameters in the System object using the values stored in hostSetInfoVecs.
	host_ptr_System->generalParams.dt = dt;
	host_ptr_System->generalParams.solve_time = solve_time;
	host_ptr_System->generalParams.tau = defaultTau;
	host_ptr_System->generalParams.kT = defaultKBT;
  // if you want to add these from the data structures then you should first modify the xml parser and add these variables there. The value can then come from the data structs as with above variables. 
  host_ptr_System->generalParams.epsilon_t = 0.1;  // example value: 10% tangential strain
  host_ptr_System->generalParams.epsilon_r = 1.0;          // example value: initial radius/scaling factor


	host_ptr_System->linearSpringInfoVecs.spring_constant = defaultLinear_Const;
	host_ptr_System->areaTriangleInfoVecs.spring_constant = defaultArea_Const;
	host_ptr_System->bendingTriangleInfoVecs.spring_constant = defaultBending_Const;

	host_ptr_System->ljInfoVecs.epsilon_M = defaultLJ_Eps;
	host_ptr_System->ljInfoVecs.Rmin_M = defaultLJ_Rmin;
	host_ptr_System->ljInfoVecs.Rcutoff_M = defaultLJ_Rmax;
	host_ptr_System->ljInfoVecs.spring_constant = defaultLJ_Const;
	host_ptr_System->ljInfoVecs.LJ_PosX = defaultLJ_X;
	host_ptr_System->ljInfoVecs.LJ_PosY = defaultLJ_Y;
	host_ptr_System->ljInfoVecs.LJ_PosZ = defaultLJ_Z;

	host_ptr_System->linearSpringInfoVecs.scalar_edge_length = defaultEdgeEq;
	host_ptr_System->areaTriangleInfoVecs.initial_area = defaultAreaEq;
	host_ptr_System->bendingTriangleInfoVecs.initial_angle = defaultAngleEq;

	// Initialize the System object with the data stored in hostSetInfoVecs.
	host_ptr_System->initializeSystem(hostSetInfoVecs);

	// Return the shared pointer to the System object.
	return host_ptr_System;
}
