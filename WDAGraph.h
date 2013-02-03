/*
 * WDAGraph.h
 *
 *	This is the header file for the WDAGraph object. WDAGraph is an implementation
 *	of a weighted directed acyclic graph.  It is implmented using adjacency lists
 *  as the typical expected use is a sparsely connected graph (e.g. a sequence
 *  graph).
 *
 *  ************* WARNING ***********************************
 *  *  There is no error handling in place for this object! *
 *  *  This means that if the file does not exist or is     *
 *  *  formatted incorrectly, you will get an error and     *
 *  *  will not be able to use this object.                 *
 *  *                                                       *
 *  *  If this code gets moved to a production setting      *
 *  *  appropriate error handling should be implemented!    *
 *  *********************************************************
 *
 *  Typical use for the file would be to use the WDAGraph(graphFileName)
 *  constructor to create the object.  This will automatically open the
 *  graph describer file specified by graphfileName, and read its contents
 *  storing them in the vertices vector and edges map.
 *
 *  The file is expected to be formatted as follows:
 *
 *	  1. A list of vertices, with each vertex on a separate line.  The vertices
 *		 are in depth order (parents preced children), and this order will be
 *		 used in the findHieghestWeightPath() method.  Each line for a vertex 
 *		 should have the following format
 *
 *			V label <START or END>
 *
 *		  - The 'V" char indicates the line is for a vertex.
 *		  - The label for the vertex should be unique, i.e. different for
 *			different vertices.
 *		  - The string "START" if the path is to be constrained to start 
 *			at this vertex
 *		  - The string "END" if the path is to be constrained to end at this
 *			vertex
 *
 *		 At most one vertex should be designated the START and at most one
 *		 vertex should be designated the END. If none are, the path is assumed
 *		 to be unconstrained.
 *
 *	  2. A list of edges, with each edge on a separate line.  The line should have
 *		 the following format:
 *
 *			E label start_vertex end_vertex weight
 *
 *		  - The 'E" char indicates the line is for an edge.
 *		  - The label for the edge does not have to be unique (i.e. different edges
 *			can have the same label)
 *		  - start_vertex is the label of the edge's beginning vertex
 *		  - end_verex is the label of the edge's ending vertex
 *		  - weight is the numerical weight attached to the edge
 *
 *  After creating the object, typical use would be to call the findHighestWeightPath()
 *  which will find the path with the highest weight using dynamic programming.
 *
 *  Finally one would typically call the resultString() method to get a formatted set
 *	of results indicating the path with the highest weight.
 *
 *  Created on: 1-21-13
 *      Author: tomkolar
 */

#ifndef WDAGraph_H
#define WDAGraph_H
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <climits>
using namespace std;

class WDAGraph {

public:

	// Constuctors
	// ==============================================
	WDAGraph();
	WDAGraph(string& aGraphFileName);

	// Destructor
	// =============================================
	~WDAGraph();

	// Public Methods
	// =============================================

	// findHighestWeightPath()
	//  Purpose: 
	//		Uses dynamic programming to find the highest weight path
	//		through the graph.  The essential idea is to iterate 
	//		through the vertices in a depth first order and calculate
	//		the highest weight for each vertex.  If the current vertex
	//		weight is higher than any found so far, than it becomes
	//		the end point of the path.
	//
	//		The highest weight is determined by taking the max weight
	//		of:
	//			1. the trivial path of starting at the current vertex (weight = 0)
	//		and 2. for each edge that ends at the vertex
	//				  - the weight of the edges start vertex + the edge weight
	//
	//		The algorithm changes slightly if it is constrained to start or end on 
	//		a particular node.
	//
	//			Start constrained - the trivial path is not considered unless
	//								the vertex is the start vertex
	//			End Constrained - only the end vertex can be set to the highest
	//							  weight path
	//
	//  Postconditions:
	//		- vertex objects in the vertices collection will have their weight 
	//		  and edgeForHWPAth set
	//		- highestWeightPath attribute will be set
	void findHighestWeightPath();

	// string resultString()
	//  Purpose:
	//		Returns an XML formatted string representing the results of the
	//		findHighestWeightPath() function.
	//
	//		format:
	//			<results type="part?" file=" <<graphFileName>> ">";
	//			  <result type="edge_weights"> <<weights for each edge label>> </result>
	//			  <result type="edge_histogram">  << frequencies for each edge label>> </result>
	//			  <result type="score"> <<highest weight path score>> </result>
	//			  <result type="beginning_vertex"> <<start vertex for path>> </result>
	//			  <result type="ending_vertex"> <<end vertex for path>> </result>
	//			  <result type="path"> << list of path edge labels in order>> </result>
	//			</results>
	//  Preconditions:
	//		findHighestWeightPath() has been run
	string resultString();

private:

	// Attributes
	// =============================================

	// Forward delcaration of edge so it can be used as an attribute in the
	// Vertex object
	struct Edge;

	// An object for holding informatoin related to a vertex
	struct Vertex {
		string label;  // name of vertex (must be unique)
		double weight; // highest path weight to get to the vertex
		Edge* edgeForHWPath; // edge that was used for the highest weight path
	};

	// An object for holding information related to an edge
	struct Edge {
		string label;  // name of edge
		double weight;  // cost to have path us this edge
		Vertex* start;  // start vertex of the edge
		Vertex* end;  // end vertex of the edge
	};

	string graphFileName;  // name of file defining the graph
	vector<Vertex*> vertices; // depth odering of vertices
	map<string, Vertex*> verticeMap;  // map of vertex by name, used for quick lookup of vertices
	map<string, vector<Edge*>> edges; // Collection of edges
	Vertex* startNode; // start node designated in graph file (if any)
	Vertex* endNode; // end node designated in graph file (if any)
	Vertex* highestWeightNode; // Ending node of the highest weight path
	map<string, double> edgeWeights;  // map of weights for each edge label
	map<string, int> edgeFrequencies;  // map of frequencies for each edge label

	// Private Methods
	// =============================================

	// string buildGraph()
	//  Purpose:
	//		Build the graph from the information contained in the graphFile.
	//		The graph is built by iterating through the lines of the file.
	//		If the first char is a V then a vertex is added to the vertices
	//		collection, if the the first char is an E then an edge is added
	//		to the edges collection. See class header for description of file
	//		contents. 
	//  Postconditions:
	//		The following attirbutes will be populated:
	//			vertices, verticeMap, startNode, endNode,
	//			edges, edgeWeights, edgeFrequencies
	void buildGraph();

	// addVertex(vector<string>& tokens)
	//  Purpose:
	//		Adds a new vertex from the line read in from the graph file.
	//  Postconditions:
	//		vertices - vertex added
	//		veritceMap - vertex added
	//		edges - entry for vertex initialized
	void addVertex(vector<string>& tokens);

	// addEdge(vector<string>& tokens)
	//  Purpose:
	//		Adds a new edge from the line read in from the graph file.
	//  Postconditions:
	//		edges - edge added to start and end vertex vectors
	//		edgeWeights - entry added for edge (first time label encountered)
	//		edgeFrequencies - entry for edge label incremented by 1
	void addEdge(vector<string>& tokens);

	// bool isStartConstrained()
	//  Purpose:
	//		Returns true if a start vertex is designated in the graph file
	bool isStartConstrained();

	// bool isEndConstrained()
	//  Purpose:
	//		Returns true if an end vertex is designated in the graph file
	bool isEndConstrained();

	// string getEdgeWeights()
	//  Purpose:
	//		Returns a comma delimited string describing each of the edge labels
	//		and its correpsonding weight.
	//
	//		Format:		<edge label> = <edge weight>
	string getEdgeWeights();

	// string getEdgeFrequenciess()
	//  Purpose:
	//		Returns a comma delimited string describing each of the edge labels
	//		and its correpsonding frequency (the number of edges in the graph
	//		that use that label).
	//
	//		Format:		<edge label> = <edge frequency>
	string getEdgeFrequencies();

	// string getPathStartNodeLabel()
	//  Purpose:
	//		Returns the label for the start node from the highest weight path.
	string getPathStartNodeLabel();

	// string getPath()
	//  Purpose:
	//		Returns a string representing the edge labels for the highest weight
	//		path.  The edge labels are listed from the start to the end of the path.
	string getPath();

};
#endif
