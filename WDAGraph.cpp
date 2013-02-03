/*
 * WDAGraph.cpp
 *
 *	This is the cpp file for the WDAGraph object. WDAGraph is an implementation
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

#include "WDAGraph.h"
#include "StringUtilities.h"
#include <limits>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cstddef>
#include <sstream>
using namespace std;


// Constuctors
// ==============================================
WDAGraph::WDAGraph() {
	
}

WDAGraph::WDAGraph(string& aGraphFileName) {

	// Initialize pointers
	startNode = NULL;
	endNode = NULL;
	highestWeightNode = NULL;

	//  Set file name
	graphFileName = aGraphFileName;
	
	// Build the graph
	buildGraph();
}

// Destructor
// =============================================
WDAGraph::~WDAGraph(){
	// Remove start,end and highest
	startNode = NULL;
	endNode = NULL;
	highestWeightNode = NULL;

	// Remove edges from vertex
	for (Vertex* vertex : vertices) {
		vertex->edgeForHWPath = NULL;
	}

	// Remove vertices from edges
	for (auto theEdges : edges) {
		vector<Edge*>& vertexEdges = theEdges.second;
		for (Edge* edge : vertexEdges) {
			edge->start = NULL;
			edge->end = NULL;
		}
	}

	// Remove vertices
	for (Vertex* vertex : vertices) {
		vertex = NULL;
	}
	for (auto vertexEntry : verticeMap) {
		vertexEntry.second = NULL;
	}

	// Remove  edges
	for (auto theEdges : edges) {
		vector<Edge*>& vertexEdges = theEdges.second;
		for (Edge* edge : vertexEdges) {
			edge = NULL;
		}
	}


}

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
void WDAGraph::findHighestWeightPath() {

	bool startFound = false;

	// Iterate through the vertices
	for (Vertex* vertex : vertices) {
		// Check to for start constraints
		if (isStartConstrained()) {
			if (!startFound ) {
				if (vertex->label == startNode->label) {
					startFound = true;
					vertex->weight = 0;
					highestWeightNode = vertex;
				}
				else // start not found yet
					continue;
			}
		}
		else {
			// Start not constrained - so consider the trivial path of starting here
			vertex->weight = 0;
		}


		// Find the path with the highest weight to this vertex
		vector<Edge*> vertexEdges = edges.find(vertex->label)->second;
		for (Edge* edge : vertexEdges) {
			// Find edges where this vertex is the end node
			if (edge->end->label == vertex->label) {
				// If start constrained make sure edge start is from a valid path
				if  (isStartConstrained()) {
					if (edge->start->weight == INT_MIN)
						continue;
				}
				// Calculate weight (parent node weight plus edge weight)
				double pathWeight = edge->start->weight + edge->weight;

				// If path weight bigger than any other found so far then
				// it becomes the new weight for the vertex
				if (pathWeight > vertex->weight) {
					vertex->weight = pathWeight;
					vertex->edgeForHWPath = edge;
				}
			}
		}

		// Check for end constraint
		if (isEndConstrained()) {
			// Is this the end node?
			if (vertex->label == endNode->label) {
				// Found end node, so set it as highest and exit vertex for loop
				highestWeightNode = vertex;
				break;
			}
			// Continue looking for the end vertex as the highest weight path must end there
			continue;
		}

		// Set highestWeightNode (for non end constrained paths)
		if (highestWeightNode == NULL || vertex->weight > highestWeightNode->weight)
			highestWeightNode = vertex;

	}  // end vertices for loop
	
}

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
string WDAGraph::resultString() {
	stringstream ss;
	// Results header
	ss << "  <results type=\"part?\" file=\"" << graphFileName << "\">\n";

	// Edge Info (Weights and Histogram)
	ss << StringUtilities::xmlResult("edge_weights", getEdgeWeights());
	ss << StringUtilities::xmlResult("edge_histogram", getEdgeFrequencies());

	// Path Info
	if (highestWeightNode == NULL)
		ss << StringUtilities::xmlResult("path", "No Path Found!");
	else {
		ss
			<< StringUtilities::xmlResult("score",  highestWeightNode->weight, 6)
			<< StringUtilities::xmlResult("beginning_vertex",  getPathStartNodeLabel())
			<< StringUtilities::xmlResult("end_vertex", highestWeightNode->label)
			<< StringUtilities::xmlResult("path", getPath());
	}

	// Results footer
	ss << "  </results>\n";


	return ss.str();
}

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
void WDAGraph::buildGraph() {

	ifstream graphFile(graphFileName);
	string line;

	while(getline(graphFile, line)) {
		vector<string> tokens;
		StringUtilities::split(line, ' ', tokens);

		// Add vertices
		if (tokens.at(0) == "V") 
			addVertex(tokens);
		// Add Edges
		else if (tokens.at(0) == "E") 
			addEdge(tokens);
	}

	graphFile.close();

}

// addVertex(vector<string>& tokens)
//  Purpose:
//		Adds a new vertex from the line read in from the graph file.
//  Postconditions:
//		vertices - vertex added
//		veritceMap - vertex added
//		edges - entry for vertex initialized
void WDAGraph::addVertex(vector<string>& tokens) {

	// Create vetex and intialize 
	Vertex* vertex =  new Vertex();
	vertex->weight = INT_MIN;
	vertex->edgeForHWPath = NULL;

	// Populate from tokens
	vertex->label = tokens.at(1);
	if (tokens.size() > 2) {
		if (tokens.at(2) == "START")
			startNode = vertex;
		else if (tokens.at(2) == "END")
			endNode = vertex;
	}

	// Add to collections
	vertices.push_back(vertex);
	verticeMap[vertex->label] = vertex;

	// initialize edges entry for this vertex
	edges[vertex->label] = vector<Edge*>();

}

// addEdge(vector<string>& tokens)
//  Purpose:
//		Adds a new edge from the line read in from the graph file.
//  Postconditions:
//		edges - edge added to start and end vertex vectors
//		edgeWeights - entry added for edge (first time label encountered)
//		edgeFrequencies - entry for edge label incremented by 1
void WDAGraph::addEdge(vector<string>& tokens) {

	// Create Edge and populate from tokens
	Edge* edge = new Edge();
	edge->label = tokens.at(1);
	edge->start = verticeMap.find(tokens.at(2))->second;
	edge->end = verticeMap.find(tokens.at(3))->second;
	edge->weight = atof(tokens.at(4).c_str());


	// Add to edges collection
	vector<Edge*>& startEdges = edges.find((edge->start)->label)->second;
	startEdges.push_back(edge);
	vector<Edge*>& endEdges = edges.find((edge->end)->label)->second;
	endEdges.push_back(edge);

	// Add to edgeWeights if not already in map ]
	map<string, double>::iterator weightsIter = edgeWeights.find(edge->label);
	if (weightsIter == edgeWeights.end()) 
		edgeWeights[edge->label] = edge->weight;

	// Increment edgeFrequencies
	map<string, int>::iterator freqIter = edgeFrequencies.find(edge->label);
	if (freqIter != edgeFrequencies.end()) 
		edgeFrequencies[edge->label] = freqIter->second++;
	else
		edgeFrequencies[edge->label] = 1;

}

// bool isStartConstrained()
//  Purpose:
//		Returns true if a start vertex is designated in the graph file
bool WDAGraph::isStartConstrained() {
	return startNode != NULL;
}

// bool isEndConstrained()
//  Purpose:
//		Returns true if an end vertex is designated in the graph file
bool WDAGraph::isEndConstrained() {
	return endNode != NULL;
}

// string getEdgeWeights()
//  Purpose:
//		Returns a comma delimited string describing each of the edge labels
//		and its correpsonding weight.
//
//		Format:		<edge label> = <edge weight>
string WDAGraph::getEdgeWeights() {
	stringstream ss;
	ss.precision(3);

	// Iterate throght edge weights map
	for (auto edgeWeight : edgeWeights) {
		ss << edgeWeight.first << "=" << edgeWeight.second << ", ";
	}

	// Return all but last 2 char (as there will be an extra ,<space> at end)
	string temp = ss.str();
	return temp.substr(0, temp.length() -2);
}

// string getEdgeFrequenciess()
//  Purpose:
//		Returns a comma delimited string describing each of the edge labels
//		and its correpsonding frequency (the number of edges in the graph
//		that use that label).
//
//		Format:		<edge label> = <edge frequency>
string WDAGraph::getEdgeFrequencies() {
	stringstream ss;

	// Iterate throght edge weights map
	for (auto edgeFrequency : edgeFrequencies) {
		ss << edgeFrequency.first << "=" << edgeFrequency.second << ", ";
	}

	// Return all but last 2 char (as there will be an extra ,<space> at end)
	string temp = ss.str();
	return temp.substr(0, temp.length() -2);
}

// string getPathStartNodeLabel()
//  Purpose:
//		Returns the label for the start node from the highest weight path.
string WDAGraph::getPathStartNodeLabel() {

	// Base Case
	if (highestWeightNode == NULL)
		return "";

	// Search for start node (previous is NULL)
	Vertex* aNode = highestWeightNode;
	while (aNode->edgeForHWPath != NULL) {
		// Walk backwards until find the start node
		Edge* anEdge = aNode->edgeForHWPath;
		aNode = anEdge->start;
	}

	return aNode->label;
}

// string getPath()
//  Purpose:
//		Returns a string representing the edge labels for the highest weight
//		path.  The edge labels are listed from the start to the end of the path.
string WDAGraph::getPath() {

	// Base Case
	if (highestWeightNode == NULL)
		return "";

	stringstream ss;

	// Walk path backwards and build string
	Vertex* aNode = highestWeightNode;
	while (aNode->edgeForHWPath != NULL) {
		// Add the label from the edge to the string stream
		Edge* anEdge = aNode->edgeForHWPath;
		ss << anEdge->label << "\n";

		// Walk backwards one node
		aNode = anEdge->start;
	}

	// Return reverse of stringstream (since we built the string backwards)
	string reversePath = ss.str();
	return string(reversePath.rbegin(), reversePath.rend());

}
