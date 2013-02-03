/*
 * driver.cpp
 *
 *	This is the driver file for creating an alignment for 3 fasta
 *  files using a weighted directed acyclic edit graph.
 *
 *	Typical use:
 *		align fastaFile1 fastaFile2 fastaFile3
 *
 *  Created on: 1-29-13
 *      Author: tomkolar
 */
#include "FastaFile.h"
#include "WDAGraph.h"
#include "WDAGraphFileBuilder.h"
#include <string>
#include <sstream>
#include <iostream>
using namespace std;

int main( int argc, char *argv[] ) {

	// Check that file name was  entered as argument
	if (argc < 4) {
		cout << "Invalid # of arguments\n";
		cout << "usage: align fastaFile1 fastaFile2 fastaFile3\n";
		return -1;
	}

	cout << "Starting\n";

	// Get Fasta File names
	string fastaFileName1 = argv[1];
	string fastaFileName2 = argv[2];
	string fastaFileName3 = argv[3];
/*
	// Set Fasta File names
	string fastaFileName1 = "rna_test1.fna";
	string fastaFileName2 = "rna_test2.fna";
	string fastaFileName3 = "rna_test3.fna";

	// Set Fasta File names
	string fastaFileName1 = "testfile1.fna";
	string fastaFileName2 = "testfile2.fna";
	string fastaFileName3 = "testfile3.fna";
*/

	// Create graph file name
	stringstream ss;
	ss 
		<< fastaFileName1 << "_"
		<< fastaFileName2 << "_"
		<< fastaFileName3 << ".graph.txt";

	string graphFileName = ss.str();

	// Create the fasta file objects
	FastaFile* fastaFile1 = new FastaFile(fastaFileName1, false);
	FastaFile* fastaFile2 = new FastaFile(fastaFileName2, false);
	FastaFile* fastaFile3 = new FastaFile(fastaFileName3, false);

	cout << "Fasta's done\n";

	// Create the graph file
	WDAGraphFileBuilder builder;
	builder.buildGraphFile(fastaFile1, fastaFile2, fastaFile3, graphFileName);

	cout << "Graph File built\n";

	// Create the WDAGraph and find the highest weight path
	WDAGraph* aGraph =  new WDAGraph(graphFileName);

	cout << "Graph built\n";

	aGraph->findHighestWeightPath();

	// Print out the result string for the highest weight path
	cout << aGraph->resultString();

	aGraph = NULL;
	fastaFile1 = NULL;
	fastaFile2 = NULL;
	fastaFile3 = NULL;
}
