/*
 * WDAGraphFileBuilder.cpp
 *
 *	This is the cpp file for the WDAGraphFileBuilder object. The
 *  WDAGraphFileBuilder is a utility object that will build a file 
 *  representation of a Weighted Directed Acyclic Graph that can be
 *  read in by the WDAGraph object.
 *
 *  Currently the builder will only build a graph file for three fasta
 *  files.  Other methods may be added in the future to support other
 *  combinations of files.
 *
 *  Typical Use:
 *		buildGraphFile(fast1, fasta2, fasta3, graphFileName)
 *
 *  Created on: 1-29-13
 *      Author: tomkolar
 */

#include "WDAGraphFileBuilder.h"
#include "Blosum62.h"
#include <sstream>
#include <iostream>
#include <fstream>

// Constuctors
// ==============================================
WDAGraphFileBuilder::WDAGraphFileBuilder() {
	gapChar = '-';
}

// Destructor
// =============================================
WDAGraphFileBuilder::~WDAGraphFileBuilder() {
}

// Public Methods
// =============================================

// buildGraphFile(FastaFile*  fasta1, FastaFile* fasta2, FastaFile*  fasta3, string& graphFileName)
//  Purpose: 
//		Builds a graph file representing an edit graph for the 3 fasta files. 
//		The vertices in the graph are triples (i,j,k) where i is the starting
//		position of the sequence in fasta1, j - fasta2, and k - fasta3.  Edges
//		are labeled with the appropriate residue (or gap char) for the starting
//		position.  Edge weights are the sum of pairs for the 3 residues/gaps using
//		the Blosum62 scoring matrix.
//
//		Further details on graph construction:
//			0 <= i <= n1, 0 <= j <= n2, and 0 <= k <= n3
//			where n1, n2 and n3 are the lengths of the three sequences.
//
//			There is an edge from (i,j,k) to (i',j',k')
//			whenever i' = i or i+1, j' = j or j+1, and k' = k or k+1,
//		    and at least one of the equalities i'=i, j'=j, and k'=k is false.
//
//			The label attached to an edge is the corresponding column of aligned
//			residues & gap characters.  For example, the label associated to the
//			edge from (10,37,5) to (11,37,6) would be V-C if V is the 11th residue
//			in the first sequence and C is the 6th residue in the third sequence.
//
//		The file will be created such that all verticies are listed first
//		followed by all edges.
//
//		Vertex format: 
//
//			V <fasta1StartLoc,fasta2StartLoc,fasta3StartLoc>
//
//		Edge format:
//
//			E <fasta1Residue/gap fasta2Residue/gap fasta3Residue/gap> <start vertex id> <end vertex id> <weight>
//
//  Preconditions:
//		Fasta Files have been read and sequence has been populated
//  Postconditions:
//		File named aGraphFileName will be populated with the edit graph 
//		associated with the sequences from the fasta files.
void WDAGraphFileBuilder::buildGraphFile(FastaFile*  fasta1, FastaFile* fasta2, FastaFile*  fasta3, string& graphFileName) {

	// Create the vertices and edges
	stringstream vertices;
	stringstream edges;

	string& seq1 = fasta1->getSequence();
	string& seq2 = fasta2->getSequence();
	string& seq3 = fasta3->getSequence();

	int seq1Length = fasta1->getSequenceLength();
	int seq2Length = fasta2->getSequenceLength();
	int seq3Length = fasta3->getSequenceLength();

	for (int seq1Loc = 0; seq1Loc <= seq1Length; seq1Loc++) {
		char residue1;
		if (seq1Loc < seq1Length)
			residue1 = seq1.at(seq1Loc);
		for (int seq2Loc = 0; seq2Loc <= seq2Length; seq2Loc++) {
			char residue2;
			if (seq2Loc < seq2Length)
				residue2 = seq2.at(seq2Loc);
			for (int seq3Loc = 0; seq3Loc <= seq3Length; seq3Loc++) {
				char residue3;
				if (seq3Loc < seq3Length)
					residue3 = seq3.at(seq3Loc);

				// Add vertex info to vertices stream
				vertices << "V " << seq1Loc << "," << seq2Loc << "," << seq3Loc;

/*				// Insert start or end if needed
				if (seq1Loc == 0 && seq2Loc == 0 && seq3Loc == 0)
					vertices << " START";

				if (seq1Loc == seq1Length && seq2Loc == seq2Length && seq3Loc == seq3Length)
					vertices << " END";
*/				
				// Add new line char for vertex
				vertices << "\n";


				// Add edge info to edges stream
				// If not at end of any of the sequences then add all edges
				if (seq1Loc < seq1Length && seq2Loc < seq2Length && seq3Loc < seq3Length) {

					// Edges for single residue change
					addEdge(edges, residue1, gapChar, gapChar, seq1Loc, seq2Loc, seq3Loc);
					addEdge(edges, gapChar, residue2, gapChar, seq1Loc, seq2Loc, seq3Loc);
					addEdge(edges, gapChar, gapChar, residue3, seq1Loc, seq2Loc, seq3Loc);

					// Edges for two resiedue change
					addEdge(edges, residue1, residue2, gapChar, seq1Loc, seq2Loc, seq3Loc);
					addEdge(edges, gapChar, residue2, residue3, seq1Loc, seq2Loc, seq3Loc);
					addEdge(edges, residue1, gapChar, residue3, seq1Loc, seq2Loc, seq3Loc);

					// Edges for three residue change
					addEdge(edges, residue1, residue2, residue3, seq1Loc, seq2Loc, seq3Loc);
				}
				// End of first sequence
				else if (seq1Loc == seq1Length) {
					if (seq2Loc < seq2Length)
						addEdge(edges, gapChar, residue2, gapChar, seq1Loc, seq2Loc, seq3Loc);
					if (seq3Loc < seq3Length)
						addEdge(edges, gapChar, gapChar, residue3, seq1Loc, seq2Loc, seq3Loc);
					if (seq2Loc < seq2Length && seq3Loc < seq3Length)
						addEdge(edges, gapChar, residue2, residue3, seq1Loc, seq2Loc, seq3Loc);
				}
				// End of second sequence
				else if (seq2Loc == seq2Length) {
					if (seq1Loc < seq1Length)
						addEdge(edges, residue1, gapChar, gapChar, seq1Loc, seq2Loc, seq3Loc);
					if (seq3Loc < seq3Length)
						addEdge(edges, gapChar, gapChar, residue3, seq1Loc, seq2Loc, seq3Loc);
					if (seq1Loc < seq1Length && seq3Loc < seq3Length)
						addEdge(edges, residue1, gapChar, residue3, seq1Loc, seq2Loc, seq3Loc);
				}
				// End of third sequence
				else if (seq3Loc == seq3Length) {
					if (seq1Loc < seq1Length)
						addEdge(edges, residue1, gapChar, gapChar, seq1Loc, seq2Loc, seq3Loc);
					if (seq2Loc < seq2Length)
						addEdge(edges, gapChar, residue2, gapChar, seq1Loc, seq2Loc, seq3Loc);
					if (seq1Loc < seq1Length && seq2Loc < seq2Length)
						addEdge(edges, residue1, residue2, gapChar, seq1Loc, seq2Loc, seq3Loc);
				}

			} // seq3Loc - fasta3
		} // seq2Loc - fasta2
	} // seq1Loc - fasta1

	// Write vertices and edges to file
	ofstream graphFile(graphFileName);
	graphFile << vertices.rdbuf();
	graphFile << edges.rdbuf();

	graphFile.close();
}

// Private Methods
// =============================================

// addEdge(stringstream& edges, char residue1, char residue2, char residu3)
//          int residue1StartLoc, int residue2StartLoc, int residue3StartLoc)
//  Purpose: 
//	  Add an edge to the edges string stream.
//
//	  Format:
//
//			E <fasta1Residue/gap fasta2Residue/gap fasta3Residue/gap> <start vertex id> <end vertex id> <weight>
//
//  Postconditions:
//		edges string stream will have the edge information for the residues appended
void WDAGraphFileBuilder::addEdge(stringstream& edges, char residue1, char residue2, char residue3,
			 int residue1StartLoc, int residue2StartLoc, int residue3StartLoc) {

		// Edge Identifier
		edges << "E ";

		// Label
		edges << residue1 << residue2 << residue3 << " ";

		// Starting Vertex
		edges << residue1StartLoc << "," << residue2StartLoc << "," << residue3StartLoc << " ";

		// Ending Vertex  -- Add 1 to start location if there is an actual residue
		edges << ((residue1==gapChar)?residue1StartLoc:residue1StartLoc+1) << ",";
		edges << ((residue2==gapChar)?residue2StartLoc:residue2StartLoc+1) << ",";
		edges << ((residue3==gapChar)?residue3StartLoc:residue3StartLoc+1) << " ";

		// Weight
		edges << Blosum62::sumOfPairsWeight(residue1, residue2, residue3) << "\n";

}
