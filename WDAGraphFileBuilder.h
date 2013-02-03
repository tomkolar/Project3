/*
 * WDAGraphFileBuilder.h
 *
 *	This is the header file for the WDAGraphFileBuilder object. The
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
 *
 *  Created on: 1-29-13
 *      Author: tomkolar
 */

#ifndef WDAGRAPHFILEBUILDER_H
#define WDAGRAPHFILEBUILDER_H

#include "FastaFile.h"
#include <string>
#include <vector>
using namespace std;

class WDAGraphFileBuilder
{
public:

	// Constuctors
	// ==============================================
	WDAGraphFileBuilder();

	// Destructor
	// =============================================
	virtual ~WDAGraphFileBuilder();

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
	void buildGraphFile(FastaFile*  fasta1, FastaFile* fasta2, FastaFile*  fasta3, string& graphFileName);

private:

	// Attributes
	// =============================================
	char gapChar;

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
	void addEdge(stringstream& edges, char residue1, char residue2, char residu3,
		int residue1StartLoc, int residue2StartLoc, int residue3StartLoc);
};

#endif // WDAGRAPHFILEBUILDER_H

