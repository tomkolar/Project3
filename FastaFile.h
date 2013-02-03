/*
 * FastaFile.h
 *
 *	This is the header file for the FastaFile object. The FastaFile object
 *  is a utility object designed to read in a Fasta File and keep the 
 *  information for the file in memory.
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
 * Typical use for the file would be to use the FastaFile(pathName, fileName)
 * constructor to create the object.  This will automatically open the
 * Fasta File specified by the pathName and fileName, and read its contents
 * storing them in the firstLine, and dnaSequence attributes.
 *
 * buildGraphFile(graphFileName, wieghtFileName) is a convenience method that
 * will create a sequence graph file for the dnaSequence.  See the method for
 * more details on what is created.
 *
 *  Created on: 1-10-13
 *	  Modified: 1-26-13
 *      Author: tomkolar
 */

#ifndef FASTAFILE_H
#define FASTAFILE_H

#include <string>
#include <vector>
using namespace std;

class FastaFile {

public:

	// Constuctors
	// ==============================================
	FastaFile(); 
	FastaFile(string filePath, string fileName);  
	FastaFile(string filePath, string fileName, bool isDna);  
	FastaFile(string fileName, bool isDna);  

	// Destructor
	// =============================================
	virtual ~FastaFile();

	// Public Methods
	// =============================================

	// buildGraphFile(string& graphFileName, string& weightFileName)
	//  Purpose: 
	//		Build a sequence graph file from this fasta file.  The graph defined
	//		will essentially be a linked list with vertecies specified inbetween
	//		each nucleotide in the dnaSequence and edges being the nucleotide
	//		with a weight as specified in the weightFileName.
	//
	//		The file will be created such that all verticies are listed first
	//		followed by all edges.
	//
	//		Vertex format: 
	//
	//			V <sequential number as identifier>
	//
	//		Edge format:
	//
	//			E <nucleotide> <start vertex id> <end vertex id> <weight>
	//
	//  Preconditions:
	//		Fasta File has been read and dnaSequence has been populated
	//  Postconditions:
	//		File named aGraphFileName will be populated with the sequence graph 
	//		associated with the dnaSequence from the fasta file.
	void buildGraphFile(string& graphFileName, string& weightFileName);

	// string firstLineResultString()
	//  Purpose:
	//		Returns the string value of an XML element representing the first line of 
	//		the Fasta file.
	//
	//		format:
	//			<result type='first line' file='<<fileName>>' >
	//				<<firstLine>>
	//			</result>
	//  Preconditions:
	//		Fasta File has been read and firstLine has been populated
    string firstLineResultString();

	// string baseCountsResultString()
	//  Purpose:
	//		Returns the string value of an XML element representing the base counts 
	//		of the dnaSequence.
	//
	//		format:
	//			<result type='nucleotide histogram' file='<<fileName>>' >
	//				A=<<baseCountForA>>,C=<<baseCountForC>>,G=<<baseCountForG>>,
	//				A=<<baseCountForT>>,N=<<countForOtherChars>>
	//			</result>
	//  Preconditions:
	//		Fasta File has been read and dnaSequence has been populated
    string baseCountsResultString();

	// bool isDNA()
	//  Purpose:
	//		Returns true if the sequence is a DNA sequence
    bool isDNA();

	// Public Accessors
	// =============================================
	const int getSequenceLength();  // length of dnaSequence
	string& getFileName();
	string& getSequence();

private:
	// Attributes
	// =============================================
	string filePath;
    string fileName;
    string firstLine;
    string sequence;
	string reverseComplement;
	bool dna; // set to true if the sequence is a dna sequence

	// Private Methods
	// =============================================

	// populate()
	//  Purpose:
	//		Reads in the Fasta File specified by filePath and fileName and populates
	//		the object with its contents
	//	Preconditions:
	//		fileName and filePath have been set
	//  Postconditions:
	//		firstLine - populated with first line from file
	//		dnaSequence - populated with dnaSequence from file
	//		reverseComplement - populated with reverse complement of dnaSequence
    void populate();

	// createReverseComplment()
	//  Purpose:
	//		populates the reverseComplement attribute with the reverse comlpement
	//		of the dnaSequence
	//	Preconditions:
	//		dnaSequence has been set
	//  Postconditions:
	//		reverseComplement - populated with reverse complement of dnaSequence
    void createReverseComplement();

	// char complement(char aChar)
	//  Purpose:  returns the dna complement of aChar
    char complement(char aChar);

	// countBases(int counts[])
	//  Purpose:
	//		populates the counts array with the counts for base occurrences
	//		in dnaSequence.  The array is populated with the folllowing 
	//		scheme:
	//			counts[0] = counts for A
	//			counts[1] = counts for C
	//			counts[2] = counts for G
	//			counts[3] = counts for T
	//			counts[4] = counts for other characters encountered
	void countBases(int counts[]);

};

#endif /* FASTAFILE_H */