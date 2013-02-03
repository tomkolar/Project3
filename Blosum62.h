/*
 * Blosum62.h
 *
 *	This is the header file for the Blosum62 object. This object
 *  is a representation of the Blosum62 scoring matrix for RNA 
 *  sequence.  Note that this class is implemented statically, so
 *  there is no need to instantiate it.
 *
 *  Typical Use:
 *		getScore(residue1, residue2)
 *		 - returns score for aligning the two residues
 *
 *		sumOfPairsWeight(residue1, residue2, residue3)
 *		 - returns the sum of pairs weights for the three residues
 *
 *		gapCost()
 *		 - returns the cost of aligning a gap
 *
 * 	Blosum62 scoring matrix
 *	#   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V 
 *	 A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 
 *	 R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 
 *	 N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3 
 *	 D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3 
 *	 C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 
 *	 Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2 
 *	 E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2 
 *	 G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 
 *	 H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3 
 *	 I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 
 *	 L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 
 *	 K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2 
 *	 M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 
 *	 F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 
 *	 P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 
 *	 S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2 
 *	 T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 
 *	 W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 
 *	 Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 
 *	 V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 
 *
 *  Created on: 1-30-13
 *      Author: tomkolar
 */
#ifndef BLOSUM62_H
#define BLOSUM62_H

#include <string>
#include <vector>
#include <map>
using namespace std;

class Blosum62
{
public:
	// Constuctors
	// ==============================================
	Blosum62(void);

	// Destructor
	// =============================================
	~Blosum62(void);

	// Public CLass Methods
	// =============================================

	// getScore(char residue1, char residue2)
	//  Purpose: 
	//		Returns the score for aligning the two residues.
	//
	//		The score that is returned will be one of the following:
	//
	//		  - the corresponding score matrix entry if residue1 and residue2
	//			are both residues
	//		  - the gap penalty if one of residue1 and residue2 is a residue,
	//			and the other is a gap character
	//		  - 0 if both residue1 and residue2 are gap characters
    static int getScore(char residue1, char residue2);

	// gapCost()
	//  Purpose: 
	//		Returns the score for aligning a gap with a residue.
    static int gapCost();

	// sumOfPairsWeight(char residue1, char residue2, char residue3)
	//  Purpose: 
	//		Returns the sum of pairs score for aligning the three residues.
	//
	//		There are 3 different unordered pairs of the residues: 
	//		(residue1, residue2), (residue1, residue3) and (residue2, residue1).
	//		The score returned is the sum of the score for each of these
	//		unordered pairs.
	static int sumOfPairsWeight(char residue1, char residue2, char residue3);

private:

	// Class Attributes
	// =============================================
	static const char gapChar;
	static const int matrix[][20]; // scores for aligning residues
	static const map<char, int> residueMap;  // indexing for residue locations in the matrix

	// Private Class Methods
	// =============================================

	// map<char, int> createResidueMap()
	//  Purpose: 
	//		Creates a map of the index location for residues in the matrix
	//		array
	static map<char, int> createResidueMap();

	// int getIndex(char residue)
	//  Purpose: 
	//	  Returns the index in the matrix for the residue
	static int getIndex(char residue);

};

#endif // BLOSUM62_H