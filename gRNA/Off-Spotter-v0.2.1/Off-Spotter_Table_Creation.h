/***********************************************************************************************************/
/*   This code (© 2015 Thomas Jefferson University, All Rights Reserved) was created by Venetia Pliatsika  */
/*   and Isidore Rigoutsos and is an implementation of the Off-Spotter algorithm that appears in           */
/*   Pliatsika, V, and Rigoutsos, I (2015) "Off-Spotter: very fast and exhaustive enumeration of genomic   */
/*   lookalikes for designing CRISPR/Cas guide RNAs" Biol. Direct 10(1):4.                                 */
/*                                                                                                         */
/* Use of these codes is bound by the following terms and conditions:                                      */
/*                                                                                                         */
/* Terms of Use: This code can be freely used for research, academic and other non-profit activities.      */
/* Only one instance of the code may be used at a time, and then for only one concurrent user. You may not */
/* use the code to conduct any type of application service, service bureau or time-sharing operation or to */
/* provide any remote processing, network processing, network telecommunications or similar services to    */
/* any person, entity or organization, whether on a fee basis or otherwise. The code can be copied and     */
/* compiled on any platform for the use authorized by these terms and conditions. All copies of the code   */
/* must be accompanied by this note. The code cannot be modified without the written permission of the     */
/* Computational Medicine Center of Thomas Jefferson University https://cm.jefferson.edu                   */
/*                                                                                                         */
/* Commercial use is strictly prohibited.  If you wish to use these codes commercially please contact the  */
/* Computational Medicine Center of Thomas Jefferson University: https://cm.jefferson.edu/contact-us/      */
/*                                                                                                         */
/*                                                                                                         */
/*     THE CODE IS PROVIDED “AS IS” WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESSED    */
/*     OR IMPLIED. TO THE FULLEST EXTENT PERMISSIBLE PURSUANT TO APPLICABLE LAW. THOMAS JEFFERSON          */
/*     UNIVERSITY, AND ITS AFFILIATES, DISCLAIM ALL WARRANTIES, EXPRESS OR IMPLIED, INCLUDING, BUT NOT     */
/*     LIMITED TO, THE IMPLIED WARRANTIES OF TITLE, MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND  */
/*     NON-INFRINGEMENT.                                                                                   */
/*                                                                                                         */
/*     NEITHER THOMAS JEFFERSON UNIVERSITY NOR ITS AFFILIATES MAKE ANY REPRESENTATION AS TO THE RESULTS    */
/*     TO BE OBTAINED FROM USE OF THE CODE.                                                                */
/***********************************************************************************************************/

#ifndef _Off_Spotter_filecreation_h
#define _Off_Spotter_filecreation_h

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>

#include "Off-Spotter_General.h"

using namespace std;

//FUNCTIONS
//To create binary files that hold table A and table B we need to first count all hits, then allocate the space, populate the tables and then
//write each on a separate file. Table A is the index file and table B the data.
//The files are created for convinience. If you are working with the same genome as the one used in the file creation last time you don't
//have to re-use this function, go to loading directly.
//If you want to add more or some of the PAMs, you need to create more such functions and then edit the main of Off-Spotter_Table_Creation.cpp
//as well as the Off-Spotter_Results files. Keep in mind that the total number of hits of all PAMs together should fit in Table B.
//Those functions count the number of hits on the input genome of a specific PAM for ALL possible gRNAs.
unsigned int countNNG(const char *filename, char PAM, unsigned int *currentPos);
int countNNNNACA(const char *filename, unsigned int *currentPos);
int countNNGRRT(const char *filename, unsigned int *currentPos);
//Those functions store the details of each hit on the input genome of a specific PAM and ALL possible gRNAs
void resultsNNG(const char *filename, char PAM, unsigned int *currentPos, MYTYPE *tableValues, unsigned int *hashTable, unsigned int PAMID);
void resultsNNNNACA(const char *filename, unsigned int *currentPos, MYTYPE *tableValues, unsigned int *hashTable, unsigned int PAMID);
void resultsNNGRRT(const char *filename, unsigned int *currentPos, MYTYPE *tableValues, unsigned int *hashTable, unsigned int PAMID);

#endif
