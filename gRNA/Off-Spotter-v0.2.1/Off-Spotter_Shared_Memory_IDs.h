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


// Header used by Load memory, Results and Dettach memory. 

#ifndef _Off_Spotter_Shared_Memory_IDs_h
#define _Off_Spotter_Shared_Memory_IDs_h

// Those Ids are used by the shared memory functions to identify the shared memory segments
// across the programs.

/**********************************************HUMAN********************************************************/
// HASHKEY is the ID for Table A, the index table.
#define HASHKEYHUMAN 4231
// DATAKEY is the ID for Table B, the table that holds the results.
#define DATAKEYHUMAN 8675
//SIZESKEY holds the total number of entries of Table B.
#define SIZESKEYHUMAN 9000

/**********************************************MOUSE********************************************************/
// HASHKEY is the ID for Table A, the index table.
#define HASHKEYMOUSE 4232
// DATAKEY is the ID for Table B, the table that holds the results.
#define DATAKEYMOUSE 8676
//SIZESKEY holds the total number of entries of Table B.
#define SIZESKEYMOUSE 8999

/**********************************************HUMAN2*******************************************************/
// HASHKEY is the ID for Table A, the index table.
#define HASHKEYHUMAN2 4233
// DATAKEY is the ID for Table B, the table that holds the results.
#define DATAKEYHUMAN2 8677
//SIZESKEY holds the total number of entries of Table B.
#define SIZESKEYHUMAN2 8998

/**********************************************YEAST********************************************************/
// HASHKEY is the ID for Table A, the index table.
#define HASHKEYYEAST 4234
// DATAKEY is the ID for Table B, the table that holds the results.
#define DATAKEYYEAST 8678
//SIZESKEY holds the total number of entries of Table B.
#define SIZESKEYYEAST 8997

#endif
