/***********************************************************************************************************/
/*   This code (© 2015 Thomas Jefferson University, All Rights Reserved) was created by Venetia Pliatsika  */
/*   and Isidore Rigoutsos and is an implementation of the Off-Spotter algorithm that appears in       */
/*   Pliatsika, V, and Rigoutsos, I (2015) "Off-Spotter: very fast and exhaustive enumeration of genomic   */
/*   lookalikes for designing CRISPR/Cas guide RNAs" Biol. Direct 10(1):4.     */
/*     */
/* Use of these codes is bound by the following terms and conditions:      */
/*     */
/* Terms of Use: This code can be freely used for research, academic and other non-profit activities.      */
/* Only one instance of the code may be used at a time, and then for only one concurrent user. You may not */
/* use the code to conduct any type of application service, service bureau or time-sharing operation or to */
/* provide any remote processing, network processing, network telecommunications or similar services to    */
/* any person, entity or organization, whether on a fee basis or otherwise. The code can be copied and     */
/* compiled on any platform for the use authorized by these terms and conditions. All copies of the code   */
/* must be accompanied by this note. The code cannot be modified without the written permission of the     */
/* Computational Medicine Center of Thomas Jefferson University https://cm.jefferson.edu       */
/*     */
/* Commercial use is strictly prohibited.  If you wish to use these codes commercially please contact the  */
/* Computational Medicine Center of Thomas Jefferson University: https://cm.jefferson.edu/contact-us/      */
/*     */
/*     */
/*     THE CODE IS PROVIDED “AS IS” WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESSED    */
/*     OR IMPLIED. TO THE FULLEST EXTENT PERMISSIBLE PURSUANT TO APPLICABLE LAW. THOMAS JEFFERSON      */
/*     UNIVERSITY, AND ITS AFFILIATES, DISCLAIM ALL WARRANTIES, EXPRESS OR IMPLIED, INCLUDING, BUT NOT     */
/*     LIMITED TO, THE IMPLIED WARRANTIES OF TITLE, MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND  */
/*     NON-INFRINGEMENT.       */
/*     */
/*     NEITHER THOMAS JEFFERSON UNIVERSITY NOR ITS AFFILIATES MAKE ANY REPRESENTATION AS TO THE RESULTS    */
/*     TO BE OBTAINED FROM USE OF THE CODE.    */
/***********************************************************************************************************/

#include <sys/ipc.h>
#include <sys/shm.h>
#include <cstring>
#include <iostream>
#include <cstdlib>
#include <fstream>

#include "Off-Spotter_Shared_Memory_IDs.h"
#include "Off-Spotter_General.h"

using namespace std;

//This program should be used only when using shared memory
//This program detaches the 3 shared memory segment with the IDs
//found in Off-Spotter_Shared_Memory_IDs.h.
//Remember to use this program to dettach those segments when you
//don't use Off-Spotter.
int main (int argc, char *argv[]){
    char *usage = (char *) "\n\nUsage: %s -g <genome_name>\n\tPlease type hg38 or GRCh38 for human hg38\n\tPlease type hg19 or GRCh37 for human hg19\n\tPlease type mm10 or GRCm38 for mouse mm10\n\tPlease type w303 for yeast strain w303\n\n";
    //Parse input
    int a;
    char *genome = NULL;
    while ((a = getopt (argc, argv, "g:")) != -1){
        if(a == 'g'){
            genome=optarg;
            //Set the memory keys depending on the genome
            if (strcmp(genome,"hg19")==0 || strcmp(genome,"GRCh37")==0){
                HASHKEY=HASHKEYHUMAN;
                DATAKEY=DATAKEYHUMAN;
                SIZESKEY=SIZESKEYHUMAN;
            }
            else if (strcmp(genome,"mm10")==0 || strcmp(genome,"GRCm38")==0){
                HASHKEY=HASHKEYMOUSE;
                DATAKEY=DATAKEYMOUSE;
                SIZESKEY=SIZESKEYMOUSE;
            }
            else if (strcmp(genome,"hg38")==0 || strcmp(genome,"GRCh38")==0){
                HASHKEY=HASHKEYHUMAN2;
                DATAKEY=DATAKEYHUMAN2;
                SIZESKEY=SIZESKEYHUMAN2;
            }
            else if (strcmp(genome,"w303")==0){
                HASHKEY=HASHKEYYEAST;
                DATAKEY=DATAKEYYEAST;
                SIZESKEY=SIZESKEYYEAST;
            }
            else{
                fprintf(stderr, "\nUknown <genome_name> %s\n", genome);
                fprintf(stderr, usage, argv[0]);
                fflush(stderr);
                return(-1);
            }
        }
        else if (a == '?'){
            if (optopt == 'c')
                fprintf (stderr, "\nOption -%c requires an argument\n", optopt);
            else if (isprint (optopt))
                fprintf (stderr, "\nUnknown option `-%c'\n", optopt);
            else
                fprintf (stderr, "\nUnknown option character `\\x%x'\n", optopt);
            fprintf(stderr, usage, argv[0]);
            fflush(stderr);
            return -1;
        }
        else{
            fprintf(stderr, usage, argv[0]);
            fflush(stderr);
            return(1);
        }
    }
    if (genome == NULL){
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        return(-1);
    }

    //Get the size of table B
    int sizes_ID;
    if ((sizes_ID = shmget(SIZESKEY,sizeof(unsigned int), 0644))==-1){
        fprintf(stderr, "\nAn error occured please make sure that the tables for %s are loaded on memory\n",genome);
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        return(-1);
    }
    unsigned int *sizes;
    if ((sizes = (unsigned int *) shmat(sizes_ID,NULL,0))==  (void *)-1){
        fprintf(stderr, "\nAn error occured please make sure that the tables for %s are loaded on memory\n",genome);
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        return(-1);
    }
    
    //Dettach table A
    cout << "\nDetaching Table A" << endl;
    int ID;
    if (( ID = shmget(HASHKEY,KEY_NUM*sizeof(unsigned int), 0644))==-1){
        fprintf(stderr, "\nAn error occured please make sure that the tables for %s are loaded on memory\n",genome);
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        return(-1);
    }
    if (shmctl(ID, IPC_RMID, 0) == -1){
        fprintf(stderr, "\nAn error occured while deleting memory segment of %s\n",genome);
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        return(-1);
    }
    
    //Dettach table B
    cout << "\nDetaching Table B" << endl;
    if ((ID = shmget(DATAKEY,sizes[0]*sizeof(MYTYPE),0644))==-1){
        fprintf(stderr, "\nAn error occured please make sure that the tables for %s are loaded on memory\n",genome);
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        return(-1);
    }
    if (shmctl(ID, IPC_RMID, 0) == -1){
        fprintf(stderr, "\nAn error occured while deleting memory segment of %s\n",genome);
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        return(-1);
    }

    //Dettach the size variable
    if (shmctl(sizes_ID, IPC_RMID, 0)){
        fprintf(stderr, "\nAn error occured while deleting memory segment of %s\n",genome);
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        return(-1);
    }
    
    cout << "\nDone!" << endl;
    
    return 1;
}
