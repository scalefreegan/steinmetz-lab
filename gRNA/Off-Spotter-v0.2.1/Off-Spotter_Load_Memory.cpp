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

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <cerrno>

#include "Off-Spotter_Shared_Memory_IDs.h"
#include "Off-Spotter_General.h"

using namespace std;

//Table A, index table
unsigned int *hashTable;
//Table B, result table
MYTYPE *tableValues;
//Size of table B (the size of A is fixed and in Off-Spotter_General.h)
unsigned int *sizes;


int readTables (char *FILEINDEX, char *FILEDATA, char *genome, char *usage, char *prname){
    FILE *f;
    struct stat fileStat;
    
    //Allocate and attach shared memory for the size of table B
    int sizes_ID;
    if((sizes_ID = shmget(SIZESKEY,sizeof(unsigned int),IPC_CREAT | 0644))==-1){
        fprintf(stderr, "\nAn error occured please make sure that the tables for %s aren't already loaded on memory\n",genome);
        fprintf(stderr, usage, prname);
        fflush(stderr);
        return(-1);
    }
    if ((sizes = (unsigned int *) shmat(sizes_ID,NULL,0))==  (void *) -1){
        fprintf(stderr, "\nAn error occured please make sure that the tables for %s aren't already loaded on memory\n",genome);
        fprintf(stderr, usage, prname);
        fflush(stderr);
        return(-1);
    }
    
    //open file with table A
    cout << "\nLoading Table A, the index table. Please wait." << endl;
    f = fopen ( FILEINDEX , "rb" );
    if (f == NULL){
        fprintf(stderr, "\nAn error occured while opening index.bin of %s\n",genome);
        fprintf(stderr, usage, prname);
        fflush(stderr);
        return(-1);
    }
    //allocate shared memory to contain the read data
    int hashID;
    if(( hashID = shmget(HASHKEY,KEY_NUM*sizeof(unsigned int),IPC_CREAT | 0644))==-1){
        fprintf(stderr, "\nAn error occured please make sure that the tables for %s aren't already loaded on memory\n",genome);
        fprintf(stderr, usage, prname);
        fflush(stderr);
        return(-1);
    }
    //Attach shared memory
    if ((hashTable = (unsigned int *) shmat(hashID,NULL,0))==  (void *) -1){
        fprintf(stderr, "\nAn error occured please make sure that the tables for %s aren't already loaded on memory\n",genome);
        fprintf(stderr, usage, prname);
        fflush(stderr);
        return(-1);
    }
    //read index and table
    if (fread(hashTable,sizeof(unsigned int),KEY_NUM,f) == 0){
        fprintf(stderr, "\nAn error occured while reading index.bin of %s\n",genome);
        fprintf(stderr, usage, prname);
        fflush(stderr);
        return(-1);
    }
    fclose(f);
    
    //open file with table B
    cout << "\nLoading Table B, the data table. Please wait." << endl;
    f = fopen (FILEDATA , "rb" );
    if (f == NULL){
        fprintf(stderr, "\nAn error occured while opening data.bin of %s\n",genome);
        fprintf(stderr, usage, prname);
        fflush(stderr);
        return(-1);
    }
    //obtain file size and assign it to the size shared memory segment
    stat(FILEDATA,&fileStat);
    sizes[0] = fileStat.st_size/sizeof(MYTYPE);
    //allocate space in memory for table B
    int dataID;
    if ((dataID = shmget(DATAKEY,sizes[0]*sizeof(MYTYPE),IPC_CREAT | 0644))==-1){
        fprintf(stderr, "\nAn error occured please make sure that the tables for %s aren't already loaded on memory\n",genome);
        fprintf(stderr, usage, prname);
        fflush(stderr);
        return(-1);
    }
    //attach table B
    if ((tableValues = (MYTYPE *) shmat(dataID,NULL,0))== (void *) -1){
        fprintf(stderr, "\nAn error occured please make sure that the tables for %s aren't already loaded on memory\n",genome);
        fprintf(stderr, usage, prname);
        fflush(stderr);
        return(-1);
    }
    //read table B on memory and close file
    if (fread(tableValues,sizeof(MYTYPE),sizes[0],f)==0){
        fprintf(stderr, "\nAn error occured while reading data.bin of %s\n",genome);
        fprintf(stderr, usage, prname);
        fflush(stderr);
        return(-1);
    }
    fclose(f);
    return 1;
}

//This program should be used only when using shared memory
//It loads the shared memory segment of 1 genome, the one defined by the input parameters.
int main (int argc, char *argv[]){
    char usage [256] = "\n\nUsage: %s -t <input_folder_full_path> -g <genome_name>\n\tPlease type hg38 or GRCh38 for human hg38\n\tPlease type hg19 or GRCh37 for human hg19\n\tPlease type mm10 or GRCm38 for mouse mm10\n\tPlease type w303 for yeast strain w303\n\n";
    
    //Parse input
    int a;
    char *genome = NULL;
    char *in_file_path = NULL;
    while ((a = getopt (argc, argv, "g:t:")) != -1){
        if(a == 'g'){
            genome=optarg;
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
        else if (a == 't'){
            in_file_path=optarg;
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
            exit(1);
        }
    }

    //Input must include the path to the directory where index.bin and data.bin are located
    //(both must be on the same directory) and the name of the genome whose tables we want loaded
    if (genome == NULL || in_file_path==NULL) {
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        exit(1);
    }

    //Create the input filenames
    char index_file[strlen(in_file_path) + strlen("/index.bin") + 1], data_file[strlen(in_file_path) + strlen("/data.bin") + 1];
    strcpy(index_file,in_file_path);
    strcpy(data_file,in_file_path);
    if (in_file_path[strlen(in_file_path) - 1] == '/'){
        strcat(index_file,"index.bin");
        strcat(data_file,"data.bin");
    }
    else{
        strcat(index_file,"/index.bin");
        strcat(data_file,"/data.bin");
    }

    //Load tha tables (files) in memory
    if (readTables(index_file,data_file,genome,usage,argv[0]) < 0){
        return -1;
    }
    
    //detach
    if (shmdt(sizes) ==-1){
        fprintf(stderr, "\nAn error occured while detaching the size of the results segment\n");
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        return(-1);
    }
    if (shmdt(hashTable) == -1){
        fprintf(stderr, "\nAn error occured while detaching the size of the results segment\n");
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        return(-1);
    }
    if (shmdt(tableValues) == -1){
        fprintf(stderr, "\nAn error occured while detaching the size of the results segment\n");
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        return(-1);
    }
    
    cout << "\nTables loaded. You can use \"Results\" now.\n" << endl;
    
    return 0;
}
