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
#include <sys/mman.h>
#include <fcntl.h>
#include <sys/stat.h>

#include "Off-Spotter_Results.h"

using namespace std;

#define PATH_LEN 256

//Hash table address
unsigned int *hashTableADDR;
//Result table address
unsigned int *resultADDR;
//Size of table B (the size of A is fixed and is located in Off-Spotter_General.h)
unsigned int *sizes;

/***********************************************************************************************************/
/************************************** M E M O R Y  H A N D L I N G ***************************************/
/***********************************************************************************************************/
/* You can choose among shared memory and memory shared files. If using shared memory, use Load_Memory and */
/* detach memory too.                                                                                      */
/***********************************************************************************************************/

/***********************************************************************************************************/
/* This function attaches this program to the (already loaded by Load_Memory) shared memory                */
/***********************************************************************************************************/
int attach_shared_memory(char *genome,char *usage,char *name){
    int hash_ID,data_ID;
    
    //locate sizes table
    int sizes_ID;
    if ((sizes_ID = shmget(SIZESKEY,sizeof(unsigned int), 0644))== -1){
        fprintf(stderr, "\nAn error occured while attaching to shared memory for %s. Please make sure that shared memory is loaded properly.",genome);
        fprintf(stderr, usage, name);
        fflush(stderr);
        return(-1);
    }
    
    //attach size table
    unsigned int *sizes;
    if ((sizes = (unsigned int *) shmat(sizes_ID,NULL,0)) == (void *)-1){
        fprintf(stderr, "\nAn error occured while attaching to shared memory for %s. Please make sure that shared memory is loaded properly.",genome);
        fprintf(stderr, usage, name);
        fflush(stderr);
        return(-1);
    }
    
    //find the correct tables
    total_len = sizes[0];
    if ((hash_ID = shmget(HASHKEY,KEY_NUM*sizeof(unsigned int),0644)) ==-1 ){
        fprintf(stderr, "\nAn error occured while attaching to shared memory for %s. Please make sure that shared memory is loaded properly.",genome);
        fprintf(stderr, usage, name);
        fflush(stderr);
        return(-1);
    }
    if ((data_ID = shmget(DATAKEY,total_len*sizeof(MYTYPE), 0644)) == -1){
        fprintf(stderr, "\nAn error occured while attaching to shared memory for %s. Please make sure that shared memory is loaded properly.",genome);
        fprintf(stderr, usage, name);
        fflush(stderr);
        return(-1);
    }
    
    //release size table
    if (shmdt(sizes) == -1){
        fprintf(stderr, "\nAn error occured while attaching to shared memory for %s. Please make sure that shared memory is loaded properly.",genome);
        fprintf(stderr, usage, name);
        fflush(stderr);
        return(-1);
    }
    
    //attach H A S H
    if ((hashTable = (unsigned int *) shmat(hash_ID,NULL,0))==(void*)-1){
        fprintf(stderr, "\nAn error occured while attaching to shared memory for %s. Please make sure that shared memory is loaded properly.",genome);
        fprintf(stderr, usage, name);
        fflush(stderr);
        return(-1);
    }
    //attach D A T A
    if((tableValues = (MYTYPE *) shmat(data_ID,NULL,0)) == (void*)-1){
        fprintf(stderr, "\nAn error occured while attaching to shared memory for %s. Please make sure that shared memory is loaded properly.",genome);
        fprintf(stderr, usage, name);
        fflush(stderr);
        return(-1);
    }
    return 1;
}

/***********************************************************************************************************/
/* This function associates the files with some part of the memory (memory mapped files)                   */
/***********************************************************************************************************/
int temp1,temp2,temp3;
void read_tables (char *FILEINDEX, char *FILEDATA){
    int temp;
    struct stat fileStat;

    //Open hash table file
    temp1 = open(FILEINDEX, O_RDONLY);
    if (temp1 == -1) {
        perror("Error opening file for reading");
        exit(EXIT_FAILURE);
    }
    //Associate with memory
    hashTable = (unsigned int *) mmap(0,KEY_NUM*sizeof(unsigned int),PROT_READ,MAP_SHARED,temp1,0);
    if (hashTable == MAP_FAILED) {
        close(temp1);
        perror("Error mmapping the file");
        exit(EXIT_FAILURE);
    }

    //Open result table file
    temp2 = open(FILEDATA, O_RDONLY);
    if (temp2 == -1) {
        perror("Error opening file for reading");
        exit(EXIT_FAILURE);
    }
    //obtain file size and assign it to the size shared memory segment and store it too
    stat(FILEDATA,&fileStat);
    total_len = fileStat.st_size/sizeof(MYTYPE);
    //Associate with memory
    tableValues = (MYTYPE *) mmap(0,fileStat.st_size,PROT_READ,MAP_SHARED,temp2,0);
    if (tableValues == MAP_FAILED) {
        close(temp2);
        perror("Error mmapping the file");
        exit(EXIT_FAILURE);
    }
}

/***********************************************************************************************************/
/*************************************** H I T S   O N   G E N O M E ***************************************/
/***********************************************************************************************************/

/***********************************************************************************************************/
/* This function attaches this program to the (already loaded by Load_Memory) shared memory.               */
/* First it converts the given gRNA in its binary form. Then finds the limits of the corresponding 16mer in*/
/* table A (top address and bottom address). The results in table B are structured as NGG | NNNNACA |      */
/* NNGRRT | NAG, thus if PAM is NGG we go to the top address and read until the PAM is NNNNACA. If the PAM */
/* NNNNACA we read through all NGG results and then read the NNNNACA until we find the NNGRRT. If the PAM  */
/* NAG then we go to the bottom address and read from bottom to top until we find an NNGRRT result and if  */
/* is NNGRRT then we go to the bottom address, read from bottom to top, skip through the NAG results, and  */
/* stop when an NNNNACA result is found. For each result we check the extension (last 4 bases of the 20mer)*/
/* and if it is the same as the query's then we add to the results after we make all binary values strings.*/
/* Note that file_creation parses the file from start to end per chromosome. If the input genome file is   */
/* by chromosome then the results are sorted by chromosome and position too. Furtermore the results are    */
/* in a less-mismatch-first way so they are ordered by mismatch number too and no additional ordering is   */
/* required before the output is returned. An exception could be the results of NAG and NNGRRT however,    */
/* because they are read from bottom to top, however, file_creation, fills the NAG and NNGRRT results from */
/* bottom to top as well in Table B, thus everything is ordered!                                           */
/* I N P U T: The gRNA we are searching for (query), the original gRNA from which the previous one is      */
/*            derived (orig_query), the number of mismatches (dots), the positions of the mismatches (p1-5)*/
/*            the vector where the results are returned (results).                                         */
/* O U T P U T: Returns 0 if successfully completed and an error code otherwise                            */
/***********************************************************************************************************/
int find_results (char * query,char *orig_query, int dots, int p1, int p2,int p3, int p4, int p5, char PAM, vector<string> *results){
    unsigned int key = 0;
    MYTYPE extension=0;
    unsigned int twobit;
    char tempresult[256];
    int j,k;
    int positions[5];
    string smallquery(query);
    positions[0]=p1;
    positions[1]=p2;
    positions[2]=p3;
    positions[3]=p4;
    positions[4]=p5;
    
    //convert given string to key
    for (j=0;j<16;j++){
        if (query[j] == 'N')
            return 1;
        twobit = cTOint[(query[j])-ASCIIA];
        key = ( (key << 2) | twobit ) & keymask;
    }
    
    //Find the right extension
    twobit = 0;
    for (j=16;j<20;j++){
        if (query[j] == 'N')
            return 1;
        twobit = cTOint[(query[j])-ASCIIA];
        extension = ( (extension << 2) | twobit ) & mask;
    }
    
    // SEARCH Table A for the corresponding key (16mer)
    long int points = hashTable[key];
    long int pointsnext=0;
    if (key < KEY_NUM - 1){
        pointsnext = hashTable[key+1];
    }
    else{
        pointsnext = total_len;
    }
    //Iterate through the correct bucket, check PAMS, check extensions and add to results
    MYTYPE temp_value,temp_extension,temp_pos,tempInt;
    long int pos,m;
    int startsmall,temp_chrom;
    string temp_lastNs,temp_string="";
    char temp_strand;
    std::vector<string>::iterator it;
    
    if (points != total_len && points != pointsnext){
        //Handle NGG requests
        if (PAM == 'G'){
            pos = points;
            while (pos < pointsnext){
                temp_value = tableValues[pos];
                //If they have the same PAM
                if (intTOpamID[(temp_value & PAMmask) >> PAMSHIFTBITS] != PAM)
                    break;
                temp_extension = temp_value & extensionMask;
                //If they end in the same 4 bases  (table A is built based on the first 16)
                if (temp_extension == extension){
                    //add chromosome
                    temp_chrom = ((temp_value & chromMask) >> CHROMSHIFTBITS);
                    //add strand
                    if ((temp_value & strandMask)==0) temp_strand = '+';
                    else temp_strand = '-';
                    //add position
                    temp_pos = (temp_value & posMask) >> POSSHIFTBITS;
                    //add the small letters
                    startsmall = 0;
                    temp_string="";
                    for (k=0;k<dots;k++){
                        temp_string += smallquery.substr(startsmall,positions[k] - startsmall);
                        temp_string.append(1,cTOsmallC[query[positions[k]] - ASCIIA]);
                        startsmall = positions[k]+1;
                    }
                    if(startsmall <= 19)
                        temp_string+=smallquery.substr(startsmall,string::npos);
                    sprintf(tempresult,"%d\t%c\t%lu\t%lu\t%s\t%s-%cGG\t%d",temp_chrom,temp_strand,temp_pos+1,temp_pos+23,orig_query,temp_string.c_str(),intTOc[(temp_value & lastntMask) >> LASTNTSHIFTBITS],dots);
                    string tempstring(tempresult);
                    //Add result
                    results->push_back(tempstring);
                }
                pos++;
            }
        }
        else if (PAM == 'C'){ //Handle NNNNACA requests
            pos = points;
            //Skip the NGG results
            if (pos < pointsnext){
                temp_value = tableValues[pos];
                if (intTOpamID[(temp_value & PAMmask) >> PAMSHIFTBITS] == 'G'){
                    pos++;
                    temp_value = tableValues[pos];
                    while (pos < pointsnext && intTOpamID[(temp_value & PAMmask) >> PAMSHIFTBITS] == 'G'){
                        pos++;
                        temp_value = tableValues[pos];
                    }
                }
            }
            while (pos < pointsnext){
                temp_value = tableValues[pos];
                //If they have the same PAM
                if (intTOpamID[(temp_value & PAMmask) >> PAMSHIFTBITS] != PAM)
                    break;
                //If they end in the same 4 bases
                temp_extension = temp_value & extensionMask;
                if (temp_extension == extension){
                    //add chromosome
                    temp_chrom = ((temp_value & chromMask) >> CHROMSHIFTBITS);
                    //add strand
                    if ((temp_value & strandMask) == 0) temp_strand = '+';
                    else temp_strand = '-';
                    //add position
                    temp_pos = (temp_value & posMask) >> POSSHIFTBITS;
                    //add found 20mer
                    int startsmall =0;
                    temp_string="";
                    for (k=0;k<dots;k++){
                        temp_string += smallquery.substr(startsmall,positions[k] - startsmall);
                        temp_string.append(1,cTOsmallC[query[positions[k]] - ASCIIA]);
                        startsmall = positions[k]+1;
                    }
                    if(startsmall <= 19)
                        temp_string+=smallquery.substr(startsmall,string::npos);
                    //add pam
                    temp_lastNs = intTOc[(temp_value & lastntMask) >> LASTNTSHIFTBITS];
                    for (m=3;m>0;m--){
                        tempInt = ((temp_value & (last3ntMask << (2*(m-1)))) >> (LAST3NTSHIFTBITS + ((m-1)*2)));
                        temp_lastNs.append(1,intTOc[tempInt]);
                    }
                    sprintf(tempresult,"%d\t%c\t%lu\t%lu\t%s\t%s-%sACA\t%d",temp_chrom,temp_strand,temp_pos+1,temp_pos+27,orig_query,temp_string.c_str(),temp_lastNs.c_str(),dots);
                    string tempstring(tempresult);
                    //Add results
                    results->push_back(tempstring);
                }
                pos++;
            }
        }
        else if (PAM == 'A'){ //Handle NAG requests
            pos = pointsnext - 1;
            while (pos > points - 1){
                temp_value = tableValues[pos];
                //If the have the same PAM
                if (intTOpamID[(temp_value & PAMmask) >> PAMSHIFTBITS] != PAM)
                    break;
                temp_extension = temp_value & extensionMask;
                //If they end in the same 4 bases
                if (temp_extension == extension){
                    //add chromosome
                    temp_chrom = ((temp_value & chromMask) >> CHROMSHIFTBITS);
                    //add strand
                    if ((temp_value & strandMask) == 0) temp_strand = '+';
                    else temp_strand = '-';
                    //add position
                    temp_pos = ((temp_value & posMask) >> POSSHIFTBITS);
                    //add found 20mer
                    startsmall = 0;
                    temp_string="";
                    for (k=0;k<dots;k++){
                        temp_string += smallquery.substr(startsmall,positions[k] - startsmall);
                        temp_string.append(1,cTOsmallC[query[positions[k]] - ASCIIA]);
                        startsmall = positions[k]+1;
                    }
                    if(startsmall <= 19)
                        temp_string+=smallquery.substr(startsmall,string::npos);
                    //add pam
                    sprintf(tempresult,"%d\t%c\t%lu\t%lu\t%s\t%s-%cAG\t%d",temp_chrom,temp_strand,temp_pos+1,temp_pos+23,orig_query,temp_string.c_str(),intTOc[(temp_value & lastntMask) >> LASTNTSHIFTBITS],dots);
                    string tempstring(tempresult);
                    //add result
                    results->push_back(tempstring);
                }
                pos--;
            }
        }
        else if (PAM == 'R'){ //Handle NNGRRT requests
            pos = pointsnext - 1;
            //Skip NAG results
            if (pos > points - 1){
                temp_value = tableValues[pos];
                if (intTOpamID[(temp_value & PAMmask) >> PAMSHIFTBITS] == 'A'){
                    pos--;
                    temp_value = tableValues[pos];
                    while (pos > points - 1  && intTOpamID[(temp_value & PAMmask) >> PAMSHIFTBITS] == 'A'){
                        pos--;
                        temp_value = tableValues[pos];
                    }
                }
            }
            while (pos > points - 1){
                temp_value = tableValues[pos];
                //If they have the same PAM
                if (intTOpamID[(temp_value & PAMmask) >> PAMSHIFTBITS] != PAM){
                    break;
                }
                temp_extension = temp_value & extensionMask;
                //If the 4 last bases of the gRNA are the same
                if (temp_extension == extension){
                    //add chromosome
                    temp_chrom = ((temp_value & chromMask) >> CHROMSHIFTBITS);
                    //add strand
                    if ((temp_value & strandMask) == 0) temp_strand = '+';
                    else temp_strand = '-';
                    //add position
                    temp_pos = (temp_value & posMask) >> POSSHIFTBITS;
                    //add found 20mer
                    int startsmall =0;
                    temp_string = "";
                    for (k=0;k<dots;k++){
                        temp_string += smallquery.substr(startsmall,positions[k] - startsmall);
                        temp_string.append(1,cTOsmallC[query[positions[k]] - ASCIIA]);
                        startsmall = positions[k]+1;
                    }
                    if(startsmall <= 19)
                        temp_string += smallquery.substr(startsmall,string::npos);
                    //add pam
                    temp_lastNs = intTOc[(temp_value & lastntMask) >> LASTNTSHIFTBITS];
                    for (m=3;m>0;m--){
                        tempInt = ((temp_value & (last3ntMask << ((m-1)*2))) >> (LAST3NTSHIFTBITS + ((m-1)*2)));
                        temp_lastNs.append(1,intTOc[tempInt]);
                        if (m==3) temp_lastNs.append(1,'G');
                    }
                    sprintf(tempresult,"%d\t%c\t%lu\t%lu\t%s\t%s-%sT\t%d",temp_chrom,temp_strand,temp_pos+1,temp_pos+26,orig_query,temp_string.c_str(),temp_lastNs.c_str(),dots);
                    string tempstring(tempresult);
                    //add result
                    results->push_back(tempstring);
                }
                pos--;
            }
        }
    }
    return 0;
}


/***********************************************************************************************************/
/* This function creates ALL gRNAs that have at most 5 bases different than the input gRNA and passes them */
/* to the result retrieving function above. It "returns" the results in the vector parameter.              */
/* Each loop below is responsible for a specific number of mismatches. Since it is impossible to create the*/
/* same gRNA with different number of mismatches, each loop creates unique gRNAs. Furthermore, all gRNAs   */
/* that are generated within a loop that handles x mismatches are unique because there are always exactly x*/
/* mismatches and once x positions have been chosen (and all variations of bases created) they never get   */
/* selected again. Those facts are ensured by the limits of the position loops (i,j,k,l,m) and by the inner*/
/* character loops (c1,c2,c3,c4,c5) which enforce the base change per position. To understand the character*/
/* loops better, please look at the cTOnextC in the header files. Essentially given the current base in a  */
/* position it assigns a new one. Since this happens 3 times, all bases are used and none is used twice.   */
/* I N P U T: A gRNA (query), the user's PAM selection (just a char - see main), the vector for the output */
/*            of results, an integer 0-5 for the user's selection of the number of mismatches.             */
/* O U T P U T: None. (The result vector)                                                                  */
/***********************************************************************************************************/
void wildcard_find (char * query, char PAM,vector<string> *results, int mismatches) {
    int error,i,j,k,l,m,c1,c2,c3,c4,c5;
    
    
    //Call printing function for exact results (no mismatches)
    find_results(query,query, 0,-1,-1,-1,-1,-1,PAM,results);
    if (mismatches == 0)
        return;
    
    char temp [21];
    
    //1 mismatch
    //Creates all gRNAs that are different in one position (i) and the
    //difference is one out of 3 bases (c1)depending on the actual base
    //in position i
    for (i=0;i<20;i++){
        strcpy(temp,query);
        for(c1=0;c1<3;c1++){
            temp[i] = cTOnextC[temp[i] - ASCIIA];
            find_results (temp,query, 1,i,-1,-1,-1,-1,PAM,results);
        }
    }
    if (mismatches == 1)
        return;
    
    //2 mismatches
    //Creates all gRNAs that are different from the input in the positions
    //i and j and at those positions have characters c1 and c2.
    for (i=0;i<19;i++){
        for (j=i+1;j<20;j++){
            strcpy(temp,query);
            for(c1=0;c1<3;c1++){
                temp[i] = cTOnextC[temp[i] - ASCIIA];
                for (c2=0;c2<3;c2++){
                    temp[j] = cTOnextC[temp[j] - ASCIIA];
                    find_results (temp,query, 2,i,j,-1,-1,-1,PAM,results);
                }
                //Note that the base in position j needs to be reset because
                //all posiible combinations of bases in the position j need to be
                //created again for the new base (new c2) in position i.
                temp[j] = query[j];
            }
        }
    }
    if (mismatches == 2)
        return;
    
    //3 mismatches
    //Creates all gRNAs that are different from the input in the positions
    //i,j and k and at those positions have characters c1,c2 and c3.
    for (i=0;i<18;i++){
        for (j=i+1;j<19;j++){
            for (k=j+1;k<20;k++){
                strcpy(temp,query);
                for (c1=0;c1<3;c1++){
                    temp[i] = cTOnextC[temp[i] - ASCIIA];
                    for(c2=0;c2<3;c2++){
                        temp[j] = cTOnextC[temp[j] - ASCIIA];
                        for (c3=0;c3<3;c3++){
                            temp[k] = cTOnextC[temp[k] - ASCIIA];
                            find_results (temp,query, 3,i,j,k,-1,-1,PAM,results);
                        }
                        //As explained in 2 mismatches k needs to be reset for the new j
                        temp[k] = query[k];
                    }
                    //As explained in 2 mismatches j needs to be reset for the new i
                    temp[j] = query[j];
                }
            }
        }
    }
    if (mismatches == 3)
        return;
    
    //4 mismatches
    //Creates all gRNAs that are different from the input in the positions
    //i,j,k and l and at those positions have characters c1,c2,c3 and c4.
    for (i=0;i<17;i++){
        for (j=i+1;j<18;j++){
            for (k=j+1;k<19;k++){
                for (l=k+1;l<20;l++){
                    strcpy(temp,query);
                    for (c1=0;c1<3;c1++){
                        temp[i] = cTOnextC[temp[i] - ASCIIA];
                        for (c2=0;c2<3;c2++){
                            temp[j] = cTOnextC[temp[j] - ASCIIA];
                            for (c3=0;c3<3;c3++){
                                temp[k] = cTOnextC[temp[k] - ASCIIA];
                                for (c4=0;c4<3;c4++){
                                    temp[l] = cTOnextC[temp[l] - ASCIIA];
                                    find_results (temp,query, 4,i,j,k,l,-1,PAM,results);
                                }
                                //As explained in 2 mismatches l needs to be reset for the new k
                                temp[l] = query[l];
                            }
                            //As explained in 2 mismatches k needs to be reset for the new j
                            temp[k] = query[k];
                        }
                        //As explained in 2 mismatches j needs to be reset for the new i
                        temp[j] = query[j];
                    }
                }
            }
        }
    }
    if (mismatches == 4)
        return;
    
    //5 mismatches
    //Creates all gRNAs that are different from the input in the positions
    //i,j,k,l and m and at those positions have characters c1,c2,c3,c4 and c5.
    for (i=0;i<16;i++){
        for (j=i+1;j<17;j++){
            for (k=j+1;k<18;k++){
                for (l=k+1;l<19;l++){
                    for (m=l+1;m<20;m++){
                        strcpy(temp,query);
                        for (c1=0;c1<3;c1++){
                            temp[i] = cTOnextC[temp[i] - ASCIIA];
                            for (c2=0;c2<3;c2++){
                                temp[j] = cTOnextC[temp[j] - ASCIIA];
                                for (c3=0;c3<3;c3++){
                                    temp[k] = cTOnextC[temp[k] - ASCIIA];
                                    for (c4=0;c4<3;c4++){
                                        temp[l] = cTOnextC[temp[l] - ASCIIA];
                                        for (c5=0;c5<3;c5++){
                                            temp[m] = cTOnextC[temp[m] - ASCIIA];
                                            find_results(temp,query, 5,i,j,k,l,m,PAM,results);
                                        }
                                        //As explained in 2 mismatches m needs to be reset for the new l
                                        temp[m] = query[m];
                                    }
                                    //As explained in 2 mismatches l needs to be reset for the new k
                                    temp[l] = query[l];
                                }
                                //As explained in 2 mismatches k needs to be reset for the new j
                                temp[k] = query[k];
                            }
                            //As explained in j mismatches k needs to be reset for the new i
                            temp[j] = query[j];
                        }
                    }
                }
            }
        }
    }
}

/***********************************************************************************************************/
/************************************* A N N O T A T I O N   I N F O ***************************************/
/***********************************************************************************************************/

/***********************************************************************************************************/
/* This function adds the annotation info to the already calculated results. The results are calculated by */
/* the two calculation functions and kept in a table. This function attaches the annotation info per       */
/* annotation file per result. That means that each annotation file is opened once and searched for all    */
/* results instead of filling all info per result. You can imagine that the annotation info is filled by   */
/* column (eg. 3UTR) instead of by row.                                                                    */
/* More specificaly                                                                                        */
/*                                                                                                         */
/* Please note:                                                                                            */
/* THE ANNOTATION FILES SHOULD BE ALREADY SORTED. Also, their structure should be the following:           */
/* chromosome <tab> strand <tab> start coordinate <tab> end coordinate <tab> annotation info <newline>     */
/***********************************************************************************************************/
void attach_class_details(vector<string> *results, char* filepath, unsigned int backward_gap){
    string line,temp0,temp1,temp2,temp3,temp4,temp5, temp6,temp7;
    coords temp;
    std::vector<coords> file_coords;
    std::vector<string>::iterator it;
    unsigned int indexes [MAX_CHROM_NUM+1];
    unsigned int i,j,k,current_chrom=1,chrom,strand,start,end,indexstart,indexend,temp_start,temp_end,mid;
    bool found;
    
    //Initialize index array
    for (j=0;j<MAX_CHROM_NUM+1;j++){
        indexes[j]=0;
    }
    
    //Load the file
    ifstream myfile (filepath);
    if (myfile.is_open()){
        i=0;
        j=1;
        while (getline (myfile,line)){
            stringstream iss(line);
            getline(iss, temp0, '\t'); //chrom
            getline(iss, temp1, '\t'); //strand
            getline(iss, temp2, '\t'); //start
            getline(iss, temp3, '\t'); //end
            getline(iss, temp4, '\n');
            temp.chrom = atoi(temp0.c_str());
            temp.start = atoi(temp2.c_str());
            temp.end = atoi(temp3.c_str());
            temp.info = temp4;
            //if (temp1 == "+") temp.strand = 0;
            //else temp.strand = 1;
            //keep index info
            //the index at position
            while (temp.chrom != current_chrom){
                indexes[j]=i;
                j++;
                current_chrom = current_chrom +1;
            }
            //keep entry
            file_coords.push_back(temp);
            i++;
        }
        myfile.close();
        //keep info on last chrom + strand index. Equals filesize
        indexes[j] = i;
        
        
        //iterate over results
        for(it = results->begin(); it != results->end(); ++it) {
            //find each result in the file
            stringstream iss(*it);
            temp4="";
            getline(iss, temp0, '\t'); //chrom
            getline(iss, temp1, '\t'); //strand
            getline(iss, temp2, '\t'); //start coord
            getline(iss, temp3, '\t'); //end coord
            getline(iss, temp4, '\n'); //the rest
            chrom = atoi(temp0.c_str());
            //if (temp1 == "+") strand = 0;
            //else strand = 1;
            start = atoi(temp2.c_str());
            end = atoi(temp3.c_str());
            indexstart = indexes[chrom -1];
            i=indexstart;
            indexend = indexes[chrom];
            k=indexend;
            temp5 = "";
            found = false;
            if (file_coords[i].start < start){
                while (indexend > indexstart){
                    mid = (unsigned int)(indexend - indexstart)/(double)2 + indexstart;
                    //Problem if transcripts don't have the same length (that is why the -30 if is needed below)
                    if (file_coords[mid].start < start && file_coords[mid].end < start){
                        i = indexstart;
                        indexstart = mid + 1 ;
                    }
                    else if (file_coords[mid].start > end){
                        k = indexend;
                        indexend = mid - 1;
                    }
                    else{
                        break;
                    }
                }
            }
            //Search the range for results
            if (i > backward_gap)
                i=i-backward_gap;
            else
                i=indexes[chrom -1];
            while (i < k){
                temp_start = file_coords[i].start;
                temp_end = file_coords[i].end;
                if(chrom == file_coords[i].chrom && ((temp_start <= start && end <= temp_end) || (start <= temp_start && temp_start <= end) || (start <= temp_end && temp_end <= end))){
                    //expand entry
                    found = true;
                    temp5 += (file_coords[i].info + "~");
                }
                i++;
            }
            if (found){
                (*it) += ("\t") + temp5;
            }
            else{
                (*it) += "\t-";
            }
        }
    }
    else{
        //If filename doesn't exist
        for(it = results->begin(); it != results->end(); ++it) {
            (*it) += ("\tN/A");
        }
    }
}

/***********************************************************************************************************/
/************************************** g R N A   D I S C O V E R Y ****************************************/
/***********************************************************************************************************/

/***********************************************************************************************************/
/* */
/***********************************************************************************************************/

void searchSequence(char* s, char* PAM, char *revPAM, char *output){
    int i=0, lengths,end, offset1,offset2;
    char *ps,*temps, twentyMer[23];
    bool add20mer = true, atleast1 = false;
    
    if (strcmp(PAM,"ACA")==0){
        if(strlen(s) < 27){
            output[0] = '\0';
            //cout << "exited ACA" << endl;
            return;
        }
        offset1 = 7;
        offset2 = 24;
        temps = s + 24;
    }
    else{
        if(strlen(s) < 23){
            output[0] = '\0';
            //cout << "exited here " << strlen(s) << endl;
            return;
        }
        temps = s + 21;
        offset1 = 3;
        offset2 = 21;
    }
    
    lengths = strlen(s);
    
    //FORWARD
    while ((ps = strstr(temps,PAM)) != NULL){
        strncpy(twentyMer, ps - offset2, 20);
        twentyMer[20]='-';
        twentyMer[21]='-';
        twentyMer[22]='\0';
        i=0;
        add20mer = true;
        while(i<20){
            if (twentyMer[i] != 'A' && twentyMer[i] != 'T' && twentyMer[i] != 'G' && twentyMer[i] != 'C'){
                add20mer = false;
                cout << "found wrong char" << twentyMer[i] << endl;
                break;
            }
            i++;
        }
        if(add20mer){
            strcat(output,twentyMer);
            atleast1 = true;
        }
        temps = ps + 1;
    }
    
    //REVERSE
    temps = s;
    while((ps = strstr(temps,revPAM)) != NULL){
        if (strlen(ps) >= 20 + offset1){
            //check if it doesn't exceed length
            strncpy(twentyMer,ps+offset1,20);
            twentyMer[20] = '-';
            twentyMer[21]='-';
            twentyMer[22]='\0';
            i=0;
            add20mer = true;
            while(i<20){
                if (twentyMer[i] != 'A' && twentyMer[i] != 'T' && twentyMer[i] != 'G' && twentyMer[i] != 'C'){
                    add20mer = false;
                    cout << "found wrong char 2 " << twentyMer[i] << ". " << i << endl;
                    break;
                }
                i++;
            }
            if(add20mer){
                atleast1 = true;
                //copy reverse 20mer
                i=19;
                end = strlen(output);
                while(i>=0){
                    output[end] = cTOcREV[twentyMer[i]-ASCIIA];
                    i--;
                    end++;
                }
                output[end] = '-';
                output[end+1] = '-';
                output[end+2] = '\0';
            }
        }
        temps = ps + 1;
    }
    if(!atleast1){
        output[0] = '\0';
    }
}

void searchSequenceWithGaps(char *s, char *output){
    int i=0,j, lengths,end;
    char twentyMer[23],c1,c2;
    bool add20mer = true, atleast1 = false;
    
    lengths = strlen(s);
    //FORWARD
    j=0;
    while (j<lengths-25){
        c1 = s[j+23];
        c2 = s[j+24];
        if(s[j+22] == 'G' && s[j+25]=='T' && (c1 =='A' || c1 == 'G') && (c2 == 'A' || c2 == 'G')){
            strncpy(twentyMer, s + j, 20);
            twentyMer[20]='-';
            twentyMer[21]='-';
            twentyMer[22]='\0';
            i=0;
            add20mer = true;
            while(i<20){
                if (twentyMer[i] != 'A' && twentyMer[i] != 'T' && twentyMer[i] != 'G' && twentyMer[i] != 'C'){
                    add20mer = false;
                    cout << "found wrong char" << twentyMer[i] << endl;
                    break;
                }
                i++;
            }
            if(add20mer){
                strcat(output,twentyMer);
                atleast1 = true;
            }
        }
        j++;
    }
    //REVERSE
    j=0;
    while (j<lengths-25){
        c1 = s[j+1];
        c2 = s[j+2];
        if(s[j] == 'A' && s[j+3]=='C' && (c1 =='T' || c1 == 'C') && (c2 == 'T' || c2 == 'C')){
            strncpy(twentyMer, s + j + 6, 20);
            twentyMer[20]='-';
            twentyMer[21]='-';
            twentyMer[22]='\0';
            i=0;
            add20mer = true;
            while(i<20){
                if (twentyMer[i] != 'A' && twentyMer[i] != 'T' && twentyMer[i] != 'G' && twentyMer[i] != 'C'){
                    add20mer = false;
                    cout << "found wrong char" << twentyMer[i] << endl;
                    break;
                }
                i++;
            }
            if(add20mer){
                atleast1=true;
                //copy reverse 20mer
                i=19;
                end = strlen(output);
                while(i>=0){
                    output[end] = cTOcREV[twentyMer[i]-ASCIIA];
                    i--;
                    end++;
                }
                output[end] = '-';
                output[end+1] = '-';
                output[end+2] = '\0';
            }
        }
        j++;
    }
    if(!atleast1){
        output[0] = '\0';
    }
}

/***********************************************************************************************************/
/**************************************** M A I N   F U N C T I O N  ***************************************/
/***********************************************************************************************************/
/* ARGUMENTS:                                                                                              */
/* -i : Input gRNAs or filename                                                                            */
/* -f : flag input is filename (optional)                                                                  */
/* -p : pam                                                                                                */
/* -n : maximum mismatch number                                                                            */
/* -g : genome name                                                                                        */
/* -t : table (.bin) folder path                                                                            */
/* -o : output file name, if present output will be written on file otherwise printed (optional)           */
/* -m : if 0 shared memory, else 1 memory mapped files (optional - default 0)                              */
/* -a : annotation files full paths divided by --  (optional)                                              */
/***********************************************************************************************************/
int main (int argc, char *argv[]){
    int mismatchNUM =-1;
    char *inputfile=NULL,*input=NULL,*in_file_path=NULL,*memoryType=NULL,*output_file=NULL,*ptr,*ptr_prev=NULL;
    char c='0';
    std::vector<string> results;
    std::vector<string>::iterator it;
    unsigned int genome_backward_gap;
    
    //Type usage information and print that upon error
    char *usage = (char*)"\nUsage: %s -i <Full path of input file or input gRNAs>\n-f <file input flag>\n-p <PAM>\n-n <mismatch max number>\n-g <genome name>\n-m <0 for shared memory and 1 for memory mapped files>\n-t <full  path of *.bin files if using -m1>\n-o <file name for output file, if not given output is stdout OPTIONAL>\n-a<Filepaths with gene info (divided by --) OPTIONAL>\nWhere PAM is G for NGG, A for NAG, R for NNGRRT, and C for NNNNACA\n<Filepaths with gene info (divided by --)>: is optional\n\n\tPlease type hg38 or GRCh38 for human hg38\n\tPlease type hg19 or GRCh37 for human hg19\n\tPlease type mm10 or GRCm38 for mouse mm10\n\tPlease type w303 for yeast strain w303\n\n";
    
    //Parse arguments
    int a;
    char *genome=NULL, *PAM=NULL;
    
    while ((a = getopt (argc, argv, "i:fp:n:g:t:o:m:a:")) != -1){
        if(a == 'i'){
            input=optarg;
        }
        else if(a == 'f'){
            inputfile=input;
        }
        else if(a == 'p'){
            PAM=optarg;
            if (strcmp(PAM,"G")==0 || strcmp(PAM,"NGG")==0){ c = 'G';}
            else if (strcmp(PAM,"A")==0 || strcmp(PAM,"NAG")==0){ c = 'A';}
            else if (strcmp(PAM,"C")==0 || strcmp(PAM,"NNNNACA")==0){ c = 'C';}
            else if (strcmp(PAM,"R")==0 || strcmp(PAM,"NNGRRT")==0){ c = 'R';}
            else{
                fprintf(stderr, "\nUnknown PAM sequence. You selected %s\n", PAM);
                fprintf(stderr, usage, argv[0]);
                fflush(stderr);
                return(-1);
            }
        }
        else if(a == 't'){
            in_file_path = optarg;
        }
        else if(a == 'o'){
            output_file=optarg;
        }
        else if(a == 'n'){
            mismatchNUM = atoi(optarg);
            if (mismatchNUM < 0 || mismatchNUM > 5){
                fprintf(stderr, "\nUnsuported max mismatch number. Please choose a value [0-5]. You selected %d\n", mismatchNUM);
                fprintf(stderr, usage, argv[0]);
                fflush(stderr);
                return(-1);
            }
        }
        else if (a == 'm'){
            memoryType=optarg;
        }
        else if (a == 'a'){
            ptr_prev=optarg;
        }
        else if (a == 'g'){
            genome=optarg;
            //Set genome dependent variables
            if (strcmp(genome,"hg19")==0 || strcmp(genome,"GRCh37")==0){
                HASHKEY=HASHKEYHUMAN;
                DATAKEY=DATAKEYHUMAN;
                SIZESKEY=SIZESKEYHUMAN;
                //genome_backward_gap=MAX_BACKWARD_GAP_HG19;
                genome_backward_gap=MAX_BACKWARD_GAP;
            }
            else if (strcmp(genome,"mm10")==0 || strcmp(genome,"GRCm38")==0){
                HASHKEY=HASHKEYMOUSE;
                DATAKEY=DATAKEYMOUSE;
                SIZESKEY=SIZESKEYMOUSE;
                //genome_backward_gap=MAX_BACKWARD_GAP_MM10;
                genome_backward_gap=MAX_BACKWARD_GAP;
            }
            else if (strcmp(genome,"hg38")==0 || strcmp(genome,"GRCh38")==0){
                HASHKEY=HASHKEYHUMAN2;
                DATAKEY=DATAKEYHUMAN2;
                SIZESKEY=SIZESKEYHUMAN2;
                //genome_backward_gap=MAX_BACKWARD_GAP_HG38;
                genome_backward_gap=MAX_BACKWARD_GAP;
            }
            else if (strcmp(genome,"w303")==0){
                HASHKEY=HASHKEYYEAST;
                DATAKEY=DATAKEYYEAST;
                SIZESKEY=SIZESKEYYEAST;
                genome_backward_gap=MAX_BACKWARD_GAP;
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
            fprintf(stderr, "\nUknown error\n");
            fprintf(stderr, usage, argv[0]);
            fflush(stderr);
            return(-2);
        }
    }
    //General correctness checks
    //Check that all necessary variables are provided
    if (input==NULL || c=='0' || mismatchNUM==-1 || genome==NULL){
        fprintf(stderr, "\nOne of the following mandatory arguments is missing:\n");
        fprintf(stderr, "(-i) input. Value entered: %s\n", input);
        fprintf(stderr, "(-p) PAM. Value entered: %s\n",PAM);
        fprintf(stderr, "(-n) maximum mismatch number. Value entered: %d\n",mismatchNUM);
        fprintf(stderr, "(-g) genome name. Value entered: %s\n",genome);
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        return(1);
    }
    //Get gRNAs from command line input
    if (inputfile==NULL && strlen(input) < 20){
        fprintf(stderr, "\nNo gRNAs found on input %s, input is less than 20 characters\n", input);
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        return(1);
    }
    //If shared memory files are selected then the *.bin path is also necessary
    if (memoryType!=NULL){
        if (strcmp(memoryType,"1")==0 && in_file_path==NULL){
            fprintf(stderr, "\nTo use memory shared files you need to enter the path of the data.bin and index.bin files.");
            fprintf(stderr, usage, argv[0]);
            fflush(stderr);
            return(1);
        }
    }
    else{
    //if memoryType is not set, set it to default
        memoryType= new char [2];
        memoryType= (char *)"0";
    }

    //Attach shared memory or load shared memory file
    if (strcmp(memoryType,"1") == 0){
        //Use shared memory files
        //Create the input filenames
        cout << "Mapping files to memory" << endl;
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
        read_tables(index_file,data_file);
    }
    else{
        //Use shared memory (default if no choice is made)
        cout << "Attaching Shared Memory" << endl;
        if (attach_shared_memory(genome,usage,argv[0]) < 0){
            return(-1);
        }
    }
    
    //Input can get gRNAs or sequence from file OR get gRNAs from the command line input directly divided by --)
    //If input is a file we need to find the gRNAs first
    if (inputfile != NULL){
        cout << "Reading input file " << inputfile << endl;
        //Input was a file. Get gRNAs or sequence from file
        ifstream myfile (inputfile);
        myfile.seekg(0, ios::end);
        int llen = myfile.tellg();
        myfile.seekg(0, ios::beg);
        char clean_input[llen+1];
        input= new char [(llen*llen) + 1];
        input[0]='\0';
        bool sequence = false;
        string line;
        //Check if file contains gRNAs only or sequence
        //If only gRNAs (while each line has 20 characters) divide them by --
        if (myfile.is_open()){
            while (getline (myfile,line)){
                llen = line.size();
                if (llen == 0){
                    continue;
                }
                else if (llen != 20){
                    sequence = true;
                    break;
                }
                strcat(input,line.c_str());
                strcat(input,"--");
            }
        }
        else{
            fprintf(stderr, "\n\nAn error occured while reading input from file %s\n", inputfile);
            fprintf(stderr, usage, argv[0]);
            fflush(stderr);
            return(-1);
        }
        //If the input contains a sequence and not individual gRNAs find gRNAs for sequence and individual gRNAs accordingly
        clean_input[0]='\0';
        if (sequence == true){
            //move file pointer at start
            myfile.clear();
            myfile.seekg(0, ios::beg);
            //Put the input in a string and skip headers and empty lines
            while (getline (myfile,line)){
                if (line[0] != '>' && line.length() != 0){
                    strcat(clean_input,line.c_str());
                }
            }
            //Search for gRNAs
            if (c=='G') { searchSequence(clean_input,(char *)"GG",(char *)"CC",input); }
            else if (c=='A') { searchSequence(clean_input,(char *)"AG",(char*)"CT",input);}
            else if (c=='C') { searchSequence(clean_input,(char *)"ACA",(char*)"TGT",input);}
            else { searchSequenceWithGaps(clean_input,input);}
        }
        myfile.close();
    }
    
    //Find results for all gRNAs
    cout << "\n";
    int i=0;
    char gRNA[21];
    char *input_ptr = input;
    int len = strlen(input_ptr);
    while (i+19<len){
        strncpy(gRNA,input_ptr,20);
        gRNA[20]='\0';
        i+=22;
        input_ptr+=22;
        cout << "Calculating results of " << gRNA << endl;
        wildcard_find(gRNA,c,&results,mismatchNUM);
    }

    //Add annotations for all annotation files on input
    if (ptr_prev != NULL){
        char filepath[strlen(ptr_prev) +1];
        cout << "\nAnnotating results using:" << endl;
        while ((ptr = strstr(ptr_prev,"--")) != NULL){
            memset(filepath,'\0',PATH_LEN);
            strncpy(filepath,ptr_prev,ptr - ptr_prev);
            cout << filepath << endl;
            attach_class_details(&results,filepath,genome_backward_gap);
            ptr_prev = ptr+2;
        }
        //last one or only one
        strcpy(filepath,ptr_prev);
        cout << filepath << endl;
        attach_class_details(&results,filepath,genome_backward_gap);
    }

    //OUTPUT can be on a file or on stdout
    if (output_file == NULL){
        //print results on stdout
        for(it = results.begin(); it != results.end(); ++it) {
            cout << (*it) << endl;
        }
    }
    else{
        //print results on file
        ofstream outfile;
        outfile.open (output_file);
        for(it = results.begin(); it != results.end(); ++it) {
            outfile << (*it) << endl;
        }
        outfile.close();
        cout << "\nYour results are ready. Thank you for using Off-Spotter.\n" << endl;
    }

    if (strcmp(memoryType,"1") != 0){
        //unattach shared memory
        shmdt(sizes);
        shmdt(hashTableADDR);
        shmdt(resultADDR);
    }
    else{
        //unmap file
        munmap(hashTable,KEY_NUM*sizeof(unsigned int));
        munmap(tableValues,temp3);
        close(temp1);
        close(temp2);
    }
    return 0;
}
