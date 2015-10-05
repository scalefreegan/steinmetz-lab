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

/***********************************************************************************************************/
/************************************************* SUMMARY *************************************************/
/***********************************************************************************************************/
/* This program populates table A and table B. In other words pre-calculates ALL the hits of ALL gRNAs for */
/* ALL the PAMS handled. Then writes the results on 2 binary files. To calculate all the results it will   */
/* take some time and it needs to be run ONLY whenever you change the genome you are working on. To go     */
/* through the hits, use load memory and results and to dettach the memory use dettach memory at the end.  */
/* To do so, it uses 2 kinds of functions. The count functions and the populating functions. The first kind*/
/* of functions counts the total number of results as well as the results per gRNA. Those are used to      */
/* allocate space for table B and to populate table A. Table A has one entry per every possible 16mer that */
/* contains the total number of hits that start with that 16mer for all gRNAs and all PAMs. This number    */
/* works as a pointer to the results of each 16mer in Table B because table be is populated in consecutive */
/* addresses spaces by gRNAs in alphabetical order. See [1] for more details. The second type of functions */
/* is used to re-calculate the hits, put them in the correct format and put them on table B. The results   */
/* are placed on table B according to numbers found by the counting functions and temporary counters that  */
/* keep the current occupation of each 16mers spaces.                                                      */
/***********************************************************************************************************/
/***********************************************************************************************************/
/***********************************************************************************************************/

#include "Off-Spotter_Table_Creation.h"


using namespace std;


/***********************************************************************************************************/
/*******************************************   COUNT FUNCTIONS   *******************************************/
/***********************************************************************************************************/
/* Those functions count the amount of results per PAM per 16base key. The results of all the functions    */
/* kept in the same table (one counter per key), this table is currentPos. total_len is a counter that has */
/* the number of all results of all PAMs. Thus in the end of each of those functions we know how many hits */
/* there are in the genome per 16mer for this PAM and their total number. We need this info to allocate    */
/* enough space for table B and populate with the sum of hits of each individual 16mer for all PAMs table A*/
/***********************************************************************************************************/
/***********************************************************************************************************/
/***********************************************************************************************************/

/***********************************************************************************************************/
/* This function counts the appearances of ALL targets of ALL possible gRNAs followed by NGG or NAG on the */
/* filename argument. CurrentPos is keeping counts for each gRNA. It returns the total number.             */
/***********************************************************************************************************/
unsigned int countNNG(const char *filename, char PAM, unsigned int *currentPos){
    unsigned int i,j,end,size,key,total_num=0;
    char RevPAM;
    string chrom(MAX_CHROM_LEN,'\0');
    MYTYPE twobit,value,rest4;
    bool insert_entry;
    
    //Depending on the input PAM find reverse for 2nd N
    if (PAM=='G')
        RevPAM = 'C';
    else if (PAM == 'A')
        RevPAM = 'T';
    else
        return 0;
    
    //Open the input file for reading
    ifstream ifile (filename,std::ifstream::in);
    
    //Count appearances per 20mer
    //Forward
    while (!ifile.eof()){
        //header
        getline(ifile,chrom);
        if (ifile.eof()) break;
        //all lines up to next header
        getline(ifile,chrom);
        if (ifile.eof()) break;
        size = chrom.size();
        i=0;
        while (i < size - 22){
            if(chrom[i+21] == PAM && chrom[i+22] == 'G'){
                insert_entry = true;
                //Make the key
                end = i + 16;
                key=0;
                twobit=0;
                j=i;
                while(j < end){
                    if (chrom[j] == 'N'){
                        insert_entry = false;
                        break;
                    }
                    twobit = cTOint[(chrom[j])-ASCIIA];
                    key = ( (key << 2) | twobit ) & keymask;
                    j++;
                }
                if (insert_entry){
                    //Check for Ns
                    end = j;
                    for (j=end;j<end+4;j++)
                        if (chrom[j] == 'N'){
                            insert_entry = false;
                            break;
                        }
                    if (insert_entry){
                        currentPos[key]++;
                        total_num++;
                    }
                }
            }
            i++;
        }
    }
    
    //move file pointer to the beggining
    ifile.clear();
    ifile.seekg(0, ios::beg);
    
    //Reverse
    while (!ifile.eof()){
        getline(ifile,chrom);
        if (ifile.eof()) break;
        getline(ifile,chrom);
        if (ifile.eof()) break;
        size = chrom.size();
        i=0;
        while (i < size - 22){
            if(chrom[i] == 'C' && chrom[i+1] == RevPAM){
                insert_entry = true;
                //Make the key
                end = i + 6;
                key=0;
                twobit=0;
                j=i+22;
                while(j > end){
                    if (chrom[j] == 'N'){
                        insert_entry = false;
                        break;
                    }
                    twobit = cTOintREV[chrom[j]-ASCIIA];
                    key = ( (key << 2) | twobit ) & keymask;
                    j--;
                }
                //Check for Ns
                if (insert_entry){
                    end = j;
                    for (j=end;j>end-4;j--)
                        if (chrom[j] == 'N'){
                            insert_entry = false;
                            break;
                        }
                    if (insert_entry){
                        currentPos[key]++;
                        total_num++;
                    }
                }
            }
            i++;
        }
    }
    ifile.close();
    
    return total_num;
}

/***********************************************************************************************************/
/* This function counts the appearances of ALL targets of ALL possible gRNAs followed by NNNNACA on the    */
/* filename argument. CurrentPos is keeping counts for each gRNA. It returns the total number.             */
/***********************************************************************************************************/
int countNNNNACA(const char *filename, unsigned int *currentPos){
    unsigned int i,j,end,size,key,total_num=0;
    string chrom(MAX_CHROM_LEN,'\0');
    MYTYPE twobit,value,rest4;
    bool insert_entry;
    
    //Open the input file for reading
    ifstream ifile (filename,std::ifstream::in);
    
    //Count appearances per 20mer
    //Forward
    while (!ifile.eof()){
        getline(ifile,chrom);
        if (ifile.eof()) break;
        getline(ifile,chrom);
        if (ifile.eof()) break;
        size = chrom.size();
        i=0;
        while (i < size - 26){
            if(chrom[i+24] == 'A' && chrom[i+25] == 'C' && chrom[i+26] == 'A'){
                insert_entry = true;
                //Make the key
                end = i + 16;
                key=0;
                twobit=0;
                j=i;
                while(j < end){
                    if (chrom[j] == 'N'){
                        insert_entry = false;
                        break;
                    }
                    twobit = cTOint[(chrom[j])-ASCIIA];
                    key = ( (key << 2) | twobit ) & keymask;
                    j++;
                }
                if (insert_entry){
                    //Check for Ns
                    end = j;
                    for (j=end;j<end+4;j++)
                        if (chrom[j] == 'N'){
                            insert_entry = false;
                            break;
                        }
                    if (insert_entry){
                        currentPos[key]++;
                        total_num++;
                    }
                }
            }
            i++;
        }
    }
    
    //move file pointer to the beggining
    ifile.clear();
    ifile.seekg(0, ios::beg);
    
    //Reverse
    while (!ifile.eof()){
        getline(ifile,chrom);
        if (ifile.eof()) break;
        getline(ifile,chrom);
        if (ifile.eof()) break;
        size = chrom.size();
        i=0;
        while (i < size - 26){
            if(chrom[i] == 'T' && chrom[i+1] == 'G' && chrom[i+2] == 'T'){
                insert_entry = true;
                //Make the key
                end = i + 10;
                key=0;
                twobit=0;
                j=i+26;
                while(j > end){
                    if (chrom[j] == 'N'){
                        insert_entry = false;
                        break;
                    }
                    twobit = cTOintREV[chrom[j]-ASCIIA];
                    key = ( (key << 2) | twobit ) & keymask;
                    j--;
                }
                if (insert_entry){
                    //Check for Ns
                    end = j;
                    for (j=end;j>end-4;j--)
                        if (chrom[j] == 'N'){
                            insert_entry = false;
                            break;
                        }
                    if (insert_entry){
                        currentPos[key]++;
                        total_num++;
                    }
                }
            }
            i++;
        }
    }
    ifile.close();
    return total_num;
}

/***********************************************************************************************************/
/* This function counts the appearances of ALL targets of ALL possible gRNAs followed by NNGRRT on the     */
/* filename argument. CurrentPos is keeping counts for each gRNA. It returns the total number.             */
/***********************************************************************************************************/
int countNNGRRT(const char *filename, unsigned int *currentPos){
    unsigned int i,j,end,size,key,total_num=0;
    string chrom(MAX_CHROM_LEN,'\0');
    MYTYPE twobit,value,rest4;
    bool insert_entry;
    
    //Open the input file for reading
    ifstream ifile (filename,std::ifstream::in);
    
    //Count appearances per 20mer
    //Forward
    while (!ifile.eof()){
        getline(ifile,chrom);
        if (ifile.eof()) break;
        getline(ifile,chrom);
        if (ifile.eof()) break;
        size = chrom.size();
        i=0;
        while (i < size - 25){
            if(chrom[i+22] == 'G' && (chrom[i+23] == 'A' || chrom[i+23] == 'G') && (chrom[i+24] == 'A' || chrom[i+24] == 'G') && chrom[i+25] == 'T'){
                insert_entry = true;
                //Make the key
                end = i + 16;
                key=0;
                twobit=0;
                j=i;
                while(j < end){
                    if (chrom[j] == 'N'){
                        insert_entry = false;
                        break;
                    }
                    twobit = cTOint[(chrom[j])-ASCIIA];
                    key = ( (key << 2) | twobit ) & keymask;
                    j++;
                }
                if (insert_entry){
                    //Check for Ns
                    end = j;
                    for (j=end;j<end+4;j++)
                        if (chrom[j] == 'N'){
                            insert_entry = false;
                            break;
                        }
                    if (insert_entry){
                        currentPos[key]++;
                        total_num++;
                    }
                }
            }
            i++;
        }
    }
    
    //move file pointer to the beggining
    ifile.clear();
    ifile.seekg(0, ios::beg);
    
    //Reverse
    while (!ifile.eof()){
        getline(ifile,chrom);
        if (ifile.eof()) break;
        getline(ifile,chrom);
        if (ifile.eof()) break;
        size = chrom.size();
        i=0;
        while (i < size - 25){
            if(chrom[i] == 'A' && (chrom[i+1] == 'T' || chrom[i+1] == 'C') && (chrom[i+2] == 'T' || chrom[i+2] == 'C') && chrom[i+3] == 'C'){
                insert_entry = true;
                //Make the key
                end = i + 9;
                key=0;
                twobit=0;
                j=i+25;
                while(j > end){
                    if (chrom[j] == 'N'){
                        insert_entry = false;
                        break;
                    }
                    twobit = cTOintREV[chrom[j]-ASCIIA];
                    key = ( (key << 2) | twobit ) & keymask;
                    j--;
                }
                if (insert_entry){
                    //Check for Ns
                    end = j;
                    for (j=end;j>end-4;j--)
                        if (chrom[j] == 'N'){
                            insert_entry = false;
                            break;
                        }
                    if (insert_entry){
                        currentPos[key]++;
                        total_num++;
                    }
                }
            }
            i++;
        }
    }
    ifile.close();
    return total_num;
}


/***********************************************************************************************************/
/***********************************************************************************************************/
/***********************************************************************************************************/
/***********************************************************************************************************/
/******************************************   RESULT FUNCTIONS   *******************************************/
/***********************************************************************************************************/
/* Those functions populate the result table with the actual results. They get in by PAM and first the     */
/* forward strand and then the reverse strand.                                                             */
/***********************************************************************************************************/
/***********************************************************************************************************/
/***********************************************************************************************************/
/***********************************************************************************************************/

/***********************************************************************************************************/
/* This function populates table B. Table A is already populated by the results of the count functions.    */
/* The input is the filename of the genome file, the PAM (NGG or NAG), a table with the positions that are */
/* currently occupied, table A, table B (to be filled), and a number (0-3).                                */
/* Essentially the function looks for hits exactly like the count functions above, but everytime one is    */
/* found instead of increasing a counter, the first 16 bases of the 20mer are looked up on table A and then*/
/* the result is inserted on the position indicated by table_A[16mer] + currentPos[16mer]. This temporary  */
/* table is identical with table A in structure but instead of the total number of hits it holds the number*/
/* of the hits of that specific 16mer that we already had inserted on table B. Hence, we know the exact    */
/* position to place the hit. The number is an ID that's not hardcoded to make the addition/removal of PAMS*/
/* simpler. It is added to each entry and used to identify the PAM when presenting the results to the user.*/
/***********************************************************************************************************/
void resultsNNG(const char *filename, char PAM, unsigned int *currentPos, MYTYPE *tableValues, unsigned int *hashTable, unsigned int PAMID){
    unsigned int i,j,end,size,key,cnum;
    char RevPAM;
    string chrom(MAX_CHROM_LEN,'\0');
    string chromosome;
    MYTYPE twobit,value,rest4;
    bool insert_entry;
    
    //Depending on the input PAM find reverse for 2nd N
    if (PAM=='G')
        RevPAM = 'C';
    else if (PAM == 'A')
        RevPAM = 'T';
    else
        return;
    
    //Open the input file for reading
    ifstream ifile (filename,std::ifstream::in);
    
    //Add actual entries to the tables
    //Populating the hash table with forward entries first
    while (!ifile.eof()){
        getline(ifile,chrom);
        if (ifile.eof()) break;
        //Read Header
        if (!isdigit(chrom[2]))
            chromosome.assign(&chrom[1],1);
        else
            chromosome.assign(&chrom[1],2);
        if(chromosome == "X") cnum = 23;
        else if (chromosome == "Y") cnum = 24;
        else if (chromosome == "M") cnum = 25;
        else cnum = atoi(chromosome.c_str());
        getline(ifile,chrom);
        if (ifile.eof()) break;
        size = chrom.size();
        i=0;
        while (i < size - 22){
            if(chrom[i+21] == PAM && chrom[i+22] == 'G'){
                insert_entry = true;
                //Make the key
                end = i + 16;
                key=0;
                twobit=0;
                j=i;
                while(j < end){
                    if (chrom[j] == 'N'){
                        insert_entry = false;
                        break;
                    }
                    twobit = cTOint[(chrom[j])-ASCIIA];
                    key = ( (key << 2) | twobit ) & keymask;
                    j++;
                }
                if (insert_entry){
                    //Make the value
                    twobit = 0;
                    value=0;
                    rest4 =0;
                    //Fix PAMID
                    value = PAMID;
                    //Fix strand (shift one - already 0)
                    value = value << STRANDBITS;
                    //Fix last nt before GG
                    value = ((value << LASTNTBITS)  | cTOint[(chrom[i+20])-ASCIIA]) & mask;
                    //Fix chromosome
                    value = ((value << CHROMBITS) | cnum) & mask;
                    //Fix position
                    value = ((value << POSBITS) | i) & mask;
                    //Fix last4 nt
                    end = j;
                    for (j=end;j<end+4;j++){
                        if (chrom[j] == 'N'){
                            insert_entry = false;
                            break;
                        }
                        twobit = cTOint[(chrom[j])-ASCIIA];
                        rest4 = ((rest4 << 2) | twobit) & mask;
                    }
                    if (insert_entry){
                        value = ((value << EXTENSIONBITS) | rest4) & mask;
                        //Add newly constructed pair to the structure
                        tableValues[hashTable[key] + currentPos[key]] = value;
                        currentPos[key]++;
                    }
                }
            }
            i++;
        }
    }
    
    //move file pointer to the beggining
    ifile.clear();
    ifile.seekg(0, ios::beg);
    
    //for (c=0;c<MAX_CHROM_NUM;c++){
    while (!ifile.eof()){
        getline(ifile,chrom);
        if (ifile.eof()) break;
        if (!isdigit(chrom[2]))
            chromosome.assign(&chrom[1],1);
        else
            chromosome.assign(&chrom[1],2);
        if(chromosome == "X") cnum = 23;
        else if (chromosome == "Y") cnum = 24;
        else if (chromosome == "M") cnum = 25;
        else cnum = atoi(chromosome.c_str());
        getline(ifile,chrom);
        if (ifile.eof()) break;
        size = chrom.size();
        i=0;
        while (i < size - 22){
            if(chrom[i] == 'C' && chrom[i+1] == RevPAM){
                insert_entry = true;
                //Make the key
                end = i + 6;
                key=0;
                twobit=0;
                j=i+22;
                while(j > end){
                    if (chrom[j] == 'N'){
                        insert_entry = false;
                        break;
                    }
                    twobit = cTOintREV[(chrom[j])-ASCIIA];
                    key = ( (key << 2) | twobit ) & keymask;
                    j--;
                }
                if(insert_entry){
                    //Make the value
                    twobit = 0;
                    rest4 =0;
                    //Fix PAMID
                    value = PAMID;
                    //Fix strand
                    value = ((value << STRANDBITS) | 1) & mask;
                    //Fix last nt before GG
                    value = ((value << LASTNTBITS) | cTOintREV[(chrom[i+2])-ASCIIA]) & mask;
                    //Fix chromosome
                    value = ((value << CHROMBITS) | cnum) & mask;
                    //Fix position
                    value = ((value << POSBITS) | i) & mask;
                    //Fix last4 nt
                    end = j;
                    for (j=end;j>end-4;j--){
                        if (chrom[j] == 'N'){
                            insert_entry = false;
                            break;
                        }
                        twobit = cTOintREV[(chrom[j])-ASCIIA];
                        rest4 = ((rest4 << 2) | twobit) & mask;
                    }
                    if(insert_entry){
                        value = ((value << EXTENSIONBITS) | rest4) & mask;
                        //count newly constructed pair to the structure
                        tableValues[hashTable[key] + currentPos[key]] = value;
                        currentPos[key]++;
                    }
                }
            }
            i++;
        }
    }
    ifile.close();
}

/***********************************************************************************************************/
/* This function populates table B. Table A is already populated by the results of the count functions.    */
/* The input is the filename of the genome file, the PAM (NNNNACA), a table with the positions that are    */
/* currently occupied, table A, table B (to be filled), and a number (0-3).                                */
/* Essentially the function looks for hits exactly like the count functions above, but everytime one is    */
/* found instead of increasing a counter, the first 16 bases of the 20mer are looked up on table A and then*/
/* the result is inserted on the position indicated by table_A[16mer] + currentPos[16mer]. This temporary  */
/* table is identical with table A in structure but instead of the total number of hits it holds the number*/
/* of the hits of that specific 16mer that we already had inserted on table B. Hence, we know the exact    */
/* position to place the hit. The number is an ID that's not hardcoded to make the addition/removal of PAMS*/
/* simpler. It is added to each entry and used to identify the PAM when presenting the results to the user.*/
/***********************************************************************************************************/
void resultsNNNNACA(const char *filename, unsigned int *currentPos, MYTYPE *tableValues, unsigned int *hashTable, unsigned int PAMID){
    unsigned int i,j,end,size,key,cnum;
    string chrom(MAX_CHROM_LEN,'\0');
    string chromosome;
    MYTYPE twobit,value,rest4;
    bool insert_entry;
    
    
    //Open the input file for reading
    ifstream ifile (filename,std::ifstream::in);
    
    //for (c=0;c<MAX_CHROM_NUM;c++){
    while (!ifile.eof()){
        getline(ifile,chrom);
        if (ifile.eof()) break;
        if (!isdigit(chrom[2]))
            chromosome.assign(&chrom[1],1);
        else
            chromosome.assign(&chrom[1],2);
        if(chromosome == "X") cnum = 23;
        else if (chromosome == "Y") cnum = 24;
        else if (chromosome == "M") cnum = 25;
        else cnum = atoi(chromosome.c_str());
        getline(ifile,chrom);
        if (ifile.eof()) break;
        size = chrom.size();
        i=0;
        while (i < size - 26){
            if(chrom[i+24] == 'A' && chrom[i+25] == 'C' && chrom[i+26] == 'A'){
                insert_entry = true;
                //Make the key
                end = i + 16;
                key=0;
                twobit=0;
                j=i;
                while(j < end){
                    if (chrom[j] == 'N'){
                        insert_entry = false;
                        break;
                    }
                    twobit = cTOint[(chrom[j])-ASCIIA];
                    key = ( (key << 2) | twobit ) & keymask;
                    j++;
                }
                if (insert_entry){
                    //Make the value
                    twobit = 0;
                    value=0;
                    rest4 =0;
                    //Fix last 3 nt before ACA (4th below)
                    value = cTOint[(chrom[i+21])-ASCIIA];
                    value = ((value << LASTNTBITS) | cTOint[(chrom[i+22])-ASCIIA]);
                    value = ((value << LASTNTBITS) | cTOint[(chrom[i+23])-ASCIIA]);
                    //Fix PAMID
                    value = ((value << LASTNTBITS) | PAMID);
                    //Fix strand (shift one - already 0)
                    value = value << STRANDBITS;
                    //Fix last nt before GG
                    value = ((value << LASTNTBITS) | cTOint[(chrom[i+20])-ASCIIA]);
                    //Fix chromosome
                    value = ((value << CHROMBITS) | cnum) & mask;
                    //Fix position
                    value = ((value << POSBITS) | i) & mask;
                    //Fix last4 nt
                    end = j;
                    for (j=end;j<end+4;j++){
                        if (chrom[j] == 'N'){
                            insert_entry = false;
                            break;
                        }
                        twobit = cTOint[(chrom[j])-ASCIIA];
                        rest4 = ((rest4 << 2) | twobit) & mask;
                    }
                    if (insert_entry){
                        value = ((value << EXTENSIONBITS) | rest4) & mask;
                        //Add newly constructed pair to the structure
                        tableValues[hashTable[key] + currentPos[key]] = value;
                        currentPos[key]++;
                    }
                }
            }
            i++;
        }
    }
    
    
    //move file pointer to the beggining
    ifile.clear();
    ifile.seekg(0, ios::beg);
    
    //for (c=0;c<MAX_CHROM_NUM;c++){
    while (!ifile.eof()){
        getline(ifile,chrom);
        if (ifile.eof()) break;
        if (!isdigit(chrom[2]))
            chromosome.assign(&chrom[1],1);
        else
            chromosome.assign(&chrom[1],2);
        if(chromosome == "X") cnum = 23;
        else if (chromosome == "Y") cnum = 24;
        else if (chromosome == "M") cnum = 25;
        else cnum = atoi(chromosome.c_str());
        getline(ifile,chrom);
        if (ifile.eof()) break;
        size = chrom.size();
        i=0;
        while (i < size - 26){
            if(chrom[i] == 'T' && chrom[i+1] == 'G' && chrom[i+2] == 'T'){
                insert_entry = true;
                //Make the key
                end = i + 10;
                key=0;
                twobit=0;
                j=i+26;
                while(j > end){
                    if (chrom[j] == 'N'){
                        insert_entry = false;
                        break;
                    }
                    twobit = cTOintREV[(chrom[j])-ASCIIA];
                    key = ( (key << 2) | twobit ) & keymask;
                    j--;
                }
                if(insert_entry){
                    //Make the value
                    twobit = 0;
                    rest4 =0;
                    //Fix last 3 nt before ACA (4th below)
                    value = cTOintREV[(chrom[i+5])-ASCIIA];
                    value = ((value << LASTNTBITS) | cTOintREV[(chrom[i+4])-ASCIIA]);
                    value = ((value << LASTNTBITS) | cTOintREV[(chrom[i+3])-ASCIIA]);
                    //Fix PAMID
                    value = ((value << LASTNTBITS) | PAMID);
                    //Fix strand
                    value = ((value << STRANDBITS) | 1);
                    //Fix last nt before GG
                    value = ((value << LASTNTBITS) | cTOintREV[(chrom[i+6])-ASCIIA]) & mask;
                    //Fix chromosome
                    value = ((value << CHROMBITS) | cnum) & mask;
                    //Fix position
                    value = ((value << POSBITS) | i) & mask;
                    //Fix last4 nt
                    end = j;
                    for (j=end;j>end-4;j--){
                        if (chrom[j] == 'N'){
                            insert_entry = false;
                            break;
                        }
                        twobit = cTOintREV[(chrom[j])-ASCIIA];
                        rest4 = ((rest4 << 2) | twobit) & mask;
                    }
                    if (insert_entry){
                        value = ((value << EXTENSIONBITS) | rest4) & mask;
                        tableValues[hashTable[key] + currentPos[key]] = value;
                        currentPos[key]++;
                    }
                }
            }
            i++;
        }
    }
    ifile.close();
}

/***********************************************************************************************************/
/* This function populates table B. Table A is already populated by the results of the count functions.    */
/* The input is the filename of the genome file, the PAM (NNGRRT), a table with the positions that are     */
/* currently occupied, table A, table B (to be filled), and a number (0-3).                                */
/* Essentially the function looks for hits exactly like the count functions above, but everytime one is    */
/* found instead of increasing a counter, the first 16 bases of the 20mer are looked up on table A and then*/
/* the result is inserted on the position indicated by table_A[16mer] + currentPos[16mer]. This temporary  */
/* table is identical with table A in structure but instead of the total number of hits it holds the number*/
/* of the hits of that specific 16mer that we already had inserted on table B. Hence, we know the exact    */
/* position to place the hit. The number is an ID that's not hardcoded to make the addition/removal of PAMS*/
/* simpler. It is added to each entry and used to identify the PAM when presenting the results to the user.*/
/***********************************************************************************************************/
void resultsNNGRRT(const char *filename, unsigned int *currentPos, MYTYPE *tableValues, unsigned int *hashTable, unsigned int PAMID){
    unsigned int i,j,end,size,key,cnum;
    string chrom(MAX_CHROM_LEN,'\0');
    string chromosome;
    MYTYPE twobit,value,rest4;
    bool insert_entry;
    
    //Open the input file for reading
    ifstream ifile (filename,std::ifstream::in);

    //for (c=0;c<MAX_CHROM_NUM;c++){
    while (!ifile.eof()){
        getline(ifile,chrom);
        if (ifile.eof()) break;
        if (!isdigit(chrom[2]))
            chromosome.assign(&chrom[1],1);
        else
            chromosome.assign(&chrom[1],2);
        if(chromosome == "X") cnum = 23;
        else if (chromosome == "Y") cnum = 24;
        else if (chromosome == "M") cnum = 25;
        else cnum = atoi(chromosome.c_str());
        getline(ifile,chrom);
        if (ifile.eof()) break;
        size = chrom.size();
        i=0;
        while (i < size - 25){
            if(chrom[i+22] == 'G' && (chrom[i+23] == 'A' || chrom[i+23] == 'G') && (chrom[i+24] == 'A' || chrom[i+24] == 'G') && chrom[i+25] == 'T'){
                insert_entry = true;
                //Make the key
                end = i + 16;
                key=0;
                twobit=0;
                j=i;
                while(j < end){
                    if (chrom[j] == 'N'){
                        insert_entry = false;
                        break;
                    }
                    twobit = cTOint[(chrom[j])-ASCIIA];
                    key = ( (key << 2) | twobit ) & keymask;
                    j++;
                }
                if (insert_entry){
                    //Make the value
                    twobit = 0;
                    value=0;
                    rest4 =0;
                    //Fix last 3 nt before ACA (4th below)
                    value = cTOint[(chrom[i+21])-ASCIIA];
                    value = ((value << LASTNTBITS) | cTOint[(chrom[i+23])-ASCIIA]);
                    value = ((value << LASTNTBITS) | cTOint[(chrom[i+24])-ASCIIA]);
                    //Fix PAMID
                    value = ((value << LASTNTBITS) | PAMID);
                    //Fix strand (shift one - already 0)
                    value = value << STRANDBITS;
                    //Fix last nt before GG
                    value = ((value << LASTNTBITS) | cTOint[(chrom[i+20])-ASCIIA]);
                    //Fix chromosome
                    value = ((value << CHROMBITS) | cnum) & mask;
                    //Fix position
                    value = ((value << POSBITS) | i) & mask;
                    //Fix last4 nt
                    end = j;
                    for (j=end;j<end+4;j++){
                        if (chrom[j] == 'N'){
                            insert_entry = false;
                            break;
                        }
                        twobit = cTOint[(chrom[j])-ASCIIA];
                        rest4 = ((rest4 << 2) | twobit) & mask;
                    }
                    if (insert_entry){
                        value = ((value << EXTENSIONBITS) | rest4) & mask;
                        //Add newly constructed pair to the structure
                        tableValues[hashTable[key] + currentPos[key]] = value;
                        currentPos[key]++;
                    }
                }
            }
            i++;
        }
    }
    
    //move file pointer to the beggining
    ifile.clear();
    ifile.seekg(0, ios::beg);
    
    //for (c=0;c<MAX_CHROM_NUM;c++){
    while (!ifile.eof()){
        getline(ifile,chrom);
        if (ifile.eof()) break;
        if (!isdigit(chrom[2]))
            chromosome.assign(&chrom[1],1);
        else
            chromosome.assign(&chrom[1],2);
        if(chromosome == "X") cnum = 23;
        else if (chromosome == "Y") cnum = 24;
        else if (chromosome == "M") cnum = 25;
        else cnum = atoi(chromosome.c_str());
        getline(ifile,chrom);
        if (ifile.eof()) break;
        size = chrom.size();
        i=0;
        while (i < size - 25){
            if(chrom[i] == 'A' && (chrom[i+1] == 'T' || chrom[i+1] == 'C') && (chrom[i+2] == 'T' || chrom[i+2] == 'C') && chrom[i+3] == 'C'){
                insert_entry = true;
                //Make the key
                end = i + 9;
                key=0;
                twobit=0;
                j=i+25;
                while(j > end){
                    if (chrom[j] == 'N'){
                        insert_entry = false;
                        break;
                    }
                    twobit = cTOintREV[(chrom[j])-ASCIIA];
                    key = ( (key << 2) | twobit ) & keymask;
                    j--;
                }
                if(insert_entry){
                    //Make the value
                    twobit = 0;
                    rest4 =0;
                    //Fix last 3 nt before ACA (4th below)
                    value = cTOintREV[(chrom[i+4])-ASCIIA];
                    value = ((value << LASTNTBITS) | cTOintREV[(chrom[i+2])-ASCIIA]);
                    value = ((value << LASTNTBITS) | cTOintREV[(chrom[i+1])-ASCIIA]);
                    //Fix PAMID
                    value = ((value << LASTNTBITS) | PAMID);
                    //Fix strand
                    value = ((value << STRANDBITS) | 1);
                    //Fix last nt before GG
                    value = ((value << LASTNTBITS) | cTOintREV[(chrom[i+5])-ASCIIA]) & mask;
                    //Fix chromosome
                    value = ((value << CHROMBITS) | cnum) & mask;
                    //Fix position
                    value = ((value << POSBITS) | i) & mask;
                    //Fix last4 nt
                    end = j;
                    for (j=end;j>end-4;j--){
                        if (chrom[j] == 'N'){
                            insert_entry = false;
                            break;
                        }
                        twobit = cTOintREV[(chrom[j])-ASCIIA];
                        rest4 = ((rest4 << 2) | twobit) & mask;
                    }
                    if (insert_entry){
                        value = ((value << EXTENSIONBITS) | rest4) & mask;
                        tableValues[hashTable[key] + currentPos[key]] = value;
                        currentPos[key]++;
                    }
                }
            }
            i++;
        }
    }
    ifile.close();

}

/***********************************************************************************************************/
/* The input of this program must be a genome file which should have one line per chromosome that is       */
/* preceeded by a header line that contains only or starts with the chromosome name. The output is 2 files */
/* that are used by load memory to load all hits of all gRNAs of all PAMs into memory. The structure can be*/
/* queried efficiently by Results.o once loaded on memory.                                                 */
/***********************************************************************************************************/
int main (int argc, char *argv[]){
    char *usage = (char *) "\n\nUsage: %s -i <genome_filename> -o <output_full_path>\n\n";
    
    //Parse input
    int a;
    char *genome_file = NULL;
    char *out_file_path = NULL;
    while ((a = getopt (argc, argv, "i:o:")) != -1){
        if(a == 'i'){
            genome_file = optarg;
        }
        else if (a == 'o'){
            out_file_path = optarg;
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

    if (genome_file == NULL || out_file_path==NULL) {
        fprintf(stderr, usage, argv[0]);
        fflush(stderr);
        exit(1);
    }
    
    
    //Allocate space for the counters
    unsigned int *currentPos = (unsigned int*) calloc (KEY_NUM,sizeof(unsigned int));
    
    //Count all the results and the results per 16mer
    cout << "\nCounting results per PAM. Please wait." << endl;
    unsigned int total_len=0;
    total_len = countNNG(genome_file,'G', currentPos);
    cout << "1/4: Counted for PAM = NGG." << endl;
    total_len += countNNNNACA(genome_file, currentPos);
    cout << "2/4: Counted for PAM = NNNNACA." << endl;
    total_len += countNNGRRT(genome_file, currentPos);
    cout << "3/4: Counted for PAM = NNGRRT." << endl;
    total_len += countNNG(genome_file,'A', currentPos);
    cout << "4/4: Counted for PAM = NAG." << endl;
    
    //Allocate space equal to the total number of possible 16mers for the index table (table A)
    unsigned int *hashTable = (unsigned int*) calloc (KEY_NUM,sizeof(unsigned int));
    
    //Fill the index table A with the counters per 16mer
    hashTable[0]=0;
    unsigned long v;
    for(v=1;v<KEY_NUM;v++){
        hashTable[v] = hashTable[v-1] + currentPos[v-1];
    }
    
    //Allocate space for the result table
    MYTYPE *tableValues = new MYTYPE [total_len];
    
    //re allocating the counters (so that they'll be 0's). This will keep the
    //number of results put in table B at every step of the populating functions
    free(currentPos);
    currentPos = (unsigned int*) calloc (KEY_NUM,sizeof(unsigned int));
    
    //Fill table B with the results
    cout << "\nCalculating results per PAM. Please wait." << endl;
    resultsNNG(genome_file,'G',currentPos,tableValues,hashTable,0);
    cout << "1/4: Calculated for PAM = NGG." << endl;
    resultsNNNNACA(genome_file,currentPos,tableValues,hashTable,1);
    cout << "2/4: Calculated for PAM = NNNNACA." << endl;
    resultsNNGRRT(genome_file,currentPos,tableValues,hashTable,2);
    cout << "3/4: Calculated for PAM = NNGRRT." << endl;
    resultsNNG(genome_file,'A',currentPos,tableValues,hashTable,3);
    cout << "4/4: Calculated for PAM = NAG." << endl;
    
    //Write index and result tables on binary files
    cout << "\nWriting tables on input.bin and data.bin" << endl;
    FILE *f1,*f2;
    char index_file[strlen(out_file_path) + strlen("/index.bin") + 1], data_file[strlen(out_file_path) + strlen("/data.bin") + 1];
    strcpy(index_file,out_file_path);
    strcpy(data_file,out_file_path);
    if (out_file_path[strlen(out_file_path) - 1] == '/'){
        f1 = fopen ( strcat(index_file,"index.bin"), "wb" );
        f2 = fopen ( strcat(data_file,"data.bin"), "wb" );
    }
    else{
        f1 = fopen ( strcat(index_file,"/index.bin"), "wb" );
        f2 = fopen ( strcat(data_file,"/data.bin"), "wb" );
    }
    fwrite(hashTable,sizeof(unsigned int),KEY_NUM,f1);
    fwrite(tableValues,sizeof(MYTYPE),total_len,f2);
    fclose(f1);
    fclose(f2);
    
    cout << "\nTables created.\n" << endl;
    
    return 1;
}

