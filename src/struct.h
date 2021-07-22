//
// Created by parthajit on 13/7/20.
//

#ifndef CPPMET_SECSTRUCT_H
#define CPPMET_SECSTRUCT_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>



class Structure {
public:
    char* secseq;
    char* dbn;
    char* priseq;
    char chain[200][6];
    long int num_chain;
    long int num_res;
public:
    Structure(char* dat_file, char* dbn_file, long size);
    ~Structure();
};

Structure::Structure(char *dat_file, char *dbn_file, long size) {
    this->num_res = size;
    if(this->num_res <= 0) return;
    this->secseq = (char*) malloc ((num_res+1) * sizeof(char));
    this->dbn = (char*) malloc ((num_res+200+1) * sizeof(char)); // 200 extra to accomodate & (cnain break)
    this->priseq = (char*) malloc ((num_res+200+1) * sizeof(char)); // 200 extra to accomodate & (cnain break)
    this->num_chain = 0;
    FILE* fpdat = fopen(dat_file, "r");
    assert( fpdat != NULL);
    char line[1024];
    char sep[]="\t :\n";
    char* token;
    long int index = 0;
    while(fgets(line, 1024, fpdat) != NULL){
        if(line[0] == '>'){
            token = strtok(line, sep);
            token = strtok(NULL, sep);//second one is the chain info
            strcpy(this->chain[this->num_chain], token);
            this->num_chain++;
        }else{
            token = strtok(line, sep);
            long int len = strlen(token);
            assert(len>0);
            strcpy(this->secseq+index, token);
            index = index + len;
        }
    }
    //fprintf(stdout, "size= %ld: %s\n",this->num_res, this->secseq);
    assert(this->num_res == index);

    fclose(fpdat);
   /* FILE* fpdbn = fopen(dbn_file,"r");
    assert(fpdbn != NULL);

    char* longline = (char*) malloc((this->num_res+200) * sizeof(char));

    fgets(longline, this->num_res+200, fpdbn);

    fgets(longline, this->num_res+200, fpdbn);
    token = strtok(longline,"\n");
    strcpy(this->priseq, token);

    fgets(longline, this->num_res+200, fpdbn);
    token = strtok(longline,"\n");
    strcpy(this->dbn, token);*/

    //fprintf(stdout, "%s\n",secst->dbn);
    //fclose(fpdbn);
}


Structure::~Structure() {
    if(this->num_res > 0){
        free(this->secseq);
        free(this->dbn);
        free(this->priseq);
    }
    //printf("Sec free called\n");
}


#endif //CPPMET_SECSTRUCT_H
