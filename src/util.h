//
// Created by parthajit on 11/7/20.
//

#ifndef CPPMET_UTIL_H
#define CPPMET_UTIL_H
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <time.h>

char pro_amino[20][2][10]={{"ALA","HYPHOB"}, {"ARG","BASIC"}, {"ASN","POLNUTRL"}, {"ASP", "ACIDIC"}, {"CYS",
                                  "POLNUTRL"},{"GLN", "POLNUTRL"}, {"GLU", "ACIDIC"}, {"GLY","SPL"}, {"HIS","BASIC"}, {"ILE", "HYPHOB"},
        {"LEU", "HYPHOB"}, {"LYS", "BASIC"}, {"MET", "HYPHOB"}, {"PHE","HYPHOB"}, {"PRO","SPL"},
        {"SER", "POLNUTRL"}, {"THR", "POLNUTRL"}, {"TRP","HYPHOB"}, {"TYR", "HYPHOB"}, {"VAL", "HYPHOB"}};

int pro_amino_size = 20;

int find_amino(char* amino){
    for (int i = 0; i < pro_amino_size; ++i) {
        if(strcmp(amino, pro_amino[i][0]) == 0 ) return i;
    }
    return -1;
}

int find_nucleic(char* nucleic){
    if(strcmp(nucleic, "G") == 0) return 0;
    if(strcmp(nucleic, "A") == 0) return 0;
    if(strcmp(nucleic, "C") == 0) return 0;
    if(strcmp(nucleic, "T") == 0) return 0;
    if(strcmp(nucleic, "U") == 0) return 0;
    return -1;
}
int is_HOH(char* res){
    if(strcmp(res, "HOH") == 0) return 1;
    return -1;
}

char get_res_loc(Atom* atom, char restype){
    if(restype == 'W'){
        return 'W';
    }
    if(restype == 'P'){
        return 'P';
    }
    if(restype == 'M'){
        return 'M';
    }
    if(restype == 'N'){
        long len = strlen(atom->loc);
        char symb = atom->loc[0];
        if(len == 2){ //Then it is Nucleobase for Cor dataset
            if(symb == 'O' || symb == 'N' || symb == 'C') {
                return 'N';
            }
        }else if(len == 3 && (atom->loc[2] == '*' || atom->loc[2] == '\'')){ // sugar backbone for cor file
            if(symb == 'O'    || symb == 'C' ){
                return 'S';
            }
        }else{ // phosphate case
            if((len == 3 && symb == 'O' ) ||(len == 1 && symb == 'P')){
                return 'P';
            }
        }
    }
    //return 'X';
    fprintf(stderr, "Error.... Wrong restype found \n");
    atom->print();
    exit(EXIT_FAILURE);
}


void file_name_split(char* file_path, char* file_name, char* file_ext, char* src_file){
    char* ext = strrchr(src_file,'.');
    assert(ext != NULL && ext != src_file);
    strcpy(file_ext, ext);
    char* file_sep = strrchr(src_file, '/');

    if(file_sep == NULL){
        file_path[0] ='\0' ;
        strncpy(file_name, src_file, (size_t)(ext - src_file));
        file_name[(int)(ext - src_file)] = '\0';
    }else{
        strncpy(file_path, src_file, (size_t)(file_sep-src_file)+1);
        file_path[(int)(file_sep-src_file)+1] = '\0';
        strncpy(file_name, file_sep+1, (size_t)(ext-file_sep)-1);
        file_name[(int)(ext-file_sep)-1] = '\0';
    }
}

void file_name_join(char* joined_file_name, const char* path, const char* file_name, const char* ext){
    strcpy(joined_file_name, path);
    strcat(joined_file_name, file_name);
    strcat(joined_file_name, ext);
}

void now(char output[]){
      time_t     rawtime;
      struct tm* timeinfo;
      time( &rawtime );
      timeinfo = localtime( &rawtime );
      sprintf(output, "%d:%d:%d", timeinfo->tm_hour,timeinfo->tm_min, timeinfo->tm_sec);
}
void today(char output[]){
      time_t     rawtime;
      struct tm* timeinfo;

      time( &rawtime );
      timeinfo = localtime( &rawtime );
      sprintf(output, "%d-%d-%d", timeinfo->tm_mday,timeinfo->tm_mon, timeinfo->tm_year+1900);
}
#endif //CPPMET_UTIL_H
