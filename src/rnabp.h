//
// Created by parthajit on 11/7/20.
//

#ifndef CPPMET_RNABP_H
#define CPPMET_RNABP_H

#include <cstdio>
#include <cassert>
#include <cctype>
#include <cstring>
#include <cstdlib>
#include "atom.h"


class basepair {
public:
    long corid;
    long cifid;
    char resname[4];
    char chain[4];
    char name[5];
    char type[3];
    double eval;
    struct basepair* bp[4];
    long numbp;
public:
    void basepair_free(){

	  
	  for(int i=0; i<numbp; ++i){
		free(bp[i]);

		
		bp[i] = NULL;
	  }
    }
    void fprint_bp_short(FILE* fp, char is_target, Atom* met, Atom* water, char sec_seq, char* loc){
        char str[3] = "> ";
        if(water != NULL){
            strcpy(str,"W-");
        }
        char blank[200] = "               ";

        if(numbp <= 0){
            if(is_target == 'T'){
                    fprintf(fp, "MOTF  %6ld  %-3s %2s-%2s %c %3s-           %-3s        %6ld\n",met->resid, met->chain,
                            met->resname,
                            str,
                            sec_seq,
                            resname, loc, corid);
            }else{

                fprintf(fp,"MOTF     %s%c %3s-\n", blank, sec_seq,resname);

            }
            return;
        }else{
            if(is_target == 'T'){
                if(numbp > 1){
                    fprintf(fp, "MOTF  %6ld  %-3s %2s-%2s %c %3s-%-3s %s   %-3s  +%ld    %6ld\n",met->resid,
                            met->chain,
                            met->resname,
                            str,
                            sec_seq,
                            resname,
                            bp[0]->resname,
                            bp[0]->name, loc, numbp-1, corid);
                }else{
                    fprintf(fp, "MOTF  %6ld  %-3s %2s-%2s %c %3s-%-3s %s   %-3s        %6ld\n",met->resid, met->chain,
                            met->resname,
                            str,
                            sec_seq,
                            resname,
                            bp[0]->resname,
                            bp[0]->name, loc,
                            corid);
                }

            }else{

                fprintf(fp, "MOTF     %s%c %3s-%-3s %s \n",blank, sec_seq, resname, bp[0]->resname, bp[0]->name);

            }
        }
    }



    void fprint_bp(FILE* fp){
        if(numbp <= 0) return;
        fprintf(fp, "%6ld  %-3s",
                cifid,
                chain);
        for(int i=0; i<numbp; ++i){
            fprintf(fp, "  %6ld  %-3s     %s:%s-%s  %6.2f  ",
                    bp[i]->cifid,
                    bp[i]->chain,
                    resname,
                    bp[i]->resname,


                    bp[i]->name, bp[i]->eval);

        }
        fprintf(fp, "\n");

    }
};




class rnabp{
public:
    basepair* bp;
    long nres;
    bool atm;
    bool hetatm;
public:
    rnabp();
    rnabp(char* outfile);
    void fprint_bp(FILE* fp){

       for(long i=0; i<nres; ++i){
           char lb = ' ';
           char rb = ' ';
           char ori = ' ';
           if(bp[i].numbp > 0){
               if(bp[i].bp[0]->name[0] == 'W') lb = '-';
               else if(bp[i].bp[0]->name[0] == 'S') lb = '~';
               else if(bp[i].bp[0]->name[0] == 'H') lb = '+';
               else lb = bp[i].bp[0]->name[0] ;

               if(bp[i].bp[0]->name[3] == 'T') ori = '.';

               if(bp[i].bp[0]->name[2] == 'W') rb = '-';
               else if(bp[i].bp[0]->name[2] == 'S') rb = '~';
               else if(bp[i].bp[0]->name[2] == 'H') rb = '+';
               else rb = bp[i].bp[0]->name[2];


               fprintf(fp, "        %s%c%c%s %c", this->bp[i].resname,
                                            lb,
                                            rb,
                                            bp[i].bp[0]->resname,
                                            ori);
               if(bp[i].numbp > 1){
                   fprintf(fp, "  +%ld", bp[i].numbp-1);
               }
               fprintf(fp, "\n");
           }else{
               fprintf(fp, "        %s%c%c \n", this->bp[i].resname,
                       lb,
                       rb);
           }

       }
    }


    ~rnabp(){
        if(bp != NULL){
	      for(int i=0; i<nres; ++i){
		    bp[i].basepair_free();
	      }
            free(bp);
        }else{
            ; //printf("RNABP not freed\n");
        }
    }
};

rnabp::rnabp(char* outfile){
    FILE* fp = fopen(outfile, "r");
   if(fp == NULL){    /* Exception Handling */ 
	fprintf(stderr, "Error in function %s(). (File: %s, Line %d)... Cannot open out file.\n", __func__, __FILE__, __LINE__);
       exit(EXIT_FAILURE);
   }       
    
    char line[256];
    int flag = 0;
    while(fgets(line, 256, fp) != NULL){
        if(strncmp(line, "#HEADER   Choise of HETATM entries:", 34) == 0){
            if(strncmp(line+34,"YES",3) ==0){
                this->hetatm = true;
            }
            continue;
        }
        if(strncmp(line, "#HEADER   Cleaned number", 24) != 0) continue;
        this->nres = atoi(line+38);
        flag = 1;
//        printf("NUMBER OF BP = %ld\nHETATM =%d \n",this->nres, this->hetatm);
        break;
    }
    if(this->nres == 0){
        fprintf(stderr, "No RNA found\n");
	this->bp = NULL;
        fclose(fp);
        return;
    }
    this->bp = (basepair*) malloc(this->nres * sizeof(basepair));
    if(flag == 0){
        fprintf(stderr, "Error in .out file could'nt locate Cleaned residue number\n");
        fclose(fp);
        exit(EXIT_FAILURE);
    }


    while(fgets(line, 256, fp) != NULL){
        if(line[0] == '#') continue;
        break; /* pass the lines */
    }

    int count = 0;
    while(line != NULL){
        if(strlen(line)<=6) break;
        basepair* bp = this->bp+count;
        count ++;
        char sep[] = "\t \n";
        char *token;
        bp->numbp = 0;

        token = strtok(line, sep);
        bp->corid = atoi(token);

        token = strtok(NULL, sep);
        bp->cifid = atoi(token);

        token = strtok(NULL, sep);
        strcpy(bp->resname, token);

        token = strtok(NULL, sep); /* skip ? token */

        token = strtok(NULL, sep);
        strcpy(bp->chain, token);

        /* For base pairs */
        while((token = strtok(NULL, sep)) != NULL){
            basepair* bp1 = (basepair*) malloc(sizeof(basepair));
            bp->bp[bp->numbp] = bp1;
            bp->numbp++;
            bp1->corid = atoi(token);

            token = strtok(NULL, sep);
            bp1->cifid = atoi(token);

            token = strtok(NULL, sep);
            strcpy(bp1->resname, token);

            token = strtok(NULL, sep); /* skip ? token */

            token = strtok(NULL, sep);
            strcpy(bp1->chain, token);

            token = strtok(NULL, sep);
            strcpy(bp1->name, token);

            token = strtok(NULL, sep);
            strcpy(bp1->type, token);

            token = strtok(NULL, sep);
            bp1->eval = atof(token);
        }
        fgets(line, 256, fp);

    }
    fclose(fp);
}


rnabp::rnabp() {
    bp = NULL;
}


#endif //CPPMET_RNABP_H
