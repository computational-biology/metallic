//
// Created by parthajit on 10/7/20.
//

#ifndef CPPMET_MOLECULE_H
#define CPPMET_MOLECULE_H

#include <iostream>

#include "atom.h"
#include "residue.h"

using namespace std;


class Molecule {

public:
    long max;
    Residue *residue;
    long size;
    bool atm;
    bool hetatm;
public:
    void reset(bool atm, bool hetatm);
    Molecule(char* corfile, rnabp* rnabp);
    Molecule(long maxsize, bool atm, bool hetatm);
    void scan_cif(char* ciffile, int (*pf)(char*));
    Molecule(Molecule* mol);
    long atom_count(){
        long count = 0;
        for(long i=0; i<size; ++i){
            count = count + residue[i].size;
        }
        return count;
    }
    void get_all_atoms(Atom* atmarray[], long* atmcount){
        long count =0;
        for(long i=0; i<this->size; ++i){
            for(long j=0; j<this->residue[i].size; ++j){
                atmarray[count] = this->residue[i].atom+j;
                count ++;
            }
        }
        *atmcount = count;
    }
    ~Molecule();
};

void Molecule::scan_cif(char *ciffile,  int (*pf)(char *)) {
    //residue = new Residue[size];
    //this->residue = new  Residue[size];
    //this->max = size;
    this ->size = 0;

    
    char tokary[100][100];
    FILE* fp = fopen(ciffile, "r");
    assert(fp != NULL);
    char line[1024];
    int grouploc=-1;  // group_PDB
            int type_symbloc = -1;

    int idloc = -1;   // id
    int locloc = -1;  // label_atom_id
    int resnameloc = -1;  // label_comp_id
    int chainloc = -1; // auth_asym_id
    int residloc = -1; // auth_seq_id
    int modelloc = -1; // pdbx_PDB_model_num
    int occuloc  = -1; //occupancy
    int xloc = -1;
    int yloc = -1;
    int zloc = -1;
    int count = -1;
    while(fgets(line, 1024, fp) != NULL){
        if(strncmp(line, "_atom_site.", 11) != 0 && count == -1) continue;
        if(strncmp(line, "_atom_site.", 11) != 0 && count != -1) break;
        count ++;
        if(strncmp(line,"_atom_site.group_PDB", 20) == 0 && line[20] != '_' ){
            grouploc = count;
        }else if(strncmp(line,"_atom_site.id", 13) == 0 && line[13] != '_' ){
            idloc = count;
        }else if(strncmp(line,"_atom_site.type_symbol", 22) == 0 && line[22] != '_' ){
            type_symbloc = count;
        }else if(strncmp(line,"_atom_site.label_atom_id", 24) == 0 && line[24] != '_' ){
            locloc = count;
        }else if(strncmp(line,"_atom_site.label_comp_id", 24) == 0 && line[24] != '_' ){
            resnameloc = count;
        }else if(strncmp(line,"_atom_site.auth_asym_id", 23) == 0 && line[23] != '_' ){
            chainloc = count;
        }else if(strncmp(line,"_atom_site.auth_seq_id", 22) == 0 && line[22] != '_' ){
            residloc = count;
        }else if(strncmp(line,"_atom_site.Cartn_x", 18) == 0 && line[18] != '_' ){
            xloc = count;
        }else if(strncmp(line,"_atom_site.Cartn_y", 18) == 0 && line[18] != '_' ){
            yloc = count;
        }else if(strncmp(line,"_atom_site.Cartn_z", 18) == 0 && line[18] != '_' ){
            zloc = count;
        }else if(strncmp(line,"_atom_site.occupancy", 20) == 0 && line[20] != '_' ){
            occuloc = count;
        }else if(strncmp(line,"_atom_site.pdbx_PDB_model_num", 29) == 0 && line[29] != '_' ){
            modelloc = count;
        }else{
            ;
        }

    }
    char sep[10] = "\t \"\n";
    char* token;
    long prevresid = -999;
    char prevchain[10] = "-";
    long resindx = -1;
    int modelflag=0;
    int modelnum = -1;
    while(line[0] != '#'){
        int indx=0;
        token = strtok(line, sep);
        while(token != NULL){
//	         printf("%s, ",token);
	      strcpy(tokary[indx], token);
	      indx++;
	      token = strtok(NULL, sep);
        }

        if(tokary[type_symbloc][0] == 'H'){ // We dont consider H
            if(fgets(line, 1024, fp) == NULL) break;
            continue;
        }
        if(pf(tokary[resnameloc]) < 0 ){
            if(fgets(line, 1024, fp) == NULL) break;
            continue;
        }

        long resid = atol(tokary[residloc]);
        char chn[10];
        strcpy(chn, tokary[chainloc]);
	if(modelflag==0){
	      modelnum = atoi(tokary[modelloc]);
	      modelflag = 1;
	}
	if(modelnum != atoi(tokary[modelloc])){
	      fclose(fp);
	      return;
	}
        if(prevresid != resid || strcmp(prevchain, chn) != 0){
            resindx ++;
            this->size ++;
            assert(this->size < this->max );
            prevresid = resid;
            strcpy(prevchain, chn);
        }

        Atom atm;
        char atmhetatm[20];
        strcpy(atmhetatm, tokary[grouploc]);

        atm.id = atol(tokary[idloc]);
        atm.resid = atol(tokary[residloc]);
        strcpy(atm.loc, tokary[locloc]);
        strcpy(atm.resname, tokary[resnameloc]);
        strcpy(atm.chain, tokary[chainloc]);
        if(strncmp(tokary[grouploc], "HETATM", 5) == 0)
            atm.type = 'H';
        else
            atm.type = 'A';

        atm.x = atof(tokary[xloc]);
        atm.y = atof(tokary[yloc]);
        atm.z = atof(tokary[zloc]);
        atm.occu = atof(tokary[occuloc]);

        this->residue[resindx].add(&atm);

        if(fgets(line, 1024, fp) == NULL) break;
        //printf("%s",line);
    }
    fclose(fp);
    //if(this->size>0)
    //printf("Moltest1=%s %lf\n", this->residue[0].atom[0].resname, residue[0].atom[0].x);
    //else{
     //   printf("No Mol found\n");
    //}

}

Molecule::Molecule(char *corfile, rnabp *rnabp) {
    FILE* fp = fopen(corfile, "r");


    if(fp == NULL){    /* Exception Handling */ 
	  fprintf(stderr, "Error in function %s(). (File: %s, Line %d)... Cannot open <accn>_rna.pdb\n", __func__, __FILE__, __LINE__);
	  exit(EXIT_FAILURE);
    }

    this->size = rnabp->nres;
    this->hetatm = rnabp->hetatm;
    this->atm = rnabp->atm;

    this->residue = new Residue[this->size];
    //cout<<"SIZE="<<size<<endl;
    char line[1024];

    char sep[] = "\t \n";
    char *token;
    long int prev_resid = -999;
    long int i= -1;
    char tmptoken[2];
    Atom atom ;//= new Atom();
    while(fgets(line,1024, fp) != NULL){
        if(line[0] == '#') continue;
        long resid;
        if(strncmp(line,"END", 3)==0){
            break;
        }
        token = strtok(line, sep); /* skip "ATOM"*/

        token = strtok(NULL, sep);
        atom.id = atol(token);
        if(line[16] != ' '){
            tmptoken[0] = line[16];
            tmptoken[1] = '\0';
            line[16] = ' ';
            token = strtok(NULL, sep);
            strcpy(atom.loc, token);
            strcat(atom.loc, tmptoken);
           // printf("Found %s in %s\n", atom->loc, corfile);
        }else{

            token = strtok(NULL, sep);
            strcpy(atom.loc, token);
        }
        if(atom.loc[0] == 'H') continue; // We dont consider Hydrogen

        token = strtok(NULL, sep);
        strcpy(atom.resname, token);

        token = strtok(NULL, sep);
        resid = atol(token);
        atom.resid = rnabp->bp[resid-1].cifid;

        char token1[20];
        memset(token1, '\0', 20);
        strncpy(token1, line+30,8);
        atom.x = atof(token1);

        memset(token1, '\0', 20);
        strncpy(token1, line+38,8);
        atom.y = atof(token1);



        memset(token1, '\0', 20);
        strncpy(token1, line+46,8);
        atom.z = atof(token1);
        memset(token1, '\0', 10);
        strncpy(token1, line+54,6);
        atom.occu = atof(token1);
        strcpy(atom.chain, rnabp->bp[resid-1].chain); // cor file does not contain chain info
        atom.type ='A'; // normal atom. not hetatm
        if(prev_resid != resid) {
            prev_resid = resid;
            i++;
        }
        this->residue[i].add(&atom);
    }
    assert(i == this->size-1);
//    long last= residue[i].size-1;
//    long last1= residue[i-1].size-1;
    //printf("z=%ld ",    this->residue[i-1].atom[last1].resid);
    //printf("z=%ld ",    this->residue[i].atom[0].resid);
   // printf(" z=%ld\n", this->residue[i].atom[last].resid);
    fclose(fp);
}


Molecule::~Molecule() {
    delete [] residue;
   // printf("Free invoked for Molecule %ld\n", max);

}

Molecule::Molecule(Molecule *mol) {
    this->size = mol->size;
    this->residue = new Residue[this->size];
    this->max = mol->max;
    this->atm = mol->atm;
    this->hetatm = mol->hetatm;
    for(long i=0; i<this->size; ++i){
        this->residue[i].copy(mol->residue+i);
    }
}

Molecule::Molecule(long maxsize, bool atm, bool hetatm) {
    this->max = maxsize;
    this->size = 0;
    this->residue = new Residue[this->max];
    this->atm = atm;
    this->hetatm = hetatm;
}

void Molecule::reset(bool atm, bool hetatm) {
    for(long i=0; i<this->max; ++i){
        this->residue[i].size = 0;
    }
    this->size = 0;
    this->atm = atm;
    this->hetatm = hetatm;
}


#endif //CPPMET_MOLECULE_H
