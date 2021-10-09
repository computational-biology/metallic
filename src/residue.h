//
// Created by parthajit on 10/7/20.
//

#ifndef CPPMET_RESIDUE_H
#define CPPMET_RESIDUE_H

#include "atom.h"

class Residue {
public:
    Atom atom[40];
    long size;
public:
    Residue();
    void add(Atom* atom);
    void copy(Residue* residue);
    void fprint_charmm_format(FILE* fp, long corid);
    ~Residue();
};



void Residue::add(Atom* atm) {
    if(size>=40){
        fprintf(stderr, "Error.... too many atoms in residue %s for atom id  %ld\n",this->atom[0].resname, this->atom[0].id);
        exit(1);
    }


    this->atom[size].type = atm->type;
    this->atom[size].x = atm->x;
    this->atom[size].y = atm->y;
    this->atom[size].z = atm->z;
    this->atom[size].id = atm->id;
    this->atom[size].resid = atm->resid;
    this->atom[size].occu = atm->occu;
    strcpy(this->atom[size].chain, atm->chain);
    strcpy(this->atom[size].resname, atm->resname);
    strcpy(this->atom[size].loc, atm->loc);

    this->size ++;
}

Residue::Residue() {
    size = 0;
}


void Residue::copy(Residue *residue) {
    this->size = residue->size;
    for(long i=0; i<size; ++i){
        this->atom[i].copy(residue->atom+i);
    }
}

Residue::~Residue() {
    ;//fprintf(stderr, "Residue delete invoked\n");
}

void Residue::fprint_charmm_format(FILE *fp, long corid) {
    for(long i=0; i<size; ++i){
        this->atom[i].fprint_charmm_format(fp, corid);
    }
}


#endif //CPPMET_RESIDUE_H
