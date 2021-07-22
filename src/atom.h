//
// Created by parthajit on 10/7/20.
//

#ifndef CPPMET_ATOM_H
#define CPPMET_ATOM_H

#include <stdio.h>
#include <cmath>
#include <cstring>


class Atom {
public:
    long id;
    char loc[5];
    char resname[5];
    long resid;
    char chain[5];
    double x;
    double y;
    double z;
    double occu;
    char type; // A for Atom H for hetatm.
public:
    Atom();
    Atom(Atom* a);

    void fprint_charmm_format(FILE* fp, long corid){
        fprintf(fp, "ATOM %6ld", id);
        char atom_loc[6];
        strcpy(atom_loc, this->loc);
        long len = strlen(atom_loc);
        for(int i=0; i<len; ++i){
            if(atom_loc[i] == '*'){
                atom_loc[i] = '\'';
            }
        }
        if(strcmp(atom_loc, "OP1") == 0){
            strcpy(atom_loc, "O1P");
        }else if(strcmp(atom_loc,"OP2") == 0){
            strcpy(atom_loc, "O2P");
        }else if(strcmp(atom_loc,"OP3") == 0){
            strcpy(atom_loc, "O3P");
        }

        fprintf(fp, "  %-3s %3s %5ld    %8.3lf%8.3lf%8.3lf %5.2lf\n",
                atom_loc,
                resname,
                corid,
                x,
                y,
                z,
                occu);
    }

    double angl_deg(Atom* b, Atom* c){
        double rad = this->angle_rad(b, c);
        return rad * (180.0/3.14159);
    }
    void print(){
        printf("%6ld %5s %5s %lf %lf %lf \n", id, loc, resname, x, y, z);
    }
    double angle_rad(Atom* b, Atom* c);
    double distsqr(Atom* a);
    void copy(Atom* a);

};

double Atom::angle_rad(Atom* b, Atom* c) {
    Atom* a = this;
    double x = (a->x - b->x);
    double y = (a->y - b->y);
    double z = (a->z - b->z);
    double ab = sqrt(x*x + y*y + z*z);

    x = (b->x - c->x);
    y = (b->y - c->y);
    z = (b->z - c->z);
    double bc = sqrt(x*x + y*y + z*z);

    x = (c->x - a->x);
    y = (c->y - a->y);
    z = (c->z - a->z);
    double ac = sqrt(x*x + y*y + z*z);

    double angl =  acos((ab*ab + ac*ac - bc*bc)/(2.0* ab * ac));
    return angl;
}

Atom::Atom(Atom* a) {
    id    = a->id;
    resid = a->resid;
    type  = a->type;
    x     = a->x;
    y     = a->y;
    z     = a->z;
    occu  = a->occu;
    strcpy(this->resname, a->resname);
    strcpy(this->loc, a->loc);
    strcpy(this->chain, a->chain);
}

double Atom::distsqr(Atom * b) {
    Atom * a = this;
    return (a->x - b->x)*(a->x - b->x) + (a->y - b->y)*(a->y - b->y) + (a->z - b->z)*(a->z - b->z);
}

Atom::Atom() {
    id = -999;

}

void Atom::copy(Atom *a) {
    this->id    = a->id;
    this->resid = a->resid;
    this->type  = a->type;
    this->x     = a->x;
    this->y     = a->y;
    this->z     = a->z;
    this->occu  = a->occu;
    strcpy(this->resname, a->resname);
    strcpy(this->loc, a->loc);
    strcpy(this->chain, a->chain);
}


#endif //CPPMET_ATOM_H
