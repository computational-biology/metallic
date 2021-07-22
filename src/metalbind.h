//
// Created by parthajit on 13/7/20.
//

#ifndef CPPMET_METALBIND_H
#define CPPMET_METALBIND_H

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "atom.h"
#include "residue.h"
#include "molecule.h"
#include "struct.h"
#include "metaldefs.h"
#include "rnabp.h"
#include "util.h"

#define MAX_SIZE 90
#define MAX_PROX 400
class Map;




class Proximity{
public:
    Molecule* mol[MAX_PROX];
    long      resindx[MAX_PROX];
    char restype[MAX_PROX];
    long size;
public:
    Proximity(){
        size = 0;
    }
    void add(Molecule* mol, long resindx, char type){
        assert(size < MAX_PROX);
        this->mol[size] = mol;
        this->resindx[size] = resindx;
        this->restype[size] = type;
        size++;
    }
    Residue* get_residue(long index){
        return &mol[index]->residue[resindx[index]];
    }
    Molecule* get_molecule(long index){
        return mol[index];
    }
};


class Ligand{
public:
    Molecule* mol[MAX_SIZE];
    long resindx[MAX_SIZE];
    long offset[MAX_SIZE];
    char loc[MAX_SIZE];      // S for sugar, N for base, P fpr php etc.
    char restype[MAX_SIZE]; //p for protein n for nuc w for water M for met
    char rule[MAX_SIZE];
    double dist[MAX_SIZE];
    long size;
public:
    Ligand(){
        size = 0;
    }
    void add(Proximity* prox, long index, long offset, char loc, char rule, double dist){
        assert(size < MAX_SIZE);
        mol[size] = prox->mol[index];
        resindx[size] = prox->resindx[index];
        this->loc[size] = loc;
        restype[size] = prox->restype[index];
        this->rule[size] = rule;
        this->offset[size] = offset;
        this->dist[size] = dist;
        size++;
    }
};

class Mediated{
public:
    long waterindx[MAX_SIZE];
    long link_count[MAX_SIZE];
    double angle[MAX_SIZE];
    Ligand ligand;
    long size;
public:
    Mediated(){
        size = 0;
        for(long i=0; i<MAX_SIZE; ++i){
            link_count[i] = 0;
        }
    }
    void add(Proximity* prox,  long index, long offset, char loc, char rule, long water_indx, double dist, double
    angle){
        assert(this->size < MAX_SIZE);
        this->waterindx[size] = water_indx;
        this->angle[size] = angle;
        ligand.add(prox, index, offset, loc, rule, dist);
        long max_link = link_count[size];
        this->size++;
        for(long i=0; i<size; ++i){
            if(waterindx[i] == water_indx){
                if(this->link_count[i] > max_link){
                    max_link = link_count[i];
                }
            }
        }
        for(long i=0; i< size; ++i){
            if(waterindx[i] == water_indx){
                link_count[i] = max_link +1;
            }
        }
    }
};

class Site{
public:
    Atom* metal;
    int atomic_no;

    Proximity prox;



    Ligand ligand;
    Mediated wmed;



public:
    Site() {
        this->metal = NULL;
        this->atomic_no = -999;
        this->prox = Proximity();
        this->ligand = Ligand();
    }

    void fprint_summary(FILE* fp){
        int pcount = 0;
        int ncount = 0;
        int wcount = 0;
        int mcount = 0;
        int wmpcount = 0;
        int wmncount = 0;
        int wmwcount = 0;
        int wmmcount = 0;
        for(long i=0; i<ligand.size; ++i){
            if(ligand.restype[i] == 'N') ncount++;
            else if(ligand.restype[i] == 'P') pcount++;
            else if(ligand.restype[i] == 'W') wcount ++;
            else if(ligand.restype[i] == 'M') mcount ++;
        }
        for(long i=0; i<wmed.size; ++i){
            if(wmed.ligand.restype[i] == 'N') wmncount++;
            else if(wmed.ligand.restype[i] == 'P') wmpcount++;
            else if(wmed.ligand.restype[i] == 'W') wmwcount ++;
            else if(wmed.ligand.restype[i] == 'M') wmmcount ++;
            else{
                fprintf(stderr,"Error... Wrong type %c\n", wmed.ligand.restype[i]);
                exit(EXIT_FAILURE);
            }
        }
        fprintf(fp, "SUMM  %6ld  %-3s  %-3s  %3ld =  %3d  %3d  %3d  %3d    %3ld  =  %3d  %3d  %3d  %3d\n",
                metal->resid,
                metal->chain,
                metal->resname,
                ligand.size, ncount, pcount, wcount, mcount,
                wmed.size,
                wmncount,
                wmpcount,
                wmwcount,
                wmmcount);
    }
    void fprint_wmed_basepair(FILE *fp, rnabp *bp, int *flag){
        if(bp == NULL) return;

        //long resindx1;
        //int link = 80;



        /*for(long i=0; i<ligand.size; ++i){
            if(ligand.restype[i] == 'W'){
                resindx1 = ligand.resindx[i];
                if(bp->bp[resindx1].numbp > 0){
                    link++;
                }
            }
        }*/
        int tag  = 0;
        for(long i=0; i< wmed.size; ++i){
            long windx = wmed.waterindx[i];

            Atom* water = ligand.mol[windx]->residue[ligand.resindx[windx]].atom+ ligand.offset[windx];
            if(wmed.ligand.restype[i] != 'N') continue;

            Molecule* mol = wmed.ligand.mol[i];
            long resindx = wmed.ligand.resindx[i];
            long offset  = wmed.ligand.offset[i];
            Atom* atm = mol->residue[resindx].atom + offset;
                if(atm->resid != bp->bp[resindx].cifid){
                    fprintf(stderr, "Error... invalid mapping in RNA and Basepair...\n");
                    exit(EXIT_FAILURE);
                }
                if(bp->bp[resindx].numbp <= 0) continue;
                char location[5] = "---";
                if(wmed.ligand.loc[i] == 'N'){
                    strcpy(location,"NUC");
                }else if(wmed.ligand.loc[i] == 'P'){
                    strcpy(location, "PHP");
                }else if(wmed.ligand.loc[i] == 'S'){
                    strcpy(location, "SUG");
                }

                if(*flag == 0){
                    *flag = 1;
                    fprintf(fp, "\n\n      |  Metal Detail   | Water Detsils   |      Base Pair Details             |"
                                "  outcome      "
                                "    |\n");
        fprintf(fp,      "      "
                         "----------------------------------------------------------------------------------------------\n");
        fprintf(fp, "       resid  chn  mtl    resid  chn  res  loc  atm   resid  chn   resid  chn     pair------>> "
                    "\n");
        fprintf(fp, "      "
            "----------------------------------------------------------------------------------------------\n");


                }
                tag = 1;
                fprintf(fp, "WMBP  %6ld  %-3s  %-3s   %6ld  %-3s  %-3s  %-3s  %-3s  ",
                        metal->resid,
                        metal->chain,
                        metal->resname,
                        water->resid,
                        water->chain,
                        water->resname,
                        location,
                        atm->loc);

                bp->bp[resindx].fprint_bp(fp);
                //fprintf(fp,"\n");

            }

            //fprintf(fp, "WMBP  %6ld  %-3s  %-3s\n",
            //        water->resid,
            //        water->chain,
            //        water->resname);
        if(tag == 1)
                fprintf(fp,"\n");

    }

    void fprint_basepair(FILE* fp, rnabp* bp, int* flag){
        if(bp == NULL) return;
        long resindx;
        int link = 0;
        for(long i=0; i<ligand.size; ++i){
            if(ligand.restype[i] == 'N'){
                resindx = ligand.resindx[i];
                if(bp->bp[resindx].numbp > 0){
                    link++;
                }
            }
        }
        for(long i=0; i< this->ligand.size; ++i){
            if(ligand.restype[i] == 'N'){
                resindx = ligand.resindx[i];
                Molecule* mol = ligand.mol[i];
                Atom* atm = mol->residue[resindx].atom +ligand.offset[i];
                if(atm->resid != bp->bp[resindx].cifid){
                    fprintf(stderr, "Error... invalid mapping in RNA and Basepair...\n");
                    exit(EXIT_FAILURE);
                }
                if(bp->bp[resindx].numbp <= 0) continue;
                char location[5] = "---";
                if(ligand.loc[i] == 'N'){
                    strcpy(location,"NUC");
                }else if(ligand.loc[i] == 'P'){
                    strcpy(location, "PHP");
                }else if(ligand.loc[i] == 'S'){
                    strcpy(location, "SUG");
                }

                if(*flag == 0){
                    *flag = 1;
                    fprintf(fp, "\n\n      |  Metal Detail       |   Base Pair Detsils               |   outcome      "
                                "         "
                                "          |\n");
        fprintf(fp,      "      "
                         "----------------------------------------------------------------------------------------------\n");
        fprintf(fp, "       resid  chn  mtl  lnk  loc  atm   resid  chn   resid  chn     pair------>> \n");
        fprintf(fp, "      "
            "----------------------------------------------------------------------------------------------\n");


                }

                fprintf(fp, "BP    %6ld  %-3s  %-3s  %3d  %-3s  %-3s  ",
                        metal->resid,
                        metal->chain,
                        metal->resname,
                        link,
                        location,
                        atm->loc);

                bp->bp[resindx].fprint_bp(fp);
               // bp->bp[resindx-1].fprint_bp_short(fp,'F', this->metal->resname);
                //bp->bp[resindx].fprint_bp_short(fp, 'T', this->metal->resname);
                //bp->bp[resindx+1].fprint_bp_short(fp, 'F', this->metal->resname);
                //fprintf(fp,"\n\n");

            }
        }
        //if(*flag == 1){
        //    fprintf(fp, "#");
        //}
        //fprintf(fp, "\n");
    }

    void fprint_basepair_motifs(FILE* fp, rnabp* bp, Structure* struc, int* flag){
        if(bp == NULL) return;
        long resindx;
        int link = 0;
        for(long i=0; i<ligand.size; ++i){
            if(ligand.restype[i] == 'N'){
                resindx = ligand.resindx[i];
                if(bp->bp[resindx].numbp > 0){
                    link++;
                }
            }
        }
        for(long i=0; i< this->ligand.size; ++i){
            if(ligand.restype[i] == 'N'){
                resindx = ligand.resindx[i];
                Molecule* mol = ligand.mol[i];
                Atom* atm = mol->residue[resindx].atom +ligand.offset[i];
                if(atm->resid != bp->bp[resindx].cifid){
                    fprintf(stderr, "Error... invalid mapping in RNA and Basepair...\n");
                    exit(EXIT_FAILURE);
                }
                //if(bp->bp[resindx].numbp <= 0) continue;
                char location[5] = "---";
                if(ligand.loc[i] == 'N'){
                    strcpy(location,"NUC");
                }else if(ligand.loc[i] == 'P'){
                    strcpy(location, "PHP");
                }else if(ligand.loc[i] == 'S'){
                    strcpy(location, "SUG");
                }

                if(*flag == 0){
                    *flag = 1;
                    fprintf(fp, "\n\n      |  Metal Detail  |   Base Pair Detsils                                   "
                                "         "
                                "           |\n");
                    fprintf(fp,      "      "
                                     "----------------------------------------------------------------------------------------------\n");
                    fprintf(fp, "       resid  chn  mtl sec  bp    type   atm         resid \n");
                    fprintf(fp, "      "
                                "----------------------------------------------------------------------------------------------\n");


                }

                /*fprintf(fp, "BP    %6ld  %-3s  %-3s  %3d  %-3s  %-3s  ",
                        metal->resid,
                        metal->chain,
                        metal->resname,
                        link,
                        location,
                        atm->loc);*/

                //bp->bp[resindx].fprint_bp(fp);
                if(resindx !=0){
                    bp->bp[resindx-1].fprint_bp_short(fp,'F', metal, NULL, struc->secseq[resindx-1], atm->loc);
                }
                bp->bp[resindx].fprint_bp_short(fp, 'T', metal, NULL, struc->secseq[resindx], atm->loc);
                if(resindx != bp->nres-1){
                    bp->bp[resindx+1].fprint_bp_short(fp, 'F', metal, NULL, struc->secseq[resindx+1], atm->loc);
                }
                fprintf(fp,"\n\n");

            }
        }
        //if(*flag == 1){
        //    fprintf(fp, "#");
        //}
        //fprintf(fp, "\n");
    }
    void fprint_wmed_basepair_motifs(FILE* fp, rnabp* bp, Structure* struc, int* flag){
        if(bp == NULL) return;
        long resindx;
        int link = 0;
        for(long i=0; i<ligand.size; ++i){
            if(wmed.ligand.restype[i] == 'N'){
                resindx = wmed.ligand.resindx[i];
                if(bp->bp[resindx].numbp > 0){
                    link++;
                }
            }
        }
        for(long i=0; i< this->wmed.ligand.size; ++i){
            if(wmed.ligand.restype[i] == 'N'){
                resindx = wmed.ligand.resindx[i];
                Molecule* mol = wmed.ligand.mol[i];
                Atom* atm = mol->residue[resindx].atom +wmed.ligand.offset[i];
                if(atm->resid != bp->bp[resindx].cifid){
                    fprintf(stderr, "Error... invalid mapping in RNA and Basepair...\n");
                    exit(EXIT_FAILURE);
                }
                //if(bp->bp[resindx].numbp <= 0) continue;
                char location[5] = "---";
                if(wmed.ligand.loc[i] == 'N'){
                    strcpy(location,"NUC");
                }else if(wmed.ligand.loc[i] == 'P'){
                    strcpy(location, "PHP");
                }else if(wmed.ligand.loc[i] == 'S'){
                    strcpy(location, "SUG");
                }

                if(*flag == 0){
                    *flag = 1;
                    fprintf(fp, "\n\n      |  Metal Detail       |   Base Pair Detsils               |   outcome      "
                                "         "
                                "          |\n");
                    fprintf(fp,      "      "
                                     "----------------------------------------------------------------------------------------------\n");
                    fprintf(fp, "       resid  chn  mtl  lnk  loc  atm   resid  chn   resid  chn     pair------>> \n");
                    fprintf(fp, "      "
                                "----------------------------------------------------------------------------------------------\n");


                }

                /*fprintf(fp, "BP    %6ld  %-3s  %-3s  %3d  %-3s  %-3s  ",
                        metal->resid,
                        metal->chain,
                        metal->resname,
                        link,
                        location,
                        atm->loc);*/

                //bp->bp[resindx].fprint_bp(fp);
                long waterindx = wmed.waterindx[i];
                Atom* water = ligand.mol[waterindx]->residue[ligand.resindx[waterindx]].atom+ligand.offset[waterindx];
                if(resindx != 0){
                    bp->bp[resindx-1].fprint_bp_short(fp,'F', metal, water,  struc->secseq[resindx-1], atm->loc);
                }

                bp->bp[resindx].fprint_bp_short(fp, 'T', metal, water, struc->secseq[resindx], atm->loc);
                if(resindx != bp->nres-1){
                    bp->bp[resindx+1].fprint_bp_short(fp, 'F', metal,water, struc->secseq[resindx+1], atm->loc);
                }
                fprintf(fp,"\n\n");

            }
        }
        //if(*flag == 1){
        //    fprintf(fp, "#");
        //}
        //fprintf(fp, "\n");
    }

    void fprint_wmed(FILE* fp){
        if(this->wmed.size <= 0) return;
        fprintf(fp, "      |  Metal Detail  |   Water atom detail    |  Mediated atom detail  |      outcome   "
                    "         |\n");
        fprintf(fp,      "      "
                         "----------------------------------------------------------------------------------------------\n");
        fprintf(fp, "       resid  chn  mtl   resid  chn  res  lnk      resid  chn  res  atm     loc    dist    "
                    "angle\n");
        fprintf(fp, "      "
            "----------------------------------------------------------------------------------------------\n");
        for(long i=0; i<wmed.size; ++i){
            char location[5] = "---";
            char restype1 = wmed.ligand.restype[i];
            char loc1 = wmed.ligand.loc[i];
            if(restype1 == 'P'){
                strcpy(location, "PRO");
            }else if(restype1 == 'W'){
                strcpy(location, "H2O");
            }else if(restype1 == 'M'){
                strcpy(location, "MET");
            }else if(restype1 == 'N'){
                if(loc1 == 'N'){
                    strcpy(location, "NUC");
                }else if(loc1 == 'P'){
                    strcpy(location, "PHP");
                }else if(loc1 == 'S'){
                    strcpy(location, "SUG");
                }
            }
            long k = wmed.waterindx[i];
            Atom* water = ligand.mol[k]->residue[ligand.resindx[k]].atom + ligand.offset[k];
            assert(strcmp(water->resname, "HOH")==0);
            Residue* res1 = wmed.ligand.mol[i]->residue+wmed.ligand.resindx[i];
            Atom* medatm = res1->atom + wmed.ligand.offset[i];
            double dist = sqrt(wmed.ligand.dist[i]);
            double angle = wmed.angle[i];
            long link = wmed.link_count[i];

            fprintf(fp, "WMED  %6ld  %-3s  %-3s  %6ld  %-3s  %-3s  %3ld     %6ld  %-3s  %-3s  %-3s     %-3s  %6.3lf   "
                        "%6.3lf\n",
                    metal->resid,
                    metal->chain,
                    metal->resname,
                    water->resid,
                    water->chain,
                    water->resname,
                    link,
                    medatm->resid,
                    medatm->chain,
                    medatm->resname,
                    medatm->loc,
                    location,
                    dist,
                    angle

                    );
        }
        fprintf(fp, "#\n");
    }



    void fprint_angle(FILE* fp){
        if(ligand.size < 2) return;
        fprintf(fp, "      |  Metal Detail  |  C-1 atom detail      |   C-2 atom detail      |      outcome   "
                    "         |\n");
        fprintf(fp,      "      "
                         "---------------------------------------------------------------------------------------------\n");
        fprintf(fp, "       resid  chn  mtl   resid  chn  res  atm      resid  chn  res  atm     angle   energy\n");
        fprintf(fp, "      "
                    "---------------------------------------------------------------------------------------------\n");
        for(long i=0; i<ligand.size; ++i){
            for(long j=i+1; j<ligand.size; ++j){
                Atom* a = ligand.mol[i]->residue[ligand.resindx[i]].atom+ligand.offset[i];
                Atom* c = ligand.mol[j]->residue[ligand.resindx[j]].atom+ligand.offset[j];
                Atom* b = this->metal;
                double angdeg = b->angl_deg(a, c);
                fprintf(fp, "ANGL  %6ld  %-3s  %-3s  %6ld  %-3s  %-3s  %-3s  :  %6ld  %-3s  %-3s  %-3s    %6.2lf\n",
                       metal->resid,
                       metal->chain,
                       metal->resname,
                       a->resid,
                       a->chain,
                       a->resname,
                       a->loc,
                       c->resid,
                       c->chain,
                       c->resname,
                       c->loc,
                       angdeg);

            }
        }
        fprintf(fp, "#\n");
    }

    void fprint_inligand(FILE *fp){
        if(this->ligand.size == 0) return;
        fprintf(fp,"\n\n");
        fprintf(fp, "        +---------------------S I T E    S U M M A R Y---------------------+\n\n");
        fprintf(fp, "      |  Metal Detail  |      Coordination atom detail    |\n");
        fprintf(fp, "      -----------------------------------------------------\n");
        fprintf(fp, "       resid  chn  mtl   loc  res  atm   resid  chn   dist\n");
        fprintf(fp, "      -----------------------------------------------------\n");
        for(long i=0; i<this->ligand.size; ++i){
            char loca[6];
            if(ligand.restype[i] == 'N'){
                if(ligand.loc[i] == 'N'){
                    strcpy(loca, "NUC");
                }else if(ligand.loc[i] == 'P'){
                    strcpy(loca, "PHP");
                }else if(ligand.loc[i] == 'S'){
                    strcpy(loca, "SUG");
                }else{
                    fprintf(stderr,"Error... Wrong type given\n");
                    exit(EXIT_FAILURE);
                }
            }else if(ligand.restype[i] == 'P'){
                strcpy(loca, "PRO");
            }else if(ligand.restype[i] == 'W'){
                strcpy(loca, "H2O");
            }else{
                fprintf(stderr,"Error... Wrong restype given\n");
                exit(EXIT_FAILURE);
            }
            Atom* a = this->ligand.mol[i]->residue[ligand.resindx[i]].atom+ ligand.offset[i];
            fprintf(fp, "BIND  %6ld  %-3s  %-3s   %-3s  %-3s  %-3s  %6ld  %-3s %6.3lf\n",
                    this->metal->resid,
                    this->metal->chain,
                    this->metal->resname,
                    loca,
                    a->resname,
                    a->loc,
                    a->resid,
                    a->chain,
                    ligand.dist[i]);

        }
        fprintf(fp,"#\n");

        //fprintf(fp, "TOTAL LIGANDS FOR %3s: %4ld",this->metal->resname, ligand.size);
    }


    void add_metal(Atom* met){
        this->metal = met;
        this->atomic_no = get_atomic_no(met->resname);
        if(this->atomic_no <0){
            fprintf(stderr, "Error... not  valid metal: %s\n", met->resname);
        }
    }

    void fill_proximity(Molecule *mol, char moltype);
    void comp_nucleic(char rule,
                 Current_params* mprm,
    Map* map);
    void comp_protein(char rule,
                 Current_params* mprm,
    Map* map);
    void comp_water(char rule, Current_params* mprm);

    void comp_water_mediated(char rule, Current_params* mprm);

    void site_populate(Molecule* rna,
                       Molecule* pro,
                       Molecule* hoh,
                       Molecule* met,
                       char rule,
                       Paremeters* prm,
                       Map* rnamap,
                       Map* promap);

};
class Map{
public:
    Molecule* mol;
    Site** metal1;
    Site** metal2;
    Site** metal3;
public:
    Map(Molecule* mol);
    void add_site(Site* site, Ligand* lig, long ligindx);
    void gen_verna(char *file){
        //FILE* fp = fopen(file, "w");
        //assert(fp != NULL);
        //fprintf(fp, "<applet  code=\"VARNA.class\" codebase=\"/usr/local/\" archive=\"VARNAv3-93-src.jar\" "
          //          "height=6000 width=4000>\n");
        int count = 0;
        for(long i=0; i< mol->size; ++i){
            if(this->metal1[i] != NULL){
                   //fprintf(fp, " WORKING\n");
                count ++;
            }

        }
        //fclose(fp);
        if(count >0){
	      ;
        }
    }
};



void Site::site_populate(
                   Molecule* rna,
                   Molecule* pro,
                   Molecule* hoh,
                   Molecule* met,
                   char rule,
                   Paremeters* prm,
                   Map* rnamap,
                   Map* promap){

    //site_init(site);

    Current_params mprm;
    //printf("Atomic %d\n", this->atomic_no);
    mprm.adjst(prm, this->atomic_no,rule, 0.5);


    //metal_params_adjust(&mprm, prm, this->atomic_no, 0.5);



    this->fill_proximity(rna, 'N');
    this->comp_nucleic(rule, &mprm, rnamap);
    this->fill_proximity(pro, 'P');
    this->comp_protein(rule, &mprm, promap);
    this->fill_proximity(hoh, 'W');
    this->comp_water(rule, &mprm);

    this->comp_water_mediated(rule, &mprm);

    //site_addproxatm(site, hoh, 'W');
    //site_addhoh(site, hoh, rule, &mprm);

}


void comp_metal_sites(Molecule* met,
                      Molecule* hoh,
                      Molecule* rna,
                      rnabp* rnabp,
                      Molecule* pro,
                      Paremeters* prm,
                      Structure* sec,
                      char rule,
                      char* file_path,
                      char* file_name,
                      FILE* fp) {
    long nsites = met->atom_count();
    Atom** met_atoms = new Atom*[nsites];
    long count; // should be same as nsites;
    met->get_all_atoms(met_atoms, &count);
    Site* sites = new Site[nsites];//(struct Site*) malloc(met->size * sizeof(struct Site));
    //long nsites = 10*met->size;
//	for(long i=0; i<nsites; ++i){
//		sites[i] = (struct Site*) malloc(sizeof(struct Site));
//	}
    Map rnamap = Map(rna);
    Map promap = Map(pro);
    //molecule_metal_init(&rnamap, rna);
    for(long i=0; i<nsites; ++i){


            Atom* metal = met_atoms[i];
            sites[i].add_metal(metal);
            //site_init(sites+i, met->atom+met->resbeg[i], get_atomic_no(metal->loc));
            sites[i].site_populate(rna, pro, hoh, met, rule, prm, &rnamap, &promap); //site_populate(sites+i, rna, pro,

        // hoh,
        // met, rule,
            //  prm, &rnamap);

        //printf("Total atoms in site %s %ld\n",sites[i].metal->resname, sites[i].ligand.size);

    }

    for(long i=0; i<nsites; ++i){
        if(i == 0){
            fprintf(fp, "\n\n\n\n            +---------------------- S U M M A R Y  R E P O R T "
                        "---------------------+\n\n\n\n");
            fprintf(fp, "      |  Metal Detail  |   Inner Coordination     |     Water Mediated "
                        "         |\n");

        fprintf(fp, "      ---------------------------------------------------------------------------\n");
        fprintf(fp, "       resid  chn  mtl  cord   NUC  PRO  H2O  MET    cnt     NUC  PRO  H2O  MET \n");
        fprintf(fp, "      ---------------------------------------------------------------------------\n");
        }
        sites[i].fprint_summary(fp);
    }
    fprintf(fp, "#\n");
            fprintf(fp, "\n\n\n\n        +--------------------- D E T A I L  R E P O R T "
                        "---------------------+\n\n\n\n");
    int flag=0;
    for(long i=0; i<nsites; ++i){
        sites[i].fprint_basepair(fp, rnabp, &flag);
    }
    fprintf(fp,"#\n");

    flag=0;
    for(long i=0; i<nsites; ++i){
        sites[i].fprint_basepair_motifs(fp, rnabp, sec, &flag);
        sites[i].fprint_wmed_basepair_motifs(fp, rnabp, sec, &flag);
    }
    fprintf(fp,"#\n");

    flag=0;
    for(long i=0; i<nsites; ++i){
        sites[i].fprint_wmed_basepair(fp, rnabp, &flag);
    }
    //fprintf(fp,"#\n");

    //rnabp->fprint_bp(fp);
    for(long i=0; i<nsites; ++i){
        sites[i].fprint_inligand(fp);
        sites[i].fprint_angle(fp);
        sites[i].fprint_wmed(fp);
    }
    fprintf(fp,"TER\n");

    //gen_verna(rna,rnabp, &rnamap, sec, sites, nsites);
    char nuc_file[512];

    file_name_join(nuc_file, file_path,file_name, ".met");

    rnamap.gen_verna(nuc_file);

    delete [] sites;
    delete [] met_atoms;
}



void Site::fill_proximity(Molecule *mol, char moltype) {
    for(long i=0; i< mol->size; ++i){
        Atom* a = mol->residue[i].atom+0;
        double val = this->metal->distsqr(a);
        double radius = 20.0*20.0;
        if(moltype == 'W'){
            radius = 8.5 * 8.5;
        }
        if(val< radius){  // skip resi which are far than 20A
            this->prox.add(mol, i, moltype);
        }
    }
}



Map::Map(Molecule* mol) {
this->mol = mol;
long nres = mol->size;
this->metal1 = (Site**) malloc(nres * sizeof(Site*));
this->metal2 = (Site**) malloc(nres * sizeof(struct Site*));
this->metal3 = (Site**) malloc(nres * sizeof(struct Site*));
for(long i=0; i<nres; ++i){
    this->metal1[i] = NULL;
    this->metal2[i] = NULL;
    this->metal3[i] = NULL;
}
}

void Map::add_site(Site* site, Ligand* lig, long ligindx) {
assert(this->mol == lig->mol[ligindx]);
long resindx = lig->resindx[ligindx];
if(this->metal1[resindx] == NULL){
    this->metal1[resindx] = site;
}else if(this->metal2[resindx] == NULL){
    this->metal2[resindx] = site;
}else if(this->metal3[resindx] == NULL){
    this->metal3[resindx] = site;
}else {
    fprintf(stderr, "Error... too many sites bind to same molecule\n");
    //exit(1);
}
}


void Site::comp_nucleic(char rule,
                     Current_params* mprm,
                     Map* map){
        Atom* currmetal = this->metal;

//mprm->print();
        double dst2;
        for(long resindx=0; resindx<prox.size; ++resindx){
           // if(prox.restype[resindx] != 'N') continue;


            if(prox.restype[resindx] != 'N') continue; // Not necleic acid

            //long start = rna->resbeg[this->proxres[resindx]];
            //long end   = rna->resend[this->proxres[resindx]];
            Residue* resi = this->prox.get_residue(resindx); //  ->residue+resindx;
            for(long offset=0; offset<resi->size; ++offset){
                Atom* curatom = resi->atom+offset;
                //currmetal->print();
                //curatom->print();
                dst2 = currmetal->distsqr(curatom)+0.00001;
                //printf("dst=%lf\n", dst2);

                long len = strlen(curatom->loc);
                char symb = curatom->loc[0];
                if(len == 2){ //Then it is Nucleobase for Cor dataset
                    if((symb == 'O' && dst2 <= mprm->o_dst2)
                       ||( symb == 'N' && dst2 <= mprm->n_dst2)
                       ||( symb == 'C' && dst2 <= mprm->c_dst2)){

                        ligand.add(&prox, resindx, offset,
                                   'N',
                                   rule,
                                   sqrt(dst2));
                        map->add_site(this, &this->ligand, ligand.size-1);
                        //printf("Added NUC\n");
                    }
                }else if(len == 3 && (curatom->loc[2] == '*' || curatom->loc[2] == '\'')){ // sugar backbone for cor
                    // file
                    if((symb == 'O' && dst2 <= mprm->o_dst2)
                       ||( symb == 'C' && dst2 <= mprm->c_dst2)){

                        ligand.add(&prox, resindx, offset, 'S', rule, sqrt(dst2));
                        map->add_site(this, &this->ligand, ligand.size-1);
                        //printf("Added SUG\n");
                        //molecule_metal_add_site(map, rna, site, this->proxres[resindx]);
                    }
                }else{ // phosphate case
                    if((len == 3 && symb == 'O' && dst2 <= mprm->o_dst2)
                       ||(len == 3 && symb == 'P' && dst2 <= mprm->p_dst2)){

                        ligand.add(&prox, resindx, offset, 'P', rule, sqrt(dst2));
                        map->add_site(this, &this->ligand, ligand.size-1);
                        //printf("Added PHP\n");
                        //molecule_metal_add_site(map, rna, site, this->proxres[resindx]);
                    }
                }
            }
        }
    }
void Site::comp_protein(char rule,
                     Current_params* mprm,
                     Map* map){
        Atom* currmetal = this->metal;

//mprm->print();
        double dst2;
        for(long resindx=0; resindx<prox.size; ++resindx){
           // if(prox.restype[resindx] != 'N') continue;


            if(prox.restype[resindx] != 'P') continue; // Not Protein

            //long start = rna->resbeg[this->proxres[resindx]];
            //long end   = rna->resend[this->proxres[resindx]];
            Residue* resi = this->prox.get_residue(resindx); //  ->residue+resindx;
            for(long offset=0; offset<resi->size; ++offset){
                Atom* curatom = resi->atom+offset;
                //currmetal->print();
                //curatom->print();
                dst2 = currmetal->distsqr(curatom)+0.00001;
                //printf("dst=%lf\n", dst2);

                char symb = curatom->loc[0];
                //if(len == 2){ //Then it is Nucleobase for Cor dataset
                    if((symb == 'O' && dst2 <= mprm->o_dst2)
                       ||( symb == 'N' && dst2 <= mprm->n_dst2)
                       ||( symb == 'C' && dst2 <= mprm->c_dst2)){

                        ligand.add(&prox, resindx, offset, 'P', rule, sqrt(dst2));
                        map->add_site(this, &this->ligand, ligand.size-1);
                        //printf("Added NUC\n");
                    }
                //}
            }
        }
    }

void Site::comp_water(char rule, Current_params *mprm) {
Atom* currmetal = this->metal;
        double dst2;
        for(long resindx=0; resindx<prox.size; ++resindx){
            if(prox.restype[resindx] != 'W') continue; // Not Water acid
            Residue* resi = this->prox.get_residue(resindx); //  ->residue+resindx;
            for(long offset=0; offset<resi->size; ++offset){
                Atom* curatom = resi->atom+offset;
                dst2 = currmetal->distsqr(curatom)+0.00001;
                    if(strcmp(curatom->resname, "HOH") == 0 && dst2 <= mprm->hoh_dst2){
                        ligand.add(&prox, resindx, offset, 'W', rule, sqrt(dst2));
                    }
            }
        }
}

void Site::comp_water_mediated(char rule, Current_params *mprm) {
    double metdist2=0.0, wdist2=0.0, includist2=0.0;
    includist2 = mprm->hb_dst + mprm->o_dst2 + 1.0;
    includist2 *= includist2;
    double hbdst2 = (mprm->hb_dst +0.0001) * (mprm->hb_dst + 0.0001);
    for(long i=0; i<ligand.size; ++i){
        if(ligand.restype[i] != 'W') continue;
        Residue* wres = ligand.mol[i]->residue + ligand.resindx[i];

        Atom* water = wres->atom + ligand.offset[i];
        for(long j=0; j<this->prox.size; ++j){
            Residue* res = prox.mol[j]->residue+ prox.resindx[j];
            for(long k=0; k<res->size; ++k){
                Atom* curratm = res->atom + k;
                metdist2 = metal->distsqr(curratm);
                wdist2 = water->distsqr(curratm);
                int flag = 1;
                if(rule != '0'){
                    if(curratm->loc[0] == 'N' || curratm->loc[0] == 'O')
                        flag = 1;
                    else
                        flag = 0;
                }
                double angle = water->angl_deg(metal, curratm);
                if(metdist2 <= includist2 &&
                        metdist2 > hbdst2 &&
                        wdist2 < hbdst2   &&
                        angle >=-520.0 &&
                        flag == 1
                        ){
                    char loc= get_res_loc(curratm, prox.restype[j]);
                    this->wmed.add(&prox, j, k, loc, rule, i, wdist2, angle);
                }
            }
        }
    }
}

#endif //CPPMET_METALBIND_H
