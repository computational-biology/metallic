//
// Created by parthajit on 13/7/20.
//

#ifndef CPPMET_METALDEFS_H
#define CPPMET_METALDEFS_H

#include <cstdio>
#include <cassert>
#include <cstring>
#include <cstdlib>


#define NUM_METAL 128
#define PI 3.14159

char _global_metal[NUM_METAL][5];
int _global_metal_atomic[NUM_METAL];
long _global_metal_size=0;

class Paremeters{
public:
    int n;
    char name[NUM_METAL][5];
    int atomic_no[NUM_METAL];
    int coord_no[NUM_METAL];
    double hb_dst[NUM_METAL];
    double vdw_rad[NUM_METAL];
    double ion_rad[NUM_METAL];
    double o_dst[NUM_METAL];
    double n_dst[NUM_METAL];
    double hoh_dst[NUM_METAL];
    double c_dst[NUM_METAL];
    double p_dst[NUM_METAL];
    double met_dst[NUM_METAL];
public:
    Paremeters(){
        n = 0;
    }
    Paremeters(const char* param_file);
};

Paremeters::Paremeters(const char *param_file) {
   FILE* fp = fopen(param_file,"r");
   if(fp == NULL){    /* Exception Handling */ 
	 fprintf(stderr, "Error in function %s. (File: %s, Line %d)... Cannot open metal param from %s. make sure that it exists and readable.\n", __func__, __FILE__, __LINE__, param_file);
	 exit(EXIT_FAILURE);
   }
   const char* sep="\t \n";
   char line[1024];
   char* token;
   char* token1;
   while(fgets(line, 1024, fp) != NULL){
      if(strncmp(line, "_metal_dist.metal_symb", 22) == 0){
         while(fgets(line, 1024, fp) != NULL){
            if(line[0] != '_'){
               break;
            }
         }
         break;
      }
   }
   for(int i=0; i< NUM_METAL; ++i){
      strcpy(this->name[i],"-");
   }
   this->n = 0;
   _global_metal_size = 0;
   do{

      if(line[0] == '#') break;

      token1 = strtok(line, sep);


      token = strtok(NULL, sep); // atomic_no; // atomic_no is the index;
      int atomic_no = atoi(token);
      strcpy(this->name[atomic_no], token1);
      strcpy(_global_metal[_global_metal_size], token1);
      _global_metal_atomic[_global_metal_size] = atomic_no;
      _global_metal_size++;

      this->atomic_no[atomic_no] = atomic_no;

      this->hb_dst[atomic_no] = 3.8;

      token = strtok(NULL, sep);
      this->coord_no[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atoi(token);

      token = strtok(NULL, sep);
      this->vdw_rad[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);



      token = strtok(NULL, sep);
      this->ion_rad[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);




      token = strtok(NULL, sep);
      this->o_dst[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);



      token = strtok(NULL, sep);
      this->n_dst[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);


      token = strtok(NULL, sep);
      this->hoh_dst[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);



      token = strtok(NULL, sep);
      this->c_dst[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);



      token = strtok(NULL, sep);
      this->p_dst[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);




      token = strtok(NULL, sep);
      this->met_dst[atomic_no] = (strcmp(token, "?") == 0) ? -9.9 : atof(token);

      this->n ++;
   }while(fgets(line, 1024, fp) != NULL);
   fclose(fp);
}

int is_metal(char* metal){
    for(long i=0; i<_global_metal_size; ++i){
        if(strcmp(metal, _global_metal[i]) == 0) return 0;
    }
    return -1;
}

int is_metal1(char* metal){
   if(strcmp(metal, "MG") == 0) return 0;
   if(strcmp(metal, "CA") == 0) return 0;
   if(strcmp(metal, "ZN") == 0) return 0;
   if(strcmp(metal, "MN") == 0) return 0;
   return -1;

}

int get_atomic_no(char* metal){
    for(int i=0; i<_global_metal_size; ++i){
        if(strcmp(metal, _global_metal[i]) == 0){
            return _global_metal_atomic[i];
        }
    }
    return -1;
}

double calc_angle_energy(char* metal, double angle){
   double ktheta;

   if(strcmp(metal, "MG") == 0){
      ktheta = 50.0;
      if(angle > (3.0*PI)/4.0){
         return ktheta * (PI - angle) * (PI - angle);
      }else{
         return ktheta * ((PI/2.0) - angle) * ((PI/2.0) - angle);
      }
   }
   return -1.0;
}
class Current_params{
public:
    char name[5];
    int atomic_no;
    int coord_no;
    double hb_dst;
    double vdw_rad;
    double ion_rad;
    double o_dst2;
    double n_dst2;
    double hoh_dst2;
    double c_dst2;
    double p_dst2;
    double met_dst2;
public:
    void print(){
        printf("o2= %lf, n2=%lf, c2=%lf, p2=%lf\n", o_dst2,
        n_dst2, c_dst2, p_dst2);
    }
    void adjst(Paremeters* prm, int atomic_no, char rule, double adjust){
        //printf("%d %d\n", prm->atomic_no[atomic_no], atomic_no);
        assert(prm->atomic_no[atomic_no] == atomic_no);
        strcpy(this->name, prm->name[atomic_no]);
        this->atomic_no = prm->atomic_no[atomic_no];
        this->coord_no = prm->coord_no[atomic_no];
        this->hb_dst = prm->hb_dst[atomic_no];
        this->vdw_rad = prm->vdw_rad[atomic_no];
        this->ion_rad = prm->vdw_rad[atomic_no];
        this->o_dst2 = prm->o_dst[atomic_no];
        this->n_dst2 = prm->n_dst[atomic_no];
        this->hoh_dst2 = prm->hoh_dst[atomic_no];
        this->c_dst2 = prm->c_dst[atomic_no];
        this->p_dst2 = prm->p_dst[atomic_no];
        this->met_dst2 = prm->met_dst[atomic_no];
        if(rule != '0'){
            this->c_dst2 = -9.9;
            this->p_dst2 = -9.9;
        }

        if(this->o_dst2 >= 0.0)
            this->o_dst2 = (this->o_dst2 + adjust) * (this->o_dst2 + adjust);
        if(this->n_dst2 >= 0.0)
            this->n_dst2 = (this->n_dst2 + adjust) * (this->n_dst2 + adjust);
        if(this->hoh_dst2 >= 0.0)
            this->hoh_dst2 = (this->hoh_dst2 + adjust) * (this->hoh_dst2 + adjust);
        if(this->c_dst2 >= 0.0)
            this->c_dst2 = (this->c_dst2 + adjust) * (this->c_dst2 + adjust);
        if(this->p_dst2 >= 0.0)
            this->p_dst2 = (this->p_dst2 + adjust) * (this->p_dst2 + adjust);
        if(this->met_dst2 >= 0.0)
            this->met_dst2 = (this->met_dst2 + adjust) * (this->met_dst2 + adjust);
    }
};


#endif //CPPMET_METALDEFS_H
