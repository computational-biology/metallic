//
// Created by parthajit on 23/2/20.
//

#ifndef CPPMET_SYSPARAMS_H
#define CPPMET_SYSPARAMS_H
#include <string>
using namespace std;
class sysparams{
      public:


	    char cifparam[512];
	    char accnparam[512];
	    char htparam[512];
	    char hdparam[512];
	    char hdvalparam[512];
	    char chainparam[512];
	    char chainvalparam[512];
	    char angparam[512];
	    char angvalparam[512];
	    char chparam[512];
	    char sgparam[512];
	    char evaltypeparam[512];
	    char corparam[50];
	    string type;
	    int cleaned_res;
	    string file_dir;
	    string accn;
	    string ext;
	    sysparams(){
		  strcpy(cifparam, "-dummyval");
		  strcpy(accnparam, "-dummyval");
		  strcpy(htparam, "-HT");
		  strcpy(hdparam, "-dummyval");
		  strcpy(hdvalparam, "-dummyval");
		  strcpy(chainparam, "-dummyval");
		  strcpy(chainvalparam, "-dummyval");
		  strcpy(angparam, "-dummyval");
		  strcpy(angvalparam, "-dummyval");
		  strcpy(chparam, "-dummyval");
		  strcpy(sgparam, "-dummyval");
		  strcpy(evaltypeparam, "-dummyval");
		  strcpy(corparam, "-dummyval");
	    }
	    void print_params(FILE* fp)
	    {
		  fprintf(fp, "\n---------------------------P A R A M S    O P T E D ---------------------\n");
		        if(strcmp(htparam, "-dummyval") == 0){
			      fprintf(fp,"PARAM   HETATM NOT REQUESTED\n");
			}else{
			      fprintf(fp,"PARAM   HETATM REQUESTED\n");
			}
		        if(strcmp(hdparam, "-dummyval") == 0){
			      fprintf(fp,"PARAM   HYDROGEN BOND DIST(FOR BASE PAIR)    3.8A (DEFAULT)\n");
			}else{
			      fprintf(fp,"PARAM   HYDROGEN BOND DIST(FOR BASE PAIR)    %sA (REQUESTED)\n",hdvalparam);
			}
		        if(strcmp(angparam, "-dummyval") == 0){
			      fprintf(fp,"PARAM   ANGLE IN DEGREE(FOR BASE PAIR)    120 (DEFAULT)\n");
			}else{
			      fprintf(fp,"PARAM   ANGLE IN DEGREE(FOR BASE PAIR)    %s (REQUESTED)\n",angvalparam);
			}
		        if(strcmp(chparam, "-dummyval") == 0){
			      fprintf(fp,"PARAM   C-H...O/N MEDIATED BASE PAIR   REQUESTED\n");
			}else{
			      fprintf(fp,"PARAM   C-H...O/N MEDIATED BASE PAIR   NOT REQUESTED\n");
			}
		        if(strcmp(sgparam, "-dummyval") == 0){
			      fprintf(fp,"PARAM   SUGAR O2' MEDIATED BASE PAIR   REQUESTED\n");
			}else{
			      fprintf(fp,"PARAM   SUGAR O2' MEDIATED BASE PAIR   NOT REQUESTED\n");
			}
			
	    }
};


#endif //CPPMET_SYSPARAMS_H
