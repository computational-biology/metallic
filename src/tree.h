//
// Created by parthajit on 16/7/20.
//

#ifndef CPPMET_TREE_H
#define CPPMET_TREE_H
#include<cstdio>
#include<cstring>
#include <cstdlib>
#define SNIPPET_SIZE 10
struct counting_tree_t{
    char text[SNIPPET_SIZE];
    int txtlen;
    long freq;
    struct counting_tree_t* left;
    struct counting_tree_t* right;
};


struct counting_tree_t* counting_tree_getnode(char* line){
    struct counting_tree_t* node = (struct counting_tree_t*) malloc(sizeof(struct counting_tree_t));

    node->txtlen = strlen(line);
    node->freq = 1;
    strcpy(node->text,line);
    node->text[node->txtlen] = '\0';
    node->left = NULL;
    node->right = NULL;
    return node;
}


void counting_tree_insert(struct counting_tree_t* root, char* line){
    int cmp = strcmp(line, root->text);
    if(cmp < 0){
        if(root->left == NULL){
            root->left = counting_tree_getnode(line);
        }else{
            counting_tree_insert(root->left, line);
        }
    }else if(cmp > 0){
        if(root->right == NULL){
            root->right = counting_tree_getnode(line);
        }else{
            counting_tree_insert(root->right, line);
        }
    }else{

        root->freq ++;
    }
    return;
}
long counting_tree_node_count(struct counting_tree_t* t){
    if( t != NULL){
        long l = counting_tree_node_count(t->left);
        long r = counting_tree_node_count(t->right);
        return l + r + 1;
    }
    return 0;
}

void counting_tree_fprintf(FILE* fp, struct counting_tree_t* t, int* count, char newline){
    if(t != NULL){
        counting_tree_fprintf(fp, t->left, count, newline);
        if((*count)%5 == 0 && newline == 'F'){
            fprintf(fp, "\n            ");

        }
        (*count) ++;
        if(newline == 'T'){
            fprintf(fp, "              %-3s : %-5ld\n",t->text, t->freq);
        }else{
            fprintf(fp, "%-3s : %-5ld   ",t->text, t->freq);
        }
        counting_tree_fprintf(fp, t->right, count, newline);
    }
    return;
}
void counting_tree_free(struct counting_tree_t* t){
    if(t != NULL){
        counting_tree_free( t->left);
        counting_tree_free(t->right);
        free(t);
    }
    return;
}


int counting_tree_search(struct counting_tree_t* t, char* word){
    if(t == NULL){
        return 0;
    }
    int test = strcmp(word, t->text);
    if(test == 0){
        return 1;
    }else if(test < 0){
        return counting_tree_search(t->left, word);
    }else{
        return counting_tree_search(t->right, word);
    }
}
#endif //CPPMET_TREE_H
