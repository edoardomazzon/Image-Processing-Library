/*
 Created by Sebastiano Vascon on 23/03/20.
*/

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"
#include <math.h>
#include <time.h>
#include <assert.h>

void ip_mat_show(ip_mat * t){
    unsigned int i,l,j;
    printf("Matrix of size %d x %d x %d (hxwxk)\n",t->h,t->w,t->k);
    for (l = 0; l < t->k; l++) {
        printf("Slice %d\n", l);
        for(i=0;i<t->h;i++) {
            for (j = 0; j < t->w; j++) {
                printf("%f ", get_val(t,i,j,l));
            }
            printf("\n");
        }
        printf("\n");
    }
}

void ip_mat_show_stats(ip_mat * t){
    unsigned int k;
    
    compute_stats(t);
    
    for(k=0;k<t->k;k++){
        printf("Channel %d:\n", k);
        printf("\t Min: %f\n", t->stat[k].min);
        printf("\t Max: %f\n", t->stat[k].max);
        printf("\t Mean: %f\n", t->stat[k].mean);
    }
}

ip_mat * bitmap_to_ip_mat(Bitmap * img){
    unsigned int i=0,j=0;

    unsigned char R,G,B;

    unsigned int h = img->h;
    unsigned int w = img->w;

    ip_mat * out = ip_mat_create(h, w,3,0);

    for (i = 0; i < h; i++)              /* rows */
    {
        for (j = 0; j < w; j++)          /* columns */
        {
            bm_get_pixel(img, j,i,&R, &G, &B);
            set_val(out,i,j,0,(float) R);
            set_val(out,i,j,1,(float) G);
            set_val(out,i,j,2,(float) B);
        }
    }

    return out;
}

Bitmap * ip_mat_to_bitmap(ip_mat * t){

    Bitmap *b = bm_create(t->w,t->h);

    unsigned int i, j;
    for (i = 0; i < t->h; i++)              /* rows */
    {
        for (j = 0; j < t->w; j++)          /* columns */
        {
            bm_set_pixel(b, j,i, (unsigned char) get_val(t,i,j,0),
                    (unsigned char) get_val(t,i,j,1),
                    (unsigned char) get_val(t,i,j,2));
        }
    }
    return b;
}

float get_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k){
    if(i<a->h && j<a->w &&k<a->k){  /* j>=0 and k>=0 and i>=0 is non sense*/
        return a->data[i][j][k];
    }else{
        printf("Errore get_val!!!");
        exit(1);
    }
}

void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
    if(i<a->h && j<a->w &&k<a->k){
        a->data[i][j][k]=v;
    }else{
        printf("Errore set_val!!!");
        exit(1);
    }
}

float get_normal_random(){
    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    return cos(2*PI*y2)*sqrt(-2.*log(y1));
}

ip_mat * ip_mat_create(unsigned int h, unsigned int w,unsigned  int k, float v) {    
    ip_mat *miamatrice = (ip_mat*)malloc(sizeof(ip_mat)+sizeof(stats)*k);

    /* Inizializzo la matrice con dimensioni h w e k */
    miamatrice->h = h;
    miamatrice->w = w;
    miamatrice->k = k;

    int i,j,l;

    miamatrice->data=(float ***)malloc(sizeof(float **)*h);
    for(i=0;i<h;i++){
        miamatrice->data[i]=(float **)malloc(sizeof(float*)*w);
        for(j=0;j<w;j++){
            miamatrice->data[i][j]=(float *)malloc(sizeof(float)*k);
        }
    }
    for(i=0;i<h;i++){
        for(j=0;j<w;j++){
            for(l=0;l<k;l++){
                miamatrice->data[i][j][l] = v;
            }
        }
    }

    /* Creo un vettore di stats per contenere le statische sui singoli canali */
    
    miamatrice->stat = (stats*)malloc(sizeof(stats)*k);
    if(miamatrice==NULL){
        printf("Errore malloc create stat");
        exit(1);
    }
    /* Riempio i k-esimi stat con le statistiche iniziali, cioè v */
    for(int x = 0; x < k; x++){
        (miamatrice->stat+x)->min = v;
        (miamatrice->stat+x)->max = v;
        (miamatrice->stat+x)->mean = v;
    }
    return miamatrice;
}

void ip_mat_free(ip_mat *a) {

    free(a->stat);
    
    /* Libero t->data */
    for(int i = 0; i< a->h; i++)
    {
        for (int j = 0; j< a->w; j++)
            free(a->data[i][j]);

        free(a->data[i]);
    }
    free(a); 
    /* Libero l'intera struct ip_mat */
}

void compute_stats(ip_mat * t) {    
    for (int x=0; x < t->k; x++) {
        float sum = 0;
        int nums = 0;
        /* Dichiara qui max e min e ad essi associa il primo valore della matrice */
        float currentmax = t->data[0][0][x];
        float currentmin = t->data[0][0][x];
        for(int i = 0; i < t->h; i++){
            for(int j = 0; j < t->w; j++){
                /* Confronta min e max e aumento i vari contatori per la media*/
                nums++;
                sum += t->data[i][j][x];
                if(t->data[i][j][x] > currentmax){
                    currentmax = t->data[i][j][x];
                }
                if(t->data[i][j][x] < currentmin){
                    currentmin = t->data[i][j][x];
                }
            }
        }
        /* Associa massimo, minimo e media alla x-esima, cioè t->k-esima, stat */    
        (t->stat+x)->max = currentmax;
        (t->stat+x)->min = currentmin;
        (t->stat+x)->mean = sum/nums;
    }    
}

/* Ausiliare per ip_mat_init_random */
double rand_normal(double mean, double stddev) {
    double n2 = 0.0;
    int n2_cached = 0;
    if (!n2_cached) {
        double x, y, r;
        do {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}
/* Inizializza una ip_mat con dimensioni w h e k.
 * Ogni elemento è generato da una gaussiana con media mean e varianza var */
/* Da correggere con mean e var */
void ip_mat_init_random(ip_mat * t, float mean, float var) { 
    srand(time(NULL));
    for(int i=0; i<t->h; i++){
        for(int j=0; j<t->w; j++){
            for(int l=0; l<t->k; l++){
                t->data[i][j][l] = rand_normal(mean, sqrt(var));/* variazione standard = sqrt(varianza) */
            }
        }
    }
}

ip_mat * ip_mat_copy(ip_mat * in) {
    ip_mat *newmat = ip_mat_create(in->h, in->w, in->k, 0.0);

    for(int i=0; i < newmat->h; i++){
        for(int j=0; j < newmat->w; j++){
            for(int l=0; l < newmat->k; l++){
                newmat->data[i][j][l] = in->data[i][j][l];
            }
        }
    }
    
    compute_stats(newmat);
    return newmat;
}

ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end) {
    unsigned int row = row_end-row_start;
    unsigned int col = col_end-col_start;
    
    ip_mat *submat = ip_mat_create(row, col, t->k, 0.0);

    int i, j, l;
    int o = 0, p = 0, q = 0;
    /* Riempio la sottomatrice da row_start a row_end e da col_start a col_end */
    for(i=row_start;i<row_end;i++){
        p=0;
        for(j=col_start;j<col_end;j++){
            q=0;
            for(l=0;l<t->k;l++){
                submat->data[o][p][q] = t->data[i][j][l];
                q++;
            }
            p++;
        }
        o++;
    }

    compute_stats(submat);
    return submat;
}

ip_mat * ip_mat_concat(ip_mat * a, ip_mat * b, int dimensione);