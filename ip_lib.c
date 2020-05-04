/*
 Created by Sebastiano Vascon on 23/03/20.
*/

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"
#include <math.h>
#include <time.h>

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
    ip_mat *miamatrice = (ip_mat*)malloc(sizeof(ip_mat));

    // Inizializzo la matrice con dimensioni h w e k
    miamatrice->h = h;
    miamatrice->w = w;
    miamatrice->h = h;

    // Ogni elemento è inizializzato a v.
    miamatrice->data = (float***)malloc(sizeof(float**));

    int i,j,l;

    miamatrice->data=(float ***)malloc(sizeof(float ***)*h);
    for(i=0;i<h;i++)
    {
        miamatrice->data[i]=(float **)malloc(sizeof(int*)*w);
        for(j=0;j<w;j++)
        {
            miamatrice->data[i][j]=(float *)malloc(sizeof(float)*k);
        }
    }
    for(i=0;i<h;i++)
    {
        for(j=0;j<w;j++)
        {
            for(l=0;l<k;l++)
            {
                miamatrice->data[i][j][l] = v;
            }
        }
    }

    // creo un vettore di stats per contenere le statische sui singoli canali
    miamatrice->stat = (stats*)malloc(sizeof(stats));
    /*miamatrice->stat->min = v;
    miamatrice->stat->max = v;
    miamatrice->stat->mean = v;*/

    return miamatrice;
}

void ip_mat_free(ip_mat *a) {
    // libero la stat
    for(int x=0; x<a->h; x++) {
        free(a->stat+x);
    }
    // libero data
    for(int i = 0; i< a->h; i++)
    {
        for (int j = 0; j< a->w; j++)
            free(a->data[i][j]);

        free(a->data[i]);
    }
    free(a); // libero la struttura
}

/* Restituisce il valore in posizione i,j,k */
float get_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k) {
    return a->data[i][j][k];
}

/* Setta il valore in posizione i,j,k a v*/
void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v) {
    a->data[i][j][k];
}

/* Da completare */
void compute_stats(ip_mat * t) {
    int n_can = t->k;
    /* Creo statarray come array di k stat, quindi una stat per colore */
    stats *statarray = (stats*)malloc(sizeof(stats)*n_can);

    float currentmin;
    float currentmax;
    for (int x=0; x<n_can; x++) {
        float sum = 0;
        int nums = 0;
        for(int i = 0; i < t->h;i++){
            for(int j = 0; j < t->w;j++){
                /* Confronta min e max */
               sum += t->data[i][j][x];
               nums++;
            }
        }
        (statarray+x)->max = currentmax;
    }

    t->stat = statarray;
}

/* Ausiliare per ip_mat_init_random */
double rand_normal(double mean, double stddev) {
    srand(time(NULL));
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
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
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}
/* Inizializza una ip_mat con dimensioni w h e k.
 * Ogni elemento è generato da una gaussiana con media mean e varianza var */
void ip_mat_init_random(ip_mat * t, float mean, float var) { /* converti var in stddev (stddev=sqrt(var))*/
    for(int i=0; i<t->h; i++){
        for(int j=0; j<t->w; j++){
            for(int l=0; l<t->k; l++){
                t->data[i][j][l] = rand_normal(mean, sqrt(var));
            }
        }
    }
}
