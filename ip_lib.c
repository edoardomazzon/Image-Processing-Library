/*
 Created by Sebastiano Vascon on 23/03/20.
*/

#include <stdio.h>
#include <stdlib.h>
#include "ip_lib.h"
#include "bmp.h"

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

    ip_mat * out = ip_mat_create(h, w,3,0.0);

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

/**** PARTE 1: TIPO DI DATI ip_mat E MEMORIA ****/

ip_mat * ip_mat_create(unsigned int h, unsigned int w,unsigned  int k, float v) {    
    unsigned int i;
    unsigned int j;
    unsigned int l;    
    unsigned int x;

    ip_mat *miamatrice = (ip_mat*)malloc(sizeof(ip_mat));

    /* Inizializzo la matrice con dimensioni h w e k */
    miamatrice->h = h;
    miamatrice->w = w;
    miamatrice->k = k;

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
    for(x = 0; x < k; x++){
        (miamatrice->stat+x)->min = v;
        (miamatrice->stat+x)->max = v;
        (miamatrice->stat+x)->mean = v;
    }
    return miamatrice;
}

void ip_mat_free(ip_mat *a) {
    unsigned int i;
    unsigned int j;
    /* Libero t->data */
    if(a){
        for(i = 0; i< a->h; i++) {
            for (j = 0; j< a->w; j++) {
                free(a->data[i][j]);
            }
            free(a->data[i]);
        }
        free(a->data);
        free(a->stat);
        free(a);
    }
}

void compute_stats(ip_mat * t) {   
    unsigned int x; 
    int nums;
    float sum;
    float currentmax;
    float currentmin;
    unsigned int i;
    unsigned int j;
    for (x=0; x < t->k; x++) {        
        sum = 0;
        nums = 0;
        /* Dichiara qui max e min e ad essi associa il primo valore della matrice */
        currentmax = t->data[0][0][x];
        currentmin = t->data[0][0][x];
        for(i = 0; i < t->h; i++){
            for(j = 0; j < t->w; j++){
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

void ip_mat_init_random(ip_mat * t, float mean, float var) { 
    unsigned int i;
    unsigned int j;
    unsigned int l;
    for(i=0; i<t->h; i++){
        for(j=0; j<t->w; j++){
            for(l=0; l<t->k; l++){
                t->data[i][j][l] = (get_normal_random()-mean)/var;
            }
        }
    }
}

ip_mat * ip_mat_copy(ip_mat * in) {
    ip_mat *newmat;
    unsigned int i;
    unsigned int j;
    unsigned int l;
    newmat = ip_mat_create(in->h, in->w, in->k, 0.0);
    for(i=0; i < newmat->h; i++){
        for(j=0; j < newmat->w; j++){
            for(l=0; l < newmat->k; l++){
                newmat->data[i][j][l] = in->data[i][j][l];
            }
        }
    }
    compute_stats(newmat);
    return newmat;
}

ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end) {
    /* Crea una nuova mat grande hxw dove h è la differenza tra le righe, w la differenza tra le colonne */
    ip_mat *submat;
    unsigned int i;
    unsigned int j;
    unsigned int l; 
    unsigned int o;
    unsigned int p;
    unsigned int q;/* per scorrere la nuova matrice */
    submat = ip_mat_create(row_end-row_start, col_end-col_start, t->k, 0.0);
    o = 0;
    p = 0;
    q = 0;
    /* Riempio la sottomatrice da row_start a row_end esclusa e da col_start a col_end esclusa */
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
    /* Ricalcolo le statistiche della nuova sottomatrice */
    compute_stats(submat);
    return submat;
}

ip_mat * ip_mat_concat(ip_mat * a, ip_mat * b, int dimensione) {
    ip_mat *out;
    unsigned int i;
    unsigned int j;
    unsigned int l;
    unsigned int o;
    unsigned int p;
    unsigned int q;
    if(dimensione==0){
        if(a->k == b->k && a->w == b->w){
            out = ip_mat_create((a->h+b->h), a->w, a->k, 0.0);
            for(i = 0; i < a->h; i++){
                for(j = 0; j < a->w; j++){
                    for(l=0; l< a->k; l++){
                        out->data[i][j][l]=a->data[i][j][l];
                    }
                }
            }
            o = 0;
            for(i = a->h; i<(a->h + b->h); i++){
                for(j = 0; j < b->w; j++){
                    for(l=0; l < b->k; l++){
                        out->data[i][j][l]=b->data[o][j][l];
                    }
                }
                o++;
            }
            compute_stats(out);
            return out;
        }
        else{
            printf("Errore in ip_mat_concat : le due matrici non hanno lo stesso numero di colonne o di canali\n");
            exit(1);
        }
    }
    else if(dimensione==1){
        if(a->k == b->k && a->h == b->h) {
            out = ip_mat_create(a->h, (a->w+b->w), a->k, 0.0);
            for(i = 0; i < a->h; i++){
                for(j = 0; j < a->w; j++){
                    for(l=0; l< a->k; l++){
                        out->data[i][j][l]=a->data[i][j][l];
                    }
                }
            }
            p = 0;
            for(i = 0; i < b->h; i++){
                p=0;
                for(j = a->w; j < (a->w + b->w); j++){
                    for(l=0; l < b->k; l++){
                        out->data[i][j][l]=b->data[i][p][l];
                    }
                    p++;
                }
            }
            compute_stats(out);
            return out;
        }
        else {
            printf("Errore in ip_mat_concat : le due matrici non hanno lo stesso numero di righe o canali\n");
            exit(1);
        }
    }
    else if(dimensione==2){
        if(a->h == b->h && a->w == b->w){
            out = ip_mat_create(a->h, a->w, (a->k+b->k), 0.0);
            for(i = 0; i < a->h; i++){
                for(j = 0; j < a->w; j++){
                    for(l=0; l< a->k; l++){
                        out->data[i][j][l]=a->data[i][j][l];
                    }
                }
            }
            q = 0;
            for(i = 0; i < b->h; i++){
                for(j = 0; j < b->w; j++){
                    q=0;
                    for(l=a->k; l < (a->k + b->k); l++){
                        out->data[i][j][l]=b->data[i][j][q];
                        q++;
                    }
                }
            }
            compute_stats(out);
            return out;
        }
        else{
            printf("Errore in ip_mat_concat : le due matrici non hanno lo stesso numero di colonne o di righe\n");
            exit(1);
        }
    }
    else{
        printf("Errore in ip_mat_concat : 'dimensione' deve valere 0, 1 o 2");
        exit(1);
    }
}


/**** PARTE 1: OPERAZIONI MATEMATICHE FRA IP_MAT ****/

ip_mat * ip_mat_sum(ip_mat * a, ip_mat * b) {
    ip_mat *sum;
    unsigned int i;
    unsigned int j;
    unsigned int l;
    if(a->h == b->h && a->w==b->w && a->k==b->k){
        sum = ip_mat_create(a->h, a->w, a->k, 0.0);
        for(i = 0; i < a->h; i++){
            for(j = 0; j < a->w; j++){
                for(l=0; l< a->k; l++){
                    sum->data[i][j][l] = a->data[i][j][l] + b->data[i][j][l];
                }
            }
        }
        return sum;
    }
    else{
        printf("Errore in ip_mat_sum : le due matrici devono avere le stesse dimensioni \n");
        exit(1);
    }
}

ip_mat * ip_mat_sub(ip_mat * a, ip_mat * b) {
    ip_mat *sub;
    unsigned int i;
    unsigned int j;
    unsigned int l;
    if(a->h == b->h && a->w==b->w && a->k==b->k){
        sub = ip_mat_create(a->h, a->w, a->k, 0.0);
        for(i = 0; i < a->h; i++){
            for(j = 0; j < a->w; j++){
                for(l=0; l< a->k; l++){
                    sub->data[i][j][l] = a->data[i][j][l] - b->data[i][j][l];
                }
            }
        }
        return sub;
    }
    else{
        printf("Errore in ip_mat_sub : le due matrici devono avere le stesse dimensioni \n");
        exit(1);
    }
}

ip_mat * ip_mat_mul_scalar(ip_mat *a, float c){
    ip_mat *mult;
    unsigned int i;
    unsigned int j;
    unsigned int l;
    mult = ip_mat_create(a->h, a->w, a->k, 0.0);
    for(i = 0; i < a->h; i++){
        for(j = 0; j < a->w; j++){
            for(l=0; l< a->k; l++){
                mult->data[i][j][l] = a->data[i][j][l]*c;
            }
        }
    }
    return mult;
}

ip_mat *  ip_mat_add_scalar(ip_mat *a, float c) {
    ip_mat *add;
    unsigned int i;
    unsigned int j;
    unsigned int l;
    add = ip_mat_create(a->h, a->w, a->k, 0.0);
    for(i = 0; i < a->h; i++){
        for(j = 0; j < a->w; j++){
            for(l=0; l< a->k; l++){
                add->data[i][j][l] = a->data[i][j][l]+c;
            }
        }
    }
    return add;
}

ip_mat * ip_mat_mean(ip_mat * a, ip_mat * b) {
    ip_mat *mean;
    unsigned int i;
    unsigned int j;
    unsigned int l;
    if(a->h==b->h && a->w==b->w && a->k==b->k){
        mean = ip_mat_create(a->h, a->w, a->k, 0.0);
        for(i = 0; i < a->h; i++){
            for(j = 0; j < a->w; j++){
                for(l=0; l< a->k; l++){
                    mean->data[i][j][l] = (a->data[i][j][l]+b->data[i][j][l])/2;
                }
            }
        }
        return mean;
    }
    else{
        printf("Errore in ip_mat_mean : le due matrici hanno dimensioni diverse\n");
        exit(1);
    }    
}


/**** PARTE 2: SEMPLICI OPERAZIONI SU IMMAGINI ****/

ip_mat * ip_mat_to_gray_scale(ip_mat * in) {
    ip_mat *bw;
    unsigned int i;
    unsigned int j;
    unsigned int l;
    float sum;
    float cont;
    bw = ip_mat_create(in->h, in->w, in->k, 0.0);
    for(i = 0; i < in->h; i++){
        for(j = 0; j < in->w; j++){
            sum = 0;
            cont = 0;
            for(l=0; l< in->k; l++){
                cont++;
                sum += in->data[i][j][l];
            }
            for(l=0; l< in->k; l++){
                bw->data[i][j][l]=sum/cont;
            }
        }
    }
    return bw;
}

ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha) {
    ip_mat *c;
    unsigned int i;
    unsigned int j;
    unsigned int l;
    if(a->h == b->h && a->w == b->w && a->k == b->k){
        c = ip_mat_create(a->h, a->w, a->k, 0.0);
        for(i = 0; i < a->h; i++){
            for(j = 0; j < a->w; j++){
                for(l=0; l< a->k; l++) {
                    c->data[i][j][l] = alpha*(a->data[i][j][l]) + (1-alpha)*(b->data[i][j][l]);
                }
            }
        }
        return c;
    }
    else{
        printf("Errore in ip_mat_blend : le due immagini devono avere la stessa dimensione\n");
        exit(1);
    }
}

ip_mat * ip_mat_brighten(ip_mat * a, float bright) {
    ip_mat *br;
    br = ip_mat_add_scalar(a, bright);
    return br;
}

ip_mat * ip_mat_corrupt(ip_mat * a, float amount) {
    ip_mat *out;
    float gauss_noise;
    unsigned int i;
    unsigned int j;
    unsigned int l;
    out = ip_mat_copy(a);
    for(i = 0; i < a->h; i++){
        for(j = 0; j < a->w; j++){
            for(l=0; l< a->k; l++) {
                gauss_noise = get_normal_random();
                if(a->data[i][j][l] + (gauss_noise*amount) > 255.0){
                    out->data[i][j][l] = 255.0;
                }
                if(a->data[i][j][l] + (gauss_noise*amount) < 0.0){
                    out->data[i][j][l]=0.0;
                }
                else{
                    out->data[i][j][l] = a->data[i][j][l] + (gauss_noise*amount);
                }
            }
        }
    }
    return out;
}


/**** PARTE 3: CONVOLUZIONE E FILTRI *****/

/*
* Dato un filtro, una sottomatrice delle stesse dimensioni del filtro e un canale,
* applica il calcolo della convoluzione sul pixel della matrice e il corrispondente
* pixel del filtro. Ritorna la somma di tutte le moltiplicazioni.
*/

float get_convolved_value(ip_mat *filtro, ip_mat *sottomatrice, int canale){
    float ris;
    unsigned int i;
    unsigned int j;
    ris=0.0;
    if(filtro->h == sottomatrice->h && filtro->w == sottomatrice->w){
        for(i = 0; i < filtro->h; i++){
            for(j = 0; j < filtro->w; j++){
                ris += (sottomatrice->data[i][j][canale]) * (filtro->data[i][j][0]);
            }
        }
        return ris;
    }
    else{
        printf("Errore in get_convolved_value : Il filtro e la sottomatrice hanno dimensioni h e w diverse \n");
        exit(1);
    }
}


/* Effettua la convoluzione di un ip_mat "a" con un ip_mat "f".
 * La funzione restituisce un ip_mat delle stesse dimensioni di "a".
 * */
ip_mat * ip_mat_convolve(ip_mat * a, ip_mat * f){
    int pad;
    unsigned int o;
    unsigned int p;
    unsigned int i;
    unsigned int j;
    unsigned int l;
    ip_mat *convolvedaux;
    ip_mat *convolved;
    ip_mat *auxsubset;
    if(a->h >= f->h || a->w >= f->w || a->k >= f->k){
        /* Calcolo il padding in base alle dimensioni del filtro */
        pad = ((f->h)-1)/2; 
        /* Calcolo il padding in base alle dimensioni del filtro */
        /* 
        * Le dimensioni della convolved sono:
        * 
        * PRIMA DEL PADDING
        * Altezza H convolved = Altezza H di a - (altezza H del filtro - 1)
        * Larghezza W convolved = Larghezza W di a - (Larghezza W del filtro - 1)
        *
        * DOPO IL PADDING
        * Altezza H convolved + pad
        * Larghezza W convolved + pad
        */
        convolved = ip_mat_create((a->h - (f->h - 1)), (a-> w - (f->w - 1)), a->k, 0.0);
        o = 0;
        p = 0;
        /*
        * In questi cicli innestati scorriamo il filtro sull'intera matrice a
        * e assegnamo ad ogni canale di ogni pixel di a il valore calcolato
        * dalla get_convolved_value, che prende in ingresso il filtro, la sottomatrice
        * (ottenuta da ip_mat_subset) di dimensioni f->h x f->w e il canale corrente, 
        * indicizzato da l.
        */
        for(i = 0; i <= (a->h - f->h); i++){
            p=0;
            for(j = 0; j <= (a->w - f->w); j++){
                auxsubset = ip_mat_subset(a, i, i+(f->h), j, j+(f->w));              
                for(l = 0; l < a->k; l++){     
                    /* Calcola la convoluzione del pixel [i][j] al canale l */               
                    convolved->data[o][p][l] = get_convolved_value(f, auxsubset, l);
                }
                p++;
                /* Libera il subset ausiliare */
                ip_mat_free(auxsubset);
            }
            o++;
        }
        /* Effettua il padding dopo aver convoluto l'immagine */
        convolvedaux = ip_mat_padding(convolved, pad, pad);
        ip_mat_free(convolved);
        return convolvedaux;
    }
    else{
        printf("Errore in ip_mat_convolve : il filtro è più grande dell'immagine \n");
        exit(1);
    }    
}

ip_mat * ip_mat_padding(ip_mat * a, int pad_h, int pad_w){
    unsigned int i;
    unsigned int j;
    unsigned int l;
    unsigned int o;
    unsigned int p;
    unsigned int q;

    o = 0;
    p = 0;
    q = 0;
    ip_mat *out;
    out = ip_mat_create((pad_h * 2)+ a->h, (pad_w * 2)+a->w, a->k, 0.0);
    for(i=pad_h; i<=a->h; i++) {
        p=0;
        for(j=pad_w; j<=a->w; j++) {
            q=0;
            for(l=0; l<a->k; l++) {
                out->data[i][j][l] = a->data[o][p][q];
                q++;
            }
            p++;
        }
        o++;
    }
    return out;   
}

ip_mat * create_sharpen_filter(){
    ip_mat *sharpen_filter;
    sharpen_filter = ip_mat_create(3, 3, 1, 0.0);
    set_val(sharpen_filter, 0, 1, 0, -1.0);
    set_val(sharpen_filter, 1, 0, 0, -1.0);
    set_val(sharpen_filter, 1, 1, 0, 5.0);
    set_val(sharpen_filter, 1, 2, 0, -1.0);
    set_val(sharpen_filter, 2, 1, 0, -1.0);
    return sharpen_filter;
}

ip_mat * create_edge_filter(){
    ip_mat *edge_filter;
    edge_filter = ip_mat_create(3, 3, 1, -1.0);
    set_val(edge_filter, 1, 1, 0, 8.0);
    return edge_filter;
}

ip_mat * create_emboss_filter(){
    ip_mat *emboss_filter;
    emboss_filter = ip_mat_create(3, 3, 1, 1.0);
    set_val(emboss_filter, 0, 0, 0, -2.0);
    set_val(emboss_filter, 0, 1, 0, -1.0);
    set_val(emboss_filter, 0, 2, 0, 0.0);
    set_val(emboss_filter, 1, 0, 0, -1.0);
    set_val(emboss_filter, 2, 0, 0, 0.0);
    set_val(emboss_filter, 2, 2, 0, 2.0);
    return emboss_filter;
}

ip_mat * create_average_filter(int w, int h, int k){
    float c;
    ip_mat *average_filter;
    c = 1.0/(w*h);
    average_filter = ip_mat_create(3, 3, k, c);
    return average_filter;
}

ip_mat * create_gaussian_filter(int w, int h, int k, float sigma){
    int cx;
    int cy;
    int x;
    int y;
    int i;
    int j;
    float sum;
    ip_mat *gaussian_filter;
    if(w%2 != 0 && h%2 != 0){
        cx = w/2;
        cy = h/2;
    }
    else{
        printf("Errore in create_gaussian_filter: il filtro ha dimensioni pari\n");
        exit(1);
    }
    gaussian_filter = ip_mat_create(w, h, k, 0.0);
    for(i=0; i < w; i++){
        x=i-cx;
        for(j=0; j < h; j++){
            y=j-cy;
            gaussian_filter->data[i][j][0]= (1/(2*PI*(sigma*sigma)))*exp(-((x*x)+(y*y))/(2*(sigma*sigma)));
        }
    }
    sum = 0.0;
    for(i = 0; i < w; i++){
        for(j = 0; j < h; j++){
            sum += gaussian_filter->data[i][j][0];
        }
    }
    for(i = 0; i < w; i++){
        for(j = 0; j < h; j++){
            /* set_val(gaussian_filter, i, j, 0, ((gaussian_filter->data[i][j][0])/sum)); */
            gaussian_filter->data[i][j][0] = gaussian_filter->data[i][j][0] / sum;
        }
    }
    return gaussian_filter;
}

void rescale(ip_mat * t, float new_max) {
    unsigned int i;
    unsigned int j;
    unsigned int l;
    float minimo;
    float massimo;
    compute_stats(t);
    for(i = 0; i < t->h; i++){
        for(j = 0; j < t->w; j++){
            for(l = 0; l < t->k; l++){      
                minimo = (t->stat+l)->min;
                massimo = (t->stat+l)->max;          
                t->data[i][j][l] = ((t->data[i][j][l] -  minimo)/(massimo - minimo))*new_max;
                /* if(t->data[i][j][l] < 0){
                    printf("Negativo\n");
                    // Non viene mai stampato
                } */
                
            }
        }
    }
}

/* Nell'operazione di clamping i valori <low si convertono in low e i valori >high in high.*/
void clamp(ip_mat * t, float low, float high){
    unsigned int i;
    unsigned int j;
    unsigned int l;
    for(i = 0; i < t->h; i++){
        for(j = 0; j < t->w; j++){
            for( l = 0; l < t->k; l++){
                /* if(t->data[i][j][l] > 200){
                    t->data[i][j][l]+=50;
                }
                if(t->data[i][j][l] < 200){
                    t->data[i][j][l]-=120;
                } */
                if(t->data[i][j][l] < low){
                    t->data[i][j][l] = low;
                }
                if(t->data[i][j][l] > high){
                    t->data[i][j][l] = high;
                }
            }
        }
    }
}