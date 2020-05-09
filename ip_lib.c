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

/* Inizializza una ip_mat con dimensioni h w e k. Ogni elemento è inizializzato a v.
 * Inoltre crea un vettore di stats per contenere le statische sui singoli canali. */
ip_mat * ip_mat_create(unsigned int h, unsigned int w,unsigned  int k, float v) {    
    /* Definizione indici */
    unsigned int i;
    unsigned int j;
    unsigned int l;    
    unsigned int x;
    
    /* Allochiamo spazio per la struct ip_mat */
    ip_mat *miamatrice = (ip_mat*)malloc(sizeof(ip_mat));

    /* Inizializziamo la matrice con dimensioni h w e k */
    miamatrice->h = h;
    miamatrice->w = w;
    miamatrice->k = k;

    /* Allochiamo memoria per la struttura data per ognuna delle sue dimensioni */
    miamatrice->data=(float ***)malloc(sizeof(float **)*h);
    for(i=0;i<h;i++){
        miamatrice->data[i]=(float **)malloc(sizeof(float*)*w);
        for(j=0;j<w;j++){
            miamatrice->data[i][j]=(float *)malloc(sizeof(float)*k);
        }
    }
    
    /* Inizializziamo ogni elemento a v. */
    for(i=0;i<h;i++){
        for(j=0;j<w;j++){
            for(l=0;l<k;l++){
                miamatrice->data[i][j][l] = v;
            }
        }
    }

    /* Allochiamo memoria per il vettore stats che contiene le statische sui singoli k canali */
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

/* Libera la memoria (data, stat e la struttura) */
void ip_mat_free(ip_mat *a) {
    /* Definizione indici */
    unsigned int i;
    unsigned int j;
    
    /* Controlliamo che a non punti a NULL */
    if(a){
        /* Liberiamo la matrice data */
        for(i = 0; i< a->h; i++) {
            for (j = 0; j< a->w; j++) {
                free(a->data[i][j]);
            }
            free(a->data[i]);
        }
        free(a->data);
        /* Liberiamo a->stat e a stesso */
        free(a->stat);
        free(a);
    }
}

/* Calcola il valore minimo, il massimo e la media per ogni canale
 * e li salva dentro la struttura ip_mat stats */
void compute_stats(ip_mat * t) {   
    /* Definizione indici */
    unsigned int x; 
    unsigned int i;
    unsigned int j;
    
    /* Definizione variabili per calcolare il massimo, minimo e media*/
    int nums;
    float sum;
    float currentmax;
    float currentmin;

    /* Iteriamo sui k canali di t */
    for (x=0; x < t->k; x++) {       
        /* Inizializziamo le variabili */ 
        sum = 0;
        nums = 0;
        
        /* Inizializziamo max e min al primo valore della matrice */
        currentmax = t->data[0][0][x];
        currentmin = t->data[0][0][x];
        
        /* Iteriamo sui pixel*/
        for(i = 0; i < t->h; i++){
            for(j = 0; j < t->w; j++){
                /* Confrontiamo e riassegnamo min e max e aumentiamo i vari contatori 
                 * per la media */
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
        
        /* Associamo massimo, minimo e media alla x-esima, cioè t->k-esima, di stat */    
        (t->stat+x)->max = currentmax;
        (t->stat+x)->min = currentmin;
        (t->stat+x)->mean = sum/nums;
    }
}

/* Inizializza una ip_mat con dimensioni w h e k.
 * Ogni elemento è generato da una gaussiana con media mean e varianza var */
void ip_mat_init_random(ip_mat * t, float mean, float var) {
    /* Definizione indici */
    unsigned int i;
    unsigned int j;
    unsigned int l;

    /* Iteriamo su ogni canale di ogni pixel */
    for(i=0; i<t->h; i++){
        for(j=0; j<t->w; j++){
            for(l=0; l<t->k; l++){
                /* Sottraiamo la media alla distribuzione normale generata da
                 * get_normal_random e dividiamo il risultato per la nuova
                 * varianza "var".
                 */
                t->data[i][j][l] = (get_normal_random()-mean)/var;
            }
        }
    }
    
    /* Aggiorniamo il vettore stats */
    compute_stats(t);
}

/* Crea una copia di una ip_mat e lo restituisce in output */
ip_mat * ip_mat_copy(ip_mat * in) {
    /* Definizione indici */
    unsigned int i;
    unsigned int j;
    unsigned int l;
    
    /* Definizione e allocamento matrice di ritorno */
    ip_mat *newmat;
    newmat = ip_mat_create(in->h, in->w, in->k, 0.0);

    /* Iteriamo su ogni canale di ogni pixel */
    for(i=0; i < newmat->h; i++){
        for(j=0; j < newmat->w; j++){
            for(l=0; l < newmat->k; l++){
                /* Copiamo i valori della matrice vecchia in quella nuova*/
                newmat->data[i][j][l] = in->data[i][j][l];
            }
        }
    }
    
    /* Aggiorniamo il vettore stats */
    compute_stats(newmat);
    
    return newmat;
}

/* Restituisce una sotto-matrice, ovvero la porzione individuata da:
 * t->data[row_start...row_end][col_start...col_end][0...k]
 * La terza dimensione la riportiamo per intero, stiamo in sostanza prendendo un sottoinsieme
 * delle righe e delle colonne. */
ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end) {
    /* Definizione indici per scorrere t */
    unsigned int i;
    unsigned int j;
    unsigned int l; 
    /* Definizione indici per scorrere la nuova matrice */
    unsigned int o;
    unsigned int p;
    unsigned int q;

    /* Creiamo una nuova mat grande hxw dove h è la differenza tra le righe, w la differenza tra le colonne */
    /* Definizione matrice di ritorno */
    ip_mat *submat;
    submat = ip_mat_create(row_end-row_start, col_end-col_start, t->k, 0.0);
    
    /* Inizializziamo gli indici */
    o = 0;
    p = 0;
    q = 0;

    /* Riempiamo la sottomatrice da row_start a row_end esclusa e da col_start a col_end esclusa */
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
    
    /* Aggiorniamo il vettore stats */
    compute_stats(submat);
    
    return submat;
}

/* Concatena due ip_mat su una certa dimensione.*/
ip_mat * ip_mat_concat(ip_mat * a, ip_mat * b, int dimensione) {
    /* Definizione indici */
    unsigned int i;
    unsigned int j;
    unsigned int l;

    /* Definizione indici per scorrere la nuova matrice */
    unsigned int o;
    unsigned int p;
    unsigned int q;

    /* Definizione matrice di ritorno */
    ip_mat *out;

    if(dimensione==0){
        /* Controlliamo che le dimensioni w e k delle due matrici siano uguali */
        if(a->k == b->k && a->w == b->w){
            /* La matrice di ritorno ha h = somma delle h delle due matrici concatenate */
            out = ip_mat_create((a->h+b->h), a->w, a->k, 0.0);
            
            for(i = 0; i < a->h; i++){
                for(j = 0; j < a->w; j++){
                    for(l=0; l< a->k; l++){
                        out->data[i][j][l]=a->data[i][j][l];
                    }
                }
            }
            
            /* Inizializziamo l'indice o che scorrerà l'altezza di b */
            o = 0;

            /* Iteriamo su ogni canale di ogni pixel partendo dall'altezza a->h
             * e iterando fino alla somma delle altezze a->h + b->h */      
            for(i = a->h; i<(a->h + b->h); i++){
                for(j = 0; j < b->w; j++){
                    for(l=0; l < b->k; l++){
                        out->data[i][j][l]=b->data[o][j][l];
                    }
                }
                o++;
            }

            /* Aggiorniamo il vettore stats */
            compute_stats(out);

            return out;
        }
        else {
            /* Errore nel caso in cui i parametri passati in input come parametro non siano legali per il metodo richiamato */
            printf("Errore in ip_mat_concat : le due matrici non hanno lo stesso numero di colonne o di canali\n");
            exit(1);
        }
    }
    else if(dimensione==1){
        /* Controlliamo che le dimensioni h e k delle due matrici siano uguali */
        if(a->k == b->k && a->h == b->h) {
            /* La matrice di ritorno ha w = somma delle w delle due matrici concatenate */
            out = ip_mat_create(a->h, (a->w+b->w), a->k, 0.0);
            
            for(i = 0; i < a->h; i++){
                for(j = 0; j < a->w; j++){
                    for(l=0; l< a->k; l++){
                        out->data[i][j][l]=a->data[i][j][l];
                    }
                }
            }

            /* Inizializziamo l'indice p che scorrerà la larghezza di b */
            p = 0;

            /* Iteriamo su ogni canale di ogni pixel partendo dalla larghezza a->w
             * e iterando fino alla somma delle larghezze a->w + b->w */
            for(i = 0; i < b->h; i++){
                p=0;
                for(j = a->w; j < (a->w + b->w); j++){
                    for(l=0; l < b->k; l++){
                        out->data[i][j][l]=b->data[i][p][l];
                    }
                    p++;
                }
            }
            
            /* Aggiorniamo il vettore stats */
            compute_stats(out);
            
            return out;
        }
        else {
            /* Errore nel caso in cui i parametri passati in input come parametro non siano legali per il metodo richiamato */
            printf("Errore in ip_mat_concat : le due matrici non hanno lo stesso numero di righe o canali\n");
            exit(1);
        }
    }
    else if(dimensione==2){
        /* Controlliamo che le dimensioni h e w delle due matrici siano uguali */
        if(a->h == b->h && a->w == b->w){
            /* La matrice di ritorno ha k = somma delle k delle due matrici concatenate */
            out = ip_mat_create(a->h, a->w, (a->k+b->k), 0.0);

            for(i = 0; i < a->h; i++){
                for(j = 0; j < a->w; j++){
                    for(l=0; l< a->k; l++){
                        out->data[i][j][l]=a->data[i][j][l];
                    }
                }
            }
            
            /* Inizializziamo l'indice q che scorrerà i canali di b */
            q = 0;

            /* Iteriamo su ogni canale di ogni pixel partendo dal canale a->k
             * e iterando fino alla somma deli canali a->k + b->k */
            for(i = 0; i < b->h; i++){
                for(j = 0; j < b->w; j++){
                    q=0;
                    for(l=a->k; l < (a->k + b->k); l++){
                        out->data[i][j][l]=b->data[i][j][q];
                        q++;
                    }
                }
            }
            
            /* Aggiorniamo il vettore stats */
            compute_stats(out);
            
            return out;
        }
        else{
            /* Errore nel caso in cui i parametri passati in input come parametro non siano legali per il metodo richiamato */
            printf("Errore in ip_mat_concat : le due matrici non hanno lo stesso numero di colonne o di righe\n");
            exit(1);
        }
    }
    else{
        /* Errore nel caso in cui i parametri passati in input come parametro non siano legali per il metodo richiamato */
        printf("Errore in ip_mat_concat : 'dimensione' deve valere 0, 1 o 2");
        exit(1);
    }
}



/**** PARTE 1: OPERAZIONI MATEMATICHE FRA IP_MAT ****/


/* Esegue la somma di due ip_mat (tutte le dimensioni devono essere identiche)
 * e la restituisce in output. */
ip_mat * ip_mat_sum(ip_mat * a, ip_mat * b) {
    /* Definizione indici */
    unsigned int i;
    unsigned int j;
    unsigned int l;
    
    /* Definizione matrice di ritorno */
    ip_mat *sum;

    /* Controlliamo che le dimensioni h, w e k delle due ip_mat siano uguali */
    if(a->h == b->h && a->w==b->w && a->k==b->k){
        /* Inizializziamo la matrice di ritorno a h, w e k */
        sum = ip_mat_create(a->h, a->w, a->k, 0.0);
        
        for(i = 0; i < a->h; i++){
            for(j = 0; j < a->w; j++){
                for(l=0; l< a->k; l++){
                    /* Effettuiamo la somma */
                    sum->data[i][j][l] = a->data[i][j][l] + b->data[i][j][l];
                }
            }
        }
        
        /* Aggiorniamo il vettore stats */
        compute_stats(sum);

        return sum;
    }
    else{
        /* Errore nel caso in cui i parametri passati in input come parametro non siano legali per il metodo richiamato */
        printf("Errore in ip_mat_sum : le due matrici devono avere le stesse dimensioni \n");
        exit(1);
    }
}

/* Esegue la sottrazione di due ip_mat (tutte le dimensioni devono essere identiche)
 * e la restituisce in output. */
ip_mat * ip_mat_sub(ip_mat * a, ip_mat * b) {
    /* Definizione indici */
    unsigned int i;
    unsigned int j;
    unsigned int l;
    
    /* Definizione matrice di ritorno */
    ip_mat *sub;
    
    /* Controlliamo che le dimensioni h, w e k delle due ip_mat siano uguali */
    if(a->h == b->h && a->w==b->w && a->k==b->k){
        /* Inizializziamo la matrice di ritorno a h, w e k */
        sub = ip_mat_create(a->h, a->w, a->k, 0.0);
        
        for(i = 0; i < a->h; i++){
            for(j = 0; j < a->w; j++){
                for(l=0; l< a->k; l++){
                    /* Effettuiamo la sottrazione */
                    sub->data[i][j][l] = a->data[i][j][l] - b->data[i][j][l];
                }
            }
        }
        
        /* Aggiorniamo il vettore stats */
        compute_stats(sub);
        
        return sub;
    }
    else{
        /* Errore nel caso in cui i parametri passati in input come parametro non siano legali per il metodo richiamato */
        printf("Errore in ip_mat_sub : le due matrici devono avere le stesse dimensioni \n");
        exit(1);
    }
}

/* Moltiplica un ip_mat per uno scalare c. Si moltiplica c per tutti gli elementi di "a"
 * e si salva il risultato in un nuovo tensore in output. */
ip_mat * ip_mat_mul_scalar(ip_mat *a, float c){
    /* Definizione indici */
    unsigned int i;
    unsigned int j;
    unsigned int l;

    /* Definiamo e inizializziamo la matrice di ritorno a h, w e k */
    ip_mat *mult;
    mult = ip_mat_create(a->h, a->w, a->k, 0.0);

    for(i = 0; i < a->h; i++){
        for(j = 0; j < a->w; j++){
            for(l=0; l< a->k; l++){
                /* Effettuiamo la moltiplicazione per lo scalare */
                mult->data[i][j][l] = a->data[i][j][l]*c;
            }
        }
    }

    /* Aggiorniamo il vettore stats */
    compute_stats(mult);
    
    return mult;
}

/* Aggiunge ad un ip_mat uno scalare c e lo restituisce in un nuovo tensore in output. */
ip_mat *  ip_mat_add_scalar(ip_mat *a, float c) {
    /* Definizione indici */
    unsigned int i;
    unsigned int j;
    unsigned int l;
    
    /* Definiamo e inizializziamo la matrice di ritorno a h, w e k */
    ip_mat *add;
    add = ip_mat_create(a->h, a->w, a->k, 0.0);

    for(i = 0; i < a->h; i++){
        for(j = 0; j < a->w; j++){
            for(l=0; l< a->k; l++){
                /* Effettuiamo la somma con uno scalare */
                add->data[i][j][l] = a->data[i][j][l]+c;
            }
        }
    }

    /* Aggiorniamo il vettore stats */
    compute_stats(add);
    
    return add;
}

/* Calcola la media di due ip_mat a e b e la restituisce in output.*/
ip_mat * ip_mat_mean(ip_mat * a, ip_mat * b) {
    /* Definizione indici */
    unsigned int i;
    unsigned int j;
    unsigned int l;

    /* Definizione matrice di ritorno */
    ip_mat *mean;
    
    /* Controlliamo che le dimensioni h, w e k delle due ip_mat siano uguali */
    if(a->h==b->h && a->w==b->w && a->k==b->k){
        /* Inizializziamo la matrice di ritorno a h, w e k */
        mean = ip_mat_create(a->h, a->w, a->k, 0.0);

        for(i = 0; i < a->h; i++){
            for(j = 0; j < a->w; j++){
                for(l=0; l< a->k; l++){
                    /* Calcoliamo la media tra i due canali dei due pixel corrispondenti delle due matrici */
                    mean->data[i][j][l] = (a->data[i][j][l]+b->data[i][j][l])/2;
                }
            }
        }

        /* Aggiorniamo il vettore stats */
        compute_stats(mean);
        
        return mean;
    }
    else{
        /* Errore nel caso in cui i parametri passati in input come parametro non siano legali per il metodo richiamato */
        printf("Errore in ip_mat_mean : le due matrici hanno dimensioni diverse\n");
        exit(1);
    }    
}



/**** PARTE 2: SEMPLICI OPERAZIONI SU IMMAGINI ****/


/* Converte un'immagine RGB ad una immagine a scala di grigio.
 * Quest'operazione viene fatta calcolando la media per ogni pixel sui 3 canali
 * e creando una nuova immagine avente per valore di un pixel su ogni canale la media appena calcolata.
 * Avremo quindi che tutti i canali saranno uguali. */
ip_mat * ip_mat_to_gray_scale(ip_mat * in) {
    /* Definizione indici */
    unsigned int i;
    unsigned int j;
    unsigned int l;

    /* Definizione variabile per calcolare la media */
    float sum;
    float cont;
    
    /* Definiamo e inizializziamo la matrice di ritorno a h, w e k */
    ip_mat *bw;
    bw = ip_mat_create(in->h, in->w, in->k, 0.0);
    
    for(i = 0; i < in->h; i++){
        for(j = 0; j < in->w; j++){
            /* Inizializziamo le variabili */ 
            sum = 0;
            cont = 0;
            
            /* Per ogni canale calcoliamo la somma dei valori e poi la media tra di essi */
            for(l=0; l< in->k; l++){
                cont++;
                sum += in->data[i][j][l];
            }

            /* Riassegnamo ad ogni canale la media ottenuta */
            for(l=0; l< in->k; l++){
                bw->data[i][j][l]=sum/cont;
            }
        }
    }
    
    /* Aggiorniamo il vettore stats */
    compute_stats(bw);
    
    return bw;
}

/* Effettua la fusione (combinazione convessa) di due immagini */
ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha) {
    /* Definizione indici */
    unsigned int i;
    unsigned int j;
    unsigned int l;

    /* Definizione matrice di ritorno */
    ip_mat *c;
    
    /* Controlliamo che le dimensioni h, w e k delle due ip_mat siano uguali */
    if(a->h == b->h && a->w == b->w && a->k == b->k){
        /* Inizializziamo la matrice di ritorno a h, w e k */
        c = ip_mat_create(a->h, a->w, a->k, 0.0);
        
        for(i = 0; i < a->h; i++){
            for(j = 0; j < a->w; j++){
                for(l=0; l< a->k; l++) {
                    /* Applichiamo ad ogni canale di ogni pixel la formula:  
                     * (alpha*valore del canale del pixel di a) + (1-alpha)*(valore del canale del pixel di)
                     */
                    c->data[i][j][l] = alpha*(a->data[i][j][l]) + (1-alpha)*(b->data[i][j][l]);

                    /* Controlliamo che nel caso il valore del canale del pixel calcolato sia maggiore di 255 
                     * esso venga forzato a 255 */
                    if(c->data[i][j][l] > 255.0){
                        c->data[i][j][l] = 255.0;
                    }
                    
                    /* Controlliamo che nel caso il valore del canale del pixel calcolato sia minore di 0 
                     * esso venga forzato a 0 */
                    if(c->data[i][j][l] < 0.0){
                        c->data[i][j][l] = 0.0;
                    }
                }
            }
        }
        
        /* Aggiorniamo il vettore stats */
        compute_stats(c);
        
        return c;
    }
    else{
        /* Errore nel caso in cui i parametri passati in input come parametro non siano legali per il metodo richiamato */
        printf("Errore in ip_mat_blend : le due immagini devono avere la stessa dimensione\n");
        exit(1);
    }
}

/* Operazione di brightening: aumenta la luminosità dell'immagine
 * aggiunge ad ogni pixel un certo valore */
ip_mat * ip_mat_brighten(ip_mat * a, float bright) {
    /* Definiamo e inizializziamo la matrice di ritorno a h, w e k */
    ip_mat *br;
    br = ip_mat_add_scalar(a, bright);
    
    /* Ristabiliamo i valori ottenuti dalla somma appena effettuata nell'intervallo [0, 255] */
    clamp(br, 0.0, 255.0);
    
    /* Aggiorniamo il vettore stats */
    compute_stats(br);
    
    return br;
}

/* Operazione di corruzione con rumore gaussiano:
 * Aggiunge del rumore gaussiano all'immagine, il rumore viene enfatizzato
 * per mezzo della variabile amount.
 * out = a + gauss_noise*amount  */
ip_mat * ip_mat_corrupt(ip_mat * a, float amount) {    
    /* Definizione indici */
    unsigned int i;
    unsigned int j;
    unsigned int l;

    /* Definizione variabile per il rumore gaussiano */
    float gauss_noise;

    /* Definizione matrice di ritorno */
    ip_mat *out;
    out = ip_mat_copy(a);
    
    for(i = 0; i < a->h; i++){
        for(j = 0; j < a->w; j++){
            for(l=0; l< a->k; l++) {
                /* Assegnamo a gauss_noise un valore ottenuto dalla distribuzione normale di get_normal_random */
                gauss_noise = get_normal_random();
                
                /* Controlliamo che nel caso il valore del canale del pixel calcolato sia maggiore di 255 
                 * esso venga forzato a 255 */
                if(a->data[i][j][l] + (gauss_noise*amount) > 255.0) {
                    out->data[i][j][l] = 255.0;
                }

                /* Controlliamo che nel caso il valore del canale del pixel calcolato sia minore di 0 
                 * esso venga forzato a 0 */
                if(a->data[i][j][l] + (gauss_noise*amount) < 0.0) {
                    out->data[i][j][l]=0.0;
                }
                
                /* Caso in cui il risultato dell'operazione sia legalmente compreso nell'intervallo [0,255] */
                else {
                    out->data[i][j][l] = a->data[i][j][l] + (gauss_noise*amount);
                }
            }
        }
    }
    
    /* Aggiorniamo il vettore stats */
    compute_stats(out);
    
    return out;
}



/**** PARTE 3: CONVOLUZIONE E FILTRI *****/


/* Dato un filtro, una sottomatrice delle stesse dimensioni del filtro e un canale,
* applica il calcolo della convoluzione sul pixel della matrice e il corrispondente
* pixel del filtro. Ritorna la somma di tutte le moltiplicazioni. */
float get_convolved_value(ip_mat *filtro, ip_mat *sottomatrice, int canale){
    /* Definizione indici */
    unsigned int i;
    unsigned int j;

    /* Definiamo e inizializziamo il valore di ritorno */
    float ris;
    ris=0.0;

    /* Controlliamo che le dimensioni h, w e k delle due ip_mat siano uguali */
    if(filtro->h == sottomatrice->h && filtro->w == sottomatrice->w){
        for(i = 0; i < filtro->h; i++){
            for(j = 0; j < filtro->w; j++){
                /* Sommiamo al valore di ritorno la moltiplicazione tra le corrispondenti celle
                 * della sottomatrice e del filtro */
                ris += (sottomatrice->data[i][j][canale]) * (filtro->data[i][j][0]);
            }
        }
        return ris;
    }
    else{
        /* Errore nel caso in cui i parametri passati in input come parametro non siano legali per il metodo richiamato */
        printf("Errore in get_convolved_value : Il filtro e la sottomatrice hanno dimensioni h e w diverse \n");
        exit(1);
    }
}

/* Effettua la convoluzione di un ip_mat "a" con un ip_mat "f".
 * La funzione restituisce un ip_mat delle stesse dimensioni di "a". */
ip_mat * ip_mat_convolve(ip_mat * a, ip_mat * f){
    /* Definizione indici */
    unsigned int i;
    unsigned int j;
    unsigned int l;

    /* Definizione indici per scorrere la nuova matrice */
    unsigned int o;
    unsigned int p;
    
    /* Definizione variabili per il padding */
    int pad;
    
    /* Definizamo la ip_mat di supporto per evitare perdite di memoria (vedere riga 900) */
    ip_mat *convolvedaux;
    
    /* Definizione di ip_mat ausiliari e di ritorno*/
    ip_mat *convolved;
    ip_mat *auxsubset;
    
    /* Controlliamo che le dimensioni della matrice "a" siano maggiori o uguali a quelle del filtro */
    if(a->h >= f->h || a->w >= f->w || a->k >= f->k){
        /* Calcolo il padding in base alle dimensioni del filtro */
        pad = ((f->h)-1)/2;
        
        /* Le dimensioni della convolved sono:
        * 
        * PRIMA DEL PADDING
        * Altezza H convolved = Altezza H di a - (altezza H del filtro - 1)
        * Larghezza W convolved = Larghezza W di a - (Larghezza W del filtro - 1)
        *
        * DOPO IL PADDING
        * Altezza H convolved + pad
        * Larghezza W convolved + pad
        */
        
        /* Inizializziamo la matrice di supporto */
        convolvedaux = ip_mat_create((a->h - (f->h - 1)), (a-> w - (f->w - 1)), a->k, 0.0);
        
        /* Inizializziamo gli indici */
        o = 0;
        p = 0;
        
        /* In questi cicli innestati scorriamo il filtro sull'intera matrice "a"
        * e assegnamo ad ogni canale di ogni pixel di "a" il valore calcolato
        * dalla get_convolved_value, che prende in ingresso il filtro, la sottomatrice
        * (ottenuta da ip_mat_subset) di dimensioni f->h x f->w e il canale corrente, 
        * indicizzato da l. */
        for(i = 0; i <= (a->h - f->h); i++){
            p=0;
            for(j = 0; j <= (a->w - f->w); j++){
                /* Inizializziamo la sottomatrice auxsubset con le dimensioni del filtro
                 * a partire dai valori degli indici i e j appena iterati su "a" */
                auxsubset = ip_mat_subset(a, i, i+(f->h), j, j+(f->w));              
                
                for(l = 0; l < a->k; l++){     
                    /* Calcola la convoluzione del pixel [i][j] al canale l */               
                    convolvedaux->data[o][p][l] = get_convolved_value(f, auxsubset, l);
                }
                p++;
                
                /* Libera il subset ausiliare */
                ip_mat_free(auxsubset);
            }
            o++;
        }
        
        /* Effettua il padding dopo aver convoluto l'immagine */
        convolved = ip_mat_padding(convolvedaux, pad, pad);
        /* Liberiamo convolvedaux */
        ip_mat_free(convolvedaux);
        
        /* Aggiorniamo il vettore stats */
        compute_stats(convolved);
        
        return convolved;
    }
    else{
        /* Errore nel caso in cui i parametri passati in input come parametro non siano legali per il metodo richiamato */
        printf("Errore in ip_mat_convolve : il filtro è più grande dell'immagine \n");
        exit(1);
    }    
}

/* Aggiunge un padding all'immagine. Il padding verticale è pad_h mentre quello
 * orizzontale è pad_w. */
ip_mat * ip_mat_padding(ip_mat * a, int pad_h, int pad_w){
    /* Definizione indici */
    unsigned int i;
    unsigned int j;
    unsigned int l;
    
    /* Definizione indici per scorrere la nuova matrice */
    unsigned int o;
    unsigned int p;
    unsigned int q;

    /* Inizializziamo gli indici */
    o = 0;
    p = 0;
    q = 0;
    
    /* Definizione e allocamento matrice di ritorno.
     * Le dimensioni h e w saranno entrambe sommate a 2*grandezza del padding 
     * Inizializziamo la matrice di ritorno a 0.0 e ne riempiremo la sottomatrice
     * centrale con i valori della matrice originale su cui effettuare il padding. */
    ip_mat *out;
    out = ip_mat_create((pad_h * 2)+ a->h, (pad_w * 2)+a->w, a->k, 0.0);
    
    /* Riempiamo la sottomatrice centrale, quindi cominciamo ad iterarne l'altezza
     * a partire dalla misura del padding e arriviamo all'altezza effettiva di "a" */
    for(i=pad_h; i<=a->h; i++) {
        p=0;

        /* Riempiamo la sottomatrice centrale, quindi cominciamo ad iterarne la larghezza
         * a partire dalla misura del padding e arriviamo alla larghezza effettiva di "a" */
        for(j=pad_w; j<=a->w; j++) {
            q=0;
            for(l=0; l<a->k; l++) {
                /* Effettuiamo la copia. */
                out->data[i][j][l] = a->data[o][p][q];
                q++;
            }
            p++;
        }
        o++;
    }
    
    /* Aggiorniamo il vettore stats */
    compute_stats(out);
    
    return out;
}

/* Crea un filtro di sharpening */
ip_mat * create_sharpen_filter(){
    /* Definizione e inizializzazione matrice di ritorno secondo gli standard definiti nella consegna */
    ip_mat *sharpen_filter;
    sharpen_filter = ip_mat_create(3, 3, 1, 0.0);
    
    /* Settiamo i valori della nostra matrice a quelli del filtro di sharpening */
    set_val(sharpen_filter, 0, 1, 0, -1.0);
    set_val(sharpen_filter, 1, 0, 0, -1.0);
    set_val(sharpen_filter, 1, 1, 0, 5.0);
    set_val(sharpen_filter, 1, 2, 0, -1.0);
    set_val(sharpen_filter, 2, 1, 0, -1.0);

    /* Aggiorniamo il vettore stats */
    compute_stats(sharpen_filter);
    
    return sharpen_filter;
}

/* Crea un filtro per rilevare i bordi */
ip_mat * create_edge_filter(){
    /* Definizione e inizializzazione matrice di ritorno secondo gli standard definiti nella consegna */
    ip_mat *edge_filter;
    edge_filter = ip_mat_create(3, 3, 1, -1.0);
    
    /* Settiamo i valori della nostra matrice a quelli del filtro "edge" */
    set_val(edge_filter, 1, 1, 0, 8.0);
    
    /* Aggiorniamo il vettore stats */
    compute_stats(edge_filter);
    
    return edge_filter;
}

/* Crea un filtro per aggiungere profondità */
ip_mat * create_emboss_filter(){
    /* Definizione e inizializzazione matrice di ritorno secondo gli standard definiti nella consegna */
    ip_mat *emboss_filter;
    emboss_filter = ip_mat_create(3, 3, 1, 1.0);
    
    /* Settiamo i valori della nostra matrice uguali a quelli del filtro "emboss"*/
    set_val(emboss_filter, 0, 0, 0, -2.0);
    set_val(emboss_filter, 0, 1, 0, -1.0);
    set_val(emboss_filter, 0, 2, 0, 0.0);
    set_val(emboss_filter, 1, 0, 0, -1.0);
    set_val(emboss_filter, 2, 0, 0, 0.0);
    set_val(emboss_filter, 2, 2, 0, 2.0);
    
    /* Aggiorniamo il vettore stats */
    compute_stats(emboss_filter);
    
    return emboss_filter;
}

/* Crea un filtro medio per la rimozione del rumore */
ip_mat * create_average_filter(int w, int h, int k){
    /* Definizione variabile che conterrà il valore medio con cui popolare la matrice filtro */
    float c;
    
    /* Definizione matrice di ritorno */
    ip_mat *average_filter;

    /* Calcoliamo il valore medio con cui popolare la matrice filtro */
    c = 1.0/(w*h);
    average_filter = ip_mat_create(3, 3, k, c);
    
    /* Aggiorniamo il vettore stats */
    compute_stats(average_filter);
    
    return average_filter;
}

/* Crea un filtro gaussiano per la rimozione del rumore */
ip_mat * create_gaussian_filter(int w, int h, int k, float sigma){
    /* Definizione indici che indicano le coordinate centrali del filtro */
    int cx;
    int cy;

    /* Definizine indici per il calcolo della distanza dal centro del filtro di ciascuna cella */
    int x;
    int y;

    /* Definizine indici */
    int i;
    int j;

    /* Definizione variabile per la somma */
    float sum;

    /* Definizione matrice di ritorno */
    ip_mat *gaussian_filter;
    
    /* Controlliamo che le dimensioni del filtro siano dispari */
    if(w%2 != 0 && h%2 != 0){
        cx = w/2;
        cy = h/2;
    } 
    /*********************************************************************************************************************/
    else{
        /* Errore nel caso in cui i parametri passati in input come parametro non siano legali per il metodo richiamato */
        printf("Errore in create_gaussian_filter: il filtro ha dimensioni pari\n");
        exit(1);
    }
    
    /* Inizializzazione della matrice filtro */
    gaussian_filter = ip_mat_create(w, h, k, 0.0);
    
    for(i=0; i < w; i++){
        /* Calcoliamo la distanza x dalla cella centrale del filtro */
        x=i-cx;
        for(j=0; j < h; j++){
            /* Calcoliamo la distanza y dalla cella centrale del filtro */
            y=j-cy;
            
            /*Calcoliamo il valore del filtro in base alla formula indicata nel file .pdf che illustra il progetto*/
            gaussian_filter->data[i][j][0]= (1/(2*PI*(sigma*sigma)))*exp(-((x*x)+(y*y))/(2*(sigma*sigma)));
        }
    }
    /* Inizializziamo sum */
    sum = 0.0;

    for(i = 0; i < w; i++){
        for(j = 0; j < h; j++){
            /* Sommiamo tutti i valori delle celle del filtro */
            sum += gaussian_filter->data[i][j][0];
        }
    }

    for(i = 0; i < w; i++){
        for(j = 0; j < h; j++){
            /* Dividiamo tutte le celle del filtro per "sum" */
            gaussian_filter->data[i][j][0] = gaussian_filter->data[i][j][0] / sum;
        }
    }
    
    /* Aggiorniamo il vettore stats */
    compute_stats(gaussian_filter);
    return gaussian_filter;
}

/* Effettua una riscalatura dei dati tale che i valori siano in [0,new_max].
 * Utilizzate il metodo compute_stat per ricavarvi il min, max per ogni canale.
 *
 * I valori sono scalati tramite la formula (valore-min)/(max - min)
 *
 * Si considera ogni indice della terza dimensione indipendente, quindi l'operazione
 * di scalatura va ripetuta per ogni "fetta" della matrice 3D.
 * Successivamente moltiplichiamo per new_max gli elementi della matrice in modo da ottenere un range
 * di valori in [0,new_max]. */
void rescale(ip_mat * t, float new_max) {
    /* Definizione indici */
    unsigned int i;
    unsigned int j;
    unsigned int l;
    float minimo;
    float massimo;
    /* Aggiorniamo il vettore stats */
    compute_stats(t);
    for(i = 0; i < t->h; i++){
        for(j = 0; j < t->w; j++){
            for(l = 0; l < t->k; l++){      
                minimo = (t->stat+l)->min;
                massimo = (t->stat+l)->max; 
                /* Calcoliamo il valore di ciascun pixel in ciascun canale con la formula 
                (valore-min)/(max - min)*/         
                t->data[i][j][l] = ((t->data[i][j][l] -  minimo)/(massimo - minimo))*new_max;
            }
        }
    }
}

/* Nell'operazione di clamping i valori <low si convertono in low e i valori >high in high.*/
void clamp(ip_mat * t, float low, float high){
    /* Definizione indici */
    unsigned int i;
    unsigned int j;
    unsigned int l;
    
    for(i = 0; i < t->h; i++){
        for(j = 0; j < t->w; j++){
            for( l = 0; l < t->k; l++){
                if(t->data[i][j][l] < low){
                    t->data[i][j][l] = low;
                }
                if(t->data[i][j][l] > high){
                    t->data[i][j][l] = high;
                }
            }
        }
    }

    /* Aggiorniamo il vettore stats */
    compute_stats(t);
}
