#include "bmp.h"
#include "ip_lib.h"
#include "ip_lib.c"
#include <math.h>
#include <stdio.h>
#include <assert.h>


Bitmap * fractal(int hxres, int hyres) {
    /* Code adapted from http://www.physics.emory.edu/faculty/weeks/software/mand.html */

    Bitmap * b = bm_create(hxres,hyres);

    double x,xx,y,cx,cy;
    int iteration,hx,hy;
    int itermax = 100;
    double magnify=1.0;

    for (hy=0;hy<hyres;hy++)  {
        for (hx=0;hx<hxres;hx++)  {
            cx = (((float)hx)/((float)hxres)-0.5)/magnify*3.0-0.7;
            cy = (((float)hy)/((float)hyres)-0.5)/magnify*3.0;
            x = 0.0; y = 0.0;
            for (iteration=0;iteration<itermax;iteration++)  {
                xx = x*x-y*y+cx;
                y = 2.0*x*y+cy;
                x = xx;
                if (x*x+y*y>100.0)  iteration = 999999;
            }
            if (iteration<99999){
                bm_set(b,hx, hy, bm_rgb(0,255,255));
            }
            else bm_set(b,hx, hy, bm_rgb(180,0,0));
        }
    }
    return b;
}

int main(){
    Bitmap *immagine = bm_create(700, 464);
    immagine = bm_load("flower.bmp");

    ip_mat *img = bitmap_to_ip_mat(immagine);
    ip_mat *edge = create_edge_filter();

    ip_mat *conv = ip_mat_convolve(img, edge);
    rescale(conv, 255.0);
    clamp(conv, 0.0, 255.0);

    Bitmap *newimg = ip_mat_to_bitmap(conv);
    bm_save(newimg, "ProvaEdge.bmp");
    
    ip_mat_show(edge);
    ip_mat_free(edge);

        
    return 0;
}
