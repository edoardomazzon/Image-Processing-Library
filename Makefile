# VECCHIO MAKEFILE FUNZIONANTE

#objects = bmp.o
#eseguibile = execute
#cflags = -Wall --ansi --pedantic -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra main_iplib.c
#libs = -lm

#$(eseguibile) : $(objects)
#	gcc -Wall -c bmp.c
#	gcc $(cflags) -o $(eseguibile) $(objects) $(libs)

#clean:
#   rm $(eseguibile)  $(objects)


#NUOVO MAKEFILE
image_processer: bmp.o
	gcc -Wall --ansi --pedantic -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra main_iplib.c -o execute bmp.o -lm

bmp.o: bmp.c bmp.h
	gcc -Wall -c bmp.c
	


	

	

