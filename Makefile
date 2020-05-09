objects = bmp.o
eseguibile = execute
cflags = -Wall --ansi --pedantic -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra main_iplib.c
libs = -lm

$(eseguibile) : $(objects)
	gcc -Wall -c bmp.c
	gcc $(cflags) -o $(eseguibile) $(objects) $(libs)

clean:
	rm $(eseguibile)  $(objects)


