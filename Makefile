bin/radexi : main.o linenoise.o dialogue.o options.o moldata.o
	gcc -Wall -Wextra -pedantic -I./inc -lm bin/main.o bin/linenoise.o bin/moldata.o bin/dialogue.o bin/options.o -o bin/radexi

main.o : src/main.c inc/radexi.h
	gcc -Wall -Wextra -pedantic -I./inc -c src/main.c -o bin/main.o

options.o : src/options.c inc/radexi.h
	gcc -Wall -Wextra -pedantic -I./inc -c src/options.c -o bin/options.o

dialogue.o : src/dialogue.c inc/radexi.h inc/linenoise.h
	gcc -Wall -Wextra -pedantic -I./inc -lm -c src/dialogue.c -o bin/dialogue.o

moldata.o : src/moldata.c inc/radexi.h inc/linenoise.h
	gcc -Wall -Wextra -pedantic -I./inc -lm -c src/moldata.c -o bin/moldata.o

linenoise.o : src/linenoise.c inc/linenoise.h
	gcc -Wall -Wextra -pedantic -I./inc -c src/linenoise.c -o bin/linenoise.o
