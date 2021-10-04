bin/radexi : main.o linenoise.o dialogue.o options.o readdata.o moldata.o bgfield.o
	gcc -Wall -Wextra -pedantic -I./inc -lm bin/main.o bin/linenoise.o bin/moldata.o bin/dialogue.o bin/readdata.o bin/bgfield.o bin/options.o -o bin/radexi

main.o : src/main.c inc/radexi.h
	gcc -Wall -Wextra -pedantic -I./inc -c src/main.c -o bin/main.o

options.o : src/options.c inc/radexi.h
	gcc -Wall -Wextra -pedantic -I./inc -c src/options.c -o bin/options.o

readdata.o : src/readdata.c inc/radexi.h
	gcc -Wall -Wextra -pedantic -I./inc -c src/readdata.c -o bin/readdata.o

dialogue.o : src/dialogue.c inc/radexi.h inc/linenoise.h
	gcc -Wall -Wextra -pedantic -I./inc -lm -c src/dialogue.c -o bin/dialogue.o

moldata.o : src/moldata.c inc/radexi.h inc/linenoise.h
	gcc -Wall -Wextra -pedantic -I./inc -lm -c src/moldata.c -o bin/moldata.o

bgfield.o : src/bgfield.c inc/radexi.h inc/linenoise.h
	gcc -Wall -Wextra -pedantic -I./inc -lm -c src/bgfield.c -o bin/bgfield.o

linenoise.o : src/linenoise.c inc/linenoise.h
	gcc -Wall -Wextra -pedantic -I./inc -c src/linenoise.c -o bin/linenoise.o
