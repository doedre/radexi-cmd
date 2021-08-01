bin/radexi : main.o dialogue.o options.o
	gcc -Wall -Wextra -pedantic -I./inc bin/main.o bin/dialogue.o bin/options.o -o bin/radexi

main.o : src/main.c inc/radexi.h
	gcc -Wall -Wextra -pedantic -I./inc -c src/main.c -o bin/main.o

options.o : src/options.c inc/radexi.h
	gcc -Wall -Wextra -pedantic -I./inc -c src/options.c -o bin/options.o

dialogue.o : src/dialogue.c inc/radexi.h
	gcc -Wall -Wextra -pedantic -I./inc -c src/dialogue.c -o bin/dialogue.o
