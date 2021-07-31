bin/radexi : main.o
	gcc -Wall -Wextra -pedantic bin/main.o -o bin/radexi

main.o : src/main.c
	gcc -Wall -Wextra -pedantic -c src/main.c -o bin/main.o
