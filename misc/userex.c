#include <stdlib.h>

void userex_(void (*f)()){
	atexit(f);
}
