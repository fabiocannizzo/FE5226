/*
 Author: Fabio Cannizzo
 NUS FE5226 Example Program

 GetInteger function head
*/

// We just state the head of the function, to tell the compiler that the function is in some other CPP file
// At link stage, the linker will search across all OBJ files for this funcion
int getInteger(const char *promptMsg);
