/*
Author: Fabio Cannizzo
NUS FE5226 Example Program

Prompt the user to repeadetly input positive integers numbers from the console
A negative number cause termination of the input phase
Upon termination of the input phase, all inputed numbers are output to the console
Data structure used is a dynamic array
*/

#include <iostream>

// We just state the head of the function, to tell the compiler that the function is in some other CPP file
// At link stage, the linker will search across all OBJ files for this funcion
int getInteger(const char *promptMsg);

using namespace std;

struct MyArray
{
    unsigned int m_capacity;   // storage capacity
    unsigned int m_size;       // how many numbers have been stored so far
    unsigned int *m_numbers;   // allocate memory storage

    void initialize()
    {
        m_capacity = 10;   // storage capacity
        m_size = 0;        // how many numbers have been stored so far
        m_numbers = new unsigned int[m_capacity];  // allocate memory storage
    }

    void addNewNumber(unsigned userInput)  // no need to pass size and and numbers 
    {
        if (m_size == m_capacity)
            growArray();
        m_numbers[m_size++] = userInput;
    }

    // no need to pass capacity and numbers
    // return type is now void
    void growArray()
    {
        unsigned size = m_capacity;

        // increase capacity variable
        m_capacity *= 2;

        // allocate new storage space with the new larger capacity
        // note that the pointer variable 'tmp' is local to this statement block
        unsigned int *tmp = new unsigned int[m_capacity];

        // copy from previous small storage space to new large storage space
        for (unsigned i = 0u; i < size; ++i)
            tmp[i] = m_numbers[i];

        // release old storage space
        delete[] m_numbers;

        // point the 'numbers' pointer to the new storage space
        m_numbers = tmp;
    }

    void printArray()
    {
        for (unsigned i = 0u; i < m_size; ++i)
            cout << m_numbers[i] << " ";
        cout << endl;
    }


    void release()
    {
        // release memory
        delete[] m_numbers;
    }
};





int main()
{
    MyArray v;

    v.initialize();

    // we start an infinite loop, which we will stop with 'break'
    // we keep asking the user to enter a number until a negative number is entered 
    while (true) {

        int userInput = getInteger("Enter a non-negative integer or a negative to terminate: ");

        // At this point we are sure have a valid integer input

        // if the input is negative, we terminate the input loop
        if (userInput < 0)
            break;

        // At this point we are sure have a valid integer input and we know it is positive

        // At this point we are sure have a valid integer input, we know it is positive
        // and we know we have enough storage space to add it
        // So we add the new number to the existing storage
        v.addNewNumber(userInput);
    }

    // the input loop is completed.

    // we print the numbers
    v.printArray();

    // release memory
    v.release();

    return 0;
}