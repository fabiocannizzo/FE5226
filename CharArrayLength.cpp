#include <iostream>
#include <cstring>

using namespace std;

int main()
{
    const char msg[10] = "Hello";

    // print all elements of the array, converted to integers
    cout << "Binary content: ";
    for (int i = 0; i < 10; ++i)
        cout << static_cast<int>(msg[i]) << " ";
    cout << endl;

    // print array interpreted as a null-terminated string
    cout << "String: " << msg << endl;

    // print size of the array ion bytes
    cout << "Size in bytes: " << sizeof(msg) << endl;

    // print length of array computed with the CRT function strlen
    cout << "String length: " << strlen(msg) << endl;
}