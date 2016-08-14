/*
Author: Fabio Cannizzo
NUS FE5226 Example Program

GetInteger function
*/

#include <iostream>

using namespace std;

// obtain a valid integer from the console
int getInteger(const char *promptMsg)
{
    int userInput;

    bool invalidInput;  // we assume the input will be valid
    do {  // this loop is executed at least once
        invalidInput = false;  // we assume the input will be valid

                               // we prompt the user to enter a number
                               // note there is no end of line
        cout << promptMsg;
        cin >> userInput;

        if (!cin) {  // was the input invalid?
            cout << "Invalid input\n";
            // discard the bad input and clean up buffers 
            cin.clear();
            cin.ignore(10, '\n');
            invalidInput = true; // this will cause the loop to continue
        }
    } while (invalidInput);

    return userInput;
}

