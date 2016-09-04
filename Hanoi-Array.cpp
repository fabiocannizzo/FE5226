/*
 Author: Fabio Cannizzo
 NUS FE5226 Example Program

 Solve the Hanoi tower puzzle
*/

#include <iostream>
#include <iomanip>

using namespace std;

const unsigned nMaxDiscs = 100;

class Tower
{
    char     m_towerId;           // tower identifier
    unsigned m_size;             // number of non-empty slots
    unsigned m_slots[nMaxDiscs]; // slots for discs

public:
    Tower(char id) : m_towerId(id), m_size(0) {}

    // push a disc at the top of the tower
    void push(unsigned discId)
    {
        m_slots[m_size++] = discId;
    }

    // pop a disc from the top of the tower
    unsigned pop()
    {
        return m_slots[--m_size];
    }

    char name() const { return m_towerId; }
};

void move(unsigned nDiscs, Tower& p1, Tower& p2, Tower& p3, unsigned& nMoves)
{
    if (nDiscs > 1) {
        // recursion
        move(nDiscs-1, p1, p3, p2, nMoves);
        move(1, p1, p2, p3, nMoves);
        move(nDiscs-1, p3, p2, p1, nMoves);
    }
    else {
        // move one disc and print message
        unsigned discid = p1.pop();
        p2.push(discid);
        cout << setw(4) << ++nMoves << ".  " 
             <<  p1.name() << " --> " << p2.name() 
             << ": " << discid << endl;
    }
}


int main()
{
    Tower p1('a');
    Tower p2('b');
    Tower p3('c');

    const unsigned nDiscs = 4;

    // initialize, all discs on tower1
    for (unsigned i = 0; i < nDiscs; ++i)
        p1.push(nDiscs-i);

    // move nDiscs discs from tower1 to tower2 using tower3
    unsigned nMoves = 0;
    move(nDiscs, p1, p2, p3, nMoves);
    cout << "Total number of moves: " << nMoves << endl; 

    return 0;
}
