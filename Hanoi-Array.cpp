/*
 Author: Fabio Cannizzo
 NUS FE5226 Example Program

 Solve the Hanoi tower puzzle
*/

#include <iostream>
#include <iomanip>

using namespace std;

const unsigned nMaxDiscs = 100;

struct Pole
{
    Pole(char id) : m_poleId(id), m_size(0) {}

    // push a disc at the top of the pole
    void push(unsigned discId)
    {
        m_slots[m_size++] = discId;
    }

    // pop a disc from the top of the pole
    unsigned pop()
    {
        return m_slots[--m_size];
    }

    char name() const { return m_poleId; }

private:
    char     m_poleId;           // pole identifier
    unsigned m_size;             // number of non-empty slots
    unsigned m_slots[nMaxDiscs]; // slots for discs
};

void move(unsigned nDiscs, Pole& p1, Pole& p2, Pole& p3, unsigned& nMoves)
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
    Pole p1('a');
    Pole p2('b');
    Pole p3('c');

    const unsigned nDiscs = 4;

    // initialize, all discs on pole1
    for (unsigned i = 0; i < nDiscs; ++i)
        p1.push(nDiscs-i);

    // move nDiscs discs from pole1 to pole2 using pole3
    unsigned nMoves = 0;
    move(nDiscs, p1, p2, p3, nMoves);
    cout << "Total number of moves: " << nMoves << endl; 

    return 0;
}
