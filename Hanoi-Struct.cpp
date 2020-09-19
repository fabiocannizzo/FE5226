#include <iostream>
using namespace std;

struct Pole {
    unsigned disk_ids[5];  // indices of the disks (5 disks maximum)
    unsigned n_disks = 0;       // number of disks in the pole
};

unsigned RemoveTopDisk(Pole& p)
{
    return p.disk_ids[--p.n_disks];
}

void AddTopDisk(Pole& p, unsigned id)
{
    p.disk_ids[p.n_disks++] = id;
}

// move n disks from one pole to the other
void moveDisks(Pole poles[3], unsigned from, unsigned to, unsigned aux, unsigned n) {
    if (n == 1) {
        // move top disk
        auto disk_id = RemoveTopDisk(poles[from]);
        AddTopDisk(poles[to], disk_id);
        cout << "move disk " << disk_id << " from pole " << from << " to pole " << to << endl;
    }
    else {
        moveDisks(poles, from, aux, to, n - 1);
        moveDisks(poles, from, to, aux, 1);   // move one disk
        moveDisks(poles, aux, to, from, n - 1);
    }
}

int main()
{
    Pole poles[3];

    const unsigned nDisks = 3;

    // init pole[0] with 3 disks
    for (unsigned i = 0; i < nDisks; ++i)
        AddTopDisk(poles[0], nDisks - i);

    // print moves
    moveDisks(poles, 0, 2, 1, nDisks);

    return 0;
}
