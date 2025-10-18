#include <iostream>
#include <vector>

using namespace std;

void print (const char *name, const vector<int>& v)
{
    cout << name << "(size="<<v.size() <<",cap=" << v.capacity()<<"): ";
    for (auto& x : v)
        cout << x << ',';
    cout << '\n';
}

int main()
{
    vector<int> v(4);
    print("v",v);

    v.push_back(4);
    v.push_back(5);
    print("v",v);

    v.pop_back();
    print("v",v);

    v.resize(8);
    print("v",v);

    v.clear();
    print("v",v);
    return 0;
}

