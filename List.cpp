#include <iostream>
#include <list>

using namespace std;


template <typename C>
void print(const C& c)
{
    for (auto& v : c)
        cout << v << " ";
    cout << endl;
}

int main()
{
    list<int> l;

    for (size_t i = 0; i < 10; ++i)
        l.push_back(rand() % 6);
    print(l);

    for (auto i = l.crbegin(); i != l.crend(); ++i)
        cout << *i << " ";
    cout << endl;

    // 4th element
    auto it = l.begin();
    for (size_t i = 0; i < 3; ++i)
        ++it;
    cout << "4th: " << *it << endl;

    l.pop_back();
    l.pop_front();
    print(l);

    for (auto i = l.begin(); i != l.end(); i++)
        if (*i > 2) {
            l.erase(i);
            break;
        }
    print(l);

    for (auto i = l.begin(); i != l.end(); i++)
        if (*i > 2) {
            l.insert(i, 100);
            break;
        }
    print(l);

    l.sort();
    print(l);

    return 0;
}
