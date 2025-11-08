#include <iostream>
#include <limits>
#include <vector>
#include <numbers>
#include <memory>

using namespace std;

struct Shape {
    virtual double perimeter() const = 0; // pure virtual
    virtual const char *name() const = 0; // pure virtual
    virtual ~Shape() {}
};

struct Circle : Shape {
    double m_radius;
    Circle( double r) : m_radius(r) {}
    double perimeter() const override { return 2.0 * numbers::pi * m_radius; }
    const char *name() const override { return "circle";}
    virtual ~Circle() { cout << "destroying " << name() << "\n"; }
};

struct Square : Shape {
    double m_side;
    Square(double l) : m_side(l) {}
    double perimeter() const override { return 4 * m_side; }
    const char *name() const override { return "square";}
    virtual ~Square() { cout << "destroying " << name() << "\n"; }
};


int main()
{
    unique_ptr<Shape> s[] = { make_unique<Circle>(10), make_unique<Square>(5) };

    for (auto& p : s )
        cout << "object of type " << p->name() << " has perimter " << p->perimeter() << "\n";

    return 0;
}
