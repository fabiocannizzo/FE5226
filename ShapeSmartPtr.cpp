#include <iostream>
#include <algorithm>
#include <vector>
#include <iostream>
#include <memory>

using namespace std;

struct Shape {
   virtual double perimeter() const = 0;  // pure virtual method
};

struct Circle : Shape {
   double radius;
   Circle( double r) : radius(r) {}
   virtual double perimeter() const { return 2*3.14*radius; }
};

struct Square : Shape {
   double side;
   Square(double l) : side(l) {}
   virtual double perimeter() const { return 4*side; }
};

int main()
{
    typedef shared_ptr<Shape> shape_ptr;
   // create a vector of 2 pointers to shape
    vector<shape_ptr> s = { shape_ptr(new Square(3)), shape_ptr(new Circle(2)) };

   // we use Shape* to access the implementation of perimeter
   // in each specific class
   // Note the use of reference to avoid creating copy of shape_ptr objects
   for (const auto& sh : s)
      cout << sh->perimeter() << endl;

   return 0;
}
