#include <iostream>
#include <algorithm>
#include <vector>
using namespace std;
#include <iostream>
#include <memory>
#include <vector>

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
    typedef unique_ptr<Shape> shape_ptr;
   // create a vector of 2 pointers to shape
   vector<shape_ptr> s;
   s.push_back(shape_ptr(new Square(3)));
   s.push_back(shape_ptr(new Circle(2)));

   // we use Shape* to access the implementation of perimeter
   // in each specific class
   for (int i = 0; i < 2; ++i)
      cout << s[i]->perimeter() << endl;

   return 0;
}