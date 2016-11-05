#include <iostream>
#include <algorithm>
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
   // create a vector of 2 pointers to shape and initialize it
   // with a circle and a square 
   Shape *s[2] = { new Square(3), new Circle(2) };

   // we use Shape* to access the implementation of perimeter
   // in each specific class
   for (int i = 0; i < 2; ++i)
      cout << s[i]->perimeter() << endl;

   // remember we need to release the objects created with new
   for (int i = 0; i < 2; ++i)
      delete s[i];

   return 0;
}
