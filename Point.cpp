#include <iostream>
#include <cmath>
#include <array>

using namespace std;

struct Point
{
	// contructor
	Point() : m_x(0), m_y(0) {}  // set in the origin
	Point( double x, double y) : m_x(x), m_y(y) {}

	// disatance
    double distance(const Point& p) const
	{
		double dx = m_x - p.m_x;
		double dy = m_y - p.m_y;
		return sqrt(dx*dx+dy*dy);
	}

	// accessor methods
	double x() const {return m_x;}
	double y() const {return m_y;}

	// vertical and horizontal distance
	double xdistance(const Point& p) const { return abs(m_x - p.m_x); }
	double ydistance(const Point& p) const { return abs(m_y - p.m_y); }

private:
	double m_x, m_y;
};

struct Quadrilateral
{
	// constructors
	Quadrilateral() {}  // degenerate quadrilateral
	Quadrilateral(const Point& p1, const Point& p2, const Point& p3, const Point& p4) : m_vertex{ p1, p2, p3, p4 } {}

	// length of a side between points i and j
	double distance(size_t i, size_t j) const { return m_vertex[i].distance(m_vertex[j]); }

	// perimeter
	double perimeter() const
	{
		double d = distance(0, 3);
		for (size_t i = 1; i < m_vertex.size(); ++i)
			d += distance(i-1, i);
		return d;
	}

	// accessor
	const Point& p(size_t i) const { return m_vertex[i]; }

private:
	array<Point, 4> m_vertex;
};

struct Rectangle : Quadrilateral
{
	Rectangle() {}
    Rectangle(const Point& topLeft, const Point& btmRight)
		: Quadrilateral{ topLeft, {btmRight.x(), topLeft.y()}, btmRight, {topLeft.x(), btmRight.y()} }
		{}
    double base() const { return distance(0,1); }
    double height() const { return distance(1,2); }
	double area() const { return base() * height(); };
};

struct Square : Rectangle
{
	Square() {}
    Square(const Point& topLeft, double sidelength)
		: Rectangle{ topLeft, {topLeft.x() + sidelength, topLeft.y() - sidelength } }
		{}
    double side() const { return base(); }
};


struct Trapeze : Quadrilateral
{
	Trapeze() {}
    Trapeze(const Point& topLeft, const Point& btmLeft, double topBase, double btmBase)
		: Quadrilateral{ topLeft, {topLeft.x()+topBase, topLeft.y()}, {btmLeft.x() + btmBase, btmLeft.y()}, btmLeft }
		{}
	double upperBase() const { return distance(0,1); }
	double lowerBase() const { return distance(2,3); }
	double height() const { return distance(2,3); }
	double area() const { return (upperBase() + (3,2)) / 2.0 * p(1).ydistance(p(2)); }
};

int main()
{
	Quadrilateral q({0,0}, {0,1}, {1,1}, {1,0});
	cout << "perimeter quadrilateral " << q.perimeter() << "\n";

	Rectangle r{{0,1}, {4,0}};
	cout << "perimeter rectangle " << r.perimeter() << endl;
	cout << "area rectangle " << r.area() << endl;

	Square sq{{0,2}, 2};
	cout << "perimeter square " << sq.perimeter() << endl;
	cout << "area square " << sq.area() << endl;

	Trapeze tr{{1,1}, {0,0}, 1, 3};
	cout << "perimeter trapeze " << tr.perimeter() << endl;
	cout << "area trapeze " << tr.area() << endl;

	return 0;
}
