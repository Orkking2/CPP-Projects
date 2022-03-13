#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>
#include <optional>

// #include <Python.h>

struct Point{
    Point(): x(0),y(0){}
    Point(double XComponent, double YComponent){
        x = XComponent;
        y = YComponent;
    }
    double x;
    double y;

    friend std::ostream& operator << (std::ostream& os, const Point& p);
    friend bool operator < (const Point& l, const Point& r);

    Point operator + (const Point& p){
        return Point(this->x + p.x, this->y + p.y);
    }

    bool operator = (const Point& p){
        return (this->x == p.x && this->y == p.y);
    }

    bool operator != (const Point& p){
        return !(this->x == p.x && this->y == p.y);
    }

    void operator += (const Point& p){
        this->x += p.x;
        this->y += p.y;
    }

    bool inRadius(Point p, double r);
};


