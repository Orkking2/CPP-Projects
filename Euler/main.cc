#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>
#include <optional>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#define PNG_DEBUG 3
#include <png.h>



double slope(double x, double y){
    return (pow(x,2) + pow(y,2) - 2*x*y - 3);  // x^2 + y^2 - 2xy - 3
}

struct Point{
    Point(): x(0),y(0){}
    Point(double XComponent, double YComponent){
        x = XComponent;
        y = YComponent;
    }
    double x, y;
    void Scale(double scaler){
    	double r = sqrt(pow(this->x, 2) + pow(this->y, 2));
		double theta = atan(this->x*this->y);
		r *= scaler;
		this->x = r*cos(theta);
		this->y = r*sin(theta);
    }

    friend std::ostream& operator << (std::ostream& os, const Point& p);
    friend bool operator < (const Point& l, const Point& r);
    friend bool operator > (const Point& l, const Point& r);
    friend bool operator <= (const Point& l, const Point& r);
    friend bool operator >= (const Point& l, const Point& r);

    Point operator + (const Point& p){
        return Point(this->x + p.x, this->y + p.y);
    }

    bool operator = (const Point& p){
        return (this->x == p.x && this->y == p.y);
    }

    bool operator != (const Point& p){
        return !(this->x == p.x || this->y == p.y);
    }

    void operator += (const Point& p){
        this->x += p.x;
        this->y += p.y;
    }

    bool inRadius(Point p, double r){
        return sqrt(pow(p.x - this->x, 2) + pow(p.y - this->y, 2)) <= r;
    }

};

bool operator < (const Point& l, const Point& r){
    return l.x < r.x && l.y < r.y;
}
bool operator > (const Point& l, const Point& r){
    return l.x > r.x && l.y > r.y;
}
bool operator <= (const Point& l, const Point& r){
    return l.x <= r.x && l.y <= r.y;
}
bool operator >= (const Point& l, const Point& r){
    return l.x >= r.x && l.y >= r.y;
}

std::ostream& operator << (std::ostream& os, const Point& p){
    os << "(" << p.x << ", " << p.y << ")";
    return os;
}

struct Segment{
    Segment() = default;
    Segment(Point point1, Point point2){
        p1 = point1;
        p2 = point2;
    }
    Point p1, p2;
    void Scale(double scaler){
    	this->p1.Scale(scaler);
    	this->p2.Scale(scaler);
    }
    void Reverse(){
    	Segment s(this->p2, this->p1);
    	this->p1 = s.p1;
    	this->p2 = s.p2;
    }
    void Order(){
    	if(this->p2 > this->p1){
    		this->Reverse();
    	}
    }
    void Translate(Point vector){
    	this->p1 += vector;
    	this->p2 += vector;
    }
    friend bool operator < (const Segment& l, const Segment& r);
    Point operator ++ (){
    	return this->p1 + this->p2;
    }
    void operator += (const Segment& s){
    	this->p1 += s.p1;
    	this->p2 += s.p1;
    }
};

bool operator < (const Segment& l, const Segment& r){
    return l.p1 < r.p1 && l.p2 < r.p2;
}

/*
template <typename T>
class container{
    T x;
    T y;
};
*/

class Matrix{
    public:
        Matrix() = default;
        Matrix(Point pMin, Point pMax, std::vector<int> dimensions){
            xDist = pMax.x - pMin.y;
            yDist = pMax.y - pMin.y;

            startMin = pMin;
            startMax = pMax;

            d = {dimensions[0], dimensions[1], dimensions[0]*dimensions[1]};

            data.resize(d[2]);

            double xSpacing = xDist/(d[0]-1);
            double ySpacing = yDist/(d[1]-1);

            for(int i = 0; i < d[0]; i++){
                xList.push_back(i*xSpacing + pMin.x);
            }
            for(int i = 0; i < d[1]; i++){
                yList.push_back(i*ySpacing + pMin.y);
            }

            for(int i = 0; i < d[2]; i++){
                data.emplace_back(Point(xList[i % d[0]],yList[floor(i/d[0])]));
            }
        }
        const std::vector<Point>& getMatrix(){
            return data;
        }
    protected:
        std::vector<Point> data;
        double yDist = 0, xDist = 0;
        std::vector<int> d;
        Point startMin, startMax;
        std::vector<double> xList, yList;
};

class SlopeField : public Matrix{
    public:
        SlopeField(Point pMin, Point pMax, std::vector<int> dimensions) : Matrix(pMin, pMax, dimensions){
            genSegments();
        }
        const std::vector<Segment>& getSegments() const {
            return segments;
        }
        const std::vector<Point>& getBounds () const {
            return bounds;
        }
    private:
        void genSegments(){ // *Also gens bounds
            double r = std::min(xDist/(2*d[0]-1),yDist/(2*d[1]-1));
            std::vector<double> k;
            std::optional<Point> pLow, pHigh;

            for(int i = 0; i < d[2]; i++){
            	/* k explanation:
             	 	*<--->.<--->* (diagram line)
                    |  k  |  k  | -- k = dist
                    * = p1|     * = p2
                          . = data[i]           */
            	k.push_back(r/pow(pow(slope(data[i].x,data[i].y), 2) + 1, 1/2));

                // Defining points in segment
                Point pMin (data[i].x - k[i], data[i].y - slope(data[i].x,data[i].y)*k[i]);
                Point pMax (data[i].x + k[i], data[i].y + slope(data[i].x,data[i].y)*k[i]);

                // Evaluating bounds
				if(!pLow.has_value() && !pHigh.has_value()){
					pLow.value().x = pMin.x;
					pHigh.value().x = pMax.x;
					if(pMin.y < pMax.x){
						pLow.value().y = pMin.y;
						pHigh.value().y = pMax.y;
					} else {
						pLow.value().y = pMax.y;
						pHigh.value().y = pMin.y;
					}
				} else {
					for(Point pCurr : {pMin, pMax}){

						// pLow
						if(pLow.value().x > pCurr.x){
							pLow.value().x = pCurr.x;
						}
						if(pLow.value().y > pCurr.y){
							pLow.value().y = pCurr.y;
						}

						// pHigh
						if(pHigh.value().x < pCurr.x){
							pHigh.value().x = pCurr.x;
						}
						if(pHigh.value().y < pCurr.y){
							pHigh.value().y = pCurr.y;
						}
					}
				}

                // Constructing Segment list
                segments.push_back(Segment(pMin, pMax));

            }

            // Constructing bounds list
            bounds.push_back(pLow.value());
            bounds.push_back(pHigh.value());
        }

        // Members
        std::vector<Point> data;
        std::vector<Segment> segments;
        std::vector<Point> bounds;
};

std::vector<double> StartSelect(){
    double startX, startY, targetX;
    // questions
    std::cout << "Start X: ";
    std::cin >> startX;
    std::cout << "Start Y: ";
    std::cin >> startY;
    std::cout << "Target X: ";
    std::cin >> targetX;
    // Out construction
    std::vector<double> out = {startX, startY, targetX};
    // Out
    return out;
}

std::vector<double> NHSelect(double targetX, double x){
    int n;
    double h, nCap = pow(10, 7);
    char nhSelect;
    std::vector<double> out;

    std::cout << "N or H: ";
    std::cin >> nhSelect;
    if(tolower(nhSelect) == 'n'){
        std::cout << "N: ";
        std::cin >> n;
        if(n > nCap){
            std::cout << "N > " << nCap << " -- N reset to " << nCap << std::endl;
            n = nCap;
        } else if(n <= 0){
            std::cout << "N <= 0, reset to " << nCap << std::endl;
            n = nCap;
        }
        h = (targetX-x)/static_cast<double>(n);
        std::cout << "H = " << h << std::endl;
    } else if(tolower(nhSelect) == 'h'){
        std::cout << "H: ";
        std::cin >> h;
        n = (targetX-x)/h;
        if(n > nCap){
            std::cout << "N > " << nCap << " -- N reset to " << nCap << std::endl;
            n = nCap;
            h = (targetX-x)/static_cast<double>(n);
            std::cout << "H = " << h << std::endl;
        } else if(n <= 0){
            std::cout << "N <= 0, reset to " << nCap << std::endl;
            n = nCap;
            h = (targetX-x)/static_cast<double>(n);
            std::cout << "H = " << h << std::endl;
        }
    } else if(nhSelect == 'm'){
        n = pow(10, 8);
        h = (targetX-x)/static_cast<double>(n);
        std::cout << std::endl << "Max N selected " << std::endl << "N = " << n << std::endl << "H = " << h << std::endl;
    } else {
        std::cout << "Invalid char" << std::endl;
        out = NHSelect(targetX, x);
    }

    std::cout << std::endl;

    out.push_back(n);
    out.push_back(h);

    return out;
}

std::vector<Point> Approximate(std::vector<double> in, bool isTesting){
	int n;
	double h, tRadius;
	Point startPoint;

	if(isTesting){
		tRadius = 1/5000;
		n = pow(10, 7);
		h = 1/n;
	} else {
		startPoint = Point (in[0],in[1]);

		std::vector<double> k;
		double targetX = in[2], xDist = abs(targetX - startPoint.x);
		tRadius = xDist/5000;

		std::vector<double> nh = NHSelect(targetX, startPoint.x);

		n = static_cast<int> (nh[0]);
		h = nh[1];
	}

	Point iterativePoint = startPoint;
	std::vector<Point> pointList = {startPoint};

	int j = 0;

    for(int i = 0; i < n; i++){
        iterativePoint += Point(h, slope(iterativePoint.x,iterativePoint.y)*h);
        if(!pointList[j].inRadius(iterativePoint, tRadius)){
            pointList.push_back(iterativePoint);
            j++;
        }
    }

    if(pointList[pointList.size()-1] != iterativePoint){
        pointList.push_back(iterativePoint);
    }

    return pointList;
}

std::vector<std::vector<double>> GetXYList (std::vector<Point> pList){
	std::vector<double> xList, yList;
	for(Point p : pList){
		xList.push_back(p.x);
		yList.push_back(p.y);
	}
	return {xList, yList};
}

void PlotLine(std::vector<Point> pList){
	auto L = GetXYList(pList);
	// plt::plot(L[0], L[1], "g-");
}


/*
void PlotSegments(std::vector<Segment> segments){
	for(Segment s : segments){
		// plt::plot({s.p1.x, s.p2.x}, {s.p1.y, s.p2.y}, "b-");
	}
}
*/

struct Pixel{
	// Pixel() = y(0), x(0), RGB({0,0,0});
	Pixel(double xIn, double yIn, std::vector<int> RGBIn){
		x = xIn;
		y = yIn;
		for(int i : RGBIn){
			RGB.push_back(i);
		}
	}
	int x, y;
	std::vector<int> RGB;
};

class PersonalImage{
	public:
		PersonalImage(std::vector<Pixel> baseImage, double scalerIn){
			for(Pixel p : baseImage){
				pixels.push_back(p);
			}
			scaler = scalerIn;
		}
		PersonalImage(std::vector<Point> xyBounds, std::vector<int> sizeIn){
			bounds = xyBounds;
			size = sizeIn;
			double areaPoints = (abs(bounds[0].x)+abs(bounds[1].x))*(abs(bounds[0].y + abs(bounds[1].y)));
			double areaReal = size[0]*size[1];
			scaler = areaReal/areaPoints;
			for(Point p : Matrix(Point(), Point(sizeIn[0], sizeIn[1]), sizeIn).getMatrix()){
				pixels.push_back(Pixel(p.x,p.y,{0,0,0}));
			}
			Pixel centerPixel = pixels[(sizeof(pixels)-1)/2];
			TranslationVector (centerPixel.x, centerPixel.y); // New origin
		}
		std::vector<Pixel> GetImage(){
			return pixels;
		}
		std::vector<Point> GetBounds(){
			return bounds;
		}
		std::vector<int> GetSize(){
			return size;
		}
		void PlotPoint(Point p){
			p.Scale(scaler);
			p += TranslationVector;
			IncludePixel(Pixel (p.x, p.y, {255,255,255}));
		}
		void PlotSegment(Segment s){
			s.Scale(scaler);
			s.Translate(TranslationVector);
			s.Order();
			double slope = (s.p1.y - s.p2.y)/(s.p1.x - s.p2.x);
			for(Pixel &p : pixels){ // Find a better way to do this
				if(p >= s.p1 && p <= s.p2 && p.y == slope*p.x){
					p.RGB = {255,255,255};
				}
			}

		}
	private:
		void IncludePixel(Pixel p){
			pixels[abs((int)p.x) + abs((int)p.y)*size[0]].RGB = p.RGB;
		}
		std::vector<Point> bounds;
		std::vector<int> size;
		double scaler;
		std::vector<Pixel> pixels;
		Point TranslationVector;
};


int main(int /*argc*/, char ** argv)
{



    return 0;
}
