#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include <math.h>
#include <algorithm>
#include <optional>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>

#define PNG_DEBUG 3
#include <png.h>



double slope(double x, double y){
    return (pow(x, 2) + pow(y, 2) - 2*x*y - 3);  // x^2 + y^2 - 2xy - 3
}

struct Point{
    Point(): x(0),y(0){}
    Point(double XComponent, double YComponent){
        x = XComponent;
        y = YComponent;
    }
    double x, y;
    void Scale(double scaler){
    	double r = sqrt(pow(x, 2) + pow(y, 2));
		double theta = atan(x*this->y);
		r *= scaler;
		x = r*cos(theta);
		y = r*sin(theta);
    }
    void ToInt(){
    	x = (int)x;
    	y = (int)y;
    }

    friend std::ostream& operator << (std::ostream& os, const Point& p);
    friend bool operator < (const Point& l, const Point& r);
    friend bool operator > (const Point& l, const Point& r);
    friend bool operator <= (const Point& l, const Point& r);
    friend bool operator >= (const Point& l, const Point& r);

    Point operator + (const Point& p){
        return Point(x + p.x, y + p.y);
    }

    bool operator == (const Point& p){
        return (x == p.x && y == p.y);
    }

    bool operator != (const Point& p){
        return !(x == p.x || y == p.y);
    }

    void operator += (const Point& p){
        x += p.x;
        y += p.y;
    }

    bool inRadius(Point p, double r){
        return sqrt(pow(p.x - x, 2) + pow(p.y - y, 2)) <= r;
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
    Segment(Point point1 = Point(), Point point2 = Point()){
        p1 = point1;
        p2 = point2;
    }
    Point p1, p2;
    void Scale(double scaler){
    	p1.Scale(scaler);
    	p2.Scale(scaler);
    }
    void Reverse(){
    	Segment s(p2, this->p1);
    	p1 = s.p1;
    	p2 = s.p2;
    }
    void Order(){
    	if(p2 > this->p1){
    		Reverse();
    	}
    }
    void Translate(Point vector){
    	p1 += vector;
    	p2 += vector;
    }
    void ToInt(){
    	p1.ToInt();
    	p2.ToInt();
    }
    friend bool operator < (const Segment& l, const Segment& r);
    Point operator ++ (){
    	return p1 + this->p2;
    }
    void operator += (const Segment& s){
    	p1 += s.p1;
    	p2 += s.p1;
    }
};

bool operator < (const Segment& l, const Segment& r){
    return l.p1 < r.p1 && l.p2 < r.p2;
}

struct Pixel{
	Pixel(double xIn = 0, double yIn = 0, std::vector<int> RGBIn = {0,0,0}){
		x = xIn;
		y = yIn;
		RGB = RGBIn;
	}
	void ToInt(){
		x = (int)x;
		this->y = (int)y;
	}
	int x, y;
	std::vector<int> RGB;
};

double InvSigmoid(double x, double a){
	return -log(a/(x + 0.5*a) - 1);
}

struct Particle{
	Particle(Point origin = Point(), double thetaIn = 0){
		p = origin;
		theta = thetaIn;
		vector = Point (r*cos(theta), r*sin(theta));
	}
	double deltaThetaGen(int seedIn = seed){
		return 1/3*InvSigmoid(((rand() + seedIn) % 1000 - 500)/50, 10); // https://www.desmos.com/calculator/hmxc9uegiz
	}
	void UpdateVector(int seedIn = seed){
		theta += deltaThetaGen(seedIn);
		vector = Point(r*cos(theta), r*sin(theta));
	}
	void UpdatePoint(Point vector){
		p += vector;
	}
	void UpdatePoint(){
		p += vector;
	}
	void Update(){
		UpdatePoint();
		UpdateVector();
	}
	Point p, vector;
	double theta;
	static int seed = rand(), r = 0.5;
};

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
        std::vector<int> d;
        double yDist = 0, xDist = 0;
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

std::vector<std::vector<double>> GetXYList (std::vector<Point> pList){
	std::vector<double> xList, yList;
	for(Point p : pList){
		xList.push_back(p.x);
		yList.push_back(p.y);
	}
	return {xList, yList};
}

/*
void PlotLine(std::vector<Point> pList){
	auto L = GetXYList(pList);
	// plt::plot(L[0], L[1], "g-");
}

void PlotSegments(std::vector<Segment> segments){
	for(Segment s : segments){
		// plt::plot({s.p1.x, s.p2.x}, {s.p1.y, s.p2.y}, "b-");
	}
}
*/

// std::vector<int> RGBInit = {0,0,0};

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
/*
template <typename T>
class container{
    T x;
    T y;
};
*/

/*
template <typename T>
class PlusIndex{
	T i;
	int index;
};
*/

struct Token{
    enum struct Type{
        Unknown,
        Number,
		Operator,
        LeftParen,
        RightParen,
		X,
    };

    Token(Type t, const std::string& s, int prec = -1, bool ra = false)
        : type {t}, str (s), precedence {prec}, rightAssociative {ra}
    {}

    const Type type;
    const std::string str;
    const int precedence;
    const bool rightAssociative;
};

std::deque<Token> ParseExpression(std::string& string){
	std::deque<Token> tokens;

	for(const auto* c = string.c_str(); *c; c++){
		if(isblank(*c)){
			// pass
		} else if(isdigit(*c)){
			const auto* b = c;
			auto* e = c;

			while(isdigit(*e)){
				e++;
			}
			const auto s = std::string(b, e);
			tokens.push_back(Token {Token::Type::Number, s});
			c += s.size();
		} else {
			Token::Type t = Token::Type::Unknown;
			int pr = -1;
			bool ra = false;
			switch(*c){
			default:                                    break;
			case '(':   t = Token::Type::LeftParen;     break;
			case ')':   t = Token::Type::RightParen;    break;
			case '_':   t = Token::Type::Operator;      pr = 5; ra = true; break;
			case '^':   t = Token::Type::Operator;      pr = 4; ra = true; break;
			case '*':   t = Token::Type::Operator;      pr = 3; break;
			case '/':   t = Token::Type::Operator;      pr = 3; break;
			case '%':   t = Token::Type::Operator;      pr = 3; break;
			case '+':   t = Token::Type::Operator;      pr = 2; break;
			case '-':   t = Token::Type::Operator;      pr = 2; break;
			case 'x':   t = Token::Type::X;             break;
			}
			const auto s = std::string(1, *c);
			tokens.push_back(Token {t, s, pr, ra});
		}
	}

	// ShuntingYard

	std::deque<Token> queue;
	std::vector<Token> stack;

	for(Token token : tokens){
		switch(token.type){
		case Token::Type::Number:
			queue.push_back(token);
			break;

		case Token::Type::Operator:
			{
				const auto o1 = token;
				while(!stack.empty()){
					const auto o2 = stack.back();
					if((!o1.rightAssociative && o1.precedence <= o2.precedence) || (o1.rightAssociative && o1.precedence <  o2.precedence)){
						stack.pop_back();
						queue.push_back(o2);

						continue;
					}
					break;
				}
				stack.push_back(o1);
			}
			break;

		case Token::Type::LeftParen:
			stack.push_back(token);
			break;

		case Token::Type::RightParen:
			{
				bool match = false;
				while(! stack.empty() && stack.back().type != Token::Type::LeftParen){
					queue.push_back(stack.back());
					stack.pop_back();
					match = true;
				}
				stack.pop_back();

				if(!match && stack.empty()){
					printf("RightParen error (%s)\n", token.str.c_str());
					return {};
				}
			}
			break;

		default:
			printf("error (%s)\n", token.str.c_str());
			return {};
		}
	}
	while(!stack.empty()){
		if(stack.back().type == Token::Type::LeftParen){
			printf("Mismatched parentheses error\n");
			return {};
		}

		queue.push_back(std::move(stack.back()));
		stack.pop_back();
	}
	return queue;
}

double EvaluateExpression(std::deque<Token> queue, double x){
	std::vector<int> stack;
	while(!queue.empty()){
		const auto token = queue.front();
		queue.pop_front();
		switch(token.type){
		case Token::Type::Number:
			stack.push_back(std::stoi(token.str));
			break;

		case Token::Type::Operator:
			{
				if(token.str[0] == 'x'){
					stack.push_back(x);
					break;
				}
				const auto r = stack.back();
				stack.pop_back();
				if(token.str[0] == '_'){
					stack.push_back(-r);
					break;
				}
				const auto l = stack.back();
				stack.pop_back();

				switch(token.str[0]) {
				default:
					printf("Operator error [%s]\n", token.str.c_str());
					exit(0);
					break;
				case '^':
					stack.push_back((int)(pow(l, r)));
					break;
				case '*':
					stack.push_back(l * r);
					break;
				case '/':
					stack.push_back(l / r);
					break;
				case '+':
					stack.push_back(l + r);
					break;
				case '-':
					stack.push_back(l - r);
					break;
				case '%':
					stack.push_back(l % r);
				}
			}
			break;

		default:
			printf("Token error\n");
			exit(0);
		}
	}
	return stack.back();
}

class PersonalImage{
	public:
		PersonalImage(std::vector<Pixel> baseImage){
			pixels.resize(sizeof(baseImage));
			int i = 0;
			for(Pixel p : baseImage){
				pixels[i] = p;
				i++;
			}
		}
		PersonalImage(std::vector<Point> xyBounds, std::vector<int> sizeIn){
			bounds = xyBounds;
			size = sizeIn;
			double areaReal = size[0]*size[1];
			pixels.resize(areaReal);
			int i = 0;
			for(Point p : Matrix(Point(0,0), Point(sizeIn[0], sizeIn[1]), sizeIn).getMatrix()){
				pixels[i] = Pixel(p.x,p.y,{0,0,0});
				i++;
			}
			TransSlope[0] = -size[0]/(bounds[0].x - bounds[1].x);
			TransSlope[1] = -size[1]/(bounds[0].y - bounds[1].y);
			TransVector = Point (-TransSlope[0]*bounds[0].x, -TransSlope[1]*bounds[0].y); // (Min x, Min y) -> (0,0)
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
			p = TranslatePoint(p);
			IncludePixel(Pixel (p.x, p.y, {255,255,255}));
		}
		void PlotSegment(Segment s){
			s = TranslateSegment(s);
			s.ToInt();
			if(s.p2.x > s.p1.x){
				s.Reverse();
			}
			double slope = (s.p1.y - s.p2.y)/(s.p1.x - s.p2.x);
			int xDist = s.p1.x - s.p2.x;
			for(int i = 0; i <= xDist; i++){
				IncludePixel(Pixel(i + s.p1.x, slope*i + s.p1.y, {255,255,255}));
			}
		}
		void PlotParticle(Particle p){
			IncludePixel(Pixel(p.p.x, p.p.y, {255,255,255}));
			while(p.p.x > 0 && p.p.y > 0 && p.p.x <= size[0] && p.p.y <= size[1]){
				p.Update();
				IncludePixel(Pixel(p.p.x, p.p.y, {255,255,255}));
			}
		}
		void PlotEXFunction(std::string f, int bounds[2]){
			std::deque<Token> queue = ParseExpression(f);
			int xDist = (int)(bounds[1]*TransSlope[0]) - (int)(bounds[0]*TransSlope[0]);
			for(int i = 0; i <= xDist; i++){
				int x = i/TransSlope[0] + bounds[0];
				int y = (int)EvaluateExpression(queue, x);
				Point p = TranslatePoint(Point(x, y));
				IncludePixel(Pixel(p.x, p.y, {255,255,255}));
			}
		}
	private:
		void IncludePixel(Pixel p){
			p.ToInt();
			pixels[abs(p.x) + abs(p.y)*size[0]].RGB = p.RGB;
		}
		Point TranslatePoint (Point p){
			p.x *= TransSlope[0];
			p.y *= TransSlope[1];
			p += TransVector;
			return p;
		}
		Segment TranslateSegment (Segment s){
			s.p1 = TranslatePoint(s.p1);
			s.p2 = TranslatePoint(s.p2);
			return s;
		}
		std::vector<Point> bounds;
		std::vector<int> size;
		std::vector<Pixel> pixels;
		Point TransVector;
		double TransSlope[2];
};


int main(int /*argc*/, char ** argv)
{


	return 0;
}
