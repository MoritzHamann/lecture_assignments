#include <math.h>
#include <iostream>
#include <string>
#include <sstream>

struct Point{
	double x;
	double y;
	double value;

	Point(){}

	Point(const double _x, const double _y, const double _value){
		x = _x;
		y = _y;
		value = _value;
	}

	Point(const Point &p){
		x=p.x;
		y=p.y;
		value = p.value;
	}
};

struct Quad{
	Point x1;
	Point x2;
	Point x3;
	Point x4;

	Quad(const Point _x1, const Point _x2, const Point _x3, const Point _x4 ){
		x1 = _x1;
		x2 = _x2;
		x3 = _x3;
		x4 = _x4;
	}

	Quad(const Quad &q){
		x1 = q.x1;
		x2 = q.x2;
		x3 = q.x3;
		x4 = q.x4;
	}
};

struct LocaleCoordinats{
	double a1;
	double a2;

	LocaleCoordinats(const double _a1, const double _a2){
		a1 = _a1;
		a2 = _a2;
	}
};

Point interpolateCoords(LocaleCoordinats lc, Quad q) {
    Point x12 = Point(0,0,0);
    x12.x = lc.a1*q.x1.x + (1-lc.a1)*q.x2.x;
    x12.y = lc.a1*q.x1.y + (1-lc.a1)*q.x2.y;
    
    Point x34 = Point(0,0,0);
    x34.x = lc.a1*q.x4.x + (1-lc.a1)*q.x3.x;
    x34.y = lc.a1*q.x4.y + (1-lc.a1)*q.x3.y;
    
    Point x = Point(0,0,0);
    x.x = lc.a2*x12.x + (1-lc.a2)*x34.x;
    x.y = lc.a2*x12.y + (1-lc.a2)*x34.y;
    
    return x;
}

double interpolateValue(LocaleCoordinats lc, Quad q){
    double value = 0.0;
    
    value = lc.a2*(lc.a1*q.x1.value + (1-lc.a1)*q.x2.value) + (1-lc.a2)*(lc.a1*q.x4.value + (1-lc.a1)*q.x3.value);
    
    return value;
}

LocaleCoordinats newtonIteration(Point target, Quad q){
        
    //define error
    double const err = 0.2;

	//define detla alphas
    //to improve accuracy, compute these with the help of the Jacobi matrix.
	double const delta_a1 = 0.01;
    double const delta_a2 = 0.03;


    //define seed Points
    LocaleCoordinats lc = LocaleCoordinats(0.0,0.0);
    
    //start/initialize iteration
    Point phi = interpolateCoords(lc, q);
    double absErrorVector = sqrt(pow(target.x-phi.x,2) + pow(target.y-phi.y,2));
    
	while( absErrorVector > err)
    {
        lc.a1 += delta_a1;
        lc.a2 += delta_a2;
		
		//calculate point according to local Coordinates
        phi = interpolateCoords(lc, q);

		//caculate error of calculated point in respects to target point.
        absErrorVector = sqrt(pow(target.x-phi.x,2) + pow(target.y-phi.y,2));
    }
	
	return lc;	
}

int main(){
	
	Point target = Point(1,1,0);
	Quad q = Quad(Point(2,4,1),Point(0,0,2),Point(5,0,3),Point(4,3,4));
	
	//compute local coordinates.
	std::cout << "Computing Local Coordinates...\n";
	LocaleCoordinats lc = newtonIteration(target, q);

	//interpolate value with computed local coordinates
	std::cout << "Interpolating Value...\n";
	target.value = interpolateValue(lc, q);
    
    std::cout << "Target:\n  x: " << target.x << "\n  y: " << target.y << "\n  value: " << target.value <<"\n\n"; 
	
	system("Pause");

	return 0;
}



