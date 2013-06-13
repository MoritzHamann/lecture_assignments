#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
#include <list>

#include "vtkMaskPoints.h"
#include "vtkPointData.h"
#include "vtkStructuredGrid.h"
#include "vtkIdList.h"

#include "vtkGlyphSource2D.h"
#include "vtkGlyph3D.h"

#include "vtkArrayCalculator.h"
#include "vtkContourFilter.h"

#include "vtkDataSetSurfaceFilter.h"

#include "vtkGenericDataObjectReader.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"

#include "vtkSmartPointer.h"
#include "vtkCellLocator.h"
#include "vtkInformationVector.h"
#include "vtkCell.h"
#include "vtkQuad.h"

struct Point{
	double x;
	double y;
	double value;

	Point(){}

	Point(const double _x, const double _y, const double _value){
		x = _x;
		y = _y;
		value =  _value;
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
	Quad(){}
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

Point operator-(const Point &a, const Point &b){
	return Point(a.x - b.x, a.y - b.y, 0);
}

Point operator+(const Point &a, const Point &b){
	return Point(a.x + b.x, a.y + b.y, 0);
}

Point operator*(const double &a, const Point &b){
	return Point(b.x*a, b.y*a, b.value);
}

Point operator*(const Point &b, const double &a){
	return Point(b.x*a, b.y*a, b.value);
}

LocaleCoordinats operator-(const LocaleCoordinats &a, const LocaleCoordinats &b){
	return LocaleCoordinats(a.a1 - b.a1, a.a2 - b.a2);
}

LocaleCoordinats operator+(const LocaleCoordinats &a, const LocaleCoordinats &b){
	return LocaleCoordinats(a.a1 + b.a1, a.a2 + b.a2);
}


using namespace std;



//--------------------------------------Stencil-Walk-------------------------------------------------

int newtonIterationStencilWalk(Point target, Quad &q, LocaleCoordinats &lc);
int checkNewPointPosition(Quad q, Point x_old, Point x_new);
LocaleCoordinats edgeIntsection(Point x_old, Point x_new, Point x1, Point x2);
int getQuadContainingPoint(Point p, Quad & q);
int getNeighborCell(int edge_number, Point p, Quad &q);
void getQuadGridValues(Quad &q, Point p);
//--------------------------------------------------------------------------------------------------

vtkSmartPointer<vtkGenericDataObjectReader> vtk_reader;

int main(){
//----------------Read-File---------------------
	//create object for reading VTK files
	vtk_reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	//set file name
	vtk_reader->SetFileName("visit_ex_db.vtk");
	//update
	vtk_reader->Update();

	//vtk_reader->GetStructuredGridOutput()->GetPointData()->Print(std::cout);
	//set active scalar field
	vtk_reader->GetStructuredGridOutput()->GetPointData()->SetActiveScalars("v");
//----------------------------------------------



	Point target = Point(3,3,0);
	Point start  = Point(1,2,0);

	Quad q;
	LocaleCoordinats lc(0,0);


	//--------------Assignment_07--------------------
    
    
    // get inital quad
    getQuadContainingPoint(start, q);
    
    // set lc for first quad at (0.5,0.5) (we could compute the localcoordinates for start, but this is to complicated...)
    lc.a1 = 0.5;
    lc.a2 = 0.5;
    
   
    newtonIterationStencilWalk(target, q, lc);
    
    getQuadGridValues(q, target);
    // bilinearinterpolation at lc
    float result = 0;
    result = (1-lc.a1)*(1-lc.a2)*q.x2.value + (1-lc.a2)*lc.a1*q.x3.value + (1-lc.a1)*lc.a2*q.x1.value + lc.a1*lc.a2*q.x4.value;
    stc::cout << "Result: " << result << std::endl;
    
	//-----------------------------------------------
}

//--------------------------------------Stencil-Walk-------------------------------------------------

/*
Function: newtonIterationStencilWalk

Input : Point target  
Output: Quad q = Quad in which target is located
        LocaleCoordinats lc

Return value is of type int, to give back a state :  1 quad that contains target was found and the LocaleCoordinats lc are valid
													-1 quad that contains target was not found and the LocaleCoordinats lc are not valid

*/
int newtonIterationStencilWalk(Point target, Quad &q, LocaleCoordinats &lc){
    
    // jacobian and inverse jacobian
    float j11 = 0;
    float j12 = 0;
    float j21 = 0;
    float j22 = 0;
    float invj11 = 0;
    float invj12 = 0;
    float invj21 = 0;
    float invj22 = 0;
    
    Point temp = Point(0,0,0);
    Point difference = Point(0,0,0);
    bool notfound =  true;
    
    while (notfound){
        // calculate temp
        temp.x = lc.a1*lc.a2*q.x1.x + lc.a2*q.x2.x - lc.a1*lc.a2*q.x2.x + lc.a1*q.x4.x - lc.a1*lc.a2*q.x4.x + (1-lc.a1-lc.a2+lc.a1*lc.a2) * q.x3.x;
        temp.y = lc.a1*lc.a2*q.x1.y + lc.a2*q.x2.y - lc.a1*lc.a2*q.x2.y + lc.a1*q.x4.y - lc.a1*lc.a2*q.x4.y + (1-lc.a1-lc.a2+lc.a1*lc.a2) * q.x3.y;
        
        // calculate difference
        difference.x = target.x - temp.x;
        difference.y = target.y - temp.y;
        
        if (sqrt(pow(difference.x, 2) + pow(difference.y, 2)) <= 0.1){
            notfound = false;
        }
        else{
            // calcualte jacobian
            j11 = lc.a2*q.x1.x - lc.a2*q.x2.x + q.x4.x - lc.a2*q.x4.x - q.x3.x + lc.a2*q.x3.x;
            j12 = lc.a2*q.x1.y - lc.a2*q.x2.y + q.x4.y - lc.a2*q.x4.y - q.x3.y + lc.a2*q.x3.y;
            j21 = lc.a1*q.x1.x + q.x2.x - lc.a1*q.x2.x - lc.a1*q.x4.x - q.x3.x + lc.a1*q.x3.x;
            j22 = lc.a1*q.x1.y + q.x2.y - lc.a1*q.x2.y - lc.a1*q.x4.y - q.x3.y + lc.a1*q.x3.y;
            
            // invert jacobian
            invj11 = 1/(j11*j22 - j12*j21) * j22;
            invj12 = 1/(j11*j22 - j12*j21) * -j12;
            invj21 = 1/(j11*j22 - j12*j21) * -j21;
            invj22 = 1/(j11*j22 - j12*j21) * j11;
            
            // add diference in local coordinates
            lc.a1 = lc.a1 + (invj11 * difference.x + invj12 * difference.y);
            lc.a2 = lc.a2 + (invj21 * difference.x + invj22 * difference.y);
            
            // check where next cell is
            if (lc.a1 < 0 || lc.a1 > 1 || lc.a2 < 0 || lc.a2 >1){
                if (fabs(lc.a1) >= fabs(lc.a2)){
                    if (lc.a1 > 0){
                        getNeighborCell(4, temp, q);
                    }
                    else{
                        getNeighborCell(2, temp, q);
                    }
                }
                else{
                    if (lc.a2 > 0){
                        getNeighborCell(3, temp, q);
                    }
                    else{
                        getNeighborCell(1, temp, q);
                    }
                }
                // start at midpoint of new cell
                lc.a1 = 0.5;
                lc.a2 = 0.5;
            }
        }
    }
    // assume function always exits with success... ;)
    return 1;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------

/*
Function: checkNewPointPosition

Estimate whether the new Point x_new is within quad q. 

Returns -1 if this is the case. 
Returns 1 or 2 or 3 or 4, if this is not the case. This number specifies, in which neighbor cell or in which direction the new point is.
1 = left; 2 = bottom, 3 = right, 4 = top.

*/
int checkNewPointPosition(Quad q, Point x_old, Point x_new){

}

//----------------------------------------------------------------------------------------------------------------------------------------------------

/*
Function: getQuadContainingPoint

Input : Point p
Output: Quad q = quad in which p should be located

Overrides the quad q with the quad information of the new quad, which contains point p.

The return value is 1 if a new quad can be found, -1 if not.


HINT:
You can use the FindCell() function from the vtkPointSet parten class of vtkStructuredGrid to get a cell id of a quad, which contains a point.
You can use the GetPoint() function from the vtkPointSet parten class of vtkStructuredGrid to get the point associated with a point id.

Be careful with the order of the points of the new quad that should be returned!
*/
int getQuadContainingPoint(Point p, Quad & q){
	// use FindCell to find cellid, then GetCell, then GetPointIds and for each Point get the coordinates, and assign them to the right position at the quad
    
}

//----------------------------------------------------------------------------------------------------------------------------------------------------

/*
Function: getNeighborCell

Input : int edge_number
        Point p
Output: Quad q = Should contain the neighbor of the cell which contains point p

Overrides the quad q with the quad information of a neighbor, which contains point p.

edge_number specifies which border of the old quad should also be part of the new one.
1 = left; 2 = bottom, 3 = right, 4 = top.

The return value is 1 if a neighbor can be found and -1 if not.

HINT:
You can use the FindCell() function from the vtkPointSet parten class of vtkStructuredGrid to get a cell id of the quad, which contains a point.
You can use the GetCellNeighbors() function from the vtkStructuredGrid class to get neighbors of a cell.


Be careful with the order of the points of the new quad which should be returned.
*/
int getNeighborCell(int edge_number, Point p, Quad &q){
    double x[3] = {p.x, p.y, 0};
    double x2[3] = {0,0,0};
    double *weights = 0;
    vtkIdType cell_id = 0;
    int sublid = 0;
    //vtkCell *cell = new vtkCell();
    // vtkCell is abstract method an cann not be used. but with vtcCell *cell = NULL i get an segementation fault
	vtk_reader->GetStructuredGridOutput()->FindCell(x, NULL, cell_id, 0.0, sublid, x2, weights);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------

/*
Function: getQuadGridValues

Reads scalar values from the mesh and assigns them to the quad points.

HINT:
You can use the GetPointData()->GetScalars()->GetComponent() function of the vtkStructuredGrid class to get a scalar value component of a point. Select the 0 component.
*/
void getQuadGridValues(Quad &q, Point p){
	// normaly i would use the FindCell method to find the cellid, then i would use the GetCell and GetPointIds to get the ids
    // of the points for this cell. after that i would use GetComponent(pointid, 0) to get the value of thesepoints and asign them to the quad
}
//----------------------------------------------------------------------------------------------------------------------------------------------------
