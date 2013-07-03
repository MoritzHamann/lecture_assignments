// Lines.cpp 
//

#include "stdafx.h"
#include "Lines.h"



Lines::Lines(COGL4CoreAPI *Api) : RenderPlugin(Api) {
    this->myName = "Lines";
    this->myDescription = "";

	streamline_steps = 200;
	streamline_delta = 0.01;
}

Lines::~Lines() {

}
bool Lines::Resize(int width, int height) {
   if(height==0) height=1;		//no div 0
   glViewport(0, 0, width, height);	
   
   glm::mat4 ortogonal = glm::ortho(0.0f,(float)width/(float)height,0.0f,1.0f);

   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glLoadMatrixf(glm::value_ptr(ortogonal));

   return false;
}


bool Lines::Activate(void) {

	streamline_delta.Set(this,"streamline_delta", nullptr);
    streamline_delta.Register("step=0.001 min=0.001 max=0.5");

	streamline_steps.Set(this,"streamline_steps", nullptr);
    streamline_steps.Register("step=1 min=1 max=1000");
    
	this->AddManipulator("Panning", &this->trafo, Manipulator::MANIPULATOR_OBJ_PAN_ROT_SCALE_2D);
    return true;
}

bool Lines::Deactivate(void) {
    return true;
}

bool Lines::Init(void) {
    return true;
}

bool Lines::Render(void) {
    glClearColor( 0.0, 0.5, 1.0, 1.0 );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glLoadMatrixf(glm::value_ptr(this->trafo));
	//----------------------------

	glEnable(GL_BLEND);
  	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);	
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);


	glBegin(GL_LINE_LOOP);

		glColor3d(1.0, 1.0, 1.0);

		glVertex2d(0.0, 0.0);
		glVertex2d(1.0, 0.0);
		glVertex2d(1.0, 1.0);
		glVertex2d(0.0, 1.0);

	glEnd();

	this->draw_line_glyph(32, 32);

	double pos[2];

	pos[0] = 0.2;
	pos[1] = 0.2;


	glLineWidth(1.0);
	glColor3d(1.0, 0.0, 0.0);

	this->draw_streamline(pos, this->streamline_delta, this->streamline_steps, 0);


	pos[0] = 0.2;
	pos[1] = 0.2;

	glLineWidth(1.0);
	glColor3d(0.0, 1.0, 0.0);

	this->draw_streamline(pos, this->streamline_delta, this->streamline_steps, 1);


    return false;
}

//----------------------------------------

vector::vector()
{
	this->x = 0.0;
	this->y = 0.0;
}

vector vector::operator*(double val)
{
	this->x *= val;
	this->y *= val;

	return *this;
}

vector vector::operator+(vector vec)
{
	this->x += vec.x;
	this->y += vec.y;

	return *this;
}

//--------------------------------------

//generates the vector field
//returns the vector at position pos
vector Lines::get_vector(double pos[2])
{
	vector result;

	if(pos[0] < 0.0f || pos[0] > 1.0f || pos[1] < 0.0f || pos[1] > 1.0f)
	{
		result.x = result.y = 0.0;
	}
	else
	{
		double dx = 2.0 * pos[0];
		double dy = 2.0 * pos[1];

		double vx = -M_PI * sin(M_PI * dx) * cos(M_PI*dy);
		double vy = M_PI * cos(M_PI * dx) * sin(M_PI*dy);

		result.x = vx;
		result.y = vy;

	}



	return result;
}

//plug in your code for a 4th order Runke-Kutta integration step
//after one iteration, pos contains the new position
//delta is the step size
void Lines::do_rk4_step(double pos[2], double delta)
{
    double tmpPos[2];
    
    vector k1 = this->get_vector(pos);
    
    tmpPos[0] = pos[0] + k1.x*(delta/2.0);
    tmpPos[1] = pos[1] + k1.y*(delta/2.0);
    
    vector k2 = this->get_vector(tmpPos);
    
    tmpPos[0] = pos[0] + k2.x*(delta/2.0);
    tmpPos[1] = pos[1] + k2.y*(delta/2.0);
    
    vector k3 = this->get_vector(tmpPos);
    
    tmpPos[0] = pos[0] + k3.x*(delta);
    tmpPos[1] = pos[1] + k3.y*(delta);
    
    pos[0] = tmpPos[0];
    pos[1] = tmpPos[1];
}

//plug in your code for an euler integration step
//after one iteration, pos contains the new position
//delta is the step size
void Lines::do_euler_step(double pos[2], double delta)
{

    vector vec = this->get_vector(pos);
    
    pos[0] = pos[0] + delta*vec.x;
    pos[1] = pos[1] + delta*vec.y;

}


//plugin your code for computing and drawing the streamline,
//starting at position "pos" with step size "delta" and "steps" iterations
//"mode" selects the two integration schemes, euler or RK4,
// 0 - euler
// 1 - RK4
//
//you can use GL_LINES or GL_LINE_STRIP for rendering
//look at Lines::draw_line_glyph() for usage
//and the OpenGL reference
void Lines::draw_streamline(double pos[2], double delta, int steps, int mode)
{
	vector v;
	double prevPos[2];

    glLineWidth(1.0f);
	glBegin(GL_LINES);
    
    for( int i = 0; i < steps; i++ )
    {
		prevPos[0] = pos[0];
		prevPos[1] = pos[1];

        if(mode == 0)
        {
			
			glColor3d(1, 0, 0);
            this->do_euler_step(pos,delta);
        }
        else
        {
			glColor3d(0, 1.0, 0);
            this->do_rk4_step(pos,delta);
        }
        

        glVertex2d(prevPos[0], prevPos[1]);
        glVertex2d(pos[0], pos[1]);
    }
    
    glEnd();
}



void Lines::draw_line_glyph(int res_x, int res_y)
{
	vector v;
	double pos[2];
	double scale = 0.3 / (double) std::max(res_x, res_y);

	glLineWidth(3.0f);

	glBegin(GL_LINES);

	for(int y = 1; y < res_y; y++)
	{
		for(int x = 1; x < res_x; x++)
		{
			pos[0] = (double) x / (double) res_x;
			pos[1] = (double) y / (double) res_y;

			v = this->get_vector(pos);
			v.x *= scale;
			v.y *= scale;

			glColor3d(0.2, 0.2, 0.2);
			glVertex2d(pos[0], pos[1]);

			glColor3d(1.0, 1.0, 1.0);
			glVertex2d(pos[0] + v.x, pos[1] + v.y);


		}
	}

	glEnd();
}
