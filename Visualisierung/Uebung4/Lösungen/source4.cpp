#include "vtkPointData.h"

#include "vtkArrayCalculator.h"

#include "vtkPlane.h"
#include "vtkCutter.h"
#include "vtkScalarBarActor.h"

#include "vtkOutlineFilter.h"

#include "vtkGenericDataObjectReader.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"

#include "vtkSmartPointer.h"

int main()
{


		//create object for reading VTK files
		vtkSmartPointer<vtkGenericDataObjectReader> vtk_reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
		//set file name
		vtk_reader->SetFileName("field.vtk");
		//update
		vtk_reader->Update();

		

		//Modules to display outline of the dataset
	
		//Create module to display the outline of the dataset
		vtkSmartPointer<vtkOutlineFilter> outline = vtkSmartPointer<vtkOutlineFilter>::New();
		//Connect it to dataset
		outline->SetInputConnection(vtk_reader->GetOutputPort());
		outline->Update();

		//create object to map the polygons of the outline
		vtkSmartPointer<vtkOpenGLPolyDataMapper> mapper_outline = vtkSmartPointer<vtkOpenGLPolyDataMapper>::New();
		//connect it to outline module
		mapper_outline->SetInputConnection(outline->GetOutputPort());
		mapper_outline->Update();

		//create actor object - required by VTK
		vtkSmartPointer<vtkActor> actor_outline = vtkSmartPointer<vtkActor>::New();
		//connect to mapper module
  		actor_outline->SetMapper(mapper_outline);







		//vtkArrayCalculator is used to compute the magnitude from the vector field

		//create calculator object
		vtkSmartPointer<vtkArrayCalculator> calculator = vtkSmartPointer<vtkArrayCalculator>::New();
		//connect with reader
		calculator->SetInput(vtk_reader->GetOutput());
		//set array name to compute on
		calculator->AddVectorArrayName("Field");
		// change the calculation string to compute the magnitude of the vectors with "calculator->SetFunction("norm(Field)");" 
		// commands for the calculation string can be found in the documentation of vtkArrayCalculator: http://www.vtk.org/doc/release/5.8/html/a00160.html
		calculator->SetFunction("norm(Field)");
		//set output name
		calculator->SetResultArrayName("Result");
		//update module
  		calculator->Update();





		// modules to extract a slice through the data and display it with color mapping

		// define a cutting plane with vtkPlane 
		vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
		// set origin of the plane to 32,32,32
		plane->SetOrigin(32.0,32.0,32.0);
		// set normal of the plane to 1,1,1
		plane->SetNormal(1.0,1.0,1.0);
		 
	  
		// use vtkCutter to extract a slice
		// create vtkCutter module
		vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
		// use the defined cutting plane as implicit function for the cutter module
		cutter->SetCutFunction(plane);
		// connect it to the data source - the array calculator
		cutter->SetInputConnection(calculator->GetOutputPort());
		//update module
		cutter->Update();
	 

		// create object to map the polygons 
		vtkSmartPointer<vtkOpenGLPolyDataMapper> mapper = vtkSmartPointer<vtkOpenGLPolyDataMapper>::New();
		// connect it to the cutting plane
		mapper->SetInputConnection(cutter->GetOutputPort());

		// get value range of the data to adapt the color mapping
		double value_range[2];
		calculator->GetOutput()->GetPointData()->GetScalars()->GetRange(value_range, 0);

		// set value range of the color mapping
		mapper->SetScalarRange(value_range);

		//create actor object - required by VTK
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		//connect to mapper module
		actor->SetMapper(mapper);



		// create the color legend with vtkScalarBarActor
		vtkSmartPointer<vtkScalarBarActor> scalar_bar = vtkSmartPointer<vtkScalarBarActor>::New();
		// set the used color map - get the color map from the mapper
		scalar_bar->SetLookupTable(mapper->GetLookupTable());
		// set the title of the legend
		scalar_bar->SetTitle("Legend");
		// set the number of labels
		// INSERT






		//create render object
		vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
		//connect with actor objects
  		renderer->AddActor(actor_outline);
		renderer->AddActor(actor);
		renderer->AddActor(scalar_bar);
		
		//set background color, default is black
		renderer->SetBackground(0.5, 0.5, 0.5);


		//create render window object
		vtkSmartPointer<vtkRenderWindow> render_window = vtkSmartPointer<vtkRenderWindow>::New();
		//connect with renderer
		render_window->AddRenderer(renderer);
		//set size of the output window
		render_window->SetSize(800, 800);


		//create interactor object handling user input
		vtkSmartPointer<vtkRenderWindowInteractor> render_window_interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
		//connect with render window
		render_window_interactor->SetRenderWindow(render_window);

		//create object defining the input style
		vtkSmartPointer<vtkInteractorStyleTrackballCamera> interactor_style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
		//connect with interactor object
		render_window_interactor->SetInteractorStyle( interactor_style );




		//init interactor
		render_window_interactor->Initialize();

		//start interaction and event loop
		render_window_interactor->Start();


  		return 0;
}


