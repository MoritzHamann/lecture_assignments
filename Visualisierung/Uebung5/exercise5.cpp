#include "vtkMaskPoints.h"
#include "vtkPointData.h"

#include "vtkGlyphSource2D.h"
#include "vtkGlyph3D.h"

#include "vtkArrayCalculator.h"
//#include "vtkGradientFilter.h"
#include "my3dGradient.h"


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



	//create calculator object
	vtkSmartPointer<vtkArrayCalculator> calculator = vtkSmartPointer<vtkArrayCalculator>::New();
	//connect with reader
	calculator->SetInput(vtk_reader->GetOutput());
	//set array name
	calculator->AddVectorArrayName("Field");
	//define calculation string
	calculator->SetFunction("mag(Field)"); 
	//set output name
	calculator->SetResultArrayName("Result");
	//update module
  	calculator->Update();

    // insert our module
    vtkSmartPointer<my3dGradient> grad = vtkSmartPointer<my3dGradient>::New();
    grad->SetInput(calculator->GetOutput());
    grad->SetResultArrayName("Gradient");
    grad->Update();

	//create object to reduce number of data samples
	//without this module, all data samples (64x64x64 here) would be used
	//this can result in a slow and cluttered visualization
	vtkSmartPointer<vtkMaskPoints> mask = vtkSmartPointer<vtkMaskPoints>::New();
	//connect with reader object
	//mask->SetInput(vtk_reader->GetOutput());
    mask->SetInput(grad->GetOutput());
    

	//set number of sampling points
	mask->SetMaximumNumberOfPoints(10000);
	//activate random sampling
	//mask->RandomModeOn();
	mask->SetOnRatio(30);
	//update module
	mask->Update();

	//pitfall of VTK!!!
	//active vector field has to be set after vtkMaskPoints
	mask->GetOutput()->GetPointData()->SetActiveVectors("Gradient");


	//create object for generating the used glyph - a 2D arrow here
	vtkSmartPointer<vtkGlyphSource2D> arrow = vtkSmartPointer<vtkGlyphSource2D>::New();
	//set glyph type
	arrow->SetGlyphTypeToArrow();
	//set size
	arrow->SetScale(1.0);
	//set filled glyphs
	arrow->FilledOn();
	//update module
	arrow->Update();


	//create object for visualizing the vector field with glyphs
	vtkSmartPointer<vtkGlyph3D> glyph = vtkSmartPointer<vtkGlyph3D>::New();
	//connect to data source
	glyph->SetInput(mask->GetOutput());
	//connect to glyph source
	glyph->SetSource(arrow->GetOutput());
	//scale glyphs with vector magnitude
	glyph->ScalingOn();
	//orient glyphs along vectors
	glyph->OrientOn();
	//color glyphs by vector magnitude
	glyph->SetColorModeToColorByVector();
	//scale with vector magnitude
	glyph->SetScaleModeToScaleByVector();
	//set scale factor
	glyph->SetScaleFactor(30.0f);
	//update module
	glyph->Update();



	//create object to map the polygons of the glyph	
	vtkSmartPointer<vtkOpenGLPolyDataMapper> mapper = vtkSmartPointer<vtkOpenGLPolyDataMapper>::New();
	//connect with the glyph generator
	mapper->SetInputConnection( glyph->GetOutputPort() );
	//set scalar range for color mapping
	mapper->SetScalarRange(0.0, 0.2);


	//create actor object - required by VTK
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	//connect to mapper module
  	actor->SetMapper( mapper );







	//create render object
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	//connect with actor object
  	renderer->AddActor( actor );
	//set background color, default is black
	renderer->SetBackground(0.5, 0.5, 0.5);


	//create render window object
	vtkSmartPointer<vtkRenderWindow> render_window = vtkSmartPointer<vtkRenderWindow>::New();
	//connect with renderer
	render_window->AddRenderer( renderer );
	//set size of the output window
	render_window->SetSize( 300, 300 );


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


