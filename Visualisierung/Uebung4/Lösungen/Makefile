includeDir += -I /home/boy2melo/Dokumente/Visualisierung/vtk_build
includeDir += -I /home/boy2melo/Dokumente/Visualisierung/vtk-5.8.0/Filtering
includeDir += -I /home/boy2melo/Dokumente/Visualisierung/vtk-5.8.0/Common
includeDir += -I /home/boy2melo/Dokumente/Visualisierung/vtk-5.8.0/Graphics
includeDir += -I /home/boy2melo/Dokumente/Visualisierung/vtk-5.8.0/Rendering
includeDir += -I /home/boy2melo/Dokumente/Visualisierung/vtk-5.8.0/IO

libDir += -L /usr/lib/vtk

libs += -l vtkCommon
libs += -l vtkFiltering
libs += -l vtkGraphics
libs += -l vtkHybrid
libs += -l vtkImaging
libs += -l vtkIO
libs += -l vtkpng
libs += -l vtkjpeg
libs += -l vtkRendering
libs += -l vtkzlib

all:
	g++ -Wfatal-errors source4.cpp -o ex4 $(includeDir) $(libDir) $(libs)
