from pylab import *

x = linspace(-3, 3, 20)
y = linspace(-3, 3, 20)

x,y = meshgrid(x,y)

for alpha in range(-2,3):
    vx = -y
    vy = alpha*x

    quiver(x,y,vx,vy, pivot="middle", headwidth=3, headlength=4)
    savefig("vectorfield-alpha-"+str(alpha)+".png")
    show()
