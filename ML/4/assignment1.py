import os, os.path
import Image
import numpy as np

# assignment a)
print("Importing all images")
print("This may take a while...")
print("")

# allocate empty array of size (166, 77760)
# uncroped images have a size of (320*243) = 77760

dataset = np.zeros(shape=(166,77760), dtype=np.integer)

for root, _, files in os.walk("./yalefaces"):
    i = 0
    for f in files:
        if not f == "Readme.txt":
            im = Image.open(os.path.join(root,f))
            dataset[i]=np.array(im.getdata())
            i += 1

print("All files loaded, calculate mean Face")
print("-------------------------------------")
print("")

# assignment b)
mean = np.zeros(shape=(1,77760))
for i in xrange(len(dataset)):
    mean += dataset[i]
mean = mean/len(dataset)

for i in xrange(len(dataset)):
    dataset[i] - mean

print("Mean face calculated and dataset centered")
print("-----------------------------------------")
print("")
Image.fromarray(mean.reshape((243, 320), order="C")).show(title="Mean Face")

print("Singular Value decomposition")
print("----------------------------")
print("")
# assignment c)
U,d,V = np.linalg.svd(dataset, full_matrices=False)

# assignment d) we use p = 166 and restrict dataset later
Z = dataset.dot(V.T)



# assignment e)
first_image = Z[0]
axis=[1,2,3,4,5,6,7,8]
print("show images with p=" + str(axis))
print("---------------------------------")
print("")


for p in axis:
    temparray = first_image[0:p].dot(V[0:p,0:])
    Image.fromarray(temparray.reshape((243,320), order="C")).show(title='P='+str(p))
    raw_input("Press Enter to continue...")

