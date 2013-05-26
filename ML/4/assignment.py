import os, os.path
import Image
import numpy as np

print("Importing all images")
print("This may take a while...")

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



mean = np.zeros(shape=(1,77760))
for i in xrange(len(dataset)):
    mean += dataset[i]
mean = mean/len(dataset)

for i in xrange(len(dataset)):
    dataset[i] - mean

print("Mean face calculated and dataset centered")
Image.fromarray(mean.reshape((243, 320), order="C")).show(title="Mean Face")

U,d,V = np.linalg.svd(dataset)

Z = dataset.dot(V)

