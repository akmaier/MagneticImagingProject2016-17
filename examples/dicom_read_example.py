import dicom
import os
import numpy
from matplotlib import animation
from matplotlib import pyplot, cm


PathDicom = "../data/"
lstFilesDCM = []  # create an empty list
for dirName, subdirList, fileList in os.walk(PathDicom):
    for filename in fileList:
        if ".dcm" in filename.lower():  # check whether the file's DICOM
            lstFilesDCM.append(os.path.join(dirName,filename))
print (lstFilesDCM)

# Get ref file
RefDs = dicom.read_file(lstFilesDCM[0])

# Load dimensions based on the number of rows, columns, and slices (along the Z axis)
ConstPixelDims = (int(RefDs.Rows), int(RefDs.Columns), len(lstFilesDCM))

# Load spacing values (in mm)
ConstPixelSpacing = (float(RefDs.PixelSpacing[0]), float(RefDs.PixelSpacing[1]), float(RefDs.SliceThickness))

x = numpy.arange(0.0, (ConstPixelDims[0]+1)*ConstPixelSpacing[0], ConstPixelSpacing[0])
y = numpy.arange(0.0, (ConstPixelDims[1]+1)*ConstPixelSpacing[1], ConstPixelSpacing[1])
z = numpy.arange(0.0, (ConstPixelDims[2]+1)*ConstPixelSpacing[2], ConstPixelSpacing[2])

# The array is sized based on 'ConstPixelDims'
ArrayDicom = numpy.zeros(ConstPixelDims, dtype=RefDs.pixel_array.dtype)

# loop through all the DICOM files
for filenameDCM in lstFilesDCM:
    # read the file
    ds = dicom.read_file(filenameDCM)
    # store the raw image data
    ArrayDicom[:, :, lstFilesDCM.index(filenameDCM)] = ds.pixel_array  
    
    


    
fig=pyplot.figure(dpi=150)
pyplot.subplot(1,2,1)
#pyplot.axes().set_aspect('equal', 'datalim')
pyplot.set_cmap(pyplot.gray())
pyplot.pcolormesh(x, y, numpy.flipud(ArrayDicom[:, :, 17]))


def init():
    pyplot.subplot(1,2,1)
    return pyplot.pcolormesh(x, y, numpy.flipud(ArrayDicom[:, :, 0])),

def animate(i):
    pyplot.subplot(1,2,1)
    return pyplot.pcolormesh(x, y, numpy.flipud(ArrayDicom[:, :, i])),


anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=17, interval=1000, blit=False)

def onclick(event):
    print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          (event.button, event.x, event.y, event.xdata, event.ydata))
    pyplot.subplot(1,2,2)
    # note that the coordinates need to be converted from world coordinates to pixels.
    OnePoint = ArrayDicom[int(round(RefDs.Columns-(event.ydata/ConstPixelSpacing[1]))), int(round((event.xdata/ConstPixelSpacing[0]))), :]
    #print (ConstPixelSpacing[0], ConstPixelSpacing[1], round((event.xdata)), round((event.ydata)), OnePoint)
    pyplot.plot(OnePoint)
    pyplot.draw()

cid = fig.canvas.mpl_connect('button_press_event', onclick)

pyplot.subplot(1,2,2)
#Identified a good point at 120, 52
OnePoint = ArrayDicom[119, 52, :]
#print (OnePoint)
pyplot.plot(OnePoint)

pyplot.show()

