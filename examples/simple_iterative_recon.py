import dicom
import os
import numpy
from matplotlib import animation
from matplotlib import pyplot
import scipy.io as sio
import pywt

file = "../data/k-space/himbeere.mat"

raw = sio.loadmat(file)
kspace = raw.get('k_space_Himbeere')

numpy.random.seed(0)
recon = (numpy.fft.fftshift(numpy.fft.fft2(kspace)))
pattern = numpy.random.random_sample(kspace.shape)
percent = 0.98
low_values_indices = pattern <= percent  # Where values are low
pattern[low_values_indices] = 0  # All low values set to 0

low_values_indices = pattern > percent  # Where values are low
pattern[low_values_indices] = 1  # All low values set to 0

kspace = kspace * pattern

def shrink(coeff, epsilon):
	shrink_values = (coeff < epsilon) #and (wavelet > -epsilon)
	high_values = coeff >= epsilon
	low_values = coeff <= -epsilon
	coeff[shrink_values] = 0
	coeff[high_values] -= epsilon
	coeff[low_values] += epsilon

def waveletShrinkage(current, epsilon):
	cA, (cH, cV, cD)  = pywt.dwt2(current, 'Haar')
	shrink(cA, epsilon)
	shrink(cH, epsilon)
	shrink(cV, epsilon)
	shrink(cD, epsilon)
	wavelet = cA, (cH, cV, cD)
	return pywt.idwt2(wavelet, 'Haar')
	

def updateData(k_space, current, step):
	# go to k-space
	update = numpy.fft.ifft2(numpy.fft.fftshift(current))
	# compute difference
	update = k_space - update
	# return to iamge space
	update = numpy.fft.fftshift(numpy.fft.fft2(update))
	update = current + (step * update)
	return update


current = numpy.zeros(kspace.size).reshape(kspace.shape)
first = updateData(kspace, current, 1)
early = first
i = 0
while i < 5:
	current = updateData(kspace, current, 0.01)
	current = waveletShrinkage(current, 0.001)
	if (i==0):
		early = current
	i += 1

#current = updateData(kspace, current, 0.1)

fig=pyplot.figure(dpi=180)
pyplot.subplot(221)
pyplot.set_cmap(pyplot.gray())
pyplot.imshow(abs(recon))
pyplot.subplot(222)
pyplot.set_cmap(pyplot.gray())
pyplot.imshow(abs(first))
pyplot.subplot(223)
pyplot.set_cmap(pyplot.gray())
pyplot.imshow(abs(early))
pyplot.subplot(224)
pyplot.set_cmap(pyplot.gray())
pyplot.imshow(abs(current))
pyplot.show()
