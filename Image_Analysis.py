# -*- coding: utf-8 -*-
"""
Created on Wed May 26 17:08:02 2021

@author: Administrator
"""
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
import os
from collections import Counter
import math
from numpy.fft import fft, ifft, ifft2, fftfreq, fft2, fftshift, fftn, ifftn

class MedicalImage:
    def NormalizeImage(image):
        temp = image.reshape(image.shape[0], image.shape[1], -1)
        return temp
    
    def __init__(self, path:str(), name:str()):
        self.image = nib.load(path)
        self.image_array = self.image.get_fdata().astype(int)
        self.n_image = NormalizeImage(self.image_array)
        self.mid_slice = self.n_image[:,:,self.n_image.shape[2]//2]
        self.mid_pos = self.n_image.shape[2]//2
        self.name = name
        self.m_contrast = self.MichelsonContrast()
        self.rms_contrast = self.RMSContrast()
        self.entropy = self.Entropy()

    def MichelsonContrast(self):
        mintemp, maxtemp = np.min(self.image_array), np.max(self.image_array)
        if mintemp < 0:
            temp = self.image_array.copy()
            temp -= mintemp
            mintemp = 0
            maxtemp = maxtemp - mintemp
        return (maxtemp-mintemp) / (maxtemp+mintemp)

    def RMSContrast(self):
        return np.std(self.n_image)

    def Entropy(self):
        mintemp, maxtemp = np.min(self.n_image), np.max(self.n_image)
        nor_image = 255 * (self.n_image-mintemp) // (maxtemp-mintemp)
        hist = Counter(nor_image.flatten())
        total_count = nor_image.size
        entropy = 0
        for key, value in hist.items():
            entropy -= (value/total_count)*math.log(value/total_count)
        return entropy

    def SNR(self, piece, signal_patch, noise_patch):
        slice1 = self.n_image[:,:,piece]
        signal = slice1[signal_patch[0]:signal_patch[2],
                         signal_patch[1]:signal_patch[3]]
        noise = slice1[noise_patch[0]:noise_patch[2], noise_patch[1]:noise_patch[3]]
        return noise, signal.mean()/noise.std()

    def SNR_Hist(self, piece, signal_patch, noise_patch):
         slice1, snr = self.SNR(piece, signal_patch, noise_patch)
         return slice1, snr
 
    def ShowSNR(self, piece, signal_patch, noise_patch):
        snr = self.SNR(piece, signal_patch, noise_patch)
        print(snr)
        plt.subplot(1,2,1)
        plt.imshow(self.n_image[:,:,piece])
        plt.axis('off')
        plt.subplot(1,2,2)
        plt.hist(self.n_image[:,:,piece].flatten(), bins=255)
        plt.xlim(0, 255)
        plt.title("SNR: %s"%snr)
        plt.show()

    def LinearFilter(self, sigma):
        img_freqs = fftshift(fftn(self.n_image))
        sz_x, sz_y, sz_z = self.n_image.shape
        [X, Y, Z] = np.mgrid[0:sz_x, 0:sz_y, 0:sz_z]
        xpr, ypr, zpr = X-sz_x//2, Y-sz_y//2, Z-sz_z//2
    
        gaussfilt = np.exp(-((xpr**2+ypr**2+zpr**2)/(2*sigma**2))) /(2*np.pi*sigma**2)
        gaussfilt = gaussfilt / np.max(gaussfilt)
    
        filtered_freqs = img_freqs * gaussfilt
    
        filtered = np.abs(ifftn(fftshift(filtered_freqs)))
        return filtered

    def ShowConstrast(self):
        print(self.name, "Michelson Contrast: %.2f"%self.m_contrast)
        print(self.name, "RMS Contrast: %.2f"%self.rms_contrast)
        print(self.name, "Entropy: %.2f"%self.entropy)
        print()



path="C:/Users/Administrator.BF-201906221526\Desktop/Github upload\Medical_Imaging/modalities/"
files = os.listdir(path)
images ={}
medical_images={}
for file in files:
    if not os.path.isdir(file):
        if file[-3:] == '.gz':
            medical_images[file.replace(".nii.gz", "")] = MedicalImage(path+file,file.replace(".nii.gz", "") )

for i, (filename, img) in enumerate(medical_images.items()):
    plt.subplot(3,3,i+1)
    plt.title(filename)
    plt.axis('off')
    plt.imshow(img.mid_slice, cmap='jet')
plt.subplots_adjust(wspace =0.5, hspace =0.5)
plt.suptitle('Images')
plt.show()
swi = np.min(medical_images['swi'].image_array[:,:,200:300], axis=2)
tof = np.max(medical_images['tof'].image_array, axis=2)
plt.imshow(swi, cmap='jet')
plt.axis('off')
plt.title('swi mip')
plt.show()
plt.imshow(tof,cmap='jet', vmin=30, vmax=200)
plt.axis('off')
plt.title('tof mip')
plt.show()

for i , (name, medical_img) in enumerate(medical_images.items()):
    plt.subplot(3,3,i+1)
    plt.title(name+' %.2f'%medical_img.m_contrast)
    plt.axis('off')
    plt.imshow(medical_img.mid_slice, cmap='jet')
plt.subplots_adjust(wspace =0.5, hspace =0.5)
plt.suptitle('Mechelson Contrast')
plt.show()
for i , (name, medical_img) in enumerate(medical_images.items()):
    plt.subplot(3,3,i+1)
    plt.title(name+' %.2f'%medical_img.rms_contrast)
    plt.axis('off')
    plt.imshow(medical_img.mid_slice, cmap='jet')
plt.subplots_adjust(wspace =1, hspace =1)
plt.suptitle('RMS Contrast')
plt.show()
for i , (name, medical_img) in enumerate(medical_images.items()):
    plt.subplot(3,3,i+1)
    plt.title(name+' %.2f'%medical_img.entropy)
    plt.axis('off')
    plt.imshow(medical_img.mid_slice, cmap='jet')
plt.subplots_adjust(wspace =1, hspace =1)
plt.suptitle('Entropy')
plt.show()

hist_img,snr=medical_images['cardiac_axial'].SNR_Hist(medical_images['cardiac_axial'].mid_pos, [150,130,200,180], [220,240,270,290])
plt.subplot(3,3,1)
plt.title('cardiac_axial '+ '%.2f'%snr)
plt.hist(hist_img)
# no noise
hist_img,snr = medical_images['cardiac_realtime'].SNR_Hist(medical_images['cardiac_realtime'].mid_pos, [70,60,90,80], [20,100,40,120])
plt.subplot(3,3,2)
plt.title('cardiac_realtime '+ '%.2f'%snr)
plt.hist(hist_img)
# no noise
hist_img, snr=medical_images['ct'].SNR_Hist(medical_images['ct'].mid_pos,
[200,150,300,250], [450,300,500,350])
plt.subplot(3,3,3)
plt.title('ct '+ '%.2f'%snr)
plt.hist(hist_img)
# little noise
hist_img, snr=medical_images['fmri'].SNR_Hist(medical_images['fmri'].mid_pos,
[60,40,90,70], [120,25,150,55])
plt.subplot(3,3,4)
plt.title('fmri '+ '%.2f'%snr)
plt.hist(hist_img)
# little noise
hist_img, snr=medical_images['swi'].SNR_Hist(medical_images['swi'].mid_pos,[200,100,300,200], [400,0,500,100])
plt.subplot(3,3,5)
plt.title('swi '+ '%.2f'%snr)
plt.hist(hist_img)
# Gaussian noise
hist_img,snr=medical_images['T1_with_tumor'].SNR_Hist(medical_images['T1_with_tumor'].mid_pos, [60,80,100,120], [120,10,150,40])
plt.subplot(3,3,6)
plt.title('T1_with_tumor '+ '%.2f'%snr)
plt.hist(hist_img)
# few noise
hist_img,snr=medical_images['meanpet'].SNR_Hist(medical_images['meanpet'].mid_pos,[100,80,150,130], [150,0,200,50])
plt.subplot(3,3,7)
plt.title('meanpet '+ '%.2f'%snr)
plt.hist(hist_img)
# Gaussian noise
hist_img, snr=medical_images['tof'].SNR_Hist(medical_images['tof'].mid_pos,[100,100,200,200], [0,0,50,50])
plt.subplot(3,3,8)
plt.title('tof '+ '%.2f'%snr)
plt.hist(hist_img)
plt.subplots_adjust(wspace =1, hspace =1)
plt.suptitle("SNR")
plt.show()

for i, (filename, img) in enumerate(medical_images.items()):
    linear_img = img.LinearFilter(2)
    plt.subplot(3,3,i+1)
    plt.title(filename)
    plt.axis('off')
    plt.imshow(linear_img[:,:,linear_img.shape[2]//2])
plt.suptitle("Sigma=2")
plt.show()
for i, (filename, img) in enumerate(medical_images.items()):
    linear_img = img.LinearFilter(4)
    plt.subplot(3,3,i+1)
    plt.title(filename)
    plt.axis('off')
    plt.imshow(linear_img[:,:,linear_img.shape[2]//2])
plt.suptitle("Sigma=4")
plt.show()
for i, (filename, img) in enumerate(medical_images.items()):
    linear_img = img.LinearFilter(15)
    plt.subplot(3,3,i+1)
    plt.title(filename)
    plt.axis('off')
    plt.imshow(linear_img[:,:,linear_img.shape[2]//2])
plt.suptitle("Sigma=15")
plt.show()


















