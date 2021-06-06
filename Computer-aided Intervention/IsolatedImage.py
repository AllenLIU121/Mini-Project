#!/usr/bin/env python
# coding: utf-8
#SHOW IN COMMAND:
#python A3.py ./brain.mha ./brainSegmentation.mha


import SimpleITK as sitk 
import sys

Input_path = sys.argv[1]
Output_path = sys.argv[2]

img_Input = sitk.ReadImage(Input_path,sitk.sitkFloat32) 
smoothing = sitk.CurvatureAnisotropicDiffusionImageFilter() 
smoothing.SetTimeStep(0.125) 
smoothing.SetNumberOfIterations(2) 
smoothing = smoothing.Execute(img_Input)

segmentation = sitk.IsolatedConnectedImageFilter() 
segmentation.SetReplaceValue(255) 
segmentation.SetLower(90) 
segmentation.SetUpper(250)

seed1 = [111,150] 
seed2 = [265,277]

segmentation.SetSeed1(seed1) 
segmentation_left = segmentation.Execute(smoothing)

segmentation.SetSeed1(seed2) 
segmentation_right = segmentation.Execute(smoothing)

segmentation_right = sitk.GetArrayFromImage(segmentation_right) 
segmentation_left = sitk.GetArrayFromImage(segmentation_left)

for x in range(segmentation_right.shape[0]): 
    for y in range(segmentation_right.shape[1]): 
        if segmentation_right[x][y]== 255: 
            segmentation_left[x][y] = 255

segmentation = sitk.GetImageFromArray(segmentation_left)

writer = sitk.ImageFileWriter() 
writer.SetFileName(Output_path) 
writer.Execute(segmentation)





