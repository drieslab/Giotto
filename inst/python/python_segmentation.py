#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue July 2 2024
@author: Pratishtha Guckhool
"""

import numpy as np
import time, os, sys
import matplotlib.pyplot as plt
import matplotlib as mpl
#matplotlib inline
mpl.rcParams['figure.dpi'] = 300

try:
    from skimage.segmentation import find_boundaries

    SKIMAGE_ENABLED = True 
except:
    SKIMAGE_ENABLED = False

#input: files = #name of image files - can be a single image or multiple images
#channels = [[cytoplasm, nucleus]]
#diameter = None (default) - should be adjusted

def cellsegmentation_cellpose(files, img_channels, diameter):

    from cellpose import models, io, plot, utils

    # DEFINE CELLPOSE MODEL
    # model_type='cyto' or model_type='nuclei'
    model = models.Cellpose(gpu=False, model_type='cyto')

    # define CHANNELS to run segmentation on
    # grayscale=0, R=1, G=2, B=3
    # channels = [cytoplasm, nucleus]
    # if NUCLEUS channel does not exist, set the second channel to 0
    # channels = [0,0]
    # IF ALL YOUR IMAGES ARE THE SAME TYPE, you can give a list with 2 elements
    # channels = [0,0] # IF YOU HAVE GRAYSCALE
    # channels = [2,3] # IF YOU HAVE G=cytoplasm and B=nucleus
    # channels = [2,1] # IF YOU HAVE G=cytoplasm and R=nucleus

    # or if you have different types of channels in each image
    # channels = [[2,3], [0,0], [0,0]]

    # if diameter is set to None, the size of the cells is estimated on a per-image basis
    # you can set the average cell `diameter` in pixels yourself (recommended) 
    # diameter can be a list or a single number for all images

    # you can run all in a list e.g.
    # >>> imgs = [io.imread(filename) in for filename in files]
    # >>> masks, flows, styles, diams = model.eval(imgs, diameter=None, channels=channels)
    # >>> io.masks_flows_to_seg(imgs, masks, flows, diams, files, channels)
    # >>> io.save_to_png(imgs, masks, flows, files)

    # or in a loop
    img = io.imread(files)
    masks, flows, styles, diams = model.eval(img, diameter=None, channels=img_channels)

    # save results so you can load in gui
    io.masks_flows_to_seg(img, masks, flows, diams, files, img_channels)

    def outline_view(img0,maski,color=[1,0,0], mode='inner'):
        """
        Generates a red outline overlay onto image. 
        """
    #     img0 = utils.rescale(img0)
        if len(img0.shape)<3:
    #         img0 = image_to_rgb(img0) broken, transposing some images...
            img0 = np.stack([img0]*3,axis=-1)
    
        if SKIMAGE_ENABLED:
            outlines = find_boundaries(masks,mode='inner') #not using masks_to_outlines as that gives border 'outlines'
        else:
            outlines = utils.masks_to_outlines(masks) #not using masks_to_outlines as that gives border 'outlines'

        return outlines

    outline = outline_view(img, masks, color=[1,0,0], mode='inner')

    # save results as png
    # io.save_to_png(img, masks, flows, files)

    return outline, masks 

print("done")
