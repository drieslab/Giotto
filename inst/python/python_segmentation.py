#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 7 01:00:00 2022
@author: rafaelspeixoto
"""

from deepcell.applications import Mesmer
from matplotlib import pyplot as plt
import numpy as np
import glob

def python_create_mesmer_app():
    # Create model
    app = Mesmer()

    print('Training Resolution:', app.model_mpp, 'microns per pixel')

    return app

def python_segment_image(app, image_matrix, file_name):
    image_matrix = [image_matrix]
    image_matrix = np.stack([image_matrix, image_matrix], axis=-1)

    # Predict
    segmentation_predictions_nuc = app.predict(image_matrix, image_mpp=0.5, compartment='nuclear')

    # Save
    plt.imsave(file_name, segmentation_predictions_nuc[0, ..., 0])

    return 'Segmentation Concluded'