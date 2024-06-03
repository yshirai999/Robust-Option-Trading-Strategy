import numpy as np
from scipy.io import loadmat
import pandas as pd
import matlab.engine
eng = matlab.engine.start_matlab()
eng.BGDataReader(nargout = 0)
mat = loadmat('bgset512.mat')

BG = mat['A']
