openptv-python
==============

PyPTV is the Python version of OpenPTV: Python GUI and core algorithms in C. 


This is the bleeding edge development branch. If you want to try it, it's on your own risk. 


Getting started
===============
See our screencasts and tutorials:


1.  Tutorial 1: http://youtu.be/S2fY5WFsFwo
2.  Tutorial 2: http://www.youtube.com/watch?v=_JxFxwVDSt0
3.  Tutorial 3: http://www.youtube.com/watch?v=z1eqFL5JIJc



Installation
============

Prerequisites:

One of the following distributions that include: Python, Cython, Numpy, Scipy, matplotlib, Enthought Tools Suite:

1. Enthought Python Distribution http://www.enthought.com/products/epd.php (academic version at least) or
2. PythonXY distribution http://code.google.com/p/pythonxy/ or
3. WinPython http://code.google.com/p/winpython/ or
4. Anaconda https://store.continuum.io/


Instructions:

1. Download and install liboptv <https://github.com/yosefm/openptv/tree/fb_submission>
2. Download the ZIP of this repository, unzip it and install from the pyptv_gui folder by:  
    python setup.py install
3. Run using   
    python pyptv_gui 

This version is fully backward-compatible with the 3D-PTV software http://3dptv.github.com so the test folders and data 
are the same. We do not destroy the parameters, we create new copies, so no worry about destroying the experiment.








