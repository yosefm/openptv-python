OpenPTV-Python (PyPTV)
======================

**OpenPTV-Python** (a.k.a **PyPTV**) is the Python version of [OpenPTV](http://www.openptv.net). It is basically the Python Traits GUI (from Enthought Inc.) that *interfaces* the OpenPTV library that includes all the core algorithms (correspondence, tracking, calibration, etc.) written in ANSI C. 

Both PyPTV and the OpenPTV library are in the development phase and continuously refactored. Please follow the development on the community mailing list:

	openptv@googlegroups.com


## Installation

1. [Optional, recommended] Download and install the unit testing framework <http://check.sourceforge.net/> - without this one you cannot test the compilation properly, although you can install and run the PyPTV software.

2. Download the install the OpenPTV library (`liboptv`) <http://github.com/alexlib/openptv>. Please follow closely the Installation instructions. 

3. Download or clone the PyPTV repository: <http://github.com/alexlib/openptv-python>. If you want to run the tests: use the same set of `autotools`. If you don't want to test, skip this step (note that `CC="gcc -arch i386"` is only needed for Mac OS X, remove it on Linux):


		autoreconf --install
		./configure CC="gcc -arch i386"
		make
		make check
	
4. Compile the Python/Cython interface to the `liboptv` library (unrelated to the tests in step 3) and install it to the default locations:

		cd pyptv_gui
		python setup.py install
	
or if you don't want to install somewhere in the path, compile all the necessary libraries in the `../src_c` directory 

		python setup.py build_ext --inplace
	


## Getting started:

If the compilation passed, open the terminal and run:  

		python pyptv_gui/pyptv_gui.py

Follow the instructions in our **screencasts and tutorials**:
  
  *  Tutorial 1: <http://youtu.be/S2fY5WFsFwo>  
  
  *  Tutorial 2: <http://www.youtube.com/watch?v=_JxFxwVDSt0>   
  
  *  Tutorial 3: <http://www.youtube.com/watch?v=z1eqFL5JIJc>  
  
  
Ask for help on our mailing list:

	openptv@googlegroups.com

