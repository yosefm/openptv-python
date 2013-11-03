OpenPTV-Python (PyPTV)
======================

**OpenPTV-Python** (a.k.a **PyPTV**) is the Python version of [OpenPTV](http://www.openptv.net). It is basically the Python Traits GUI (from Enthought Inc.) that *interfaces* the OpenPTV library that includes all the core algorithms (correspondence, tracking, calibration, etc.) written in ANSI C. 

Both PyPTV and the OpenPTV library are in the development phase and continuously refactored. Please follow the development on the community mailing list:

	openptv@googlegroups.com


## Installation

Mac OS X users can use the [pre-compiled try-before-installing instructions](https://github.com/OpenPTV/openptv-python/wiki/Getting-started-using-OpenPTV-Python-on-Mac-OS-X)

### Prerequisites:

1. [Optional, recommended] Download and install the unit testing framework <http://check.sourceforge.net/> - without this one you cannot test the compilation properly, although you can install and run the PyPTV software.

2. Install the OpenPTV library (`liboptv`) from source: <https://github.com/openptv/openptv/>. Please follow closely the Installation instructions. 
	* If you are not interested in development, **download** the appropriate pre-compiled library files: [liboptv installers](http://goo.gl/MqDzP). 
The `include` folder content should be copied/installed to something like `/usr/local/include` and the `lib` content into `/usr/local/lib`. Otherwise
add the location of the `liboptv` files in the next step, e.g. `python setup.py install -I /usr/local/include -L /usr/local/lib`

3. Get one of the Python development environments for your platform. We recommend WinPython or PythonXY on Windows or 
Enthought Python Distribution on all the platforms. If you want to build all the prerequisities yourself (e.g. using Linux
package manager, `apt-get` or Mac OS X `homebrew` or `ports`) then you need Python 2.7, Cython, Numpy, Scipy, Matplotlib, ETS tools (TraitsUI, Chaco, Enable and many sub-packages)  
3. Download or clone the PyPTV repository: <http://github.com/openptv/openptv-python>. 
 * If you don't want to test, **skip** this step, proceed to PyPTV installation:
 * If you want to run the tests: use the same set of `autotools`. On Mac OS X you might need to add `CC="gcc -arch i386"` to the `./configure` command. 

			autoreconf --install
			./configure
			make
			make check
 

### Install PyPTV:

Install by compiling the Python/Cython interface to the `liboptv` library (unrelated to the tests in the previous step) and install it to the default locations:

		cd pyptv_gui
		python setup.py install
	
or if you don't want to install somewhere in the path, compile all the necessary libraries in the `../src_c` directory 

		python setup.py build_ext --inplace
if the setup process creates errors related to `liboptv` linking, please add the `/include` and `/lib` locations, e.g.:

		python setup.py install -I /usr/local/include -L /usr/local/lib
	


## Getting started:

If the compilation passed, open the terminal and run:  

		python pyptv_gui/pyptv_gui.py

Follow the instructions in our **screencasts and tutorials**:
  
  *  Tutorial 1: <http://youtu.be/S2fY5WFsFwo>  
  
  *  Tutorial 2: <http://www.youtube.com/watch?v=_JxFxwVDSt0>   
  
  *  Tutorial 3: <http://www.youtube.com/watch?v=z1eqFL5JIJc>  
  
  
Ask for help on our mailing list:

	openptv@googlegroups.com

