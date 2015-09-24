rPGA
=======

rPGA is a pipeline to discover  hidden  splicing  variations  by  mapping  
personal transcriptomes to personal genomes .

Requirements
------------

rPGA is intended to be used in a Unix-based environment. It has been tested
on Mac OS and Linux.

A working installation of Python is required. In addition, rPGA requires the
following python packages: . These should be
installed automatically when you install the rPGA package, if you don't
already have them.

rPGA makes use of several external software packages, these are: STAR.
These must be installed. By default, rPGA expects them to be somewhere in your
path, but you if they are not you can specify their locations when running the
configure script (see below).

Installation
------------
All of the below installation instructions assume you have write access to the
installation directory for rPGA; if you want to install into a central
location for all users, you may need to prefix these with sudo.

### Unpack ###
To begin the installation, unpack the distribution and CD into the newly created
directory.

    tar -xf rPGA-0.0.1.tar.gz
    cd rPGA-0.0.1


### Configure ###
 Now type

    ./configure

This will configure the installation. rPGA is written partially in Python.
By default, it will use whichever python interpreter is invoked when you type
'Python' on the command line. If you wish to specify a different version of
Python, you can do this when running the ``configure`` script. Additionally,
if any of the required external software packages listed above are not on your
path, you can specify their location when running ``configure``. For example,
to specify the path to python and all of the exeternal software packages,
you would type the following:

    ./configure --with-python=/some/path/to/python \
                --with-STAR=/some/path/to/STAR

You needn't include them all, just the ones that cannot be found in your path.
The configure script will tell you if any requirements cannot be found, in
which case rPGA cannot be installed.

### Install ###
After configuration, rPGA is installed by typing:

    make install

This will install Python modules for whichever version of Python you are using.

User executables and scripts will be placed into /path/to/rPGA/bin and
/path/to/rPGA/scripts respectively. These cannot be moved, as the pipeline
expects them to remain here. We suggest you modify your PATH environment
variable to include these directories, but it is not required. All of the
below instructions assume these directories are in your PATH variable; if not,
replace the names of the scripts/executable with their full path.

Usage
-----
The rPGA manual is here:

Contacts and bug reports
------------------------
Yi D. Xing
yxing@ucla.edu

Shayna Stein
sstein93@ucla.edu

Emad Bahrami-Samani
ebs@ucla.edu

If you found a bug or mistake in this project, we would like to know about it.
Before you send us the bug report though, please check the following:

1. Are you using the latest version? The bug you found may already have been
   fixed.
2. Check that your input is in the correct format and you have selected the
   correct options.
3. Please reduce your input to the smallest possible size that still produces
   the bug; we will need your input data to reproduce the problem, and the
   smaller you can make it, the easier it will be.


Copyright and License Information
---------------------------------
Copyright (C) 2015 University of California, Los Angeles (UCLA)
Shayna R. Stein, Emad Bahrami-Samani, Yi Xing

Authors: Shayna R. Stein, Emad Bahrami-Samani, Yi Xing

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see http://www.gnu.org/licenses/.
