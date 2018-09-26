FROM ubuntu:14.04


########################################################################
# Basic setup
########################################################################

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update

# Python
RUN apt-get install -y python python-setuptools python-pip
# python-dev (needed for pip installs)

########################################################################
# Travis apt packages
########################################################################

RUN apt-get install -y build-essential
RUN apt-get install -y libgeos-c1
RUN apt-get install -y libhdf5-serial-dev
RUN apt-get install -y libatlas-base-dev
RUN apt-get install -y gfortran
RUN apt-get install -y libnetcdf-dev
RUN apt-get install -y libblas3gf
RUN apt-get install -y libc6
RUN apt-get install -y libgcc1
RUN apt-get install -y libgfortran3
RUN apt-get install -y liblapack3gf
RUN apt-get install -y libstdc++6
RUN apt-get install -y libgdal1h
RUN apt-get install -y gdal-bin
RUN apt-get install -y libgdal-dev
RUN apt-get install -y python-gdal

########################################################################
# Travis install (pip install instead of conda install, pinned to conda's version)
########################################################################

# Additional steps for building python libraries
# for matplotlib
# RUN apt-get install -y libfreetype6-dev
# for basemap, as suggested by https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=741242
# RUN ln -s /usr/lib/libgeos-3.4.2.so /usr/lib/libgeos.so
# for netcdf4
# RUN pip install --upgrade setuptools
# for sqlite (which is not a python package)
# RUN apt-get install -y sqlite3 libsqlite3-dev



# TODO : move to top
RUN apt-get install -y python-dev

RUN pip install scipy==1.1.0
# begin matplotlib
    RUN apt-get install -y libfreetype6-dev
    RUN pip install matplotlib==2.2.3
# begendin matplotlib
# begin basemap
    # fix geos dep as suggested by https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=741242
    RUN ln -s /usr/lib/libgeos-3.4.2.so /usr/lib/libgeos.so
    # pip install from github as basemap is not in pypi
    RUN pip install https://github.com/matplotlib/basemap/archive/v1.0.7rel.tar.gz
# end basemap
RUN pip install shapely==1.6.4
RUN pip install nose==1.3.7
# begin netcdf4
RUN pip install --upgrade setuptools
RUN pip install netcdf4==1.4.1
# end netcdf4
RUN pip install cftime==1.0.0b1
RUN pip install coverage==4.5.1
RUN pip install coveralls==1.5.0
RUN pip install pycurl==7.43.0.2
RUN pip install pyproj==1.9.5.1
RUN pip install seaborn==0.9.0
RUN pip install simplejson==3.16.0
# begin sqlite (sqlite is not a python package)
# RUN pip install sqlite==3.24.0
RUN apt-get install -y sqlite3 libsqlite3-dev
# end sqlite
RUN pip install statsmodels==0.9.0
# begin gdal (this was already installed though apt packages)
# RUN pip install libgdal==2.3.1
# RUN pip install gdal==2.3.0
# end gdal


########################################################################
# Travis script
########################################################################

# missing an additional dep when running tests
RUN apt-get install -y python-tk

# Add the source
ADD . /home/src
WORKDIR /home/src

# Build
RUN python installer/setup.py build_ext -i

# Test (we need --exe as on Docker for windows, all mounted files are executable)
# RUN nosetests -v --with-coverage --exe

# RUN python tcevent.py -c example/yasi.ini
