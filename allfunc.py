## This file all the related functions for this project
from struct import pack, unpack, calcsize, error, Struct
import os
import sys
import time
import array
import tempfile
import warnings
import io
from datetime import date
import operator
from osgeo import gdal, gdalnumeric, ogr, osr,gdal_array,gdalconst
from PIL import Image, ImageDraw
gdal.UseExceptions()
## --------------------------------------------------------
# This function will convert the rasterized clipper shapefile
# to a mask for use within GDAL.
## --------------------------------------------------------
def imageToArray(i):
    """
    Converts a Python Imaging Library array to a
    gdalnumeric image.
    """
    a = gdalnumeric.fromstring(i.tobytes(),'b')
    a.shape=i.im.size[1], i.im.size[0]
    return a
## --------------------------------------------------------
## Function to convert Array to Image (raster)
## --------------------------------------------------------
def arrayToImage(a):
    """
    Converts a gdalnumeric array to a
    Python Imaging Library Image.
    """
    i=Image.fromstring('L',(a.shape[1],a.shape[0]),
            (a.astype('b')).tobytes())
    return i
## --------------------------------------------------------
def world2Pixel(geoMatrix, x, y):
  """
  Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
  the pixel location of a geospatial coordinate
  """
  ulX = geoMatrix[0]
  ulY = geoMatrix[3]
  xDist = geoMatrix[1]
  yDist = geoMatrix[5]
  rtnX = geoMatrix[2]
  rtnY = geoMatrix[4]
  pixel = int((x - ulX) / xDist)
  line = int((ulY - y) / xDist)
  return (pixel, line)
## --------------------------------------------------------
#
#  EDIT: this is basically an overloaded
#  version of the gdal_array.OpenArray passing in xoff, yoff explicitly
#  so we can pass these params off to CopyDatasetInfo
## --------------------------------------------------------
def OpenArray( array, prototype_ds = None, xoff=0, yoff=0 ):
    ds = gdal.Open( gdalnumeric.GetArrayFilename(array) )

    if ds is not None and prototype_ds is not None:
        if type(prototype_ds).__name__ == 'str':
            prototype_ds = gdal.Open( prototype_ds )
        if prototype_ds is not None:
            gdalnumeric.CopyDatasetInfo( prototype_ds, ds, xoff=xoff, yoff=yoff )
    return ds

## --------------------------------------------------------

def histogram(a, bins=range(0,256)):
  """
  Histogram function for multi-dimensional array.
  a = array
  bins = range of numbers to match
  """
  fa = a.flat
  n = gdalnumeric.searchsorted(gdalnumeric.sort(fa), bins)
  n = gdalnumeric.concatenate([n, [len(fa)]])
  hist = n[1:]-n[:-1]
  return hist
## --------------------------------------------------------

  def stretch(a):
    """
    Performs a histogram stretch on a gdalnumeric array image.
    """
    hist = histogram(a)
    im = arrayToImage(a)
    lut = []
    for b in range(0, len(hist), 256):
      # step size
      step = reduce(operator.add, hist[b:b+256]) / 255
      # create equalization lookup table
      n = 0
      for i in range(256):
        lut.append(n / step)
        n = n + hist[i+b]
    im = im.point(lut)
    return imageToArray(im)

'''
# Returns number of given cells
    def returnLCnumber(self,array,cl):
        return int(self.count_nonzero(array[array==cl]))

    # Returns the proportion of the labeled class in the landscape
    def returnLCproportion(self,array,cl):
        res = []
        classes = sorted(numpy.unique(array)) # get classes
        # Value 0 seems to be the default nodata-value
        for i in classes:
            arr = numpy.copy(array)
            arr[array!=i] = 0
            res.append(self.count_nonzero(arr))
        arr = numpy.copy(array)
        arr[array!=cl] = 0
        prop = self.count_nonzero(arr) / float(sum(res))
        return prop

    ## Unclassified Methods ##
    # Returns sum of clipped raster cells
    def returnArraySum(self,array):
        try:
            return numpy.sum(array[array!=self.nodata])
        except ValueError:
            return None

    # Returns mean of clipped raster cells
    def returnArrayMean(self,array):
        try:
            return numpy.mean(array[array!=self.nodata])
        except ValueError:
            return None

    # Returns standard deviation of clipped raster cells
    def returnArrayStd(self,array):
        try:
            return numpy.std(array[array!=self.nodata])
        except ValueError:
            return None

    # Returns the minimum of clipped raster cells
    def returnArrayMin(self,array):
        if numpy.size(array) != 0 and self.count_nonzero(array) != 0:
            try:
                return numpy.min(array[array!=self.nodata])
            except ValueError: # doesn't work always if no-data is a normal integer?
                return None
        else:
            return None

    # Returns the minimum of clipped raster cells
    def returnArrayMax(self,array):
        if numpy.size(array) != 0 and self.count_nonzero(array) != 0:
            try:
                return numpy.max(array[array!=self.nodata])
            except ValueError:
                return None
        else:
            return None

    # Returns the median of clipped raster cells
    def returnArrayMedi(self,array):
        if numpy.size(array) != 0 and self.count_nonzero(array) != 0:
            try:
                return numpy.median(array[array!=self.nodata])
            except ValueError:
                return None
        else:
            return None

    # Returns the weighed average of clipped raster cells
    def returnArrayLowerQuant(self,array):
        if numpy.size(array) != 0 and self.count_nonzero(array) != 0:
            try:
                return scipy.percentile(array[array!=self.nodata],25)
            except ValueError:
                return None
        else:
            return None

    # Returns the weighed average of clipped raster cells
    def returnArrayHigherQuant(self,array):
        if numpy.size(array) != 0 and self.count_nonzero(array) != 0:
            try:
                return scipy.percentile(array[array!=self.nodata],75)
            except ValueError:
                return None
        else:
            return None

    # Calculates the Shannon Index
    def f_returnShannonIndex(self,array,classes):
        sh = []
        cl_array = numpy.copy(array) # create working array
        cl_array[cl_array==int(self.nodata)] = 0
        for cl in classes:
            res = []
            for i in classes:
                arr = numpy.copy(array)
                arr[array!=i] = 0
                res.append(self.count_nonzero(arr))
            arr = numpy.copy(array)
            arr[array!=cl] = 0
            prop = self.count_nonzero(arr) / float(sum(res))
            sh.append(prop * math.log(prop))
        return sum(sh)*-1

    # Calculates the Simpson Index
    def f_returnSimpsonIndex(self,array,classes):
        si = []
        cl_array = numpy.copy(array) # create working array
        cl_array[cl_array==int(self.nodata)] = 0
        for cl in classes:
            res = []
            for i in classes:
                arr = numpy.copy(array)
                arr[array!=i] = 0
                res.append(self.count_nonzero(arr))
            arr = numpy.copy(array)
            arr[array!=cl] = 0
            prop = self.count_nonzero(arr) / float(sum(res))
            si.append(math.pow(prop,2))
        return 1-sum(si)
    # Calculates the Shannon Equitability / Eveness
    def f_returnShannonEqui(self,array,classes):
        return self.f_returnShannonIndex(array,classes) / math.log(len(classes))
'''
