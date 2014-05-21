"""
:mod:`random` -- extended version of Python's :class:`random` library
=====================================================================

This provides additional random variates beyond those in the `random` 
libray.

 - logisticvariate
 - 
 Author:
 CreationDate: Wed May 21 13:33:23 2014
 Description:
     
 Version:
 Id:

"""
import random
import math

class Random(random.Random):
    def __init__(self, x=None):
        random.Random.__init__(self, x)
    
    def logisticvariate(self, mu, sigma):
        """Logistic distribution.
        
        mu is the location parameter, sigma the scale parameter
        """
        u1 = self.random()
        return mu + sigma * math.log(u1 / (1 - u1))

    def cauchyvariate(self, x0, gamma):
        """Cauchy distribution
        
        """
        u1 = self.random()
        return x0 + gamma * math.tan(math.pi * (u1 - 0.5))