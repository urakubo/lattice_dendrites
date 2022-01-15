# 
# University of Illinois Open Source License
# Copyright 2008-2016 Luthey-Schulten Group,
# All rights reserved.
# 
# Developed by: Luthey-Schulten Group
#                           University of Illinois at Urbana-Champaign
#                           http://www.scs.uiuc.edu/~schulten
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the Software), to deal with 
# the Software without restriction, including without limitation the rights to 
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
# of the Software, and to permit persons to whom the Software is furnished to 
# do so, subject to the following conditions:
# 
# - Redistributions of source code must retain the above copyright notice, 
# this list of conditions and the following disclaimers.
# 
# - Redistributions in binary form must reproduce the above copyright notice, 
# this list of conditions and the following disclaimers in the documentation 
# and/or other materials provided with the distribution.
# 
# - Neither the names of the Luthey-Schulten Group, University of Illinois at
# Urbana-Champaign, nor the names of its contributors may be used to endorse or
# promote products derived from this Software without specific prior written
# permission.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL 
# THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
# OTHER DEALINGS WITH THE SOFTWARE.
# 
# Author(s): Michael J. Hallock and Joseph R. Peterson
# 
#

import lm

def _map(f,xs):
        return [f(x) for x in xs]

############################
# Length Wrapper Functions #
############################
## Returns a representation of a number in angstroms
# @param qty A list or singleton of a number
def angstrom(*qty):
        if(len(qty) > 1):
                return _map(angstrom, qty)
        else:
                return float(qty[0])*1e-10

## Returns a representation of a number in nanometers
# @param qty A list or singleton of a number
def nm(*qty):
        if(len(qty) > 1):
                return _map(nm, qty)
        else:
                return float(qty[0])*1e-9

## Returns a representation of a number in micrometers
# @param qty A list or singleton of a number
def micron(*qty):
        if(len(qty) > 1):
                return _map(micron, qty)
        else:
                return float(qty[0])*1e-6

## Returns a representation of a number in millimeter
# @param qty A list or singleton of a number
def mm(*qty):
        if(len(qty) > 1):
                return _map(mm, qty)
        else:
                return float(qty[0])*1e-3

## Returns a representation of a number in centimeter
# @param qty A list or singleton of a number
def cm(*qty):
        if(len(qty) > 1):
                return _map(cm, qty)
        else:
                return float(qty[0])*1e-2

##########################
# Time Wrapper Functions #
##########################
## Returns a representation of a number in nanosecond
# @param qty A list or singleton of a number
def ns(*qty):
        if(len(qty) > 1):
                return _map(ns, qty)
        else:
                return float(qty[0])*1e-9

## Returns a representation of a number in microsecond
# @param qty A list or singleton of a number
def microsecond(*qty):
        if(len(qty) > 1):
                return _map(microsecond, qty)
        else:
                return float(qty[0])*1e-6

## Returns a representation of a number in millisecond
# @param qty A list or singleton of a number
def ms(*qty):
        if(len(qty) > 1):
                return _map(ms, qty)
        else:
                return float(qty[0])*1e-3

## Returns seconds in seconds.  Seems silly, but for completeness and
#  ability to annotate the unit in code.
def second(*qty):
        if(len(qty) > 1):
                return _map(second, qty)
        else:
                return float(qty[0])

## Returns a representation of a number in minutes
# @param qty A list or singleton of a number
def minute(*qty):
        if(len(qty) > 1):
                return _map(minute, qty)
        else:
                return float(qty[0])*60.0

## Returns a representation of a number in hours
# @param qty A list or singleton of a number
def hr(*qty):
        if(len(qty) > 1):
                return _map(hr, qty)
        else:
                return float(qty[0])*3600.0

## Returns a representation of a number in days
# @param qty A list or singleton of a number
def day(*qty):
        if(len(qty) > 1):
                return _map(day, qty)
        else:
                return float(qty[0])*3600.0*24.0


