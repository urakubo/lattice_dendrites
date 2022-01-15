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

## \mainpage Description
# pyLM is a Problem Solving Environment (PSE) for biological simulations [1].  Written in Python, it wraps and extends the highly optimized multi-GPU <a href="http://www.scs.illinois.edu/schulten/lm">Lattice Microbes</a> stochastic simulation software [2,3,4].  The PSE is comprised of a base set of functionality to set up, monitor and modify simulations, as well as a set of standard post-processing routines that interface to other Python packages, including NumPy, SciPy, H5py, iGraph to name a few. See [1] for additional information as well as the user guide on the main <a href="http://www.scs.illinois.edu/schulten/lm">website</a>.  If you use pyLM in your simulations, please cite references [1] and [3] below.
#
# <h3> References </h3>
# 1. J.R. Peterson, M.J. Hallock, J.A. Cole and Z. Luthey-Schulten. A Problem Solving Environment for Stochastic Biological Simulations. <a href="http://www.scs.illinois.edu/schulten/Papers/Peterson_et_al_PyHPC_2013.pdf"><em>PyHPC 2013</em> at <em>Supercomputing 2013</em>, 2013.</a>
# 2. E. Roberts, J.E. Stone, L. Sepulveda, W.W. Hwu, and Z. Luthey-Schulten. Long time-scale simulations of in vivo diffusion using GPU hardware. In <a href="http://www.scs.illinois.edu/schulten/Papers/Roberts_et_al_IPDPS_2009.pdf"><em>Proceedings of the 2009 IEEE International Symposium on Parallel & Distributed Processing</em>, 2009.</a>
# 3. E. Roberts, J. E. Stone, and Z. Luthey-Schulten. Lattice microbes: high-performance stochastic simulation method for the reaction-diffusion master equation. <a href="http://onlinelibrary.wiley.com/doi/10.1002/jcc.23130/pdf"><em>J. Comp. Chem.</em>, 32(3), 245-55, 2013.</a>
# 4. M.J. Hallock, J.E. Stone, E. Roberts, C. Fry and Z. Luthey-Schulten. Simulation of reaction diffusion processes over biologically-relevant size and time scales using multi-GPU workstations <a href="http://www.sciencedirect.com/science/article/pii/S0167819114000398"><em>Parallel Comput.</em> 40:86-99, 2014, doi:10.1016/j.parco.2014.03.009.</a>

##
# @package pyLM The problem solving environement which interfaces to Lattice Microbes simulation software

__all__ = ['LMLogger', 'CME', 'RDME', 'units', 'ipyInterface']

