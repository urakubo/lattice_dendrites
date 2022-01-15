# -*- coding: utf-8 -*-
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
# Author(s): Joseph R. Peterson
# 
#

import lm


def getReactionString(rct, prd, rate):
	rxnStr = ""
	if isinstance(rct,str):
		if rct == '':
			rct = '∅'
		rct = [rct]
	if isinstance(prd,str):
		if prd == '':
			prd = '∅'
		prd = [prd]
	rxnStr += " + ".join(rct) 
	rxnStr += " ⟶ "
	rxnStr += " + ".join(prd)
	units = "?"
	if len(rct) == 0:
		units = "molecules/s"
	elif len(rct) == 1:
		units = "s⁻¹"
	elif len(rct) == 2:
		units = "molecule⁻¹sec⁻¹"
	return rxnStr, rate, units


def writeTable(columnNames,rows):
	# Table definitions
	cols = len(columnNames)
	headerStyle = '<td style="text-align:%s"><b>%s</b></td>'
	rowStyle    = '<td style="text-align:%s">%s</td>'
	def align(idx):
		if idx == 0:
			return "left"
		return "center"
	s  = ''
	# Write header
	s += '<table>'
	s += '<tr>%s</tr>'%("".join([headerStyle%(align(i), x) for i,x in enumerate(columnNames)]))
	for row in rows:
		s += '<tr>%s</tr>'%("".join([rowStyle%(align(i), x) for i,x in enumerate(row)]))
	# Write footer
	s += '</table>'
	return s
	
	
