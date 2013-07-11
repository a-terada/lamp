#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Table for storing P-value.
# This source is used in Fisher's exact test and Chi-square test.
# @author Terada, 16, Apr, 2013

import sys

class PvalTable():
	def __init__( self, row_size ):
		self.table = {}
		"""
		for row in xrange( 0, row_size + 1 ):
			row_list = [-1]
			for col in xrange( 0, row ):
				row_list.append(-1)
			self.table.append( row_list )
		"""
	
	def getValue( self, row, col ):
		try:
			return self.table[row][col]
		except KeyError:
			return -1
#		sys.stderr.write("row: %s, col: %s\n" % (row, col))
#		return self.table[row][col]

	def putValue( self, row, col, pval ):
		if not row in self.table:
			self.table[row] = {}
		self.table[row][col] = pval

	def hashSize( self ):
		size = 0
		for row in self.table:
			size = size + len( self.table )
		return size

	def output( self ):
		for i in xrange( 0, len(self.table) ):
			row = self.table[i]
			sys.stdout.write("[%s]" % i)
			for j in xrange( 0, len(row) ):
				sys.stdout.write(" %s:%s" % (j, row[j]))
			sys.stdout.write("\n")
	
