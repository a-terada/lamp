#!/usr/bin/env python

# Testcode for LCM-LAMP

__author__ = "Aika Terada"

import unittest, sys, datetime
import frepattern.frequentPatterns as frequentPatterns

class TestLcmLamp(unittest.TestCase):
	def setUp(self):
		self.sample_infile = "./sample/sample_item.4lcmlamp.txt"
		self.sample_n1 = 7
		self.yeast_infile = "./sample/yeast_tfsite.4lcmlamp.txt"
		self.yeast_n1 = 530
		self.egf_infile= "./sample/gse6462_item.4lcmlamp.txt"
		self.egf_n1 = 1129
		self.sig_level = 0.05
		
		self.lcm_path = "./lcm53/lcm"
		d = datetime.datetime.today()
		self.log_file = "lcmlamp_test_" + d.strftime("%Y%m%d") + "_" + d.strftime("%H%M%S") + "_log.txt"

	def getNums( self, outfile ):
		f = open( outfile, 'r' ); line = ""
		line = f.readline()
		s = line[:-1].split( ',' )
		trans_num = int( s[1].split( ' ' )[1] )
		items_num = int( s[2].split( ' ' )[1] ) - 1
		f.close()
		return trans_num, items_num

	
	def checkLCMLAMP( self, infile, n1, max_lambda, max_comb, outfile, \
					  true_items, true_trans, true_minsup ):
		outlog = open( self.log_file, "a+" )
		fre_pattern = frequentPatterns.LCM( self.lcm_path, max_lambda, outlog )
		lam = fre_pattern.runLCMLAMP( infile, max_comb, n1, self.sig_level )
		outlog.close()
		
		sys.stderr.write("check # of transactions and # of items...\n")
		trans_num, items_num = self.getNums( outfile )
		self.assertAlmostEqual( trans_num, true_trans )
		self.assertAlmostEqual( items_num, true_items )		
		sys.stderr.write("check minimum support...\n")
		self.assertAlmostEqual( lam, true_minsup )
	
	def testSample(self):
		sys.stderr.write( "\n\n#######################################\n")
		sys.stderr.write( "  Test LCM-LAMP with sample file \n" )
		sys.stderr.write( "#######################################\n")
		max_lambda = 7; true_items = 4; true_trans = 15
		
		sys.stderr.write( "--- without arity limit ---\n" )
		max_comb = -1; true_minsup = 5; true_k = 5
		outfile = self.sample_infile + ".results.lcm/sample_item.4lcmlamp.txt.lcmlamp.closed"
		self.checkLCMLAMP( self.sample_infile, self.sample_n1, max_lambda, max_comb, outfile, true_items, true_trans, true_minsup )

		sys.stderr.write( "\n--- arity limit = 2 ---\n" )
		max_comb = 2; true_minsup = 5; true_k = 7
		outfile = self.sample_infile + ".results.lcm/sample_item.4lcmlamp.txt.lcmlamp.aritylim" + str( max_comb )
		self.checkLCMLAMP( self.sample_infile, self.sample_n1, max_lambda, max_comb, outfile, true_items, true_trans, true_minsup )
		
	
	def testYeast(self):
		sys.stderr.write( "\n\n#######################################\n")
		sys.stderr.write( "  Test LCM-LAMP with yeast file \n" )
		sys.stderr.write( "#######################################\n")
		max_lambda = 226; true_items = 102; true_trans = 6074
		
		sys.stderr.write( "--- without arity limit ---\n" )
		max_comb = -1; true_minsup = 4; true_k = 303
		outfile = self.yeast_infile + ".results.lcm/yeast_tfsite.4lcmlamp.txt.lcmlamp.closed"
		self.checkLCMLAMP( self.yeast_infile, self.yeast_n1, max_lambda, max_comb, outfile, true_items, true_trans, true_minsup )

		sys.stderr.write( "\n--- arity limit = 1 ---\n" )
		max_comb = 1; true_minsup = 4; true_k = 82
		outfile = self.yeast_infile + ".results.lcm/yeast_tfsite.4lcmlamp.txt.lcmlamp.aritylim" + str( max_comb )
		self.checkLCMLAMP( self.yeast_infile, self.yeast_n1, max_lambda, max_comb, outfile, true_items, true_trans, true_minsup )

		sys.stderr.write( "\n--- arity limit = 2 ---\n" )
		max_comb = 2; true_minsup = 4; true_k = 262
		outfile = self.yeast_infile + ".results.lcm/yeast_tfsite.4lcmlamp.txt.lcmlamp.aritylim" + str( max_comb )
		self.checkLCMLAMP( self.yeast_infile, self.yeast_n1, max_lambda, max_comb, outfile, true_items, true_trans, true_minsup )

		sys.stderr.write( "\n--- arity limit = 3 ---\n" )
		max_comb = 3; true_minsup = 4; true_k = 362
		outfile = self.yeast_infile + ".results.lcm/yeast_tfsite.4lcmlamp.txt.lcmlamp.aritylim" + str( max_comb )
		self.checkLCMLAMP( self.yeast_infile, self.yeast_n1, max_lambda, max_comb, outfile, true_items, true_trans, true_minsup )

		sys.stderr.write( "\n--- arity limit = 4 ---\n" )
		max_comb = 4; true_minsup = 4; true_k = 394
		outfile = self.yeast_infile + ".results.lcm/yeast_tfsite.4lcmlamp.txt.lcmlamp.aritylim" + str( max_comb )
		self.checkLCMLAMP( self.yeast_infile, self.yeast_n1, max_lambda, max_comb, outfile, true_items, true_trans, true_minsup )

		sys.stderr.write( "\n--- arity limit = 5 ---\n" )
		max_comb = 5; true_minsup = 4; true_k = 398
		outfile = self.yeast_infile + ".results.lcm/yeast_tfsite.4lcmlamp.txt.lcmlamp.aritylim" + str( max_comb )
		self.checkLCMLAMP( self.yeast_infile, self.yeast_n1, max_lambda, max_comb, outfile, true_items, true_trans, true_minsup )


	def testEGF(self):
		sys.stderr.write( "\n\n#######################################\n")
		sys.stderr.write( "  Test LCM-LAMP with breast cancer file \n" )
		sys.stderr.write( "#######################################\n")
		max_lambda = 1129; true_items = 397; true_trans = 12773
		
		sys.stderr.write( "--- without arity limit ---\n" )
		max_comb = -1; true_minsup = 8; true_k = 3750336
		outfile = self.egf_infile + ".results.lcm/gse6462_item.4lcmlamp.txt.lcmlamp.closed"
		self.checkLCMLAMP( self.egf_infile, self.egf_n1, max_lambda, max_comb, outfile, true_items, true_trans, true_minsup )

		sys.stderr.write( "\n--- arity limit = 1 ---\n" )
		max_comb = 1; true_minsup = 4; true_k = 397
		outfile = self.egf_infile + ".results.lcm/gse6462_item.4lcmlamp.txt.lcmlamp.aritylim" + str( max_comb )
		self.checkLCMLAMP( self.egf_infile, self.egf_n1, max_lambda, max_comb, outfile, true_items, true_trans, true_minsup )
		
		sys.stderr.write( "\n--- arity limit = 2 ---\n" )
		max_comb = 2; true_minsup = 6; true_k = 46209
		outfile = self.egf_infile + ".results.lcm/gse6462_item.4lcmlamp.txt.lcmlamp.aritylim" + str( max_comb )
		self.checkLCMLAMP( self.egf_infile, self.egf_n1, max_lambda, max_comb, outfile, true_items, true_trans, true_minsup )
		
		sys.stderr.write( "\n--- arity limit = 3 ---\n" )
		max_comb = 3; true_minsup = 7; true_k = 534124
		outfile = self.egf_infile + ".results.lcm/gse6462_item.4lcmlamp.txt.lcmlamp.aritylim" + str( max_comb )
		self.checkLCMLAMP( self.egf_infile, self.egf_n1, max_lambda, max_comb, outfile, true_items, true_trans, true_minsup )
		
		
		sys.stderr.write( "\n--- arity limit = 4 ---\n" )
		max_comb = 4; true_minsup = 8; true_k = 1727712
		outfile = self.egf_infile + ".results.lcm/gse6462_item.4lcmlamp.txt.lcmlamp.aritylim" + str( max_comb )
		self.checkLCMLAMP( self.egf_infile, self.egf_n1, max_lambda, max_comb, outfile, true_items, true_trans, true_minsup )

		
if __name__ == '__main__':
	unittest.main()
