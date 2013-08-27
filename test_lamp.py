#!/usr/bin/env python

# Testcode for LAMP

__author__ = "Aika Terada"

import unittest, sys, datetime
import lamp

class TestLamp(unittest.TestCase):
	def setUp(self):
		self.csv_file = "sample/sample_item.csv"
		self.flag_file = "sample/sample_expression_over1.csv"
		self.value_file = "sample/sample_expression_value.csv"
		self.sig_level = 0.05

	def testFisher(self):
		d = datetime.datetime.today()
		log_file = "lamp_testfisher_log_" + d.strftime("%Y%m%d") + "_" + d.strftime("%H%M%S") + ".txt"
		true_k = 5; true_p = 0.00699300699301; true_comb = set(["TF1", "TF2", "TF3"])
		enrich_lst, k, columnid2name \
					= lamp.run( self.csv_file, self.flag_file, self.sig_level, "fisher", None, -1, log_file, ',' )
		sys.stderr.write("check correction factor...")
		self.assertAlmostEqual(k, true_k)
		sys.stderr.write("\n")

		sys.stderr.write("check the significance combinations...\n")
		for comb in enrich_lst:
			sys.stderr.write("   check combination name...")
			detect_set = set()
			for i in comb[0]:
				detect_set.add( columnid2name[i-1] )
			flag = detect_set.issubset( true_comb ) and true_comb.issubset( detect_set )
			self.assertTrue( flag )
			sys.stderr.write("\n")
			sys.stderr.write("   check P-value..")
			self.assertAlmostEqual( comb[1], true_p )

	def testUTest(self):
		d = datetime.datetime.today()
		log_file = "lamp_testutest_log_" + d.strftime("%Y%m%d") + "_" + d.strftime("%H%M%S") + ".txt"
		true_k = 5; true_p = 0.00602414187918; true_zscore = 2.510727
		true_comb = set(["TF1", "TF2", "TF3"])
		enrich_lst, k, columnid2name \
					= lamp.run( self.csv_file, self.value_file, self.sig_level, "u_test", None, -1, log_file, ',' )
		
		sys.stderr.write("check correction factor...")
		self.assertAlmostEqual(k, true_k)
		sys.stderr.write("\n")

		sys.stderr.write("check the significance combinations...\n")
		for comb in enrich_lst:
			sys.stderr.write("   check combination name...")
			detect_set = set()
			for i in comb[0]:
				detect_set.add( columnid2name[i-1] )
			flag = detect_set.issubset( true_comb ) and true_comb.issubset( detect_set )
			self.assertTrue( flag )
			sys.stderr.write("\n")
			sys.stderr.write("   check U-value..")
			self.assertAlmostEqual( comb[3], true_zscore )
			sys.stderr.write("\n")
			sys.stderr.write("   check P-value..")
			self.assertAlmostEqual( comb[1], true_p )

if __name__ == '__main__':
	unittest.main()
