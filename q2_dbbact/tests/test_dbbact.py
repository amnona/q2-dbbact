import unittest
import logging
from os.path import join, dirname, abspath

import pandas as pd
from pandas.testing import assert_frame_equal
import qiime2
import biom
import numpy as np

from q2_dbbact.methods import diff_abundance, enrichment


class Test_DBBact(unittest.TestCase):
	def setUp(self):
		self.test_data_dir = join(dirname(abspath(__file__)), 'data')

	def test_diff_abund(self):
		# load the data
		table = qiime2.Artifact.load(join(self.test_data_dir, 'cfs-merged.qza')).view(biom.Table)
		metadata = qiime2.Metadata.load(join(self.test_data_dir, 'map.cfs.txt'))
		table_no_repseqs = qiime2.Artifact.load(join(self.test_data_dir, 'cfs-table.qza')).view(biom.Table)
		rep_seqs = join(self.test_data_dir, 'dna-sequences.fasta')

		res = diff_abundance(table=table, metadata=metadata, field='Subject', random_seed=2021)
		# did we get the correct number of results?
		self.assertEqual(len(res), 53)
		# and correct number of positive?
		self.assertEqual(np.sum(res['dir'] == 0), 16)

		# integrate representative seqs
		res2 = diff_abundance(table=table_no_repseqs, metadata=metadata, field='Subject', random_seed=2021, repseqs=rep_seqs)
		# did we get the correct number of results?
		self.assertEqual(len(res), 53)
		# and correct number of positive?
		self.assertEqual(np.sum(res['dir'] == 0), 16)

		# and make sure results are the same
		assert_frame_equal(res, res2)

	def test_enrichment(self):
		# to generate the songbird results we use:
		# QIIME 2 Plugin 'songbird' version 1.0.4 (from package 'songbird' version 1.0.4)
		# qiime songbird multinomial --i-table cfs-merged.qza --m-metadata-file map.cfs.txt --p-formula Subject --o-differentials diff-cfs-songbird.qza --o-regression-stats sb-stats.qza --o-regression-biplot sb-plot
		sb_table = qiime2.Artifact.load(join(self.test_data_dir, 'diff-cfs-songbird.qza')).view(pd.DataFrame)
		res = enrichment(diff=sb_table, source='songbird', maxid=6279, random_seed=2021)
		self.assertEqual(len(res), 439)
		self.assertEqual(np.sum(res['odif'] > 0), 97)
		self.assertEqual(res.index.values[0], 'rural community')
		self.assertEqual(res.index.values[-1], 'LOWER IN male')

		# correlation enrichment (i.e. order instead of 2 groups)
		res = enrichment(diff=sb_table, method='correlation', source='songbird', maxid=6279, random_seed=2021)
		print(res)
		print(res[:30])
		print(res[-30:])
		print(len(res))
		print(np.sum(res['odif'] > 0))
		print(res.iloc[0]['odif'])
		print(res.index.values[0])
		self.assertEqual(len(res), 495)
		self.assertEqual(np.sum(res['odif'] > 0), 182)
		self.assertEqual(res.index.values[4], '-city')
		self.assertEqual(res.index.values[-1], '-thailand')
