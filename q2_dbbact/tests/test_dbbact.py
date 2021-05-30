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
		rep_seqs = join(self.test_data_dir, 'cfs-repseqs.fasta')

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
		repseqs = join(self.test_data_dir, 'cfs-repseqs.fasta')
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
		self.assertEqual(len(res), 495)
		self.assertEqual(np.sum(res['odif'] > 0), 182)
		self.assertEqual(res.index.values[4], '-city')
		self.assertEqual(res.index.values[-1], '-thailand')

		# to generate aldex2 results we use:
		# QIIME 2 Plugin 'aldex2' version 1.14.1 (from package 'q2-aldex2' version 0+untagged.36.ga035a04.dirty)
		#  qiime aldex2 aldex2 --i-table cfs-table.qza --m-metadata-file map.cfs.txt --m-metadata-column Subject --o-differentials diff-cfs-aldex2.qza
		ad2_table = qiime2.Artifact.load(join(self.test_data_dir, 'diff-cfs-aldex2.qza')).view(pd.DataFrame)
		# we test with sig_threshold=0.5 and do correlation method since no significant results at 0.1
		res = enrichment(diff=ad2_table, repseqs=repseqs, source='aldex2', method='correlation', sig_threshold=0.5, maxid=6279, random_seed=2021)
		self.assertEqual(len(res), 365)
		self.assertEqual(np.sum(res['odif'] > 0), 101)
		self.assertEqual(res.index.values[0], 'small village')
		self.assertEqual(res.index.values[-1], '-physical activity')

		# to generate ancom results:
		# QIIME 2 Plugin 'composition' version 2020.6.0 (from package 'q2-composition' version 2020.6.0)
		# qiime composition add-pseudocount --i-table cfs-merged.qza --o-composition-table cfs-pseudocount.qza
		# qiime composition ancom --i-table cfs-pseudocount.qza --m-metadata-file map.cfs.txt --m-metadata-column Subject --output-dir diff-cfs-ancom --o-visualization diff-cfs-ancom.qzv
		# then need to save the ancom statistical results tsv from the ancom output qzv (named ancom.tsv)
		res = enrichment(diff_tsv=join(self.test_data_dir, 'diff-cfs-ancom-export.tsv'), source='ancom', maxid=6279, random_seed=2021)
		print(len(res))
		print(np.sum(res['odif'] > 0))
		print(res.iloc[0]['odif'])
		print(res.index.values[0])
		self.assertEqual(len(res), 365)
		self.assertEqual(np.sum(res['odif'] > 0), 101)
		self.assertEqual(res.index.values[0], 'small village')
		self.assertEqual(res.index.values[-1], '-physical activity')
