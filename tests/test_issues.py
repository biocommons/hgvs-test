# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import logging
import pprint
import re
import sys
import unittest

from nose.plugins.attrib import attr

from hgvs.exceptions import HGVSError, HGVSDataNotAvailableError, HGVSParseError, HGVSInvalidVariantError, HGVSInvalidVariantError
import hgvs.dataproviders.uta
import hgvs.normalizer
import hgvs.parser
import hgvs.transcriptmapper
import hgvs.validator
import hgvs.variant
import hgvs.variantmapper



@attr(tags=["issues"])
class Test_VariantMapper(unittest.TestCase):
    def setUp(self):
        self.hdp = hgvs.dataproviders.uta.connect()
        self.hm = hgvs.variantmapper.VariantMapper(self.hdp)
        self.hp = hgvs.parser.Parser()
        self.hv = hgvs.validator.IntrinsicValidator()
        self.evm37 = hgvs.variantmapper.EasyVariantMapper(self.hdp,
                                                          replace_reference=True,
                                                          assembly_name='GRCh37',
                                                          alt_aln_method='splign')
        self.evm38 = hgvs.variantmapper.EasyVariantMapper(self.hdp,
                                                          replace_reference=True,
                                                          assembly_name='GRCh38',
                                                          alt_aln_method='splign')
        self.vn = hgvs.normalizer.Normalizer(self.hdp, shuffle_direction=3, cross_boundaries=True)


    def test_260_raise_exception_when_mapping_bogus_variant(self):
        v = self.hp.parse_hgvs_variant("NM_000059.3:c.7790delAAG")
        with self.assertRaises(HGVSInvalidVariantError):
            self.evm37.c_to_p(v)


    def test_285_partial_palindrome_inversion(self):
        # https://bitbucket.org/biocommons/hgvs/issues/285/
        # Inversion mapping is not palindrome aware
        v = self.hp.parse_hgvs_variant("NM_000088.3:c.589_600inv")
        vn = self.vn.normalize(v)
        self.assertEqual(str(vn), "NM_000088.3:c.590_599inv")


    def test_293_parser_attribute_assignment_error(self):
        # https://bitbucket.org/biocommons/hgvs/issues/293/        
        var = self.hp.parse_hgvs_variant('NG_029146.1:g.6494delG')
        self.vn.normalize(var)  # per issue, should raise error
        

    def test_307_validator_rejects_valid_interval(self):
        # https://bitbucket.org/biocommons/hgvs/issues/307/
        # before fix, raises error; after fix, should pass 
        self.hv.validate(self.hp.parse_hgvs_variant("NM_002858.3:c.1903-573_*1108del"))


    def test_314_parsing_identity_variant(self):
        v = self.hp.parse_hgvs_variant("NM_206933.2:c.6317=")
        self.assertEqual(str(v), "NM_206933.2:c.6317=")

        v = self.hp.parse_hgvs_variant("NM_206933.2:c.6317C=")
        self.assertEqual(str(v), "NM_206933.2:c.6317C=")
                                                   

    def test_322_raise_exception_when_mapping_bogus_variant(self):
        v = self.hp.parse_hgvs_variant("chrX:g.71684476delTGGAGinsAC")
        with self.assertRaises(HGVSInvalidVariantError):
            self.evm37.g_to_c(v, "NM_018486.2")


    def test_324_error_normalizing_simple_inversion(self):
        v = self.hp.parse_hgvs_variant("NM_000535.5:c.1673_1674inv")
        vn = self.vn.normalize(v)
        self.assertEqual(str(vn), "NM_000535.5:c.1673_1674inv")  # no change

        vg = self.evm37.c_to_g(v)
        v2 = self.evm37.g_to_c(vg, tx_ac = v.ac)
        self.assertEqual(str(v), str(v2))  # no change after roundtrip


    def test_326_handle_variants_in_par(self):
        "For a variant in the PAR, allow user user choose preference"
        hgvs_c = "NM_000451.3:c.584G>A"
        var_c = self.hp.parse_hgvs_variant(hgvs_c)

        self.evm37.in_par_assume = None
        with self.assertRaises(HGVSError):
            var_g = self.evm37.c_to_g(var_c)

        self.evm37.in_par_assume = 'X'
        self.assertEqual(self.evm37.c_to_g(var_c).ac, "NC_000023.10")

        self.evm37.in_par_assume = 'Y'
        self.assertEqual(self.evm37.c_to_g(var_c).ac, "NC_000024.9")

        self.evm37.in_par_assume = None


    def test_330_incorrect_end_datum_post_ter(self):
        # https://bitbucket.org/biocommons/hgvs/issues/330/
        # In a variant like NM_004006.2:c.*87_91del, the interval
        # start and end are parsed independently. The * binds to
        # start but not end, causing end to have an incorrect datum.
        v = self.hp.parse_hgvs_variant("NM_004006.2:c.*87_91del")
        self.assertEqual(v.posedit.pos.start.datum, hgvs.location.CDS_END)
        self.assertEqual(v.posedit.pos.end.datum,   hgvs.location.CDS_END)
        

    def test_334_delins_normalization(self):
        # also tests 335 re: inv including sequence (e.g., NC_000009.11:g.36233991_36233992invCA)
        v = self.hp.parse_hgvs_variant("NC_000009.11:g.36233991_36233992delCAinsAC")
        vn = self.vn.normalize(v)
        self.assertEqual(str(vn), "NC_000009.11:g.36233991_36233992delCAinsAC")

        v = self.hp.parse_hgvs_variant("NM_000535.5:c.1673_1674delCCinsGG")
        vn = self.vn.normalize(v)
        self.assertEqual(str(vn), "NM_000535.5:c.1673_1674inv")
        
        v = self.hp.parse_hgvs_variant("NM_000535.5:c.1_3delATGinsCAT")
        vn = self.vn.normalize(v)
        self.assertEqual(str(vn), "NM_000535.5:c.1_3inv")


    def test_346_reject_partial_alignments(self):
        # hgvs-346: verify that alignment data covers full-length transcript
        with self.assertRaises(HGVSDataNotAvailableError):
            hgvs.transcriptmapper.TranscriptMapper(self.hdp,
                                                   tx_ac="NM_001290223.1",
                                                   alt_ac="NC_000010.10",
                                                   alt_aln_method="splign")
