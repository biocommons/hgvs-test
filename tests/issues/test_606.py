import unittest

from hgvs.exceptions import HGVSInvalidIntervalError
import hgvs


class TestIssue606(unittest.TestCase):
    def test_606(self):
        """
            https://github.com/biocommons/hgvs/issues/606

            Occasionally, an IndexError is thrown by the _get_altered_sequence method. This seems to occur
            when there is either inconsistent data for a transcript in UTA or an invalid variant input.
            There could be other root causes, but these are two that we have observed and are illustrated in the
            examples below.

            Example 1 (NC_000001.10:g.16890441C>G):
            Transcript NM_017940.4 has inconsistent data. 25 exons in the 'transcript' alignment to itself,
            but 28 exons in the 'splign' alignment to the reference sequence.
            Ends up with an index error in _get_altered_sequence.

            Example 2 (NM_001291722.1:c.283-3C>T):
            This is invalid input -- according to the transcript in UTA, exon 5 starts at position 286 of the transcript,
            not position 283. Ends up with an index error in _get_altered_sequence.
        """

        from hgvs.easy import am37, parser

        # Example 1
        hgvs.global_config.mapping.strict_bounds = False
        hgvs_g = 'NC_000001.10:g.16890441C>G'
        var_g = parser.parse_hgvs_variant(hgvs_g)
        with self.assertRaises(HGVSInvalidIntervalError):
            var_c = am37.g_to_c(var_g, 'NM_017940.4')
        hgvs.global_config.mapping.strict_bounds = True

        # Example 2
        hgvs_c = 'NM_001291722.1:c.283-3C>T'
        var_c = parser.parse_hgvs_variant(hgvs_c)
        with self.assertRaises(HGVSInvalidIntervalError):
            var_g = am37.c_to_g(var_c)
