"""translate between HGVS and other formats"""

import bioutils.assemblies

import hgvs.normalizer


def _as_interbase(posedit):
    if posedit.edit.type == "ins":
        # ins coordinates (only) exclude left position
        start_i = posedit.pos.start.base
        end_i = posedit.pos.end.base - 1
    else:
        start_i = posedit.pos.start.base - 1
        end_i = posedit.pos.end.base
    return (start_i, end_i)


class Babelfish:
    def __init__(self, hdp, assembly_name):
        self.hdp = hdp
        self.hn = hgvs.normalizer.Normalizer(hdp,
                                             cross_boundaries=False,
                                             shuffle_direction=5,
                                             validate=False)
        self.ac_to_chr_name_map = {
            sr["refseq_ac"]: sr["name"]
            for sr in bioutils.assemblies.get_assembly(assembly_name)["sequences"]}



    def hgvs_to_vcf(self, var_g):
        """**EXPERIMENTAL**

        converts a single hgvs allele to (chr, pos, ref, alt) using
        the given assembly_name. The chr name uses official chromosome
        name (i.e., without a "chr" prefix).

        Returns None for non-variation (e.g., NC_000006.12:g.49949407=)

        """

        if var_g.type != "g":
            raise RuntimeError("Expected g. variant, got {var_g}".format(var_g=var_g))

        vleft = self.hn.normalize(var_g)

        (start_i, end_i) = _as_interbase(vleft.posedit)

        chr = self.ac_to_chr_name_map[vleft.ac]

        typ = vleft.posedit.edit.type

        if typ == "dup":
            start_i -= 1
            alt = self.hdp.seqfetcher.fetch_seq(vleft.ac, start_i, end_i)
            ref = alt[0]
            end_i = start_i
            return (chr, start_i + 1, ref, alt, typ)

        if vleft.posedit.edit.ref == vleft.posedit.edit.alt:
            return None

        alt = vleft.posedit.edit.alt or ""

        if end_i - start_i == 1 and vleft.posedit.length_change == 0:
            # SNVs aren't left anchored
            ref = vleft.posedit.edit.ref

        else:
            # everything else is left-anchored
            start_i -= 1
            ref = self.hdp.seqfetcher.fetch_seq(vleft.ac, start_i, end_i)
            alt = ref[0] + alt

        return (chr, start_i + 1, ref, alt, typ)



if __name__ == "__main__":
    """
      49949___  400       410       420
                  |123456789|123456789|
    NC_000006.12  GACCAGAAAGAAAAATAAAAC

    """

    import hgvs.easy
    import hgvs.normalizer
    from hgvs.extras.babelfish import Babelfish

    babelfish38 = Babelfish(hgvs.easy.hdp, assembly_name="GRCh38")
    hnl = hgvs.normalizer.Normalizer(hgvs.easy.hdp,
                                     cross_boundaries=False, shuffle_direction=5, validate=False)

    def _h2v (h):
        return babelfish38.hgvs_to_vcf(hgvs.easy.parser.parse(h))
    
    def _vp(h):
        v = hgvs.easy.parser.parse(h)
        vl = hnl.normalize(v)
        return (v,vl)

    for h in (
            # Non-variation
            "NC_000006.12:g.49949407=",
            
            # SNV
            "NC_000006.12:g.49949407A>T",

            # delins
            "NC_000006.12:g.49949413_49949414delinsCC",

            # del
            "NC_000006.12:g.49949415del",
            "NC_000006.12:g.49949413del",
            "NC_000006.12:g.49949414del",
            "NC_000006.12:g.49949413_49949414del",

            # ins
            "NC_000006.12:g.49949413_49949414insC",
            "NC_000006.12:g.49949414_49949415insC",
            "NC_000006.12:g.49949414_49949415insCC",

            # ins (dup)
            "NC_000006.12:g.49949413_49949414insA",
            "NC_000006.12:g.49949414_49949415insA",
            "NC_000006.12:g.49949414_49949415insAA"
            ):
        print('assert _h2v("{h}") == {res}'.format(res=str(_h2v(h)), h=h))
