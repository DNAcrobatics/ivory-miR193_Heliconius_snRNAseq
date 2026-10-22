"""Tiny synthetic fixture only; no biological analysis."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import pysam
from scsplice import cigar_segments, ribbon_bins

ROOT = Path(__file__).resolve().parent


class TestSplice(unittest.TestCase):
    def test_cigar_coordinates(self):
        self.assertEqual(cigar_segments(101, '5S10M100N5M2D5M3I5M'),
                         ([(100,110),(210,215),(217,222),(222,227)], [(110,210)]))

    def test_historical_bins(self):
        reads = [{'pos':p} for p in [10,10,20,40,50,50]]
        self.assertEqual(ribbon_bins(reads,'legacy',2,500,1,100),
                         [(10,510,3),(40,540,3)])
        self.assertEqual(ribbon_bins(reads,'fixed',2,25,1,100),
                         [(51,76,0),(76,101,0),(1,26,3),(26,51,3)])

    def test_bam_to_figures(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            bam_path = d/'test.bam'
            with pysam.AlignmentFile(bam_path,'wb',header={
                    'HD': {'VN':'1.6','SO':'coordinate'},
                    'SQ':[{'SN':'chrTest','LN':1000}]}) as b:
                for i, (bc, reverse, cigar) in enumerate([
                        ('AAA-1',False,'10M100N10M')]*6 +
                        [('BBB-1',False,'20M')]*6 +
                        [('AAA-1',True,'10M100N10M')]*6 +
                        [(None,False,'10M100N10M')]):
                    r=pysam.AlignedSegment()
                    r.query_name=f'r{i}'
                    r.query_sequence='A'*20
                    r.flag=16 if reverse else 0
                    r.reference_id=0
                    r.reference_start=100+i
                    r.mapping_quality=60
                    r.cigarstring=cigar
                    if bc: r.set_tag('CB',bc)
                    b.write(r)
            pysam.index(str(bam_path))
            meta=d/'cells.tsv'
            # Reused barcode in another sample must never leak into the join.
            meta.write_text('sample\tbarcode\tcelltype\nS1\tAAA\tScale\nS1\tBBB\tSocket\nS2\tAAA\tWRONG\n')
            def run(*args, ok=True):
                r=subprocess.run([sys.executable,str(ROOT/'scsplice.py'),*map(str,args)],
                                 text=True,capture_output=True)
                self.assertEqual(r.returncode,0,msg=r.stderr) if ok else self.assertNotEqual(r.returncode,0)
                return r
            run('extract','--bam',bam_path,'--metadata',meta,'--sample','S1',
                '--region','chrTest:90-300','--strand','forward','--spliced-only',
                '--strip-barcode-suffix','--output',d/'out')
            with (d/'out.reads.tsv').open() as f:
                rs=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(len(rs),6)
            self.assertEqual({r['celltype'] for r in rs},{'Scale'})
            self.assertEqual(json.loads((d/'out.extract.json').read_text())['counts']['missing_barcode_tag'],1)
            for method in ['legacy','fixed']:
                run('plot','--reads',d/'out.reads.tsv','--region','chrTest:90-300',
                    '--binning',method,'--log-limits','0','2','--output',d/method)
                for suffix in ['.ribbons.pdf','.cells.pdf']:
                    self.assertTrue((d/(method+suffix)).read_bytes().startswith(b'%PDF'))
            run('plot','--reads',d/'out.reads.tsv','--region','chrTest:90-300',
                '--output',d/'legacy',ok=False)
            # Exact matching deliberately fails for bare metadata vs suffixed BAM tags.
            run('extract','--bam',bam_path,'--metadata',meta,'--sample','S1',
                '--region','chrTest:90-300','--output',d/'mismatch',ok=False)
            # Save a small illustration only when explicitly requested by the caller.
            import os
            if os.environ.get('SCSPLICE_TEST_PREVIEW'):
                import shutil
                for suffix in ['.ribbons.pdf','.cells.pdf']:
                    shutil.copy(d/('legacy'+suffix), os.environ['SCSPLICE_TEST_PREVIEW']+suffix)


if __name__ == '__main__':
    unittest.main()
