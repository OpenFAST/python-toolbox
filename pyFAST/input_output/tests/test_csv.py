import unittest
import glob
import pyFAST
from pyFAST.input_output import CSVFile
import os
import numpy as np
MyDir=os.path.join(os.path.dirname(__file__),'example_files')

class Test(unittest.TestCase):

    def test_001_read_all(self, DEBUG=True):
        nError=0
        for f in glob.glob(os.path.join(MyDir,'CSV*.*')):
            if os.path.splitext(f)[-1] in ['.py','.pyc'] or f.find('_TMP')>0:
                continue
            try:
                obj = CSVFile(f)
                s=type(obj).__name__.replace('file','')[:20]
                if DEBUG:
                    print('[ OK ] {:30s}\t{:20s}'.format(os.path.basename(f)[:30],s))
            except:
                nError += 1
                if DEBUG:
                    print('[FAIL] {:30s}\tException occurred'.format(os.path.basename(f)[:30]))
                raise 
        if nError>0:
            raise Exception('Some tests failed')

    def DF(self,FN):
        """ Reads a file and return a dataframe """ 
        return CSVFile(os.path.join(MyDir,FN)).toDataFrame()
 
    def test_CSV(self):
        self.assertEqual(self.DF('CSVAutoCommentChar.txt').shape,(11,6))
 
        DF=self.DF('CSVColInHeader.csv')
        self.assertTrue(all(DF.columns.values==['ColA','ColB','ColC']))
        self.assertEqual(DF.shape,(2,3))
 
        DF=self.DF('CSVColInHeader2.csv')
        self.assertTrue(all(DF.columns.values==['ColA','ColB','ColC']))
        self.assertEqual(DF.shape,(2,3))
 
        DF=self.DF('CSVColInHeader3.csv')
        self.assertTrue(all(DF.columns.values==['ColA','ColB','ColC']))
        self.assertEqual(DF.shape,(2,3))

        #DF=self.DF('CSVComma_UTF16.csv') # TODO encoding
        #self.assertEqual(DF.shape,(4,3))
 
        self.assertEqual(self.DF('CSVComma.csv').shape,(4,2))
        self.assertEqual(self.DF('CSVDateNaN.csv').shape,(11,2))
        self.assertEqual(self.DF('CSVNoHeader.csv').shape,(4,2))
        self.assertEqual(self.DF('CSVSemi.csv').shape,(3,2))
        self.assertEqual(self.DF('CSVSpace_ExtraCol.csv').shape,(5,4))
        self.assertEqual(self.DF('CSVTab.csv').shape,(5,2))
 
        DF = self.DF('CSVTwoLinesHeaders.txt')
        self.assertEqual(DF.columns.values[-1],'GenTq_(kN m)')
        self.assertEqual(DF.shape,(9,6))

if __name__ == '__main__':
    #Test().test_CSV()
    unittest.main()
