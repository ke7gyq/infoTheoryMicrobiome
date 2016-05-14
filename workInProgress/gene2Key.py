# gene2Key.py
# Manage associations between genes ('gi|...), bacteria name (k_xxx|p_xxx ...) and
# key numbers.
#

import os,  sys, re
import numpy as np
import multiprocessing

import pickle

import pdb

class g2k:
    def __init__ (self):
        self.testTop = os.environ['TESTTOP']
        self.dbName = '/markers2clades_DB.txt'
        self.initialized = False
        self.byGene = None
        self.byBacteria = None
        self.byIndex = None

        self.dataFrame = None
        self.giColumns = None
        self.nGiColums = None

        # After data has been imported.
        self.headers = None
        self.genes = None

        
    # Create lookup tables where we can associate gene name, bacteria names
    # and index (integer) to a column

    # Attach a dataframe that contains genes to this object.
    # Initialize the data frame fields to iterate over to form values.
    #
    def initialize ( self , dataFrame ) :
        if self.initialized: return
        self.initialized = 1
        fName = self.testTop + self.dbName
        lines = [ line.rstrip('\n').split('\t') for line in open( fName )]
        self.byGene  = { v[0]: (v[1],idx) for idx, v in enumerate(lines)}  
        #self.byBacteria  = { v[1]: (v[0],idx) for idx, v in enumerate(lines)}

        pdb.set_trace()
        bacteria = np.unique([l[1] for l in lines ])
        self.byBacteria = dict().fromkeys(bacteria, [])
        for idx, v in enumerate ( lines ):
            self.byBacteria[v[1]].append ( (v[0],idx)) 

        
        self.byIndex = {idx: (v[0], v[1]) for idx, v in enumerate(lines)}

        self.dataFrame = dataFrame
        sKey = re.compile ('gi\|.*')
        self.giColumns  = [ c for c in dataFrame.columns if re.match(sKey,c) ]
        self.nGiColumns  = [ c for c in dataFrame.columns if not re.match(sKey,c)]
        
        
    # makeCompactGeneArray
    # Given a row into the data frame, copy over attribute information ( metadata )
    # and find the genes in the data.
    
    def makeCompactGeneArray ( self, row ) :
        attributes = np.array([ row.loc[a] for a in self.nGiColumns ])
        # Do all columns so that all values can be cast and compared at once.
        genes = np.array ([ row.loc[a] for  a in self.giColumns]).astype('float')
        gValue = np.array( [ (idx, v ) for idx,v in enumerate ( genes ) if v > 0.0] )
        return (attributes, gValue)

    # Self.dataframe must be set.
    # Note that this will be spawned as a multiprocessor command.
    def importData ( self , start, end , processNumber, return_dict ):
        if not self.initialized:
            raise ValueError ("Initialize failure, Data frame not set ")
        hdrData , geneData = ([],[] )
        for  rowNumber in range ( start, end+1 ):
            print "Processor %d Translating row %d " % (processNumber, rowNumber)  
            hdr, gene = self.makeCompactGeneArray( self.dataFrame.ix[rowNumber] )
            hdrData.append( hdr)
            geneData.append(gene)
        return_dict [ processNumber] = (hdrData, geneData)
 
    # Start multithreading on gene data.
    # 
    def startImport ( self, nProcess = 10 , n=None ) :
        if not n :
            n = self.dataFrame.shape[0]

        stride , start, end , pn , jobs  = (n // (nProcess-1), 0 , 0, -1 , [])
        manager = multiprocessing.Manager()
        rd = manager.dict()
        start = time.time()
        while start < n  :
            pn += 1 
            end = start + stride if start+stride < n else n-1
            p = multiprocessing.Process(target= g2k.importData, args = ( self, start , end, pn, rd ) )
            start = start+stride
            jobs.append(p)
            p.start()

        for p in jobs:
            p.join()

        hdr, gene = ( [], [] )
        for k in rd.keys():
            hdr.extend(rd[k][0])
            gene.extend(rd[k][1])

        # create data frame for header data.
        self.headers = pd.DataFrame ( hdr, columns=self.nGiColumns)        
        self.genes = gene
        end = time.time()        
        print "ImportData took %d seconds " % (end-start)
               
 

    # We don't want to save off the origional data. Just save the
    # compression artifacts.
    def save( self, fileName ):
        self.dataFrame = None
        self.initialized = False
        start = time.time()
        pickle.dump ( self, open ( fileName, 'wb'))
        end = time.time()
        print "Save took %d seconds" % (end-start)
    # Factory function.
    # Create an object. Load in the compression atrifacts.
    @staticmethod
    def load ( fileName ):
        start=time.time()
        self = pickle.load ( open ( fileName, 'rb'))
        end = time.time ()
        print "Load took %d seconds" % (end-start)
        return self

# Testing.

if False:
    # %pdb
    # To Reload:
    del sys.modules['gene2Key']

    import gene2Key

    # Read data object.
    if False:
        g2k = gene2Key.g2k()
        g2k.initialize ( ma )
        g2k.startImport ( 35 )

    # Try to save class.
    if False:
        g2k.save ( g2k.testTop +'/compressedKeys.pkl')

    # Reload class.
    if False:
        g2k = gene2Key.g2k.load(g2k.testTop +'/compressedKeys.pkl')
