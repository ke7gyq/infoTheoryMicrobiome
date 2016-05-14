# gene2Key.py
# Manage associations between genes ('gi|...), bacteria name (k_xxx|p_xxx ...) and
# key numbers.
#

import os,  sys, re
import numpy as np, pandas as pd
import multiprocessing, pickle, time
import pdb


# Note that this was origionally in 'loadData'
# Return a data frame where column name matches a query.
#
def queryMatch ( df, column , query ) :
    try:
        idxCol = df.columns==column
        idxRow = (df.loc[:,idxCol].values == query).flatten()
        return df.loc[idxRow, :]
    except:
        print ("Query of %s on %s failed" %( column, query))
        return None

# Allow user to pass in rules for query.
def queryFilter ( df, column, rule ):
    try:
        idxCol= df.columns==column
        values = pd.to_numeric( df.loc[:,idxCol].values.flatten() ,errors='coherce')
        idxRow = np.array(map( rule, values ))
        return df.loc[idxRow,:]
    except:
        print ("Query Lambda failed column:  %s " %( column))
        return None

# Return data frame queried by list of genes. 
# Normally get the gene names by g2k.backteriaToGenes ( bacteria)
def queryByGeneList ( df, geneList ):
    idxRow = np.zeros ( df.shape[0] , dtype=bool)
    for g in geneList:
        selected = np.array(df.loc[:,g].values).astype('float') > 0.0
        idxRow = [ a or  b for a,b in zip(idxRow, selected)]
    return df.loc[idxRow,:]

    



# Read in pickle file. Give a time count for the amount of time needed.
def readPickleFile ( fileName ) :
    print "Starting reading file : %s" % fileName
    start = time.time()
    df =  pd.read_pickle( fileName )
    end = time.time()
    print "Completed read file : %s %d seconds " % (fileName,(end-start))
    return df



# Read in text file. 
def readFile (fileName):
    try:
        print "Starting reading file : %s" % fileName
        start = time.time()
        df=pd.read_csv( fileName, sep='\t', header=1, engine='python', 
                        na_values=['nd'], index_col=0).transpose()
        end = time.time()
        print "Completed readfile  : %s %d seconds " % (fileName,(end-start))
        return df
    except:
        print "Could not read %s "% fileName
        os.exit(-1)

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
    # Create a list of genes that can be indexed via bacteria.
    # 
    def initialize ( self , dataFrame ) :
        if self.initialized: return
        self.initialized = 1
        fName = self.testTop + self.dbName
        lines = [ line.rstrip('\n').split('\t') for line in open( fName )]
        self.byGene  = { v[0]: (v[1],idx) for idx, v in enumerate(lines)}  

        bacteria = np.unique([l[1] for l in lines ])
        self.byBacteria = dict()
        for b in bacteria:
            self.byBacteria[b] = list()

        for idx, v in enumerate ( lines ):
            self.byBacteria[v[1]].append ( (v[0],idx)) 
        
       
        self.byIndex = {idx: (v[0], v[1]) for idx, v in enumerate(lines)}

        self.dataFrame = dataFrame
        sKey = re.compile ('gi\|.*')
        self.giColumns  = [ c for c in dataFrame.columns if re.match(sKey,c) ]
        self.nGiColumns  = [ c for c in dataFrame.columns if not re.match(sKey,c)]
        
    # Given a bacteria, find list of genes that correspond to the bacteria.
    # Return 2 lists, first list is the gene name, second is the gene number.
    #
    def bacteriaToGenes(self, bacteria ):
        gb = np.array(self.byBacteria[bacteria])
        return (gb[:,0], gb[:,1].astype(int))
        
    # given a gene number ( from bacteria to gene) return the bacteria and 
    # gene name.
    def getGeneFromGeneNumber( self, geneNumber):
        g , b = self.byIndex [ geneNumber] 
        return (g,b)
        
    # Return the bacteria that is accociated with this gene.
    # Note that this returns both the bacteria and the gene number 
    # associated with this bacteria.
    def gene2bacteria(self, gene ):
        return self.byGene[gene]


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

    # Run code snippet in debugger.
    #pdb.run ('gene2Key.queryByGeneList(ma,geneList )')
