# loadData.py
import numpy as np
import pandas as pd
import os , time , re

# Debugging.
import sys, pdb

import gene2Key 

# This file is used as the testing stub for gene2Key.py.
# Handle the initial loads of data, and  then drive various queries.



# Testtop is a environment variable that points to where the data
# is located.


testTop=os.environ['TESTTOP']
if not testTop :
    print "Must define TESTTOP"
    os.exit(-1)


# cols2File.
# Write out column names to a output file
def cols2File ( fileName, df ) :
    with open( fileName,'w') as fp:
        for c in df.columns:
            print >> f ,  "%s" %c 

    
# # Read in pickle file. Give a time count for the amount of time needed.
# def local_read ( fileName ) :
#     print "Starting reading file : %s" % fileName
#     start = time.time()
#     df =  pd.read_pickle( fileName )
#     end = time.time()
#     print "Completed read file : %s %d seconds " % (fileName,(end-start))
#     return df

# # Read in text file. 
# def readFile (fileName):
#     try:
#         print "Starting reading file : %s" % fileName
#         start = time.time()
#         df=pd.read_csv( fileName, sep='\t', header=1, engine='python', 
#                         na_values=['nd'], index_col=0).transpose()
#         end = time.time()
#         print "Completed readfile  : %s %d seconds " % (fileName,(end-start))
#         return df
#     except:
#         print "Could not read %s "% fileName
#         os.exit(-1)

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


if __name__ == "__main__":
    abundance=testTop+'/abundance.txt'
    markerAbundance=testTop+'/marker_abundance.txt'
    markerPresence=testTop+'/marker_presence.txt'

    # Do we need to create,or can they be read.
    if False:
        ab = g2k.readFile(abundance)             # Relatively small
        ab.to_pickle(testTop+'/abundance.pkl')
        ma = g2k.readFile(markerAbundance)        # Really big. Hits swap.
        ma.to_pickle(testTop+'/markerAbundance.pkl')
        mp = g2k.readFile(markerPresence)
        mp.to_pickle(testTop+'/markerPresence.pkl')
 
        # Write out some header data info...
        cols2File (testTop+'/abundanceColumns.txt', ab)
        cols2File (testTop+'/markerAbundanceColumns.txt', ma)
        cols2File (testTop+'/markerPresenceColumn.txt', mp)
  
    else:
        #ab = gene2Key.readPickleFile(testTop+'/abundance.pkl')
        ma = gene2Key.readPickleFile(testTop+'/markerAbundance.pkl')
        #mp = gene2Key.readPickleFile(testTop+'/markerPresence.pkl')
       
  
    # Test the loading / construction of the gene/ key data bases.
    if True:
        g2k = gene2Key.g2k()
        g2k.initialize(ma)
        #g2k.startImport(35)
        #g2k.save ( testTop +'/compressedKeys.pkl')
    #else:
    #    g2k = g2k.k2k.load(  testTop +'/compressedKeys.pkl')

    # Test some queries.
    qCancer = gene2Key.queryMatch( ma, 'disease', 'cancer')
    qAge    = gene2Key.queryFilter( qCancer, 'age', lambda x: x>= 50)
    print "Query Age %s " % qAge.age.values


    # Get some genes associated with a bacteria.
    # Note that first argument to query could be from other queries.
    # 

    gene=ma.columns[300]
    bacteria,_=g2k.gene2bacteria(gene)
    geneList, gNumber = g2k.bacteriaToGenes ( bacteria )
    geneQuery = gene2Key.queryByGeneList ( ma, geneList ) 



# # xx=df.iloc[:,1].values
# # np.count_nonzero ( xx == 'stool')
