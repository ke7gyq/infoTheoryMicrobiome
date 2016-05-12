# loadData.py

import numpy as np
import pandas as pd
import os , time , re

testTop=os.environ['TESTTOP']
if testTop is None:
    print "Testtop not defined"
    os.exit(-1)


# cols2File.
# Write out column names to a output file
def cols2File ( fileName, df ) :
    cols = df.columns
    f = open (fileName, 'w')
    for c in cols:
        print >> f , "%s" %c 
    f.close()
    
# Read in pickle file. Give a time count for the amount of time needed.
def local_read ( fileName ) :
    print "Starting reading file : %s" % fileName
    start = time.time()
    df =  pd.read_pickle( fileName )
    end = time.time()
    print "Completed read file : %s %d seconds " % (fileName,(end-start))
    return df


# Reaad in the marker data base and change the marker presense fields into
# the bacteria names.
def changeMarkerPresenseColumnNames (mp) :
    fileName = testTop+'/markers2clades_DB.txt'
    d={}
    lines = [line.rstrip ('\n').split('\t') for line in open ( fileName) ]
    for key,value in lines :
        d[key] = value


    columns = mp.columns
    l = []
    for c in columns :
        m = re.search ('gi\|', c )
        if m :
            l.append ( d[m.string] )
        else:
            l.append ( c )
    return l

 


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
    markers2claudes=testTop+'/markers2clades_DB.txt'

    # Do we need to create,or can they be read.
    if False:
        ab = readFile(abundance)             # Relatively small
        ab.to_pickle(testTop+'/abundance.pkl')
        ma = readFile(markerAbundance)        # Really big. Hits swap.
        ma.to_pickle(testTop+'/markerAbundance.pkl')
        mp = readFile(markerPresence)
        mp.to_pickle(testTop+'/markerPresence.pkl')
 
        # Write out some header data info...
        cols2File (testTop+'/abundanceColumns.txt', ab)
        cols2File (testTop+'/markerAbundanceColumns.txt', ma)
        cols2File (testTop+'/markerPresenceColumn.txt', mp)
  

        
    else:
        ab = local_read(testTop+'/abundance.pkl')
        ma = local_read(testTop+'/markerAbundance.pkl')
        mp = local_read(testTop+'/markerPresence.pkl')
       
  
        




        
# xx=queryMatch( 'bodysite', 'stool')

# Find disease = cancer, age > 50
# xx=queryFilter (df, 'age', lambda x: x >=50) 
# yy=queryMatch ( xx, 'disease', 'cancer' )
# yy.age.values


# # Get a column.
# # yy=xx.loc[:,'k__Bacteria|p__Firmicutes|c__Negativicutes|o__Selenomonadales|f__Veillonellaceae|g__Megasphaera|s__Megasphaera_sp_BV3C16_1']


# # xx=df.iloc[:,1].values
# # np.count_nonzero ( xx == 'stool')
