# gatherData.py
# Start looking at genomic data .
#

import numpy as np
import pandas as pd
from pandas.io import sql
import os
import re


import MySQLdb

try:
    testTop = os.environ ['TESTTOP']
except:
    testTop= '../data'


class BNames:
    def  __init__(self):
        self.allCols = None
        self.bacteria = None
        self.bDict    = None
       

    def parseColumns ( self, columns ) :
        self.allCols = columns
        # Filter list of colums to only include k__..
        self.bacteria = sorted( [ k  for k in columns if k.startswith('k__') ])

        # Create dictionary tree with column names.
        self.bDict = dict()
        
        for b in self.bacteria:
            d = self.bDict;
            fields = b.split('|')
            for f in fields [:-1]:
                d = d[f]
            d[fields[-1]] = {}        
        pass






    
# readTextData.
# Reads raw data from .txt file and stroes it into a panda data structure
# Each line ( terminated by \n) contains all data for the field
# Field name is first entry on the line.
# Number of fields is the number of lines in the data set.
#


def readTextData( inFileName, outFileName ):
    dData = {}    
    inFile = open (inFileName, 'r')
    for line in inFile:
        line,_ = line.rsplit('\n')
        data = line.split('\t')
        fieldName = data[0]
        fieldData = data[1:]
        #fieldData = data[1:10]
        dData [fieldName] =fieldData    
    pandaData = pd.DataFrame (dData)
    pandaData.to_pickle( outFileName )
    return pandaData



def fixupCols ( inCols) :
    l =[]
    cnt =0
    for n in inCols :
        lbl = "Lbl%03d" % cnt
        cnt +=1;
        l.append (lbl)

    
    # for n in inCols :
    #     v = ''
    #     if n[0] == '#':
    #         v=n[1:]
    #     elif  n[0] == '%':
    #         v=n[1:]
    #     else: v=n
        
    #     if v[0] == '_':
    #         v=v[1:]
    #     v=v.replace ( '(', "1x1")
    #     v=v.replace ( ')', "1x1")
    #     v=v.replace ( '_', '1' )
    #     v=v.replace ( '?', '1' )
    #     v=v.replace ('-', '1')
    #     v=v.replace ('/', '1')
    #     v=v.replace ('#', '1')
    #     v=v.replace ('|', '2')
    #     v=v.replace ('<', '3')
    #     v=v.replace ('>', '4')
    #     l.append(v)
    return l



if __name__ == '__main__':
    inFileName=testTop+'/abundance.txt'
    outFileName= testTop+'/abundance.pkl'
    # pData = readTextData( inFileName, outFileName)

    pData = pd.read_pickle( outFileName )

    # Rename columns to someting sane
    cols = pData.columns.values
    bNames = BNames ()
    bNames.parseColumns( pData.columns )

    # After substitution.
    pData.columns = fixupCols(cols);

    # Note that this doesn't work if using anaconda...
    #
    db= MySQLdb.connect( host='localhost',
            user='doug', passwd='ding-ding', db='testDb')

    pData.to_sql( con=db, name='abundance',
                     if_exists='replace', flavor='mysql')


