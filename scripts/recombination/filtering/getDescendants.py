#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 2/1/2018
# compareDatabases.py

import sys
import os
import datetime
import numpy
from numpy import random
import gzip
import math
import re

##########################
##### MAIN FUNCTIONS #####
##########################

def getNClosest():
    nodeToDescendants = {}
    nodeToDescendantsPlusOne = {}
    nodeToDescendantsPlusTwo = {}
    nodeToDescendantsPlusThree = {}
    with open('filtering/data/allRelevantNodeNames.txt') as f:
        for line in f:
            #print('('+str(line.strip()[5:])+')')
            #print(str(line.strip()[5:]))
            #exit()
            nodeToDescendants['('+str(line.strip()[5:])+')'] = {}
            nodeToDescendantsPlusOne['('+str(line.strip()[5:])+')'] = {}
            nodeToDescendantsPlusTwo['('+str(line.strip()[5:])+')'] = {}
            nodeToDescendantsPlusThree['('+str(line.strip()[5:])+')'] = {}
    #count = 0
    with open('filtering/data/sample_paths.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0] == 'sample_id':
                mySplit = (splitLine[1]).split('>')
                #print("MY SPLIT: ", mySplit)

                myParent = mySplit[-2].strip().split()[-1]
                #print("PARENT: ", myParent)
                #print("MY PARENT: ", myParent)
                if myParent in nodeToDescendants:
                    #print("PARENT IN nodeToDescendants")
                    #print(myParent)
                    nodeToDescendants[myParent][splitLine[0]] = True
                    #count +=1
                if len(mySplit) >= 3:
                    myParent1 = mySplit[-3].strip().split()[-1]
                    if myParent1 in nodeToDescendantsPlusOne:
                        nodeToDescendantsPlusOne[myParent1][splitLine[0]] = True
                    if len(mySplit) >= 4:
                        myParent2 = mySplit[-4].strip().split()[-1]
                        if myParent2 in nodeToDescendantsPlusTwo:
                            nodeToDescendantsPlusTwo[myParent2][splitLine[0]] = True
                        if len(mySplit) >= 5:
                            myParent3 = mySplit[-5].strip().split()[-1]
                            if myParent3 in nodeToDescendantsPlusThree:
                                nodeToDescendantsPlusTwo[myParent3][splitLine[0]] = True

    #print("REACHED HERE")
    #print("COUNT: ", count)
    #exit()
    myOutString = ''
    allDescendants = ''
    #print("nodeToDescendants LIST LEN: ", len(nodeToDescendants))
    #exit()
    #print(nodeToDescendants.)
    #exit()
    for n in nodeToDescendants:
        #print("FIRST n: ",n)
        #print(n[1:len(n)-1])
        #exit()
        #print(type(n))
        #exit()
        #n = n[6:len(n)-1]
        #n = n[1:len(n)-1]
        myDescendantsList = list(nodeToDescendants[n].keys())
        myKey = 0
        #print("LEN OF myDescendantsList: ", len(myDescendantsList))
        #exit()

        if len(myDescendantsList) > 10:
            myList = random.choice(myDescendantsList, 10, replace=False)
        else:
            myList = myDescendantsList

        #print("LIST LEN: ", len(myList))
        #exit()
        if len(myList) < 1:
            if len(list(nodeToDescendantsPlusOne[n].keys())) > 10:
                myList = random.choice(list(nodeToDescendantsPlusOne[n].keys()), 10, replace=False)
            else:
                myList = list(nodeToDescendantsPlusOne[n].keys())
            myKey = 1

            if len(myList) < 1:
                if len(list(nodeToDescendantsPlusTwo[n].keys())) > 10:
                    myList = random.choice(list(nodeToDescendantsPlusTwo[n].keys()), 10, replace=False)
                else:
                    myList = list(nodeToDescendantsPlusTwo[n].keys())
                myKey = 2


        myOutString += n+'\t'+joinerC(myList)+'\t'+str(myKey)+'\n'
        allDescendants += joinerN(myList)+'\n'
    open('filtering/data/allRelevantNodesToDescendants.txt','w').write(myOutString)
    open('filtering/data/allDescendants.txt','w').write(allDescendants)


##########################
#### HELPER FUNCTIONS ####
##########################

def getPos(myInds, intLineNumToPos):
    myReturn = []
    for k in myInds:
        myReturn.append(intLineNumToPos[k])
    return(myReturn)

def joiner(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return('\t'.join(newList))

def joinerC(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return(','.join(newList))

def joinerN(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return('\n'.join(newList))

#########################
##### FUNCTION CALL #####
#########################

def main():
    getNClosest()

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit

