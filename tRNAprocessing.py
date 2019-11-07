#!/usr/bin/env python3
# Name: Kelsey McDaid (kmcdaid)
# Group Members: None. class FastAreader provided by Instructor David Bernick.

"""
Find that which is unique to and essential in identifying each tRNA sequence, for use in analysis such as shotput sequencing.

Input: strings containing headers and sequences from a fastA file.
Output: for every header, sequence within the fastA file, prints (on separate lines) the header, sequence and each of the short sequences that are unique to and essential for identifying the tRNA sequence.
In the output, alignment of the short sequences within the tRNA sequence is denoted by '.' for every base it is offset from the start of the tRNA sequence.

Powerset = a set of all potential subsets
"""

class Compare():

    """
    Compare the instantiated powerset of an individual sequence to the subsequences of all tRNAs to determine which are unique and essential.
    Return the header, sequence and aligned essential subsequences.
    """
    masterSet = set()

    def __init__(self, sendSet):
        """Initialize the sequence, header and powerset"""
        x,y,z = sendSet #as sendSet is a set
        self.header = x
        self.sequence = y
        self.powerset = z #not frozen

    def findUnique(self):
        """
        For each powerset instantiated, determine the set difference between the powerset and the masterSet.
        Uses the powerset and masterset generated in main() from the fasta file
        Output: a set that contains the set difference of the powerset and the masterset.
        """
        uniqueSeqs = self.powerset - Compare.masterSet
        return uniqueSeqs

    def findEssential(self):
        """
        Determine whether the listing of the unique values are essential as well, building a set of only the essential.
        Input: the set of unique sequences
        Output: a set of the essential sequences, essentialSeqs
        """
        unique = self.findUnique()
        essentialSeqs = set()
        for seq in unique:
            if ((seq[1:] not in unique) and (seq[:-1] not in unique)): #if not present, then there is no smaller substring
                essentialSeqs.add(seq)
        return essentialSeqs

    def essPlace(self):
        """
        Build a dictionary of string:list for each occurence in fullSequence
        Dictionary key-sequence : value-starting positions(list)
        Input: a set from findEssential
        Output: a dictionary of the essential sequence(key) :  its starting position (value)
        """
        essSet = self.findEssential()
        startDict = {}
        for a in essSet:
            startList = [] #the list of the starting positions of each essential sequence
            index = 0
            while index < len(self.sequence): #stackoverflow
                index = self.sequence.find(a, index)
                if index == -1: #if at end
                    break
                startList.append(index)
                index += len(a)
            startDict[a] = startList
        return startDict

        #https://stackoverflow.com/questions/3873361/finding-multiple-occurrences-of-a-string-within-a-string-in-python

    def formatEssential(self):
        """
        For every value, for every starting position within the list, create the alignments for printing in accordance with how far down the sequence each starting point is.
        Input: a dictionary constructed in essPlace
        Outputs: a sorted list of strands with a count of '.' that is equivalent to the number of positions are before it in the full sequence.
        """
        starts = self.essPlace()
        printList = []
        for key, value in starts.items():
            for a in value:
                newPrint = ('.'*a + key) #concenating together a string of length a + key
                printSet = (a, newPrint)
            printList.append(printSet)

        sortList = sorted(printList, key=lambda x: x[0]) #stackoverflow
        return sortList

        #https://stackoverflow.com/questions/10695139/sort-a-list-of-tuples-by-2nd-item-integer-value

    def workflow(self):
        """
        Prints the desired output of header, sequence, aligned substrings.
        Inputs: the list created in formatEssential that orders and aligns the essential sequences based on their placement in the full sequence.
        Output: header, sequence, formatted essential sequences -- all for each strand of tRNA
        """
        formats = self.formatEssential()
        print(self.header)
        print(self.sequence)
        for length, unique in formats:
            print(unique)

#########################################
"""
This module (FastAreader) contains the instructor-provided FastA and related code.
Instructor: David Bernick
"""

import sys
class FastAreader :

    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen (self):
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta (self):

        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence

#########################################
#Main
def main():
    """
    Construct dictionary, list of tuples for use by the class Compare.

    Construction:
    Construct class attribute masterSet, which contains the intersections of powersets.
    Constructs powerset for each seq in the fastA file, which is a set of all contiguous subsets.
    Send tuple of (header, seq, powerset) through each instantiation of Compare().

    Assumes the only extraneous characters are "-" and "_".

    masterSet is a set that contains the intersections of all the powersets
    """
    myReader = FastAreader() # make sure to change this to use stdin
    Compare.masterSet = set() #make a class attribute in Compare
    powerSet = set()
    sendList = []
    newerSet = set()

    for head, seq in myReader.readFasta():
        seq = seq.replace('-','').replace('_','')
        power = set() #uses recursion
        for begin in range(len(seq)): #stackoverflow
            for end in range(begin,len(seq)):
                power.add(seq[begin:end+1])
        powerSet.add(frozenset(power))
        sendList.append((head, seq, power))

    for a in powerSet: #convert back to set
        #print(len(a)) # 22 elements
        for b in powerSet:
            if a != b:
                if a.isdisjoint(b): #case of if nothing in common
                    continue
                else:
                    Compare.masterSet.update(a.intersection(b))

    for x in sendList: #now a list
        sendList = sorted(sendList, key=lambda x: str(x[0])) #sort by header to be in alphabetical order
    #print (sendList)
    for element in sendList:
        Ess = Compare(element)
        Ess.workflow()

main()
#https://stackoverflow.com/questions/41706038/how-do-i-get-all-consecutive-substrings-of-a-string-with-recursion
