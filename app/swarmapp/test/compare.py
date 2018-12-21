import csv
import argparse

epsilon = 1.e-12

parser = argparse.ArgumentParser(description='compare two csv files for agreement')
parser.add_argument('files', nargs=2, help='files to compare')
parser.add_argument('distribution', help='distribution: 0 for random, 1 for grid')

args = parser.parse_args()

def compare(goldfile, testfile, distribution):
  goldcoords = []
  testcoords = []

  goldvals = []
  testvals = []
  
  goldreader = csv.reader(open(goldfile,"r"))
  testreader = csv.reader(open(testfile,"r"))
  goldreader.next()
  testreader.next()
  
  for line in goldreader:
    goldcoords.append([float(line[0]),float(line[1])])
    goldvals.append(float(line[2]))
  
  for line in testreader:
    testcoords.append([float(line[0]),float(line[1])])
    testvals.append(float(line[2]))
 
  if len(goldvals) != len(testvals):
    print "length of gold and test data not equal"
  
  isgood = True

  if int(distribution)==1:
    bothcoords = zip(goldcoords, testcoords, range(len(testcoords)))
    for x in bothcoords:
      isgood2 = abs(x[0][0]-x[1][0])<epsilon
      isgood = isgood and isgood2
      isgood3 = abs(x[0][1]-x[1][1])<epsilon
      isgood = isgood and isgood3
      if (not isgood2) or (not isgood3):
        print "disagree: ", x[0], x[1], "at line ", x[2]+2
  
  bothvals = zip(goldvals, testvals, range(len(testvals)))
  isgood = True
  for x in bothvals:
    isgood2 =  abs(x[0]-x[1])<epsilon
    if not isgood2:
      print "disagree: ", x[0], x[1], ' at line ', x[2]+2
    isgood = isgood and isgood2

  return isgood

print 'comparing ', args.files[0], ' and ', args.files[1]
result = compare(args.files[0], args.files[1], args.distribution)
if result: 
    print 'files agree'
exit(not result)


