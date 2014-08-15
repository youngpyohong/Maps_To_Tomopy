import sys

try:
      x=sys.argv[1]
      print x
except IndexError:
      print "no"
