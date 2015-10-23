#! /bin/bash

cdref = 6.18533
clref = 0.009401

cd1 = 6.1834911
cd2 = 6.1849822

cl1 = 0.008939246
cl2 = 0.009413395

print 'error cd midpoint {} % '.format(100*abs(cd1-cdref)/cdref)
print 'error cd circle {} %'.format(100*abs(cd2-cdref)/cdref)

print 'error cl midpoint {} %'.format(100*abs(cl1-clref)/clref)
print 'error cl circle {} %'.format(100*abs(cl2-clref)/clref)

