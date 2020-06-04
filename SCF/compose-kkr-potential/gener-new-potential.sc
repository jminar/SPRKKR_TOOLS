#!/bin/csh -f


# Take the old potential  orig-template.pot  and replace potentials for some atomic types
#    by potentials picked from file  potential-to-be-picked.pot
#
#    E.g., potential for IT=4, defined on a radial mesh with RWS=1.60409206025255
#          is going to be put at potential for type IT=5
#          and interpolated into  radial mesh with RWS=1.63973253839201



./compose-kkr-potential.x   orig-template.pot   newly-assembled-potential.pot   1.29487625186768  <<EOF
#
potential-to-be-picked.pot      1    1.90365576669073    1     1.90365576669073
potential-to-be-picked.pot      2    1.63309959129205    2     1.63309959129205
potential-to-be-picked.pot      3    2.11800255848969    3     2.11800255848969
potential-to-be-picked.pot      4    1.60409206025255    4     1.62228801607648
potential-to-be-picked.pot      4    1.60409206025255    5     1.63973253839201
potential-to-be-picked.pot      5    2.31030797838875    6     2.31030797838875
potential-to-be-picked.pot      6    1.99056239209031    7     1.99056239209031
potential-to-be-picked.pot      7    2.63000000000000    8     2.63000000000000
#
EOF
#

exit
