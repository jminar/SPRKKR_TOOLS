#!/bin/csh -f

magnet-split-kkr-potential.x   sample-input.pot   sample-output.pot  <<EOF
   10   +0.2
   11   +0.05
   13   +-0.1
   16   *3.0
   22   *-1.0
EOF

exit




# Description:
# Potential in file sample-input.pot will be manipulated in the following way:
#   a) BT part of the potential for IT=10 will be increased by 0.2 
#   b) BT part of the potential for IT=11 will be increased by 0.05
#   c) BT part of the potential for IT=13 will be decreased by 0.1 (it is negative!)
#   d) BT part of the potential for IT=16 will be multiplied by 3
#   e) BT part of the potential for IT=22 will be multiplied by -1
#
# All other potentials will be left intact.  The new potential is in file
#    sample-output.pot

