
# This is a demo for segycut 
# ChangZhimiao, 2023/11/16
# changzm23@mails.jlu.edu.cn

# This program is used to cut the segy file
# Convert the segy file to data head file and shot file

# input
input="10-shot.sgy"
dataum=0
isdepth=0
ibmieee=1
# output
output="shot.dat"
elev="elev.dat"
headRecx="headRecx.dat"
headRecy="headRecy.dat"
headShotx="headShotx.dat"
headShoty="headShoty.dat"

# segycut
../../bin/segycut input=$input output=$output dataum=$dataum isdepth=$isdepth ibmieee=$ibmieee \
    elev=$elev headRecx=$headRecx headRecy=$headRecy headShotx=$headShotx headShoty=$headShoty




