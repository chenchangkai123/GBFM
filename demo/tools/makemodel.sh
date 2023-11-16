
# This is a demo for segycut 
# ChangZhimiao, 2023/11/16
# changzm23@mails.jlu.edu.cn

# This program is used to cut the segy file
# Make a homogenous model

output=model.dat
nx=100
nz=100
velocity=1500

# segycut
../../bin/makehead output=$output nx=$nx nz=$nz velocity=$velocity




