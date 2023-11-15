
# This is a demo for migration (based on Gaussian beam method)
# ChangZhimiao, 2023/11/13
# changzm23@mails.jlu.edu.cn

# Before you use migration program, you should generate the velocity file, shot file, head file
# You can use the forward.sh in demo/forward to generate
# Or you can use the data in data folder

# head file parameters
ntr=100
dtr=12.5
ftr=0
dftr=50

nshot=150
dshot=50
fshot=50

# forward parameters
naper=150
nt=750
nx=737
nz=750
dt=0.004
dx=12.5
dz=4
amax=50
vfile=../../data/mar_vel_737_750_sun.sgy  # input
seifile=../../data/shot.dat # input
headInfo=../../data/headInfo # input
imgfile=mig.dat # output


# migration
mpirun -n 1 ../../bin/migration naper=$naper nt=$nt nx=$nx nz=$nz dt=$dt dx=$dx dz=$dz amax=$amax nshot=$nshot\
     seifile=$seifile vfile=$vfile imgfile=$imgfile headInfo=$headInfo

# SU ximage (be sure you have SU installed)
# ximage < $imgfile n1=$nz
