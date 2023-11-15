
# This is a demo for LSM
# ChangZhimiao, 2023/11/13
# changzm23@mails.jlu.edu.cn

# Before you use LSM program, you should define the velocity file and head file
# This demo gives a program to generate the head file (program: makehead)
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
vfile=../../data/mar_vel_737_750_sun.sgy # input
reftfile=reflectivity.dat # output
orgfile=shot_0.dat # output

# LSM parameters
R1name=resMig_0.dat
iters=1

../../bin/generate nx=$nx nz=$nz name=$R1name

# generate the head file
../../bin/makehead ntr=$ntr dtr=$dtr ftr=$ftr dftr=$dftr \
    nshot=$nshot dshot=$dshot fshot=$fshot 

# generate the seismic file
../../bin/forward naper=$naper nt=$nt nx=$nx nz=$nz dt=$dt dx=$dx dz=$dz amax=$amax \
    reftfile=$reftfile orgfile=$orgfile vfile=$vfile useOuterreflectivity=0

# copy the seismic file
cp $orgfile orgseis.dat
orgseis=orgseis.dat

# LSM
for ((i=1;i<=$iters;i++))
do
    echo "\033[32mLSM iter: $i\033[0m"

    R=mig_$((i-1)).dat
    seis=shot_$((i-1)).dat
    newseis=shot_${i}.dat

    resMig=resMig_$((i-1)).dat
    newresMig=resMig_${i}.dat
    

    # migration
    mpirun -n 4 ../../bin/migration naper=$naper nt=$nt nx=$nx nz=$nz dt=$dt dx=$dx dz=$dz amax=$amax \
            seifile=$seis vfile=$vfile imgfile=$R
    
    # forward
    ../../bin/forward naper=$naper nt=$nt nx=$nx nz=$nz dt=$dt dx=$dx dz=$dz amax=$amax \
            reftfile=$R orgfile=$newseis vfile=$vfile useOuterreflectivity=1
    
    # LSM
    ../../bin/lsm seisname=$newseis migname=$R resMigname=$resMig newresMIGname=$newresMig

    # forward
    ../../bin/forward naper=$naper nt=$nt nx=$nx nz=$nz dt=$dt dx=$dx dz=$dz amax=$amax \
            reftfile=$newresMig orgfile=$seis vfile=$vfile useOuterreflectivity=1

    # subtract
    ../../bin/subtract nx=$nx nz=$nz originShotname=$orgseis newShotname=$seis subtractDataname=$newseis

done
