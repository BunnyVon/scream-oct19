#!/bin/bash

# Set run options
#=============================================

resolution=AT_ne1024pg2_r0125_oRRS18to6v3
compset=ICRUELM
checkout_date=20220926  #the date you *checked out* the code
branch=ad0e203          #actual git hash of branch to check out
run_descriptor=SCREAMv0 #will be SCREAMv0 for production run
repo=scream
machine=pm-cpu
compiler=gnu
project=m4513
stop_option="nmonths"
stop_n="36"
rest_n="1"
walltime="6:00:00"
queue="regular"
debug_compile='FALSE'
date_string=`date +"%Y%m%d%H%M"`

# Setup processor layout
nnodes_atm=40
nnodes_ocn=40
nthreads=1
mpi_tasks_per_node=128
ntasks_atm=$(expr ${nnodes_atm} \* ${mpi_tasks_per_node})
ntasks_ocn=$(expr ${nnodes_ocn} \* ${mpi_tasks_per_node})
total_tasks_per_node=$(expr ${mpi_tasks_per_node} \* ${nthreads})
if [ ${ntasks_ocn} -ne ${ntasks_atm} ]; then
    nnodes=$(expr ${nnodes_atm} + ${nnodes_ocn})
else
    nnodes=${nnodes_atm}
fi
pelayout=${nnodes}x${mpi_tasks_per_node}x${nthreads}

stridc=16
npcpl=$(expr ${ntasks_atm} \/ ${stridc})
echo "npcpl=${npcpl}"

# Who to send email updates on run status to
email_address="sourav.taraphdar@pnnl.gov"

# Set flags specific to running this script
do_download=false
do_newcase=true
do_setup=true
do_build=true
do_submit=true

case_name=SCREAMv0.${compset}.${resolution}.${date_string}

# Set paths
#code_root=/global/cfs/cdirs/m1867/souravt/scream-11apr23
code_root=/global/cfs/projectdirs/m1867/souravt/scream-oct19
case_root=$SCRATCH/E3SM/cases/${case_name}

######################################################################################
### USERS PROBABLY DON'T NEED TO CHANGE ANYTHING BELOW HERE EXCEPT user_nl_* FILES ###
######################################################################################

# Make directories created by this script world-readable:
#=============================================
umask 022

# Download code
#=============================================
if [ "${do_download}" == "true" ]; then

    echo "Cloning repository repo = $repo into branch = $branch under code_root = $code_root"
    cdir=`pwd`
    mkdir -p $code_root/

    # This will put repository, with all code, in directory $tag_name
    git clone git@github.com:E3SM-Project/${repo}.git $code_root
    
    # Setup git hooks
    rm -rf $code_root/.git/hooks
    git clone git@github.com:E3SM-Project/E3SM-Hooks.git $code_root/.git/hooks
    cd $code_root
    git config commit.template $code_root/.git/hooks/commit.template

    # Bring in all submodule components
    git submodule update --init --recursive

    # Check out desired branch
    git checkout ${branch}
    cd ${cdir}

fi

# Create case
#=============================================
if [ "${do_newcase}" == "true" ]; then
	    ${code_root}/cime/scripts/create_newcase \
        --case ${case_root} --compset ${compset} --res ${resolution} --machine ${machine} --compiler ${compiler} \
        --project ${project} --queue ${queue} --walltime ${walltime}
fi

# Copy this script to case directory
#=============================================
cp -v `basename $0` ${case_root}/

# Setup
#=============================================
if [ "${do_setup}" == "true" ]; then
    cd ${case_root}

    # Set run length
    ./xmlchange STOP_OPTION=${stop_option},STOP_N=${stop_n}
    ./xmlchange REST_N=${rest_n}

    # Set processor layout
    if [ ${ntasks_ocn} -ne ${ntasks_atm} ]; then
        ./xmlchange NTASKS=${ntasks_ocn}
        ./xmlchange NTASKS_ATM=${ntasks_atm}
        ./xmlchange ROOTPE_ATM=${ntasks_ocn}
    else
        ./xmlchange NTASKS=${ntasks_atm}
    fi
    ./xmlchange NTHRDS_ATM=${nthreads}
    ./xmlchange MAX_MPITASKS_PER_NODE=${mpi_tasks_per_node}
    ./xmlchange NTASKS_CPL=${npcpl} # change tasks in CPL if using strid
    ./xmlchange PSTRID_CPL=${stridc}
    ./xmlchange NTHRDS=1 # set all to 1 first
    ./xmlchange MAX_TASKS_PER_NODE=${total_tasks_per_node}

    # Flag for debug compile
    ./xmlchange --id DEBUG --val ${debug_compile}
    
    ./xmlchange PIO_NETCDF_FORMAT="64bit_data"
    ./xmlchange PIO_VERSION="2"
    ./xmlchange PIO_BUFFER_SIZE_LIMIT=134217728

    ./xmlchange DATM_MODE=CLMCRUNCEPv7
    ./xmlchange DATM_CLMNCEP_YR_START=2017,DATM_CLMNCEP_YR_END=2017,DATM_CLMNCEP_YR_ALIGN=2017
    ./xmlchange RUN_STARTDATE=2014-08-24
    ./xmlchange START_TOD=0

    cat <<EOF >> user_nl_elm
    fsurdat = '/pscratch/sd/s/souravt/RRM_setup_tutorial/RRM_generation_inputs/surfdata_ATL_ne128x8pg2_simyr2010_c231026.nc'
    finidat = '/global/cfs/cdirs/m1867/souravt/r0125.AT.elm.r.0051-01-01-00000.nc'
EOF

    # Finally, run setup
    ./case.setup

    # The run location is determined in the bowels of CIME
    # Symlink to that location from user-chosen $case_root (=current dir)
    ln -s `./xmlquery -value RUNDIR` run


    # This disables the logic that sets tprof_n and tprof_options internally.
    #./xmlchange --file env_run.xml TPROF_TOTAL=-1
    #echo "tprof_n = 1" >> user_nl_cpl
    #echo "tprof_option = 'nsteps'" >> user_nl_cpl
fi

# Build
#=============================================
if [ "${do_build}" == "true" ]; then
    cd ${case_root}
    ./case.build
fi

# change atmospheric forcing path for ERA5
cd ${case_root}/CaseDocs
cp datm.streams.txt.CLMCRUNCEPv7.Precip ../user_datm.streams.txt.CLMCRUNCEPv7.Precip
cp datm.streams.txt.CLMCRUNCEPv7.Solar ../user_datm.streams.txt.CLMCRUNCEPv7.Solar
cp datm.streams.txt.CLMCRUNCEPv7.TPQW ../user_datm.streams.txt.CLMCRUNCEPv7.TPQW

cd ${case_root}
sed -i 's#/global/cfs/cdirs/e3sm/inputdata/atm/datm7/atm_forcing.datm7.cruncep_qianFill.0.5d.v7.c160715/#/global/cfs/projectdirs/m1867/souravt/atmospheric_forcing/#g' user_datm.streams.txt.CLMCRUNCEPv7*

./preview_namelists
# Run
#=============================================
if [ "${do_submit}" == "true" ]; then
    cd ${case_root}

    # Set file striping on run dir for writing large files

#    ./xmlchange  -file env_batch.xml -id JOB_WALLCLOCK_TIME -val '0:30:00'
#    ./xmlchange  -file env_batch.xml -id JOB_QUEUE          -val 'debub'
    ./case.submit --batch-args="--mail-type=ALL --mail-user=${email_address}"
fi

echo "Done working in ${case_root}"
