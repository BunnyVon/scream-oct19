#!/bin/bash

# Set run options
#=============================================

resolution=WNP_ne1024pg2_r0125_oRRS18to6v3
compset=F2010-SCREAM-HR-DYAMOND2
checkout_date=20220926  #the date you *checked out* the code
branch=ad0e203          #actual git hash of branch to check out
run_descriptor=SCREAMv0 #will be SCREAMv0 for production run
repo=scream
machine=pm-cpu
compiler=gnu
project=m4393
stop_option="ndays"
stop_n="8"
rest_n="1"
walltime="08:00:00"
queue="regular"
debug_compile='FALSE'
date_string=`date +"%Y%m%d%H%M"`

# Setup processor layout
nnodes_atm=80
nnodes_ocn=80
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
    
    # Set PIO format, use PIO version 2, and increase PIO buffer size 
    ./xmlchange PIO_NETCDF_FORMAT="64bit_data"
    ./xmlchange PIO_VERSION="2"
    ./xmlchange PIO_BUFFER_SIZE_LIMIT=134217728
    ./xmlchange ATM_NCPL=864

    ./xmlchange -file env_build.xml -id CAM_DYCORE -val se
    ./xmlchange EPS_AGRID=1e-9

    ./xmlchange CAM_TARGET=theta-l
    ./xmlchange SSTICE_DATA_FILENAME=/pscratch/sd/s/souravt/2023-08-11_17-07-43/HICCUP.sst_noaa.2009-01-03.nc
    ./xmlchange SSTICE_YEAR_START=2009,SSTICE_YEAR_END=2009,SSTICE_YEAR_ALIGN=2009
    ./xmlchange GLC_AVG_PERIOD=glc_coupling_period
    ./xmlchange DOUT_S=FALSE
    ./xmlchange RUN_STARTDATE=2009-01-03
    ./xmlchange START_TOD=0

    # Edit CAM namelist to set dycore options for new grid
    cat <<EOF >> user_nl_eam

    ! Users should add all user specific namelist changes below in the form of
    ! namelist_var = new_namelist_value

    !theta_hydrostatic_mode=.true.
    !tstep_type=5

    theta_hydrostatic_mode=.false.
    tstep_type=9


    se_ne=0
    se_ne_x=0
    se_ne_y=0

    nu_top=1e4
    se_tstep=8.3333333333333
    dt_tracer_factor = 6
    hypervis_subcycle_q=6

    iradsw = 3
    iradlw = 3

    empty_htapes=.TRUE.
!If using dtime=100s, 9 means 15 min.
    nhtfrq=9,36
    mfilt=1,1
    avgflag_pertape='I','A'
    fincl1='PS','CAPE','CIN','PSL','PRECL','PRECC','PRECT','SHFLX','LHFLX','TS','SST','CLDLOW','CLDMED','CLDHGH', 'CLDTOT','TMCLDLIQ', 'TMCLDICE', 'TMRAINQM', 'TMCLDRIM', 'TMQ', 'TREFHT', 'QREFHT', 'FSNTOA','FLNT', 'FLNTC', 'FSNTOAC', 'FSNS', 'FSDS', 'FLNS', 'TMNUMLIQ','TMNUMICE','TMNUMRAI','TGCLDLWP','FSUTOA' ,'SWCF','LWCF','TGCLDIWP','TUQ','TVQ','FSDSC','FSNSC','FLNSC','FLDS','FLDSC', 'FSNT','FSNTC','FLUT','FLUTC','T','U','V','Q','Z3','OMEGA','RELHUM'
    fincl2='PS','CAPE','CIN','PSL','PRECL','PRECC','PRECT','SHFLX','LHFLX','TS','SST','CLDLOW','CLDMED','CLDHGH', 'CLDTOT','TMCLDLIQ', 'TMCLDICE', 'TMRAINQM', 'TMCLDRIM', 'TMQ', 'TREFHT', 'QREFHT', 'FSNTOA','FLNT', 'FLNTC', 'FSNTOAC', 'FSNS', 'FSDS', 'FLNS', 'TMNUMLIQ','TMNUMICE','TMNUMRAI','TGCLDLWP','FSUTOA' ,'SWCF','LWCF','TGCLDIWP','TUQ','TVQ','FSDSC','FSNSC','FLNSC','FLDS','FLDSC', 'FSNT','FSNTC','FLUT','FLUTC'
  
    ncdata='/pscratch/sd/s/souravt/2023-08-11_17-07-43/HICCUP.atm_era5.2009-01-03.WNP_ne128x8pg2.L128.nc'

    !.......................................................
    ! nudging
    !.......................................................
    Nudge_Model          = .true.
    Nudge_Path           = '/global/cfs/projectdirs/m1867/souravt/ERA5_nudging/WNP_ne128x8pg2_200901/'
    Nudge_File_Template  = 'ERA5_WNP_ne128x8pg2_L128.%y-%m-%d-%s.nc'
    Nudge_Times_Per_Day  = 24  !! nudging input data frequency
    Model_Times_Per_Day  = 864 !! should not be larger than 48 if dtime = 1800s
    Nudge_Uprof          = 2
    Nudge_Ucoef          = 1.
    Nudge_Vprof          = 2
    Nudge_Vcoef          = 1.
    Nudge_Tprof          = 2
    Nudge_Tcoef          = 1.
    Nudge_Qprof          = 2
    Nudge_Qcoef          = 1.
    Nudge_PSprof         = 0
    Nudge_PScoef         = 0.
    Nudge_Beg_Year       = 0001
    Nudge_Beg_Month      = 1
    Nudge_Beg_Day        = 1
    Nudge_End_Year       = 9999
    Nudge_End_Month      = 1
    Nudge_End_Day        = 1
    Nudge_Hwin_lo        = 1.0
    Nudge_Hwin_hi        = 0.0
    Nudge_Hwin_lat0      = 48.0
    Nudge_Hwin_latWidth  = 16.
    Nudge_Hwin_latDelta  = 0.5
    Nudge_Hwin_lon0      = 236.
    Nudge_Hwin_lonWidth  = 24.
    Nudge_Hwin_lonDelta  = 0.5
    Nudge_Vwin_Lindex    = 0.
    Nudge_Vwin_Hindex    = 129.
    Nudge_Vwin_Ldelta    = 0.001
    Nudge_Vwin_Hdelta    = 0.001
    Nudge_Vwin_lo        = 0.
    Nudge_Vwin_hi        = 1.
    Nudge_Method         = 'Linear'
    Nudge_Loc_PhysOut    = .True.
    Nudge_Tau            = 3.       !! relaxation time scale, unit: 6h
    Nudge_CurrentStep    = .False.
    Nudge_File_Ntime     = 1
EOF

    cat <<EOF >> user_nl_elm
        hist_dov2xy = .true.,.true.
        hist_fincl2 = 'H2OSNO', 'QRUNOFF', 'QSNOMELT', 'SNORDSL', 'WIND', 'U10', 'VA', 'TBOT', 'THBOT', 'QBOT', 'PBOT', 'ZBOT', 'PSurf', 'RAIN','RH2M','Q2M', 'LWdown', 'TSA', 'QSNWCPICE', 'H2OROF', 'QH2OSFC', 'H2OSOI', 'TWS', 'TSOI', 'TSOI_10CM', 'TSOI_ICE', 'QSOIL', 'QOVER', 'INT_SNOW'
        hist_mfilt = 1,1
        hist_nhtfrq = 0,-1
        hist_avgflag_pertape = 'A','A'

     finidat = '/pscratch/sd/s/souravt/e3sm_scratch/pm-cpu/SCREAMv0.ICRUELM.WNP_ne1024pg2_r0125_oRRS18to6v3.202308141428/run/SCREAMv0.ICRUELM.WNP_ne1024pg2_r0125_oRRS18to6v3.202308141428.elm.r.2009-01-03-00000.nc'
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
