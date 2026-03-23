#!/bin/bash

# execute.sh
# This is a shell script that executes the program SvKolROA. More information about how to execute it and which outputs are given is in the README.
# Procedure: Get input arguments. Execute LOD.py, then Lod_plot.py for non existing output files, then ClusterandPlot.R. 



# Make a usage function that gives feedback on usage to user
usage() {
    echo -e "Usage: $0 -i [input filepath] -o [output filename template]\n\nOptional Parameters:\n-p [population config filepath]\n-t [threshold config filepath]\n-w [window size]\n-r [plotting parameters]\n\toptions:\n\tp -- plot by population\n\tc -- plot by chromosome\n\ti -- plot by individual"
}

# Set default values for parameters
w=60
r="pci"
p="nofile"
t="nofile"
r="s"

peaks="nopeaks"
populations=("running_one_population_only")
all_processed="True"


# Run programs using paths relative to this scripts folder -- allows output paths to be relative to where user is
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

# Save arguments from options as variables
while getopts "i:o:p:t:w:r:h" args; do
    case "${args}" in        
        i)
            i=${OPTARG}
            ;;
        o)
            o=${OPTARG}
            ;;
        p)
            p=${OPTARG}
            ;;
        t)
            t=${OPTARG}
            ;;
        w)
            w=${OPTARG}
            ;;
        r)
            r=${OPTARG}
            ;;
        h)
            usage
            exit 0
            ;;
        *)
            echo "incorrect usage, try -h for help"
            exit 1
            ;;
    esac
done
shift $((OPTIND-1))

# Test for if input file exist
if [[ ! -f ${i} ]]; then
    echo "input file doesn't exist"
    exit 1
fi

## ----- step 1: run LOD.py -------

# Check that the population config file exist, exiting with feedback if it doesn't
if [[ ! -f "${p}" ]] && [[ "${p}" != "nofile" ]]; then
    echo "population config file doesn't exist"
    exit 1
fi

# Calculate the lod scores, if it hasn't already been done
if [[ -f ${o}.lod ]]; then
   echo "File ${o}.lod exists. Skip calculation of LOD scores"
elif [[ ${p} != "nofile" ]]; then
    python $parent_path/LOD.py ${i} ${o}.lod ${w} ${p} 
else
   python $parent_path/LOD.py ${i} ${o}.lod ${w} 
fi
   
# Give Feedback
echo "Step 1 done"

## ------- step 2: run Lod_plot.py ---------

# Get the populations into an array
if [[ -f "${p}" ]]; then
    populations=($( cut -f2 ${p}| sort -u | tr '\n' ' '))
fi


# 
for pop in "${populations[@]}"; do
    # Get the peaks for this population from the threshold config file, with error message if the threshold config file doesn't exist.
    if [[ -f "${t}" ]]; then
        peaks=($(cat "${t}" | grep "$pop" | cut -f2,3))
    elif [[ "${t}" != "nofile" ]]; then
        echo "threshold config file doesn't exist"
        exit 1
    fi    

    # Proccess each population, taking the peaks from the peaks variable if they exist, and not taking it if it doesn't exist. 
    # If the output file already exists, don't run the code again.
    if [[ "${peaks}" != "nopeaks" ]] && [[ ! -f "${o}.${pop}.lod.sig" ]] && [[ "${peaks[0]}" != "skip" ]]; then
        echo "processing population $pop"
        python $parent_path/Lod_plot.py ${o}.lod ${o}.${pop}.lod.sig $pop ${p} ${peaks[0]} ${peaks[1]}
    elif [[ "${peaks[0]}" == "skip" ]]; then
        echo "skipping population $pop"
        touch ${o}.${pop}.lod.sig
    elif [[ ! -f "${o}.${pop}.lod.sig" ]]; then
        echo "processing population $pop"
        python $parent_path/Lod_plot.py ${o}.lod ${o}.${pop}.lod.sig $pop ${p}
    else
        echo -e "population $pop already processed\n file is here: ${o}.${pop}.lod.sig"
    fi 
done


# Concatenate the population files into one file
# First check so that they all are processed
for pop in "${populations[@]}"; do
    if [[ ! -f "${o}.${pop}.lod.sig" ]]; then
        all_processed=False
    fi
done

if [[ ${all_processed} == "True" ]]; then
    echo -e "#Start\tStop\tChrom\tInd\tLength\n" > ${o}.lod.sig
    for pop in "${populations[@]}"; do
        tail ${o}.${pop}.lod.sig
        echo ${o}.${pop}.lod.sig 
        grep -v "#" ${o}.${pop}.lod.sig >> ${o}.lod.sig
    done
fi

# Give feedback
echo "Step 2 done"

## ------- step 3: run ClusterandPlot.R --------

if [[ -f ${o}.lod.sig ]]; then
    Rscript $parent_path/ClusterandPlot.R ${o}.lod.sig ${p} ${o} ${r}
elif [[ all_processed == True ]]; then
    echo "Error: output files from Lod_plot.py missing"
    exit 1
fi

echo "Processing done!"