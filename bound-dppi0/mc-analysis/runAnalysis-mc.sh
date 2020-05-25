#!/bin/bash
for X in `seq 1 1 4`; do
    input="${WMC_DATA}/$1-$X.ems.bz2"
    #input="/data10/users/khreptak/SIMULATION/WMC/$1-$X.ems.bz2"
    output="${OUTPUT_MC}/MC-newcuts-AddGammaCut-$1-$X"
    if [ -e ${input} ]; then
        echo "MC number $X..."
        if [ -e ${output}.root ];then
            echo "...was already done."
        else
            scriptname="analysis-mc-$1-$X.sh"
            if [ -e ${scriptname} ]; then
                echo "...already in process."
            else
                echo "... PROCESSING MC $X ..."
                echo "#!/bin/bash" >> ${scriptname}
                echo "cd $PWD" >> ${scriptname}
                echo "./main -mode mc -fin file:${input} -n ${output} -abort" >> ${scriptname}
                echo >> ${scriptname}
                echo "rm -f $PWD/${scriptname}" >> ${scriptname}
                chmod u+x ${scriptname}
                #./${scriptname} &> ${output}.log
                qsub -q batch ${scriptname}
                echo "...done."
                sleep 2
            fi
        fi
    fi
done
