#!/usr/bin/env bash

function psredit_read(){
    datafile=$1
    param_name=$2
    
    param_value=`psredit -c $param_name $datafile 2> /dev/null -q -Q | tr -d ' '`
    echo $param_value 
}

echo "========================================================================="
echo "This script does the following operations."
#echo "1. Convert to PSRFITS"
echo "1. Time collapse"
echo "2. Make bandwidth positive"
#echo "4. Correct source coordinates"
echo "3. Correct the frequency between MJDs 59217 and 59424"
echo "========================================================================="
#echo "Ensure that:"
#echo "1. It is run on a copy of the original data. It rewrites data files."
#echo "2. It is run only once. Otherwise, the frequency will be over-corrected."
#echo "========================================================================="
read -n 1 -r -s -p $'Press ENTER to continue, Ctrl+C to abort.\n'

read -p "Enter Pulsar Name: " pulsarname
read -p "Enter Subbands for Band3-100 MHz: " b3_sub_100
read -p "Enter Subbands for Band5-100 MHz: " b5_sub_100
read -p "Enter Subbands for Band3-200 MHz: " b3_sub_200
read -p "Enter Subbands for Band5-200 MHz: " b5_sub_200

#exclude_dir="/Data/bcj/INPTA/Pinta.V6.2.Results/DR/Cycle46Results/:/Data/bcj/INPTA/Pinta.V6.2.Results/DR/Cycle*Results/*/BAND4/:/Data/bcj/INPTA/Pinta.V6.2.Results/DR/Master_SNR_log/"

echo "Copying Fits files for $pulsarname..."
find "/Data/bcj/INPTA/Pinta.V6.2.Results/DR/" \
    -type f \
    ! -path "/Data/bcj/INPTA/Pinta.V6.2.Results/DR/Cycle46Results/*" \
    ! -path "/Data/bcj/INPTA/Pinta.V6.2.Results/DR/Cycle*Results/*/BAND4/*" \
    ! -path "/Data/bcj/INPTA/Pinta.V6.2.Results/DR/Cycle*Results/*/BAND4_850/*" \
    ! -path "/Data/bcj/INPTA/Pinta.V6.2.Results/DR/Master_SNR_log/*" \
    ! -path "/Data/bcj/INPTA/Pinta.V6.2.Results/DR/Cycle*Results/*/BAND3PA/*" \
    ! -path "/Data/bcj/INPTA/Pinta.V6.2.Results/DR/Cycle35Results/*/BAND5PA/*" \
    ! -path "/Data/bcj/INPTA/Pinta.V6.2.Results/DR/Cycle33Results/*/BAND5PA/*" \
    -name "*$pulsarname*" \
    -name "*.rfiClean.fits" \
    -exec cp {} . \;

for archive_file in *;
do
    echo
    echo Processing $archive_file ...
    echo `realpath $archive_file`

    # Time collapse, convert to PSRFITS
    echo pam -T -m $archive_file
    pam -T -m $archive_file  
    
    # If bandwidth is negative, reverse channels
    bw=$(psredit_read $archive_file bw)
    #bw=`psredit -c bw $archive_file 2> /dev/null -q -Q | tr -d ' '`
    echo "BW: $bw"	    
    if test $bw -le 0
    then
        echo pam --reverse_freqs -m $archive_file
        pam --reverse_freqs -m $archive_file
    fi
    
    # Update the coordinates
    #echo update_coords.sh $archive_file
    #update_coords.sh $archive_file
    
    # Frequency correction
    data_mjd=$(psredit_read $archive_file "int[0]:mjd")
    #data_mjd=`psredit -c "int[0]:mjd" $archive_file 2> /dev/null -q -Q | tr -d ' '`
    echo "MJD: $data_mjd"
    if test `echo "$data_mjd>=59218 && $data_mjd<59424" | bc -l` == 1
    then
        freq_centre=$(psredit_read $archive_file freq)
        #freq_centre=`psredit -c freq $archive_file 2> /dev/null -q -Q | tr -d ' '`
        freq_centre_new=`echo $freq_centre + 0.01 | bc -l`
        echo pam -o $freq_centre_new -m $archive_file 2> /dev/null
        pam -o $freq_centre_new -m $archive_file 2> /dev/null
    fi
    bw1=100
    bw2=200
    bwarch=$(psredit_read $archive_file bw)
# This is for Band3 archive files
    if [[ -f "$archive_file" && "$archive_file" == *"500"* ]]; then
# If BW is 100 MHz, then set nchan to corresponding input value of 100 MHz 
        if [ "$bwarch" -eq "$bw1" ]; then
        	echo pam -m --setnchn $b3_sub_100 $archive_file
        	pam -m --setnchn $b3_sub_100 $archive_file
	fi
# If BW is 200 MHz, then set nchan to corresponding input value of 200 MHz 
	if [ "$bwarch" -eq "$bw2" ]; then
		echo pam -m --setnchn $b3_sub_200 $archive_file
		pam -m --setnchn $b3_sub_200 $archive_file
	fi 
    fi

# This is for Band5 archive files
    if [[ -f "$archive_file" && "$archive_file" == *"1460"* ]]; then
# If BW is 100 MHz, then set nchan to corresponding input value of 100 MHz 
        if [ "$bwarch" -eq "$bw1" ]; then
                echo pam -m --setnchn $b5_sub_100 $archive_file
                pam -m --setnchn $b5_sub_100 $archive_file
        fi
# If BW is 200 MHz, then set nchan to corresponding input value of 200 MHz 
        if [ "$bwarch" -eq "$bw2" ]; then
                echo pam -m --setnchn $b5_sub_200 $archive_file
                pam -m --setnchn $b5_sub_200 $archive_file
        fi
    fi
    
    mv "$archive_file" "${archive_file}.tmp"
done
