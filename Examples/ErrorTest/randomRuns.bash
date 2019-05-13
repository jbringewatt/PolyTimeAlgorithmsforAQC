#!/bin/bash
#generate numExamples
numExamples=2450
i=50
zero="0"
one="1"
while [ $i -lt $numExamples ]; do
    #generate a random number of qubits between 4 and 10
    n=$(( $RANDOM % 6 + 5 ))
    echo $n > "RandomizedTestCases_$i.well.input"
    echo "RandomizedTestCasesTB0_$i" > "RandomizedTestCases_$i.params.tb0.input"
    echo "RandomizedTestCasesTB1_$i" > "RandomizedTestCases_$i.params.tb1.input"
    echo "RandomizedTestCasesFull_$i" > "RandomizedTestCases_$i.params.full.input"
    #generate a random number of wells between 2 and 10
    numWells=$(( $RANDOM % 8 + 2 ))
    echo $numWells >> "RandomizedTestCases_$i.well.input"
    #generate wells and their properties 
    #keep track of max depth well and its degeneracy, don't worry about overlaps other than centers, error parameter should catch this
    j=0
    maxDepth1=0
    maxDepth2=0
    degeneracy=1
   
    while [ $j -lt $numWells ]; do
        k=0
        bitstring=""
        while [ $k -lt $n ]; do
            if [ $(( $RANDOM % 2 )) == "0" ]
			then
				bitstring="$bitstring$zero"
			else
				bitstring="$bitstring$one"
			fi
            let k=k+1
        done
        locations[j]=$bitstring
		l=0
		tryagain=0
		while [ $l -lt $j ]; do
			if [ ${locations[j]} == ${locations[l]} ]
			then
				let tryagain=1
				break
			fi
			let l=l+1
		done
		if [ $tryagain -eq 1 ]
		then
			continue
		fi
        echo $bitstring >> "RandomizedTestCases_$i.well.input"
        #randomly generate depth for well
        pt1=$(( $RANDOM % 5 + 1))
        pt2=$(( $RANDOM % 100 ))
        if [ $pt1 == $maxDepth1 ]
        then
            if [ $pt2 == $maxDepth2 ]
            then
                let pt2=pt2+1
            fi
        fi       
        if [ $pt1 -ge $maxDepth1 ]
        then
            let maxDepth1=pt1
            if [ $pt2 -gt $maxDepth2 ]
            then
                let maxDepth2=pt2
            fi
        fi
        echo -$pt1.$pt2 >> "RandomizedTestCases_$i.well.input"

        #get a well width in 1-3
        #echo $(( $RANDOM % 3 + 1 )) >> "RandomizedTestCases_$i.well.input"
        echo "1" >> "RandomizedTestCases_$i.well.input"
        let j=j+1        
    done
    #include wfs
    echo "1" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "1" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "1" >> "RandomizedTestCases_$i.params.full.input"
    #include tb wfs
    echo "1" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "1" >> "RandomizedTestCases_$i.params.tb1.input"
    #set tightbinding order
    echo "0" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "1" >> "RandomizedTestCases_$i.params.tb1.input"
    #set etol (currently hard coded...how to automate...)
    echo "1e-1" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "1e-1" >> "RandomizedTestCases_$i.params.tb1.input"
    #numEnergies for full solver
    echo $(( $numWells + 2 )) >> "RandomizedTestCases_$i.params.full.input"
    #set sList=1
    echo "1" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "1" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "1" >> "RandomizedTestCases_$i.params.full.input"

    echo "21" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "21" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "21" >> "RandomizedTestCases_$i.params.full.input"
    #set sValues to run at
    echo "0.0" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "0.0" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "0.0" >> "RandomizedTestCases_$i.params.full.input"

    echo "0.05" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "0.05" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "0.05" >> "RandomizedTestCases_$i.params.full.input"

    echo "0.1" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "0.1" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "0.1" >> "RandomizedTestCases_$i.params.full.input"
  
    echo "0.15" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "0.15" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "0.15" >> "RandomizedTestCases_$i.params.full.input"

    echo "0.2" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "0.2" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "0.2" >> "RandomizedTestCases_$i.params.full.input"

    echo "0.25" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "0.25" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "0.25" >> "RandomizedTestCases_$i.params.full.input"

    echo "0.30" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "0.30" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "0.30" >> "RandomizedTestCases_$i.params.full.input"

    echo "0.35" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "0.35" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "0.35" >> "RandomizedTestCases_$i.params.full.input"

    echo "0.40" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "0.40" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "0.40" >> "RandomizedTestCases_$i.params.full.input"

    echo "0.45" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "0.45" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "0.45" >> "RandomizedTestCases_$i.params.full.input"

    echo "0.50" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "0.50" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "0.50" >> "RandomizedTestCases_$i.params.full.input"

    echo "0.55" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "0.55" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "0.55" >> "RandomizedTestCases_$i.params.full.input"

    echo "0.60" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "0.60" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "0.60" >> "RandomizedTestCases_$i.params.full.input"

    echo "0.65" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "0.65" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "0.65" >> "RandomizedTestCases_$i.params.full.input"

    echo "0.70" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "0.70" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "0.70" >> "RandomizedTestCases_$i.params.full.input"

    echo "0.75" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "0.75" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "0.75" >> "RandomizedTestCases_$i.params.full.input"
    
    echo "0.80" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "0.80" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "0.80" >> "RandomizedTestCases_$i.params.full.input"

    echo "0.85" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "0.85" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "0.85" >> "RandomizedTestCases_$i.params.full.input"

    echo "0.90" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "0.90" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "0.90" >> "RandomizedTestCases_$i.params.full.input"

    echo "0.95" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "0.95" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "0.95" >> "RandomizedTestCases_$i.params.full.input"

    echo "1.0" >> "RandomizedTestCases_$i.params.tb0.input"
    echo "1.0" >> "RandomizedTestCases_$i.params.tb1.input"
    echo "1.0" >> "RandomizedTestCases_$i.params.full.input"
    #run codes
    ../../Codes/FullSpaceSolver "RandomizedTestCases_$i.well.input" "RandomizedTestCases_$i.params.full.input"
    ../../Codes/TBSolver "RandomizedTestCases_$i.well.input" "RandomizedTestCases_$i.params.tb0.input"
   
    let i=i+1
done