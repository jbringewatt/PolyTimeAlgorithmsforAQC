#!/bin/bash
n=20
guessDist=1;
zero="0"
one="1"
while [ $n -lt 80 ]; do
    guessDist=1;
    while [ $guessDist -lt 15 ]; do
        echo $n > "OnePriorGroverDepthM1_$n.$guessDist.well.input"
        echo "OnePriorGroverDepthM1TwoWell_$n.$guessDist" > "OnePriorGroverDepthM1_$n.$guessDist.params.twoWell.input"
        echo "OnePriorGroverDepthM1TB_$n.$guessDist" > "OnePriorGroverDepthM1_$n.$guessDist.params.tb.input"
        #add number of wells
        echo "2" >> "OnePriorGroverDepthM1_$n.$guessDist.well.input"
        #add marked item
        i=0
        location=""
        while [ $i -lt $n ]; do
            location="$location$zero"
            let i=i+1
        done
        echo $location >>  "OnePriorGroverDepthM1_$n.$guessDist.well.input"
        echo "-1" >>  "OnePriorGroverDepthM1_$n.$guessDist.well.input"
        echo "1" >>  "OnePriorGroverDepthM1_$n.$guessDist.well.input"
        #add guess
         i=0
        location=""
        while [ $i -lt $n ]; do
            if [ $i -lt $guessDist ]
            then
                location="$location$one"
            else
                location="$location$zero"
            fi            
            let i=i+1
        done
        echo $location >>  "OnePriorGroverDepthM1_$n.$guessDist.well.input"
        echo "-1" >>  "OnePriorGroverDepthM1_$n.$guessDist.well.input" #well depth is -1 (do more runs later with some other depth)
        echo "1" >>  "OnePriorGroverDepthM1_$n.$guessDist.well.input"


        #include wfs
        echo "1" >> "OnePriorGroverDepthM1_$n.$guessDist.params.twoWell.input"
        echo "1" >> "OnePriorGroverDepthM1_$n.$guessDist.params.tb.input"
        echo "1" >> "OnePriorGroverDepthM1_$n.$guessDist.params.tb.input"

        #tbOrder
        echo "0" >> "OnePriorGroverDepthM1_$n.$guessDist.params.tb.input"

        #etol
        echo "1e-3" >> "OnePriorGroverDepthM1_$n.$guessDist.params.tb.input"

        #sList=0 
        echo "0" >> "OnePriorGroverDepthM1_$n.$guessDist.params.twoWell.input"
        echo "0" >> "OnePriorGroverDepthM1_$n.$guessDist.params.tb.input"
        #set nunmSteps, iMax, telescopeParam
        echo "51" >> "OnePriorGroverDepthM1_$n.$guessDist.params.twoWell.input"
        echo "51" >> "OnePriorGroverDepthM1_$n.$guessDist.params.tb.input"
        echo "10" >> "OnePriorGroverDepthM1_$n.$guessDist.params.twoWell.input" 
        echo "10" >> "OnePriorGroverDepthM1_$n.$guessDist.params.tb.input"
        echo "0" >> "OnePriorGroverDepthM1_$n.$guessDist.params.twoWell.input"
        echo "0" >> "OnePriorGroverDepthM1_$n.$guessDist.params.tb.input"
        #degeneracy
        echo "0" >> "OnePriorGroverDepthM1_$n.$guessDist.params.tb.input"
        
        #run codes
        # ../../../Codes/TwoWellSolver_GroverWithPriors "OnePriorGroverDepthM1_$n.$guessDist.well.input" "OnePriorGroverDepthM1_$n.$guessDist.params.twoWell.input"
        ../../../Codes/TBSolver_GroverWithPriorsFinal "OnePriorGroverDepthM1_$n.$guessDist.well.input" "OnePriorGroverDepthM1_$n.$guessDist.params.tb.input"
        let guessDist=guessDist+1
    done
    let n=n+10    
done