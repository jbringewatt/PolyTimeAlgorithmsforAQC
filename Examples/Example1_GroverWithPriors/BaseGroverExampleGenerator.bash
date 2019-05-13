#!/bin/bash
n=10
zero="0"
one="1"
while [ $n -lt 80 ]; do
    echo $n > "BaseGrover_$n.well.input"
    echo "BaseGroverOneWell_$n" > "BaseGrover_$n.params.oneWell.input"
    #add number of wells
    echo "1" >> "BaseGrover_$n.well.input"
    i=0
    location=""
    while [ $i -lt $n ]; do
        location="$location$zero"
        let i=i+1
    done
    echo $location >> "BaseGrover_$n.well.input"
    echo "-1" >> "BaseGrover_$n.well.input"
    echo "1" >> "BaseGrover_$n.well.input"
    #include wfs
    echo "1" >> "BaseGrover_$n.params.oneWell.input" 
    echo "0" >> "BaseGrover_$n.params.oneWell.input" 
    #set nunmSteps, iMax, telescopeParam
    echo "51" >> "BaseGrover_$n.params.oneWell.input" 
    echo "10" >> "BaseGrover_$n.params.oneWell.input" 
    echo "0" >> "BaseGrover_$n.params.oneWell.input" 
    
    #run codes
    ../../../Codes/OneWellSolver "BaseGrover_$i.well.input" "BaseGrover_$i.params.oneWell.input"
    
    let n=n+1
done