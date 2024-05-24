#!/bin/bash

echo "What directory:"

path="/project/wangj/node819/Jaiman/Wulver/sp24/ram_4_24/test"

# e.g.: /project/wangj/node819/Jaiman/Wulver/sp24/ram_4_24/test3/VAPOR/merge/Visuals/

read test_number

directory="${path}${test_number}/"

# Confirm directory does not exist

if [ ! -d "$directory" ]; then
    echo "Press y to confirm new directory: $directory"
    read flag
    if [ "$flag" = "y" ]; then
        echo "Making ${directory}DATA/"
        echo "Making ${directory}VAPOR/Merge/Visuals/"
        mkdir -p "${directory}DATA" "${directory}VAPOR/Merge/Visuals/"
	cp "/home/jdp46/github_MHD/MHDCode/NAMELIST" "${directory}DATA/"
        exit 0
    else
        exit 0
    fi
else
    echo "Directory already exists, will not overwrite"
fi

