#!/usr/bin/env bash
# -*- coding: utf-8 -*-
# @Date    : 2019-09-04 17:43:58
# @Author  : lizd (lzd_hunan@126.com)
# @Link    : ---
# @Version : $Id$



usage() {
    echo "Usage:"
    echo "    bash $0 -i in_dir [-h help information]"
    echo ""
    echo "Description:"
    echo "    Merge STAR final.out.log in specific dir."
    echo ""
    echo "Options:"
    echo "    -i, STAR output directory."
    echo "    -h, show help information."
    exit 1
}

while getopts ':i:h' opt; do
    case $opt in
        i) in_dir="$OPTARG";;
        h) usage;;
        *) usage;;
    esac
done

run=0

if [ -z "${in_dir}" ]; then
    echo -e "Error: Parameter -i should not be empty!\n"
    usage
elif [ ! -d "${in_dir}" ]; then
    echo ${in_dir} "is invalid."
    exit 1
fi

set -e
set -u
set -o pipefail

run=1

if [ ${run} -eq 1 ]; then
    sra=$(find ${in_dir} -name "*Log.final.out" | sort)
    ( \
    echo Run ${sra} | sed -r -e 's/[^ ]*\///g' -e 's/ /\t/g' -e 's/[\.\_]{0,}Log.final.out//g'

    paste ${sra} | grep "|" | sed -n -e '{5,$p}' | \
    sed -e 's/^[ \t]*//g' -e 's/^%/percentage/g' -e 's/ %/ percentage/g' -e 's/ |//g' | \
    awk -F "\t" '{for(i=1;i<=NF;i++){if(i<=2){printf $i"\t"}else if(i%2==0){printf $i"\t"}} printf "\n"}' | \
    sed 's/\t$//g' \
    ) | \
    awk -F "\t" '{for(i=1;i<=NF;i=i+1){a[NR,i]=$i}}END{for(j=1;j<=NF;j++){str=a[1,j];for(i=2;i<=NR;i++){str=str "\t" a[i,j]}print str}}'
fi
