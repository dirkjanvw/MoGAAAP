#!/bin/bash

# First pass to get the headers
headers=$(head -n1 $@ | awk 'NR%3==2' | tr '\t' '\n' | sort -u | tr '\n' '\t' | sed 's/\t$/\n/')

# Second pass to print the headers and the lines in the correct order
awk -v headers="$headers" '
BEGIN {
    FS=OFS="\t";
    n = split(headers, head_arr);
    for (i=1; i<=n; i++) header_pos[head_arr[i]] = i;
    print headers;
}

FNR==1 {
    delete col_map;
    for (i=1; i<=NF; i++) col_map[i] = header_pos[$i];
    next;
}

{
    delete line_arr;
    for (i=1; i<=n; i++) line_arr[i] = "";
    for (i=1; i<=NF; i++) if (col_map[i]) line_arr[col_map[i]] = $i;
    line = "";
    for (i=1; i<=n; i++) line = line (i>1 ? OFS : "") line_arr[i];
    print line;
}
' $@
