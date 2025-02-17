#!/bin/bash
# Written by Dirk-Jan van Workum (2025) with help of claude.ai

# First pass to get all headers and identify common columns (starting with 'Name' if exists)
headers=$(awk '
BEGIN {
    FS=OFS="\t";
    delete col_map;
    total = 0;
}

FNR==1 {
    total++;
    for (i=1;i<=NF;i++) {
        if (!($i in col_map)) col_map[$i] = 1;
        else col_map[$i]++;
    }
}

END {
    if ("Name" in col_map) { # If Name column exists, print it first
        print "Name";
    }
    for (col in col_map) if (col_map[col] == total && col != "Name") { # Print common columns
        print col;
    }
    for (col in col_map) if (col_map[col] != total && col != "Name") { # Print non-common columns
        print col;
    }
}
' $@ | tr '\n' '\t' | sed 's/\t$//')

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
