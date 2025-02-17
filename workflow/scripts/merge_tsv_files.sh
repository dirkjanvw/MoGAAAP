#!/bin/bash
# Written by Dirk-Jan van Workum (2025) with help of claude.ai

# First pass to get all relevant headers if available (cfr. 13.statistics.smk)
headers=$(awk '
BEGIN {
    FS=OFS="\t";
    delete col_map;
}

FNR==1 {
    for (i=1;i<=NF;i++) if (!($i in col_map)) col_map[$i] = 1;
}

END {
    #printf "Name\\tTotal length\\t#sequences\\tN50\\t#genes (full)\\t#genes (coding)\\t#transcripts (coding)\\t#chromosomes\\tTotal length (chromosomes)\\t#unassigned sequences\\tTotal length (unassigned sequences)\\tTotal QV\\tK-mer completeness\\t#contigs\\tContig N50\\tInput data\\tAssembler\\n" > {output.tsv}
    if ("Name" in col_map) print "Name";
    if ("Total length" in col_map) print "Total length";
    if ("#sequences" in col_map) print "#sequences";
    if ("N50" in col_map) print "N50";
    if ("#genes (full)" in col_map) print "#genes (full)";
    if ("#genes (coding)" in col_map) print "#genes (coding)";
    if ("#transcripts (coding)" in col_map) print "#transcripts (coding)";
    if ("#chromosomes" in col_map) print "#chromosomes";
    if ("Total length (chromosomes)" in col_map) print "Total length (chromosomes)";
    if ("#unassigned sequences" in col_map) print "#unassigned sequences";
    if ("Total length (unassigned sequences)" in col_map) print "Total length (unassigned sequences)";
    if ("Total QV" in col_map) print "Total QV";
    if ("K-mer completeness" in col_map) print "K-mer completeness";
    if ("#contigs" in col_map) print "#contigs";
    if ("Contig N50" in col_map) print "Contig N50";
    if ("Input data" in col_map) print "Input data";
    if ("Assembler" in col_map) print "Assembler";
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
