#!/bin/bash

clines=$1

total_clines=$(cat $clines | grep -v -c "^g")
significant_cline_centers=$(cat $clines | grep -v "^g" | awk '$9 < 0.05 {print $0}' | wc -l)
non_significant_cline_centers=$(cat $clines | grep -v "^g" | awk '$9 >= 0.05 {print $0}' | wc -l)
significant_cline_widths=$(cat $clines | grep -v "^g" | awk '$5 < 0.05 {print $0}' | wc -l)
non_significant_cline_widths=$(cat $clines | grep -v "^g" | awk '$5 >= 0.05 {print $0}' | wc -l)


echo "Name of the File:        $clines"
echo "Total Clines:            $total_clines"
echo "Significant Cline Centers:      $significant_cline_centers"
echo "Not Significant Cline Centers:  $non_significant_cline_centers"
echo "Significant Cline Widths:	$significant_cline_widths"
echo "Not Significant Cline Widths:	$non_significant_cline_widths"
