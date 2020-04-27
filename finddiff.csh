#!/bin/csh -f
touch diff.txt
set basename1 = `pwd`
#set basename2 = "../gizmo-mufasa-phew/"
set basename2 = "../specexbin_phew/"
set list = `ls ./* | grep '\.c'`
@ Nfiles = $#list
@ i = 0
while ($i < $Nfiles)
    echo "File1: "$list[$i]
    echo "File2: "$basename2$list[$i]
    @ i++
    echo "File: "$list[$i] >> "diff.txt"
    echo "================================" >> "diff.txt"
    echo "" >> "diff.txt"
    diff $list[$i] $basename2$list[$i] >> "diff.txt"
end
set list = `ls ./* | grep '\.h'`
@ Nfiles = $#list
@ i = 0
while ($i < $Nfiles)
    echo "File1: "$list[$i]
    echo "File2: "$basename2$list[$i]
    @ i++
    echo "File: "$list[$i] >> "diff.txt"
    echo "================================" >> "diff.txt"
    echo "" >> "diff.txt"
    diff $list[$i] $basename2$list[$i] >> "diff.txt"
end

