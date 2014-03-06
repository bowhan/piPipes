awk '{if ($1~/>/) {OFS="_"; gsub ("/","-",$0) ; print $1,$2,$3}  else {print} }' $1
