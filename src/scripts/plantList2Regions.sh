echo '# Region file format: DS9 version 4.1'
echo 'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
echo 'icrs'
cat $1 | awk '{printf("circle(%fd,%fd,5p) # text={ %5.2f %5.2f %d %5.2f}\n", $2, $3, $6, $7, $1, $10)}'
