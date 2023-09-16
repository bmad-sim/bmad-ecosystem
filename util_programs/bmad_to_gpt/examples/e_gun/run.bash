# Create gpt files

rm lat.gpt
../../../../production/bin/bmad_to_gpt

# Convert ascii gpt field files to binary

ext='_ASCII.gpt'
newext='.gdf'
for f in *$ext
do
  base=`basename $f $ext`
  newfile=$base$newext
  asci2gdf -o  $newfile $f
  echo $newfile 'written'
done

asci2gdf -o beam0.gdf lat.gpt_particles

# Run gpt

gpt -v -o output.gdf gpt.in GPTLICENSE=$GPTLICENSE

# Get trajectories

gdftrans -o trajectories.gdf output.gdf time x Bx y By z Bz G fEx fEy fEz fBx fBy fBz

# Translate trajectories to GPT ASCII format.

gdf2a -w 16 trajectories.gdf > gpt_trajectories.txt
gdf2a output.gdf > gpt_particles.txt

# Get statistics.
# Will not work if no time slices are outputted.
# Some of these arguments may be Cornell custom stuff.

## gdfa -o avgs.gdf output.gdf time Q numpar avgG avgp avgx avgy avgz stdx stdy stdz stdBx stdBy stdBz stdG nemixrms nemiyrms nemizrms avgfBz avgfEz CSbetax CSbetay CSbetaz CSalphax CSalphay CSalphaz CSgammax CSgammay CSgammaz
## gdf2a -w 16 avgs.gdf > gpt_summary.txt


