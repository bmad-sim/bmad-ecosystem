#!/usr/bin/env ruby

# Script for automatically downloading data from
# Henke Layered Mirror Reflectivity
# http://henke.lbl.gov/optical_constants/layer2.html
# Produces a .dat file for use with Synrad3D
# Also produces a .gp gnuplot file and resulting .pdf
# Tested with ruby 1.8.7 and gnuplot 5.0 or 4.6
####
# STP 2018/6/14 : first version of full script
# STP 2018/6/14 : input file argument
####

if ARGV.length >= 1
  inputFile = ARGV[0]
else
  inputFile = "default.input"
end

puts "Loading input file: #{inputFile}"
load inputFile

####
#### Shouldn't need to edit anything below

RunName = "#{$substrate}_#{$layerThickness}nm#{$layer}"

def file(t,angle,suf)
  return "#{RunName}/#{t[:Emin]}_#{t[:Emax]}_#{'%05.2f'%angle}_#{suf}.dat"
end

####

puts "The surface is: #{$layerThickness}nm of #{$layer} on #{$substrate}"
puts "The directory is: #{RunName}"
puts "-----"
puts "Starting download."

`rm -rf #{RunName}`
`mkdir -p #{RunName}`

itab=0
$tables.each do |t|

  puts "Table #{itab+=1}: #{t[:Emin]}-#{t[:Emax]} eV:"

  npts = (t[:Emax]-t[:Emin])/t[:Estep]

  # Parameters to pass to the webpage
  params = %Q(
    Layer=#{$layer}
    Ldensity=-1
    Thick=#{$layerThickness}
    Sigma1=0
    Substrate=#{$substrate}
    Sdensity=-1
    Pol=-1
    Scan=Energy
    Min=#{t[:Emin]}
    Max=#{t[:Emax]}
    Npts=#{npts}
    Plot=Linear
    Output=Plot
  )

  #### download
  if true

  print "  Angles: "
  t[:angles].each do |angle|
  
    iparams = params.split
    iparams << "Fixed=#{angle}"
    params_str = iparams.join("&")
    
    # get the result webpage
    cmd = "wget -O laymir.pl --post-data \"#{params_str}\" henke.lbl.gov/cgi-bin/laymir.pl >& /dev/null"
    system(cmd)
    if $? != 0
       puts "Couldn't submit request to Henke webpage"
    end
    
    # extract the link to the data
    text = `grep dat laymir.pl`
    r = text.match(/.*HREF=\"(.*?)\"/)

    if !r
      puts "Bad response from website:"
      system('cat laymir.pl')
      exit
    end
    
    # download the data
    out = file(t,angle,'raw')
    cmd = "wget -O #{out} henke.lbl.gov/#{r[1]} >& /dev/null"
    system(cmd)
    if $? == 0
       #puts out
       print "#{angle} "
    else
       puts "Couldn't download data file from Henke website"
    end
    
    # don't hammer their server
    sleep 2
  end # loop
  puts
  end # if

  #### parse
  t[:angles].each do |angle|
    file = file(t,angle,'raw')
    file_out = file(t,angle,'refl')
    `tail -n +3 #{file} | awk '{printf "%.6f\\n", $2}' > #{file_out}`
  end

  # combine refl files
  file_comb = "#{RunName}/#{t[:Emin]}_#{t[:Emax]}_comb.dat"
  `paste -d ',' #{RunName}/#{t[:Emin]}_#{t[:Emax]}_*refl.dat > #{file_comb}`
  `sed -i -e 's/,/, /g' #{file_comb}`

end # loop over $tables

`rm -f laymir.pl`

#### create the file for Synrad3D

fOutName = "#{RunName}.dat"
fOut = File.open(fOutName, "w")

fOut.puts %Q(
&general
  name = "#{RunName}"
  description = "#{$substrate} with #{$layerThickness}nm #{$layer} layer"
  n_table = #{$tables.length}
  surface_roughness_rms = #{$s3d_roughness_rms}
  roughness_correlation_len = #{$s3d_roughness_correlation_len}
/
)

$tables.each do |t|

  fOut.puts %Q(
&table
  energy_min = #{t[:Emin]} ! eV
  energy_max = #{t[:Emax]}
  energy_delta = #{t[:Estep]}
  angles = 0, #{t[:angles].join(', ')}
/

)

  # format reflectivity data
  file_comb = "#{RunName}/#{t[:Emin]}_#{t[:Emax]}_comb.dat"
  iLine = 0
  File.open(file_comb).readlines.each do |lref|
    iLine+=1
    iLineStr = '% 3d' % iLine
    fOut.puts "&row ix_row = #{iLineStr}, p_reflect = 1.00, #{lref.chop} /"
  end

end

fOut.close
puts "Wrote to #{fOutName}"

#### write a gnuplot file

fOutName = "#{RunName}.gp"
fOut = File.open(fOutName, "w")

Emins = $tables.each.collect { |t| t[:Emin] }
EminsStr = Emins.join(' ')

### gnuplot code
fOut.puts %Q(
#!/usr/bin/env gnuplot

emins = "#{EminsStr}"

set xlabel 'Energy [eV]'
set ylabel 'Reflectivity' offset 1,0
set key right out reverse Left title 'Angle' samplen 0.5 box
unset colorbox
set xtics rotate by -45

outname = '#{RunName}.pdf'
set out outname
set term pdf size 20,5 fontscale 0.7
#set term png size 1500,400 fontscale 0.7
set multiplot layout 1,words(emins)

do for [emin in emins] {

fs = system('ls -1 #{RunName}/'.emin.'_*raw.dat')

# pull out the angle from the filename
ang(f) = sprintf("%g",1*f[strstrt(f,"_raw")-5:strstrt(f,"_raw")-1])

plot for [i=1:words(fs)] word(fs,i) u 1:2:(log(1*ang(word(fs,i)))) w lp palette pt 7 ps 0.5 t ang(word(fs,i))

}
unset multi
set out
print 'Wrote to '.outname
)
### end gnuplot code

fOut.close
puts "Wrote plot script to #{fOutName}"

`gnuplot #{fOutName}`
#`/home/stp44/sw/bin/gnuplot #{fOutName}`
