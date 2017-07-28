import subprocess

out_file = open('output.now', 'w')

results = subprocess.run(['../../production/bin/tao', '-noplot', '-lat', 'lat.bmad'], stdout=subprocess.PIPE).stdout.decode('utf-8')

if 'contact DCS' in results or 'FATAL' in results:
  out_file.write ('"Bookkeeper" STR  "BAD"')
  print(results)
else:
  out_file.write ('"Bookkeeper" STR  "GOOD"')
