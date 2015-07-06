Pkg.add("DataFrames")
using DataFrames

# Set current working directory
cd("C:\\Copy\\개방경제\\IO-table\\IO_raw")
function VS(year)
  iotable = readtable(join(["DOD",year,".csv"]), header=false)
  return iotable
end

table = VS(2008)
