make clean
make
./CalEff 1 2>&1 |tee logs/log
./Correc 1 2>&1 |tee logs/cLog
root -l draw.cxx
