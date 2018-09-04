g++ `root-config --cflags --libs` test.cxx -o TEST
./TEST 1 2>&1|tee logs/testLog
