ReverseadHome=$HOME/packages/reversead

for num_iter in 40 60 80 100 120 140 160 180; do
  command="g++ -std=c++11 -D NUM_DOUBLE_LAYERS=$num_iter -I$ReverseadHome/include reversead.cpp -o reversead -L$ReverseadHome/lib -lreversead"
  sh -c "$command"
  echo $command
  sh -c "./reversead"
  sh -c "rm reversead"
done
