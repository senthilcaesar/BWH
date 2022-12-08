
# helper script to be run in this folder, to make a single 'harm.txt' for NAP, etc
# this compiles that individual components

cat units.txt  >  harm.txt
cat eeg.txt    >> harm.txt
cat eog.txt    >> harm.txt
cat emg.txt    >> harm.txt
cat ecg.txt    >> harm.txt
cat resp.txt   >> harm.txt
cat flow.txt   >> harm.txt
cat oxy.txt    >> harm.txt
cat pulse.txt  >> harm.txt
cat leg.txt    >> harm.txt
cat misc.txt   >> harm.txt
cat snore.txt  >> harm.txt
cat cpap.txt   >> harm.txt


