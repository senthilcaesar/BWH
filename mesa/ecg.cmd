${s=ECG}
${a=[${s}][_ht_mag]}   

% restrict to NREM
MASK ifnot=N2,N3
RE

% convert all to mV & resample (most at SR=200 already)
mV 
RESAMPLE sr=200

% throw out any particularly unusual epochs
CHEP-MASK ep-th=4,4
CHEP epochs
RE

% find (simple/numeric) peaks

EPOCH 

HILBERT sig=${s} f=10,20 tw=0.5,5 ripple=0.01,0.01

Z-PEAKS cache=pks1 sig=${a} w=10 th=2 annot=P1

CACHE dump int=pks1

ANNOTS

TLOCK cache=pks1 w=1 same-channel=T sig=${s} seed-postfix=_ht_mag

