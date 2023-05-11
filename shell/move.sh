#!/bin/bash

while read -r case
do
	mv ${case}.hyp ${case}.eannot

done < hyps.txt
