
dirname=AllPhis
root -l $* 'PET.cpp+(7, 4, false, false, 0.0005)'

mv *.eps ${dirname}/
mv *.gif ${dirname}/
