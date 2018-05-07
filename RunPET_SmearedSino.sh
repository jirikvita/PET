
dirname=Smeared
root -l $* 'PET.cpp+(1, 4, false, true, 0.0005)'

mv *.eps ${dirname}/
mv *.gif ${dirname}/
