	data paul_correlation;
	input similarity $ age $  count;
	datalines;
	Young <30 6
	Old <30 3
	Young 30-39 4
	Old 30-39 6
	Young 40+ 3
	Old 40+ 8
	run;
ods html file='\\Client\T$\Desktop\MTroester\PaulSignature.xls';
	proc freq data = paul_correlation order =data;
	weight count;
	tables similarity*age / nocol chisq cmh scores=modridit; /* see if significance varies much */
	tables similarity*age / nocol chisq cmh;
	run;

	proc freq data = paul_correlation order =data;
	weight count;
	tables age*similarity/chisq cmh riskdiff;
	where age ne '30-39';
	run;
ods html close;
