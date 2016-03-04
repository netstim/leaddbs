% trigger subject:
options.prefs=ea_prefs('');
options.root='/Volumes/Neuro_Charite/regiongrowingsegmentation/';
options.patientname='KKI2009-38';
ea_subcorticalsegmentation(options);

% trigger mni:

options.prefs=ea_prefs('');
options.root='/Volumes/Neuro_Charite/regiongrowingsegmentation/';
options.patientname='template';
ea_subcorticalsegmentation(options);