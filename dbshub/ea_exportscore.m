function score=ea_exportscore()

% for now static sample export:

score.studyOrigin='Earlystim'; % Tag label study name
score.studyDoi='10.1056/NEJMoa1205158'; % Tag label study name
score.metric='UPDRS-III'; % clinical score label
score.dbsOff.vals=[1,2,0,0,0,0,0,0,0,0,0,0,0,0,2,3,2,2,3,3,1,1,1,2,2,2,2];
score.dbsOff.mean=(mean(score.dbsOff.vals));
score.dbsOn.vals=[0,2,0,0,0,0,0,2,2,0,0,0,0,0,3,3,1,2,1,3,0,2,1,2,2,2,2];
score.dbsOn.mean=(mean(score.dbsOn.vals));
score.patient.DOB='05-10-1951'; % date of birth
score.patient.DOS='03-02-2014'; % date of surgery
score.patient.disease='Parkinson''''s Disease';
score.patient.diseaseSubtype='Tremor Dominant';
score.version=1.0;

savejson('', score, 'FileName', 'sample.json');
