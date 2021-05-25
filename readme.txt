##Welcome to the clinical score generator app##
To use:
Add leaddbs to path using addpath('/path/to/leaddbs')
1. Before you use the app itself, please go to lead group and generate a "LEAD_groupanalysis.mat" file containing the patient list into the root directory of your subjects. 
The group analysis file must contain the patient list you want to generate the scores for - no other parameters (stimulation,etc) are required.
2. Each Tab specifies the score type you would like to generate scores for. UPDRS is for motor scores (you can have MDS or MDS-UPDRS scores).
3. Follow the example excel sheet to make your headers easily accessible. While no order is required to be maintained, a cleaner excel sheet will lead to faster computation.
4. To open the app, navigate to the patient folder (for e.g., '/Volumes/APFS_EXCHANGE/MDST_BER' is the overarching folder containing all my subjects, so I navigate to this folder.)
Run with:
clinical_score_generator('/path/to/excel_sheet.xslx','/path/to/postop_excel_sheet.xslx','sheet_name')
- The first argument is compulsory, but second argument is optional. If both pre-op and post-op files are in the same excel sheet, you can specify that as the first argument and the second argument can be left as "". 
- The Sheet Name specifies a spreadsheet name if you have many sheets, and want to specify one. Otherwise, just specify "Sheet1" or if you leave as "", the first sheet is selected.
Please ensure that the number of entries in the sheet name should correspond to the number of subject directories you have.
5. To enter the baseline range: either specify the column in the excel sheet which corresponds to the baseline start column (e.g., B) and the baseline end column (e.g., Z). You can also click on the import from excel button and all the headers in the excel sheet will be available, and you can click and drag the range of columns and click on "OK" to lock in the selection.
6. You can select upto three Postop scores (postop 6months,12months,etc) and enable them by clicking on the checkbox. You can edit the name of the Postop as required, but no special characters are allowed.
7. Choose your laterality by clicking Left, Right or Both (Default is Both). If you choose Left/Right, all the subscores will be calculated only with "left body/right body" scores.
8. You can specify a custom score requirement by writing in the text box provided. 
For e.g., if you require UPDRS22 and UPDRS 23 scores as a custom subscore, you can write it as: UPDRS22 
UPDRS23 (in two seperate lines). 
9. Click on Generate to Generate the subscores.
10.You can also export these subscores (along with the custom subscore) in a table using the export excel button.
11. To export these scores to your lead group folder, navigate to the lead group app and select the patient directory. The scores should auto populate.
UPDRS Sub-scores:
1. Bradykinesia
2. Rigidity
3. Tremor
4. Axial
5. Custom
If you only have one of these subscores, it will not affect the algorithm: the outputs for the rest of the subscores will just be "NaN" values.
BDI Sub-scores:
1. APATHY
2. BIS
3. EQ
4. BDI
5. QUIP
6. Custom..

