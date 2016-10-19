Analysis of CNEEs.

1. Run cons_element_prep.R to generate annotations for each CNEE. Need to change the bedfile and filepath variables as needed.

Output: ce annotation files in ce_final directory

2. Parse presence/absence based on psl code in ../alignment/cnee_presence. Output will be a Rlist object

Output: Rlist objects in cnee_presence directory

3. Merge together all acceleration results with process_acceleration_results.R

Output: R object saved in accelTests directory

4. Analyze CNEEs
