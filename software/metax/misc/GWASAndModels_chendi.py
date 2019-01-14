import logging
import os
import re
import pandas

from .. import Constants
from .. import PredictionModel
from ..gwas import Utilities as GWASUtilities

"""
This function: strand_switch, convert the allele to the paried base
Edited by Chendi on 2019-01-01
"""
def strand_switch(input_data):  ## TODO: Find a more effecient way to do this
    output=[]
    for base in input_data:
        if base == u'A':
            output_letter = u'T'
        if base == u'T':
            output_letter = u'A'
        if base == u'C':
            output_letter = u'G'
        if base == u'G':
            output_letter = u'C'
        output.append(output_letter)
            
    return output

# ## TEST case shown as below:
# testt = u'A'; 
# logging.info("testt: %s", testt)
# logging.info("switch testt: %s", switch(testt))  

"""
This function: align_data_to_alleles, handles 5 cases for allele coding scheme correction
Edited by Chendi on 2019-01-01
"""
def align_data_to_alleles(data, base, left_on, right_on):
    EA, NEA = Constants.EFFECT_ALLELE, Constants.NON_EFFECT_ALLELE
    EA_BASE, NEA_BASE = EA+"_BASE", NEA+"_BASE"
    merged = pandas.merge(data, base, left_on=left_on, right_on=right_on, suffixes=("", "_BASE"))
    
    logging.info("%s rows of Beta data", data.shape[0])
    logging.info("%s rows of Weight data", base.shape[0])
    logging.info("%s rows of merged data", merged.shape[0])

    """
    # Comment: Original code - could be used to delete strand switch
    alleles_1 = pandas.Series([set(e) for e in zip(merged[EA], merged[NEA])])
    alleles_2 = pandas.Series([set(e) for e in zip(merged[EA_BASE], merged[NEA_BASE])])
    eq = alleles_1 == alleles_2
    merged = merged[eq]
    if eq.shape[0] == 0:
        return merged # Nothing 
    logging.info("%s rows of merged data after deleting strand switch", merged.shape[0])
    """

    # This set information can be used for both strand switch correction, and discarding (ambiguous case + seem same case)
    alleles_1_2 = pandas.Series([set(e) for e in zip(merged[EA_BASE], merged[NEA_BASE], merged[EA], merged[NEA])])
    # logging.info("%s rows in alleles_1_2", alleles_1_2.shape[0])
    # logging.info("alleles1_2 ENtry 0: %s", alleles_1_2[0])
    # logging.info("alleles1_2 Entry 1: %s", alleles_1_2[1])

    # =============== Handle strand switch ===============
    switched = alleles_1_2 == set([u'A', u'T', u'C', u'G'])  
    # switched[0] = True # Toggle this line to TEST: single case strand switch
    # switched = alleles_1_2 == set([u'A', u'G']) # Toggle this line to TEST: multiple strand switch
    logging.info("%s rows to be switched", sum(switched))

    logging.info("EA Entry to be switched: %s", merged.loc[switched, EA]) ## TODO: change into --verbosity 1 level output
    merged.loc[switched, EA] = strand_switch(merged.loc[switched, EA])
    logging.info("EA Entry after being switched: %s", merged.loc[switched, EA])

    logging.info("NEA Entry to be switched: %s", merged.loc[switched, NEA])
    merged.loc[switched, NEA] = strand_switch(merged.loc[switched, NEA])
    logging.info("NEA Entry after being switched: %s", merged.loc[switched, NEA])
     
    # =============== Handle ambiguous cases (switch | flip) & seems exact cases (but with switch & flip) ===============   
    # This has to be behind the strand switch correction: avoid index errors also, if A T vs. T A etc. cases happen after strand correction, don't have to discard
    """ ### Test A&T or C&G cases: Toggle to test
    alleles_1_2[0] = set([u'A', u'T'])
    alleles_1_2[1] = set([u'C', u'G'])
    ### """

    """ ### if want to save detailed information about which cases
    delete_case_AT = alleles_1_2 == set([u'A', u'T'])
    logging.info("%s rows to be deleted due to A&T case", sum(delete_case_AT))
    delete_case_CG = alleles_1_2 == set([u'C', u'G'])
    logging.info("%s rows to be deleted due to C&G case", sum(delete_case_CG)) """

    keep_case = (alleles_1_2 != set([u'A', u'T'])) & (alleles_1_2 != set([u'C', u'G'])) 
    # keep_case[0] = False # Toggle this line to TEST: drop lines when there is no cases to be deleted
    logging.info("%s rows to be DELETED", sum(keep_case==False))
    logging.info("%s rows to be KEPT", sum(keep_case))

    merged = merged[keep_case]
    logging.info("%s rows of merged data after dropping the discarded cases", merged.shape[0])   

    # =============== Handle allele flip =============== 
    flipped = merged[EA] != merged[EA_BASE]

    logging.info("%s rows to be flipped", sum(flipped))
    # logging.info(flipped)

    Z = Constants.ZSCORE
    if Z in merged:
        merged.loc[flipped, Z] = - merged.loc[flipped, Z]
    B = Constants.BETA
    if B in merged:
        merged.loc[flipped, B] = - merged.loc[flipped, B]

    merged.loc[flipped, EA] = merged.loc[flipped, EA_BASE]
    merged.loc[flipped, NEA] = merged.loc[flipped, NEA_BASE]

    # =============== Handle random error: eg. AC vs. AT =============== 
    keep_case_final = (merged[EA] == merged[EA_BASE]) & (merged[NEA] == merged[NEA_BASE])
    logging.info("%s rows before final checking", merged.shape[0])
    merged = merged[keep_case_final]
    logging.info("%s rows after final checking", merged.shape[0])

    return merged

def gwas_model_intersection(args):
    gwas= GWASUtilities.load_plain_gwas_from_args(args)
    paths = PredictionModel._model_paths(args.models_folder, args.models_name_filter)
    PF = PredictionModel.WDBQF
    intersection = set()
    for db_path in sorted(paths):
        logging.log(9, "loading %s", db_path)
        model = PredictionModel.load_model(db_path)
        base = model.weights[[PF.K_RSID, PF.K_EFFECT_ALLELE, PF.K_NON_EFFECT_ALLELE]].drop_duplicates()
        b = align_data_to_alleles(gwas, base, Constants.SNP, PF.K_RSID)
        intersection.update(b[Constants.SNP])
    return intersection


