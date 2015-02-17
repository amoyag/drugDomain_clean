"""
Function:  master 
    Master script, calls all scripts in correct order to obtain mapping of
    small molecule binding to Proteins.

    To run, you need access to a  MySQL instance of ChEMBL > chembl_15. 
    --------------------
    Contact:
    Aurelio Moya-Garcia - aurelio.moya@ucl.ac.uk
    Felix Kruger - fkrueger@ebi.ac.uk
"""  
import yaml
import os
import queryDevice
import math
import numpy as np
import MySQLdb

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def mysql_query(query, params, config):
    ''' Extract the chembl_id, molecule and type for each compound in a list of 
    molregnos:

    Inputs:
    params          -- tuple containing query parameters
    default-file    -- path to mysql --defaults-file
    database        -- query database

    '''
    mysql = MySQLdb.connect(user=config['user'], passwd=config['pword'], port=config['port'], host=config['host'] )
    c = mysql.cursor()
    c.execute("use {0}".format(config['release']))
    c.execute(query, params)
    q_res = c.fetchall()
    c.close()
    return q_res

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getHumanTargets(config):
    query = """SELECT cs.accession, cs.component_id, tid
        FROM component_sequences cs 
            JOIN target_components tc 
            ON tc.component_id = cs.component_id  
        WHERE db_source IN('SWISS-PROT', 'TREMBL')
        AND ORGANISM = 'Homo sapiens'"""
    qres = mysql_query(query, (), config)
    targets= []
    for target in qres:
        targets.append(target[0])
    return targets

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

# Aurelio - I suggest we a query that is not limited to drugs - this way we 
# capture all domains capable of interactions with small molecules. Since 
# newer releases of chembl have the pchembl value, we can actually use that to
# filter, and drop the whole filterForTarget step \o/

def getLigandsForTarget(target, config):  
    query  = """SELECT DISTINCT ms.canonical_smiles, act.pchembl_value, act.molregno, act.activity_id
    FROM component_sequences cs 
    JOIN target_components tc ON cs.component_id=tc.component_id 
    JOIN target_dictionary td ON tc.tid=td.tid 
    JOIN assays ass ON ass.tid = tc.tid
    JOIN activities act  ON ass.assay_id=act.assay_id
    JOIN molecule_dictionary md ON act.molregno=md.molregno
    JOIN compound_structures ms ON act.molregno=ms.molregno 
    WHERE act.pchembl_value >= {threshold}
    AND potential_duplicate IS NULL
    AND(
        data_validity_comment IS NULL
        OR data_validity_comment = 'manually validated'
        )
    AND assay_type = 'B'
    AND relationship_type = 'D'
    AND target_type = 'SINGLE PROTEIN'
    AND standard_relation = '='
    AND accession = %s""".format(threshold=config['threshold'])
    qres = mysql_query(query, (target,), config)
    ligands = []
    for res in qres:
        ligands.append([res[0],res[1],res[2], res[3]])
    return ligands

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def mapProtein(targets, config):
    out = open('map_protein.txt', 'w')
    smallmProtein = {}
    for target in targets:
        
        ligands = getLigandsForTarget(target, config)
        for ligand in ligands:
            smiles = ligand[0]
            aff = ligand[1]
            molregno = str(ligand[2])
            actId = str(ligand[3])
            mapType = "NA"
            out.write('%s\t%s\t%s\t%s\t%s\n'%(actId, target, molregno,target, mapType))
            
    
                    
                    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def master():
    # Read config file.
    configFile = open('mpf.yaml')
    config = yaml.safe_load(configFile)
    configFile.close()
    
    ## Get all Human targets from ChEMBL with a Uniprot accession.
    targets = getHumanTargets(config)
    
    ## Map small molecules to proteins
    
    smallmProtein = mapProtein(targets, config)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    
if __name__ == '__main__':
    import sys
    if len(sys.argv) != 1:  # the program name and the two arguments
        sys.exit("All parameters specified in config file `mpf.yaml`.")
    master()     
                            