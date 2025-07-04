#!/usr/bin/env python
# coding: utf-8

# In[67]:


import requests
import json
from itertools import product
import datetime
import copy
import os
import glob
import pandas as pd
# import nbconvert

# !jupyter nbconvert --to script AlphaStream.ipynb


# # Helper Function: Protein-Sequence request

# In[4]:


"""
Function that takes the name of a protein as a *string* and gives back the aminoacid sequence through UniProts API
"""

def get_uniprot_sequence(protein_name, taxon_id="559292"):
    url = "https://rest.uniprot.org/uniprotkb/search"
    query = f'gene_exact:{protein_name} AND organism_id:{taxon_id}'

    params = {
        "query": query,
        "format": "fasta", # Wir wollen fasta weil es nices standard-format ist aus dem wir die info gut auslesen koennen
        "size": 1
    }

    response = requests.get(url, params=params)

    if response.status_code == 200 and response.text.startswith(">"): # check rightigen status-code und ob fasta text
        lines = response.text.splitlines()
        sequence = "".join(lines[1:]) # Nur die Sequenz rausholen und in einen string formattieren
        return sequence
    else:
        raise ValueError(f"No sequence found or request failed for {protein_name}")


# In[4]:


# Example usage
#get_uniprot_sequence("YB022C") # is equal to get_uniprot_sequence("Pim1")


# # Helper Function: Jason-Job-File creation
# 
# Noting right here, that these don't open a json-file just yet, they only return a job string in json format together.

# In[5]:


"""
Die Funktion soll dem modularem Aufbau von AlphaStream helfen - es ist eine unterfunktion von create_af3_job und kommt 
zu Nutzen in den compound_job Funktionen.
Input: 
protein_name : str
count : int
template : boolean
"""
def additional_sequence_json(protein_name, count = 1, template = True):
    if len(protein_name) == 0:
        raise ValueError("Enter a protein (or its sequence) to continue")
    if not isinstance(count, int) or count < 1:
        raise ValueError(f"count for {protein_name} must be a positive integer")
    if not isinstance(template, bool):
        raise ValueError(f"template parameter for {protein_name} must be a boolean")

    if len(protein_name) <= 15: 
        sequence = get_uniprot_sequence(protein_name)
    else: 
        sequence = protein_name.upper()

    the_additional_sequence = {"proteinChain" : {
        "sequence" : sequence, 
        "count" : count, 
        "usesStructureTemplate" : template}
                          }
    return the_additional_sequence, len(sequence)


"""
This function acts as a helper function in future functions.
It returns a json-File in accordance to AlphaFold Server prerequisites.
Input :
protein_name : one or multiple str
count : int or list of ints (MUST be length of *protein_name), or if int and multiple proteins -> int is used for all proteins
template : boolean
name_prefix : str
"""
def create_af3_job(*protein_name, count = 1, template = True, name_prefix="af3_job"): 
    if len(protein_name) > 1:
        if isinstance(count, int):
            count = [count] * len(protein_name)
        elif len(protein_name) != len(count):
            raise ValueError("Amount of Proteins must be equal to length of count parameter") 

    # I need to make sure the name can be inputted by the user, in case the protein_name is a sequence. WHat happens in split, what happens in ace and in rico? 
    #OK - thats an alternative, I will make it be the first 6 letters of the protein sequence, thats a 1/16384 chance of being the same as another protein           
    defined_protein_name = [p[:7] if len(p) > 15 else p for p in protein_name]
    job = {
        "name": f"{name_prefix}_{'_and_'.join(defined_protein_name)}",
        "version": 1,
        "dialect": "alphafoldserver",
        "sequences": []
    }
    tokens = 0
    for i, p_name in enumerate(protein_name):
        actual_count = count[i] if isinstance(count, list) else count
        json_seq, seq_len = additional_sequence_json(p_name, count = actual_count, template = template)
        job["sequences"].append(json_seq) #wie ist es mit count und template?
        tokens += actual_count*seq_len

    if tokens > 5120:
        raise ValueError(f"The Job for {name_prefix}_{'_and_'.join(defined_protein_name)} uses more than 5120 tokens - the compound of {list(zip(defined_protein_name, count))} is too big.")
    return job



# In[6]:


#Example usage
#create_af3_job("Pim1", "Pim1", count = 22)


# # Split Function

# In[49]:


"""
Single Protein Inference Tasker: split

Function that takes one to multiple protein names as strings and makes a list of Jason-Jobs out of it.
The Output is a Json-File for AlphaFold in list-style, so each protein gets own prediction.
You can add a *list* of counts to change the count for every protein
Input:
protein_names : one or multiple str
count : int OR list, order must be in accordance to protein name order
template : boolean
name_prefix : str
file_name : str
"""

def split(*protein_names, count = 1, template=True, file_name = "", name_prefix="split"):

    if len(protein_names) > 1:
        if isinstance(count, int):
            count = [count] * len(protein_names)
        elif len(protein_names) != len(count):
            raise ValueError("Amount of Proteins must be equal to length of count parameter")
    if len(protein_names) == 0:
        raise ValueError("Enter a protein (or its sequence) to continue")
    if len(protein_names) == 1 and not isinstance(count, int):
        raise ValueError("Amount of Proteins must be equal to length of count parameter")

    jobs = []

    for i, protein_name in enumerate(protein_names):
        actual_count = count[i] if isinstance(count, list) else count
        job = create_af3_job(protein_name, count = actual_count, template = template, name_prefix = name_prefix)
        jobs.append(job)

    x = datetime.datetime.now()

    defined_protein_names = [p[:7] if len(p) > 15 else p for p in protein_names]

    if len(file_name) == 0:
        with open(f"{'_and_'.join(defined_protein_names)}_{x.strftime("%d%b_%H_%M_%S")}.json", "w") as f:
            json.dump(jobs, f, indent=2)
    else:
        with open(f"{file_name}_{x.strftime("%d%b_%H_%M_%S")}.json", "w") as f:
            json.dump(jobs, f, indent=2)


# In[51]:


# split("Pim1", "MNQLGALAQVSRFTQNFSMENIKSEFQSLQSKLATLRTPQEFFNFKKISKPQNFGEVQSRVAYNLKYFSSNYGLIIGCLSIYTLLTNLLLLFVIVLVVAGIVGINKLKGEELVTPFGSFKTNQLYTGLVCVAVPIGFLASPISTLLWLIGASAVSVFGHASLMEKPIETVFDEETV", count = 3, file_name = "SPLIT")
# or more - Syntax Okay?
#split("SomeFakeProtein")


# # Ace Function
# Be minutious with your input.

# In[7]:


"""
Aggregate of Compound Entities: ace

Function that takes one or multiple *dictionaries* with protein names and their respective count.
spits out a json-list with as many jobs as there are dictionaries, and this time
every job can entail multiple proteins at the same time!
Input :
protein_compounds : one or multiple dictionaries:
    key : str, the protein name
    value : int, is count
template : boolean
name_prefix :  str
"""

def ace(*protein_compounds, template = True, name_prefix = "ace", file_name = "", name = "", ashelper = False): # protein_compounds are dictionaries

    jobs = []
    name = name

    if len(protein_compounds) == 0 or not protein_compounds:
        raise ValueError("Enter protein compounds to continue")

    for i, protein_compound in enumerate(protein_compounds):
        if len(protein_compound) == 0:
            raise ValueError("Don't enter empty dictionaries into ace")
        job = create_af3_job(*protein_compound.keys(), count = list(protein_compound.values()), template = template, name_prefix = f"{name_prefix}{i+1}")
        jobs.append(job)
        defined_names = [p[:7] if len(p) > 12 else p for p in protein_compound.keys()]
        if not ashelper:
            if i+1 < len(protein_compounds):
                name += "-".join(defined_names) + "_and_"
            else:
                name += "-".join(defined_names)


    x = datetime.datetime.now()


    if len(file_name) == 0:
        with open(f"{"ace_for_" if not ashelper else ""}{name}_{x.strftime("%d%b_%H_%M_%S")}.json", "w") as f:
            json.dump(jobs, f, indent=2)
    else:
         with open(f"{file_name}_{x.strftime("%d%b_%H_%M_%S")}.json", "w") as f:
            json.dump(jobs, f, indent=2)


# In[18]:


# ace({"Pim1" : 1, "YBL022C": 2}, {"mam33" : 4, "Pim1" : 1}, file_name="ACE") # Syntax zu schwierig?


# # Helper Function: Cartesian-Product for protein ranges 

# In[10]:


"""
Function that takes one single dictionary with protein names and their respective count, but this time,
the count can be a range.
It spits out a list of dictionaries with all possible combinations of the given ranges.
Its output it supposed to work as an input for the compound function: it's a *list*
Input :
protein_counts : dictionary:
    key : str, protein name
    value : int or range, the counts
"""

# CHATGPT FOR THE WIN - this gives back a list with dictionaries in all combinations of the given ranges

def combine_count_ranges(protein_counts):
    # Step 1: Normalize all values to lists
    normalized = {k: (v if isinstance(v, range) or isinstance(v, list) else [v]) for k, v in protein_counts.items()}

    # Step 2: Separate keys by how many values they have
    multi_keys = [k for k, v in normalized.items() if len(v) > 1]
    fixed_keys = {k: v[0] for k, v in normalized.items() if len(v) == 1}

    # Step 3: Create product of variable value combinations
    combinations = product(*(normalized[k] for k in multi_keys))

    # Step 4: Rebuild job dicts
    job_inputs = []
    for combo in combinations:
        job = {**fixed_keys}
        job.update(dict(zip(multi_keys, combo)))
        job_inputs.append(job)

    return job_inputs


# In[12]:


# What does its input look like: dictionaries with ranges

# input_dict = {
#     "Pim1": range(2,7),
#     "Fcyx": range(4, 9),
#     "Mrx6": range(3,7)
# }

# tryout = combine_count_ranges(input_dict)
# len(tryout)


# # Rico Function

# In[11]:


"""
Run Iterative Combination Operation: rico

This is an alternative to the previous compound function.
Function that takes *one single dictionary* with protein names and their respective count/ranges.
It uses the combine_count_ranges() function and then 
spits out a json-list with as many jobs as there are dictionaries (the output of combine function)
Input :
range_dictionary : same as for combine_count_ranges
template: boolean
name_prefix : str
"""

def rico(range_dictionary, template = True, name_prefix = "rico", file_name = ""):

    if len(range_dictionary) == 0:
        raise ValueError("Enter a dictionary with proteins as keys and their ranges as values to continue")
    for key, val in range_dictionary.items():
        if not isinstance(val, range):
            raise TypeError(f"Expected a range for '{key}', got {type(val).__name__}.")
        if len(val) == 0:
            raise ValueError(f"Range for '{key}' is empty. Enter a range and don't forget that the last number isnt included")

    protein_compounds = combine_count_ranges(range_dictionary)
    defined_name = [p[:7] if len(p) > 12 else p for p in range_dictionary.keys()]
    name = "rico_for_" + "_and_".join(defined_name)
    ace(*protein_compounds, template = template, name_prefix = name_prefix, file_name = file_name, name = name, ashelper = True)


# In[21]:


# rico(
#     {
#     "MNQLGALAQVSRFTQNFSMENIKSEFQSLQSKLATLRTPQEFFNFKKISKPQNFGEVQSRVAYNLKYFSSNYGLIIGCLSIYTLLTNLLLLFVIVLVVAGIVGINKLKGEELVTPFGSFKTNQLYTGLVCVAVPIGFLASPISTLLWLIGASAVSVFGHASLMEKPIETVFDEETV": range(1,3),
#     "Fcyx": range(1, 3)
# }, file_name="RICO")
# rico({"Pim1":[1,2], "Pim1" : [2,3]})
# rico({
#     "MNQLGALAQVSRFTQNFSMENIKSEFQSLQSKLATLRTPQEFFNFKKISKPQNFGEVQSRVAYNLKYFSSNYGLIIGCLSIYTLLTNLLLLFVIVLVVAGIVGINKLKGEELVTPFGSFKTNQLYTGLVCVAVPIGFLASPISTLLWLIGASAVSVFGHASLMEKPIETVFDEETV": range(50,53),
#     "Fcyx": range(1, 4),
#     "Mrx6": range(1,4)
# }, file_name="test")


# # Pair Function

# In[26]:


def pair(pairsubjects, *pairobjects, count = 0, file_name = "", name_prefix = "pairing"):
    if not isinstance(pairsubjects, list):
        raise ValueError("The first input (pairsubjects) needs to be a list (format: []) with one or more protein names or sequences")
    if len(pairsubjects) == 0:
        raise ValueError("Enter a non-empty list for parameter: pairsubjects")
    if len(pairobjects) == 0:
        raise ValueError("Enter parameter: pairsubjects as list with strings (format = ['Protein1', 'Protein2']) AND parameter: pairobjects as strings (format example = 'Protein3', 'Protein4')")
    if count == 0:
        count = [1] * len(pairobjects)
    else:
        if not isinstance(count, list):
            raise ValueError("Count must be a list (format = [])")
        if len(count) != len(pairobjects):
            raise ValueError("Count must be a list (format = []) with the same amount of inputs as your pairobjects")

    defined_name = [p[:7] if len(p) > 15 else p for p in pairsubjects]
    jobs = []
    job_template = create_af3_job(*pairsubjects)
    tokens = 0
    for i, p in enumerate(pairobjects):
        add_protein = additional_sequence_json(p, count[i])[0]
        paired_job = copy.deepcopy(job_template)
        paired_job["sequences"].append(add_protein)
        paired_job["name"] = f"{name_prefix + str(i+1)}_{"_and_".join(defined_name)}_and_{p[:7] if len(p) > 15 else p}"
        for j, seq in enumerate(paired_job["sequences"]):
            tokens += len(seq["proteinChain"]["sequence"]) * seq["proteinChain"]["count"]
        if tokens > 5120:
            raise ValueError(f"The job for pairing {'_and_'.join(defined_name)} with {count[i]} {p[:7] if len(p) > 15 else p} uses more than 5120 tokens.")
        else:
            tokens = 0
        jobs.append(paired_job)


    x = datetime.datetime.now()
    if len(file_name) == 0:
        with open(f"pair_{'_and_'.join(defined_name)}_{x.strftime("%d%b_%H_%M_%S")}.json", "w") as f:
            json.dump(jobs, f, indent=2)
    else:
        with open(f"{file_name}_{x.strftime("%d%b_%H_%M_%S")}.json", "w") as f:
            json.dump(jobs, f, indent=2)

pair(["Mrpl15", "fcyx"], "pim1", "mam33", "LALALALALLALALALALLALALALALALLALALALALLALALALALALLALALALLALAALLALALALLALA", file_name = "PAIR")



# # Search Function

# In[63]:


def search(parent_folder: str) -> pd.DataFrame:
    results = []
    # List all subfolders in the parent_folder
    for subfolder in os.listdir(parent_folder):
        subfolder_path = os.path.join(parent_folder, subfolder)
        #skip the md file and focus on the jobs that come in subfolders
        if not os.path.isdir(subfolder_path):
            continue
        # Find all summary_confidences JSON files
        pattern = os.path.join(subfolder_path, "fold_*_summary_confidences_*.json")
        files = glob.glob(pattern)
        best = None
        for f in files:
            with open(f, "r") as fh:
                data = json.load(fh)
            score = data.get("ranking_score", float('-inf'))
            if best is None or score > best["ranking_score"]:
                best = {
                    "folder": subfolder,
                    "iptm": data.get("iptm"),
                    "ptm": data.get("ptm"),
                    "ranking_score": score,
                    "fraction_disordered": data.get("fraction_disordered"),
                    "has_clash": data.get("has_clash")
                }
        # Add protein counts and name for their columns from job_request.json
        job_request_path = os.path.join(subfolder_path, "fold_" + subfolder + "_job_request.json")
        if os.path.exists(job_request_path):
            with open(job_request_path, "r") as fh:
                job_data = json.load(fh)
            # job_data might be a list or dict, handle both
            if isinstance(job_data, list):
                job_data = job_data[0]
            # print(job_data.get("name", []))
            col_names = job_data.get("name", []).split("_and_")
            # because of structure of names (starting with rico1_prot_and_prot_and...) you have to resplit the first one 
            col_names[0] = col_names[0].split("_")[1]
            for i, seq in enumerate(job_data.get("sequences", [])):
                sequence = seq["proteinChain"]["sequence"]
                count = seq["proteinChain"]["count"]
                col_name = col_names[i]
                best[col_name] = count
        if best:
            results.append(best)


    df = pd.DataFrame(results)
    df = df.sort_values("ranking_score", ascending=False).reset_index(drop=True)
    df = df.fillna("None")
    return df


# # Zusammenfassung und tests

# In[15]:


# # Get Sequence *helper function*
# get_uniprot_sequence() # Try mistake

# # Get Protein Jobs
# proteins_to_af3_job()
# count_proteins_to_af3_job()

# # Compound jobs
# single_compounds_af3({"Pim1" : 1, "fcyx": 2}, {"mam33" : 4, "ATP11" : 1})

# # Combine Ranges *helper function*
# input_dict = {}
# combine_count_ranges(input_dict)

# # better compound job
# iterative_compound_af3(input_dict)


# In[16]:


# need to create ligand, ion implementation, not that urgent
# need to mnake code accessible to osmanlab, urgent, but these steps are necessarz first:
# have functional code ready and then put it on the v0 website - question, can I access and work on the website easily?
# - update: they use jupyter lab as well


# In[65]:


# """
# Edge case tests for AlphaStream functions: split, ace, and rico
# """

# # SPLIT FUNCTION EDGE CASES TO TEST

# def test_split_edge_cases():
#     """Test cases for the split function"""

#     print("=== SPLIT FUNCTION EDGE CASES ===\n")

#     # 1. Empty protein names
#     try:
#         split("")  # Empty string
#         print("❌ ISSUE: Empty string should raise ValueError")
#     except ValueError as e:
#         print("✅ Empty string handled correctly:", e)

#     # 2. None values
#     try:
#         split(None)  # None value
#         print("❌ ISSUE: None should raise error")
#     except (ValueError, TypeError) as e:
#         print("✅ None value handled correctly:", e)

#     # 3. Count edge cases
#     try:
#         split("Pim1", count=0)  # Zero count
#         print("❌ ISSUE: Zero count should be invalid")
#     except ValueError as e:
#         print("✅ Zero count handled correctly:", e)

#     try:
#         split("Pim1", count=-1)  # Negative count
#         print("❌ ISSUE: Negative count should be invalid")
#     except ValueError as e:
#         print("✅ Negative count handled correctly:", e)

#     try:
#         split("Pim1", count="five")  # String count
#         print("❌ ISSUE: String count should be invalid")
#     except (ValueError, TypeError) as e:
#         print("✅ String count handled correctly:", e)

#     # 4. Count list mismatches
#     try:
#         split("Pim1", "Act1", count=[1])  # Too few counts
#         print("❌ ISSUE: Count list mismatch should raise error")
#     except ValueError as e:
#         print("✅ Count list mismatch handled correctly:", e)

#     try:
#         split("Pim1", count=[1, 2])  # Too many counts for single protein
#         print("❌ ISSUE: Too many counts should raise error")
#     except ValueError as e:
#         print("✅ Too many counts handled correctly:", e)

#     # 5. Very long sequences (edge case for naming)
#     very_long_seq = "M" * 1000
#     try:
#         split(very_long_seq)
#         print("✅ Very long sequence handled (check file naming)")
#     except Exception as e:
#         print("❌ Very long sequence failed:", e)

#     # 6. Special characters in file_name
#     try:
#         split("Pim1", file_name="test/file")  # Path separators
#         print("⚠️  File name with path separators - check if valid")
#     except Exception as e:
#         print("❌ Special characters in filename failed:", e)

# # ACE FUNCTION EDGE CASES TO TEST

# def test_ace_edge_cases():
#     """Test cases for the ace function"""

#     print("\n=== ACE FUNCTION EDGE CASES ===\n")

#     # 1. Empty dictionary
#     try:
#         ace({})  # Empty dict
#         print("❌ ISSUE: Empty dict should be handled")
#     except (ValueError, KeyError) as e:
#         print("✅ Empty dict handled correctly:", e)

#     # 2. Dictionary with empty protein names
#     try:
#         ace({"": 1})  # Empty protein name as key
#         print("❌ ISSUE: Empty protein name should raise error")
#     except ValueError as e:
#         print("✅ Empty protein name handled correctly:", e)

#     # 3. Dictionary with invalid counts
#     try:
#         ace({"Pim1": 0})  # Zero count
#         print("❌ ISSUE: Zero count should be invalid")
#     except ValueError as e:
#         print("✅ Zero count in ace handled correctly:", e)

#     try:
#         ace({"Pim1": -1})  # Negative count
#         print("❌ ISSUE: Negative count should be invalid")
#     except ValueError as e:
#         print("✅ Negative count in ace handled correctly:", e)

#     try:
#         ace({"Pim1": "two"})  # String count
#         print("❌ ISSUE: String count should be invalid")
#     except (ValueError, TypeError) as e:
#         print("✅ String count in ace handled correctly:", e)

#     # 4. None as dictionary key or value
#     try:
#         ace({None: 1})  # None as key
#         print("❌ ISSUE: None as key should be invalid")
#     except (ValueError, TypeError) as e:
#         print("✅ None as key handled correctly:", e)

#     try:
#         ace({"Pim1": None})  # None as value
#         print("❌ ISSUE: None as value should be invalid")
#     except (ValueError, TypeError) as e:
#         print("✅ None as value handled correctly:", e)

#     # 5. Very large compound (memory/performance test)
#     large_compound = {"Pim1": 1 for i in range(100)}
#     try:
#         ace(large_compound, file_name="large_test")
#         print("✅ Large compound handled (check performance)")
#     except Exception as e:
#         print("❌ Large compound failed:", e)

# # RICO FUNCTION EDGE CASES TO TEST

# def test_rico_edge_cases():
#     """Test cases for the rico function"""

#     print("\n=== RICO FUNCTION EDGE CASES ===\n")

#     # 1. Empty dictionary
#     try:
#         rico({})  # Empty dict
#         print("❌ ISSUE: Empty dict should raise ValueError")
#     except ValueError as e:
#         print("✅ Empty dict handled correctly:", e)

#     # 2. Non-range values
#     try:
#         rico({"Pim1": 5})  # Integer instead of range
#         print("❌ ISSUE: Non-range should raise TypeError")
#     except TypeError as e:
#         print("✅ Non-range value handled correctly:", e)

#     try:
#         rico({"Pim1": [1, 2, 3]})  # List instead of range
#         print("❌ ISSUE: List should raise TypeError")
#     except TypeError as e:
#         print("✅ List value handled correctly:", e)

#     # 3. Empty ranges
#     try:
#         rico({"Pim1": range(5, 5)})  # Empty range
#         print("❌ ISSUE: Empty range should raise ValueError")
#     except ValueError as e:
#         print("✅ Empty range handled correctly:", e)

#     try:
#         rico({"Pim1": range(5, 3)})  # Backwards range
#         print("❌ ISSUE: Backwards range should raise ValueError")
#     except ValueError as e:
#         print("✅ Backwards range handled correctly:", e)

#     # 4. Negative ranges
#     try:
#         rico({"Pim1": range(-2, 0)})  # Negative range
#         print("⚠️  Negative range - should this be allowed?")
#     except Exception as e:
#         print("❌ Negative range failed:", e)

#     # 5. Very large ranges (performance test)
#     try:
#         rico({"Pim1": range(1, 1000)})  # Large range
#         print("⚠️  Large range - this will create 999 jobs! Memory/time concern")
#     except Exception as e:
#         print("❌ Large range failed:", e)

#     # 6. Multiple large ranges (exponential explosion)
#     try:
#         rico({"Pim1": range(1, 10), "Act1": range(1, 10)})  # 9x9=81 combinations
#         print("⚠️  Multiple ranges - creates 81 jobs, check if intentional")
#     except Exception as e:
#         print("❌ Multiple ranges failed:", e)

# # ADDITIONAL CROSS-FUNCTION EDGE CASES

# def test_general_edge_cases():
#     """Test cases that affect multiple functions"""

#     print("\n=== GENERAL EDGE CASES ===\n")

#     # 1. Unicode/special characters in protein names
#     try:
#         split("Pim1α")  # Greek letter
#         print("✅ Unicode characters handled")
#     except Exception as e:
#         print("❌ Unicode failed:", e)

#     # 2. Very long protein names (affects file naming)
#     long_name = "VeryLongProteinNameThatExceedsNormalLimits" * 3
#     try:
#         split(long_name[:50])  # Truncate to reasonable length
#         print("✅ Long protein names handled")
#     except Exception as e:
#         print("❌ Long protein names failed:", e)

#     # 3. Template parameter edge cases
#     try:
#         split("Pim1", template="yes")  # String instead of bool
#         print("❌ ISSUE: Non-boolean template should raise error")
#     except (ValueError, TypeError) as e:
#         print("✅ Non-boolean template handled correctly:", e)

#     # 4. File system limitations
#     try:
#         split("Pim1", file_name="a" * 255)  # Very long filename
#         print("⚠️  Very long filename - OS may reject")
#     except Exception as e:
#         print("❌ Long filename failed:", e)




# test_split_edge_cases()
# test_ace_edge_cases() 
# test_rico_edge_cases()
# test_general_edge_cases()




# In[51]:


#TESTS
#get_uniprot_sequence("")  # Empty protein name (should raise ValueError)
# get_uniprot_sequence(None)  # None as protein name (should raise ValueError)
#get_uniprot_sequence("FAKEPROTEINNAME12345")  # Non-existent protein (should raise ValueError from API response)
# additional_sequence_json("ACT1", count=0)  # count is zero (should raise ValueError)
# additional_sequence_json("ACT1", count=-5)  # count is negative (should raise ValueError)
# additional_sequence_json("ACT1", template="yes")
# create_af3_job("", count=1)  # Empty protein name (should raise ValueError)
# create_af3_job("ACT1", "ACT2", count=[1])  # count list too short (should raise ValueError)
# create_af3_job("ACT1", "ACT2", count=[1, 0])  # count list contains zero (should raise ValueError)
# create_af3_job("ACT1", count="one")  # count is not int or list (should raise ValueError)
# create_af3_job("ACT1", name_prefix=123) 

