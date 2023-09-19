import re
from functools import reduce

codon_dict = {}
eval_dict = {}

def read_codons(file):
    global codon_dict
    codon_dict = {}
    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            if ":" in line:
                key, sequences = line.split(":")
                sequences = sequences.split(",")
                sequences = [re.sub(r"(\w)\{(\d+)\}", lambda x: x.group(1) * int(x.group(2)), s.strip()) for s in sequences]
                codon_dict[key.strip()] = sequences

def read_evals(file):
    global eval_dict
    eval_dict = {}  # Reset the dictionary
    
    with open(file, 'r') as f:
        lines = f.readlines()
        
        for line in lines:
            parts = line.strip().split(": ")
            if len(parts) != 2:
                continue  # Skip invalid lines

            order_name = parts[0]
            operations = parts[1].split(", ")

            if len(operations) != 2:
                continue  # Skip invalid lines

            read_order, op_order = operations
            if read_order not in ["L", "R"] or op_order not in ["PO", "PR", "I"]:
                continue  # Skip invalid lines

            eval_dict[order_name] = (read_order, op_order)

    return eval_dict

def encode(sequence: str) -> str:
    amino_acids = sequence.split(" ")  # Split the sequence string into individual amino acids
    rna_sequence = ""

    for aa in amino_acids:
        # Get the codons for the amino acid
        codons = codon_dict.get(aa, [])
        if codons:
            # Sort the codons by length and select the longest one
            longest_codon = max(codons, key=len)
            rna_sequence += longest_codon

    return rna_sequence

def decode(sequence):
    amino_acids = []
    while sequence:
        found = False
        # Sort the codons by length in descending order to prioritize longer codons
        for key, values in sorted(codon_dict.items(), key=lambda x: -max(len(v) for v in x[1])):
            for value in sorted(values, key=len, reverse=True):
                if sequence.startswith(value):
                    amino_acids.append(key)
                    sequence = sequence[len(value):]
                    found = True
                    break
            if found:
                break
        if not found:
            sequence = sequence[1:]
    return ' '.join(amino_acids)

def operate(sequence, eval_name):
    # First, we need to fetch the read and operation order from eval_dict
    if eval_name not in eval_dict:
        return None
    
    read_order, operation_order = eval_dict[eval_name]
    
    # Convert sequence to list for easy manipulations
    sequence_list = list(sequence)
    if read_order == "R":
        sequence_list.reverse()
    
    # Initialize our result
    result = []
    start_flag = False
    i = 0
    while i < len(sequence_list) - 2:  # We use -2 since we'll be checking in blocks of 3
        codon = ''.join(sequence_list[i:i+3])
        
        # Decode the codon to its amino acid representation
        amino_acid = decode(codon)
        
        # If not START or STOP and the start_flag is False, continue
        if not start_flag and amino_acid not in ["START", "STOP"]:
            i += 3
            continue
        
        # Handle different amino acids
        if amino_acid == "START":
            start_flag = True
        elif amino_acid == "STOP":
            start_flag = False
        elif amino_acid == "DEL":
            if operation_order == "PO":  # Postfix
                result.pop()
            else:  # Prefix or Infix
                i += 3
        elif amino_acid == "SWAP":
            if operation_order == "PO":  # Postfix
                if len(result) >= 2:
                    result[-1], result[-2] = result[-2], result[-1]
            else:  # Prefix or Infix
                next_codon_1 = ''.join(sequence_list[i+3:i+6])
                next_codon_2 = ''.join(sequence_list[i+6:i+9])
                result.extend([next_codon_2, next_codon_1])
                i += 6
        elif amino_acid == "EXCHANGE":
            if operation_order == "PO":  # Postfix
                # We'll just append the last codon since we don't have logic to fetch an equivalent codon for the same amino acid
                result.append(result[-1])
            else:  # Prefix or Infix
                next_codon = ''.join(sequence_list[i+3:i+6])
                # For now, we'll just append the next codon since we don't have logic to fetch an equivalent codon for the same amino acid
                result.append(next_codon)
                i += 3
        else:
            result.append(codon)
        
        i += 3
    
    # Convert the result list back to a string
    final_sequence = ''.join(result)
    if read_order == "R":
        final_sequence = final_sequence[::-1]
    
    return final_sequence
