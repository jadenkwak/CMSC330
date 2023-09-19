import os
from src.ribosome import *

curr_dir = os.path.dirname(__file__)

# Need this because windows and unix file system naming is diff
def get_path(f):
  return os.path.join(curr_dir,f)

def test_public_ValidFileReading():
  read_codons(get_path("inputs/codons5.txt"))
  # Encode the keys and compare the results
  assert (encode("START") == "AAA")
  assert (encode("STOP") == "CCC")
  assert (encode("DEL") == "GGGGGGGGGGGG")
  assert (encode("SWAP") == "UUU")
  assert (encode("EXCHANGE") == "ACG")
  assert (encode("Lysine") in ["UGA", "AUC"])
  assert (encode("Tyrosine") == "GGGGGA")
  assert (encode("Byrosine") == "UAC")
  assert (encode("CMSC") == "ACGU")
  assert (encode("LongSine") == "AAACCCGGGUUU")
  assert (encode("START") == "AAA")  # Repeated START key

  # Decode the values and compare the results
  assert (decode("AAA") == "START")
  assert (decode("CCC") == "STOP")
  assert (decode("GGGGGGGGGGGG") == "DEL")
  assert (decode("UUU") == "SWAP")
  assert (decode("ACG") == "EXCHANGE")
  assert (decode("UGA") == "Lysine")
  assert (decode("GGGGGA") == "Tyrosine")
  assert (decode("UAC") == "Byrosine")
  assert (decode("ACGU") == "CMSC")
  assert (decode("AAACCCGGGUUU") == "LongSine")
  assert (decode("AAA") == "START")  # Repeated START key

def test_public_InvalidFileReading():
  read_codons(get_path("inputs/codons6.txt"))
  # Encode the keys and compare the results
  assert (encode("START") == "AAA")
  assert (encode("Lysine5") == "")
  assert (encode("Tyro3sine") == "")
  assert (encode("Byrosine") == "")
  assert (encode("LongSine") == "AUG")
  assert (encode("Hello") == "")
  assert (encode("CMSC") == "")
  
def test_public_encode1():
  read_codons(get_path("inputs/codons1.txt"))
  assert (encode("START Alanine SWAP STOP") == "UAAGCUUGAUAG")
  assert (encode("START Alanine SWAP INVALID STOP") == "UAAGCUUGAUAG")
  assert (encode("START") == "UAA")

def test_public_encode2():
  read_codons(get_path("inputs/codons2.txt"))
  assert (encode("DEL") == "UACC")
  assert (encode("Alanine SWAP") == "GGUUUUUUUACCC")
  assert (encode("INVALD Alanine SWAP") == "GGUUUUUUUACCC")

def test_public_encode3():
  read_codons(get_path("inputs/codons3.txt"))
  assert (encode("SWAP") == "GACU")
  assert (encode("Tyrosine SWAP") == "AAAAUUUGACU")
  assert (encode("THIS IS INVALD") == "")

def test_public_encode4():
  read_codons(get_path("inputs/codons4.txt"))
  assert (encode("Methionine SWAP Invalid") == "AAAAAAACCG")
  assert (encode("Invalid START Glycine") == "GGUACCGGGUUUUUUAA")
  assert (encode("END Invalid EXCHANGE Methionine Tyros1ne") == "UACGAAAAAAACCG")

def test_public_decode1():
  read_codons(get_path("inputs/codons1.txt"))
  assert (decode("UAAGCUUGAUAG") == "START Alanine SWAP STOP")
  assert (decode("UAA") == "START")

def test_public_decode2():
  read_codons(get_path("inputs/codons2.txt"))
  assert (decode("UACC") == "DEL")
  assert (decode("GGUUUUUUUACCC") == "Alanine SWAP")

def test_public_decode3():
  read_codons(get_path("inputs/codons3.txt"))
  assert (decode("AAAAAUUUU") == "Tyrosine")
  assert (decode("CCCCCCCCCCCGGA") == "Lysine")
  assert (decode("AAAUGGGGACU") == "DEL SWAP")
  assert (decode("LAOENVLACOQM") == "")

def test_public_decode4():
  read_codons(get_path("inputs/codons4.txt"))
  assert (decode("GGUACCGGGUACG") == "START EXCHANGE")
  assert (decode("AAAAAAACCGUUUUUUAA") == "Methionine Glycine")
  assert (decode("UACGAAAAAAACCG") == "EXCHANGE Methionine")

def test_public_OperateLeftToRight():
  # START Lysine DEL Alanine Methionine SWAP Lysine Methionine EXCHANGE METHIONINE STOP
  # UAA AAA UAC GCU [AUG,GUA] UGA AAA [AUG,GUA] GGG [AUG,GUA] UAG
  read_codons(get_path("inputs/codons1.txt"))
  read_evals(get_path("inputs/order1.txt"))

  # Lysine Methionine Methionine Lysine Methionine(EXCH)
  assert (operate("UAAAAAUACGCUAUGUGAAAAAUGGGGAUGUAG","evalorder1")=="AAAAUGAUGAAAGUA") # PR

  # Lysine Lysine Methionine Methionine Methionine(EXCH)
  assert (operate("UAAAAAUACGCUAUGUGAAAAAUGGGGAUGUAG","evalorder3")=="AAAAAAAUGAUGGUA") # I

  # Methionine Alanine Lysine Methionine(EXCH) Methionine
  assert (operate("UAAAAAUACGCUAUGUGAAAAAUGGGGAUGUAG","evalorder4")=="AUGGCUAAAGUAAUG") # PO

def test_public_OperateRigthToLeft():
  # STOP Lysine DEL Alanine Methionine SWAP Lysine Methionine EXCHANGE Methionine START
  # UAA AAA UAC GCU [AUG,GUA] UGA AAA [AUG,GUA] GGG [AUG,GUA] UAG
  read_codons(get_path("inputs/codons1.txt"))
  read_evals(get_path("inputs/order1.txt"))

  # If input string processed right to left: 
  # STOP Lysine DEL Alanine Methionine SWAP Lysine Methionine EXCHANGE Methionine START
  # Output of PO operating: Lysine Methionine Methonine Lysine Methionine(EXCH) == GUAAAAAUGAUGAAA
  assert (operate("GAUAAACAUUCGGUAAGUAAAGUAGGGGUAAAU","evalorder2")=="GUAAAAAUGAUGAAA") # PO

  # If input string processed right to left: 
  # STOP Lysine DEL Alanine Methionine SWAP Lysine Methionine EXCHANGE Methionine START
  # Output of PR operating: Methionine Alanine Lysine Methionine(EXCH) Methionine == AUGGUAAAAGCUAUG
  assert (operate("GAUAAACAUUCGGUAAGUAAAGUAGGGGUAAAU","evalorder5")=="AUGGUAAAAGCUAUG") # PR

  # If input string processed right to left: 
  # STOP Lysine DEL Alanine Methionine SWAP Lysine Methionine EXCHANGE Methionine START
  # Output of I operating: Alanine Lysine Methionine Methionine(EXCH) Methionine == AUGGUAAUGAAAGCU
  assert (operate("GAUAAACAUUCGGUAAGUAAAGUAGGGGUAAAU","evalorder6")=="AUGGUAAUGAAAGCU") # I

def test_public_op1():
  read_codons(get_path("inputs/codons1.txt"))
  read_evals(get_path("inputs/order1.txt"))
  # START Lysine DEL Alanine Methionine SWAP Lysine Methionine EXCHANGE METHIONINE STOP
  # UAA AAA UAC GCU [AUG,GUA] UGA AAA [AUG,GUA] GGG [AUG,GUA] UAG
  assert (operate("UAAAAAUGAAUGGCU","evalorder1")=="AAAGCUAUG") 
  assert (operate("UAAAUGAAAGCUUACAUG","evalorder1")=="AUGAAAGCU")
  assert (operate("AAUAAACAUGCUGUAAGUAAAGUAGGGGUAUAG","evalorder2")=="")

def test_public_op2():
  read_codons(get_path("inputs/codons1.txt"))
  read_evals(get_path("inputs/order1.txt"))
  
  # "GAUAGUAAAGUAAAU" -> (STOP, SWAP, Lysine, Methionine, START)
	# Reading from right to left, post-fix & doing the swap -> Methionine, Lysine
	# Then returning the string in the order we read the codons (right to left)-> Lysine, Methionine
  assert (operate("GAUAGUAAAGUAAAU", "evalorder2") =="AAAAUG") 
  assert (operate("GCUUAAAAAAUGGCUUGAAAAUAG", "evalorder3")=="AAAAUGAAAGCU")

def test_public_op3():
  read_codons(get_path("inputs/codons2.txt"))
  read_evals(get_path("inputs/order1.txt"))

  # "UACCC GUA AAA GGUUUUUU UUAA " -> SWAP STOP START Alanine Methionine 
  assert (operate("UACCC GUA AAA GGUUUUUU UUAA ", "evalorder1")=="GGUUUUUUUUAA")
  # "UACCCGUAAAAUACCGGUUUUUUUUAA " -> SWAP STOP START DEL Alanine Methionine
  assert (operate("UACCCGUAAAAUACCGGUUUUUUUUAA ", "evalorder1")=="UUAA")
  assert (operate("AAUUCCCAUAACUUUUGGUUUUUUAAAUACCGGUUUUUU", "evalorder2")=="UUAA")

# Tests noise filtering for encode and decode
def test_public_NoiseFiltering():
  read_codons(get_path("inputs/codons1.txt"))
  read_evals(get_path("inputs/order1.txt"))

  assert (decode("UAAXUGAGFJAAAGCUAUGUAGCCC") == "START SWAP Lysine Alanine Methionine STOP")
  assert (decode("") == "")
  assert (encode("START SWAP Lysine Alanine Methionine STOP") in ["UAAUGAAAAGCUGUAUAG", "UAAUGAAAAGCUAUGUAG"])
