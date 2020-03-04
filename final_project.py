

from Bio.Seq import Seq
from Bio.Alphabet  import generic_dna, generic_protein, generic_rna
from Bio.Blast import NCBIWWW
import itertools


'''
Generator sekwencji kodujących (wej. sekw. aminokwasów, wyj. zbiór sekwencji kodujących)

[ok] napisany obiektowo,
[ok] musi zapisywać/wczytywać dane z plików,
[ok] musi w nim wystąpić własnoręcznie zdefiniowany iterator (ew. generator),
[ok] dotyczyć tematyki bioinformatycznej,
[ ! ] można korzystać z bibliotek (np. BioPythona),
[ok] musi być zaimplementowany w Pythonie.

'''


class My_Protein:
	def __init__(self):		
		self.seq = ""
		self.codons = {
		'A':['GCT','GCC','GCA','GCG'],
        'G':['GGT','GGC','GGA','GGG'],
        'V':['GTT','GTC','GTA','GTG'],
        'L':['TTA','TTG','CTT','CTC','CTA','CTG'],
        'I':['ATT','ATC','ATA'],
        'P':['CCT','CCC','CCA','CCG'],
        'S':['AGT','AGC'],
        'T':['ACT','ACC','ACA','ACG'],
        'C':['TGT','TGC'],
        'M':['ATG'],
        'F':['TTT','TTC'],
        'Y':['TAT','TAC'],
        'W':['TGG'],
        'D':['GAT','GAC'],
        'N':['AAT','AAC'],
        'E':['GAA','GAG'],
        'Q':['CAA','CAG'],
        'K':['AAA','AAG'],
        'R':['AGA','AGG','CGT','CGC','CGA','CGG'],
        'H':['CAT','CAC'],
        '_':['TAA','TAG','TGA']
		}
		self.input_path = ''
		self.output_path = ''
	def read_in_out(self):
		print("please enter a path to input file (example:'D:\\3 semestr\\SJP\\input.txt')")
		self.input_path = input()
		print("please enter a path to output file (example:'D:\\3 semestr\\SJP\\output.txt')")
		self.output_path = input()
	def reader(self):
		s = ""
		fh = open(self.input_path) #enter the pass
		for line in fh:
			s += line
		self.seq = s
		print("this is seq: ", self.seq)
		
	def how_many_RNAs(self):
		num = 1
		for aa in self.seq:
			num *= len(self.codons[aa])		
		print("Totally we will have {} RNA pre-seqs.".format(num))

	#Function for generation all possible sequences in RNA
	def generator_RNA_seqs(self, seq:str):
		h = [self.codons[aa] for aa in seq]
		for sub_sol in itertools.product(*h):
			yield"".join(sub_sol)
		print("Now i am calculating all posible pre-sequences from RNA.")

	def writer(self):
		set_of_seqs = self.generator_RNA_seqs(self.seq) #reeading->executing
		fh = open(self.output_path, "w")
		index = 1
		for line in set_of_seqs:
			fh.write("{}: {}\n".format(index, line))
			index += 1
		print("Done. Check solutions here: {}".format(path))









# class My_RNA:
# 	def __init__(self):		
# 		self.seq = ""
# 	def reader(self):
# 		s = []
# 		fh = open("D:\\3 semestr\\SJP\\output.txt") #enter the pass
# 		for line in fh:
# 			s.append(line)
# 		self.seq = s
# 		print("this is seq: ", self.seq)


# 	#Function to count GC content
# 	def compute_gc(self):
# 		return len([b for b in self if b in["G", "C"]]) / len(self)







def main():
	print("hello")
	amino = My_Protein()
	amino.read_in_out()
	amino.how_many_RNAs()
	amino.writer()






if __name__ == "__main__": main()