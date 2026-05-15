import requests
import sys
import re

CODON_TAB = {
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "AGA": "R", "AGG": "R", "AGU": "S", "AGC": "S", "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L"
}

def url_content(url):
    try:
        response = requests.get(url)

        if not response.ok:
            print(f"ERROR: Status - {response.status_code}")
            sys.exit()

        return response.text

    except requests.exceptions.RequestException as e:
        print(f"ERROR: No se pudo establecer conexion con el servidor.")
        sys.exit()


def get_dna(text:str) -> str:
    raw_matches = re.findall(r"[ATCG][ATCG\s]{10,}[ATCG]", text, re.IGNORECASE)
    if not raw_matches:
        return None

    # el tercer argumento quita espacios y saltos de lineas
    table = str.maketrans("", "", " \n\r\t")
    clean_sequences = [sequence.translate(table).upper() for sequence in raw_matches]

    # devuelve el elemento mas largo en caso de que haya mas de uno ("ruido")
    return max(clean_sequences, key=len)


def transcription(dna):
    return dna.replace("T", "U")


def translation(mrna):
    start = mrna.find("AUG")
    if start == -1:
        return None
    
    codons = [
    mrna[i:i+3]
    for i in range(start, len(mrna), 3)
    if len(mrna[i:i+3]) == 3 
    ]
    
    prot = []
    
    for i in codons:
        aa = CODON_TAB.get(i, "?")
        
        if aa == "*":
            break
        prot.append(aa)
    return "".join(prot)
    

def main():
    url = "https://raw.githubusercontent.com/Anscharius-anq/PILE/refs/heads/main/tarea_1.txt"

    content = url_content(url)

    dna = get_dna(content)
    if  dna == None:
        print("No se encontró una secuencia de ADN válida.")
        return

    mrna = transcription(dna)
    if mrna == None:
        print("No se encontró secuencia de inicio (AUG)")
        return

    protein = translation(mrna)

    print("DNA:", dna)
    print("mRNA:", mrna)
    print("Protein:", protein)

if __name__ == "__main__":
    main()