"""
Automatische uitvoering van meerdere ANOVA's
Martijn Landman en Alex Staritsky
Datum: 28-12-2016
"""

from scipy import stats

# Gene class om data op te slaan
class Gene:
    
    def __init__(self, locus, markers):
        self.__p_value = 0
        self.__f_value = 0
        self.__locus = locus
        self.__markers = markers
        
    def get_locus(self):
        return self.__locus
        
    def get_markers(self):
        return self.__markers
        
    def get_f_value(self):
        return self.__p_value
        
    def get_p_value(self):
        return self.__p_value
        
    def set_f_value(self, p_value):
        self.__p_value = p_value
        
    def set_p_value(self, f_value):
        self.__f_value = f_value

# Lees het locusbestand met markerdata
def read_file_locus(path):
    raise NotImplementedError

# Lees het traitbestand met traitnames
def read_file_trait(path):
    raise NotImplementedError

# Maak de a en de b groep
def create_groups(traits, markers):
    raise NotImplementedError

# Sla het eindbestand bestand op
def save_file(path, data, sep="\t"):
    raise NotImplementedError

# Main functie
def main():
    # Lees de bestanden; traits is een lijst met alle traitnames; genes is een lijst met gene objecten die een locus en markerdata bevat
    traits = read_file_trait("CvixLer.qua") 
    genes = read_file_locus("CvixLer.loc")
    
    # Voor elk gen (of locus)
    for gene in genes:
        # Maak de 2 groepen voor de ANOVA
        a, b = create_groups(traits, gene.get_markers())
        # Voer de ANOVA uit
        f, p = stats.f_oneway(a, b)
        # Voeg de f en p waarde van de ANOVA toe aan het gen
        gene.set_f_value(f)
        gene.set_p_value(p)
    # Sla het bestand op (waarschijnlijk tab-delimited op deze manier voor elk gen "locus \t f-value \t p-value \n") en dan geen spaties tussen
    save_file(genes)
    
main()