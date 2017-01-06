"""
Automatische uitvoering van meerdere ANOVA's
Martijn Landman en Alex Staritsky
Datum: 28-12-2016
"""

from scipy import stats

# Gene class om data op te slaan
# p_value is een float
# f_value is een float
# locus is een String
# markers is een lijst met strings die "a" of "b" zijn
class Gene:
    
    def __init__(self, locus, markers):
        self.__p_value = float(0)
        self.__f_value = float(0)
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
    genes = []
    with open(path, "r") as bestand:
        for regel in bestand:
            # Doe iets
            raise NotImplementedError
    bestand.close()
    return genes

# Lees het traitbestand met traitnames -> een lijst met strings die een kommagetal of "*" zijn
def read_file_trait(path):
    traits = []
    with open(path, "r") as bestand:
        for regel in bestand:
            # Doe iets
            raise NotImplementedError
    bestand.close()
    return traits

# Maak de a en de b groep
def create_groups(traits, markers):
    a, b = [], []
    for i in range(len(markers)):
        if markers[i] == "a" and traits[i] != "*":
            a.append(float(traits[i]))
        if markers[i] == "b" and traits[i] != "*":
            b.append(float(traits[i]))
    return a, b

# Sla het eindbestand bestand op
def save_file(path, data, sep="\t"):
    with open(path, "w") as bestand:
        for regel in bestand:
            # Doe iets
            raise NotImplementedError
    bestand.close()

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
    save_file("results_auto_anova.txt", genes)
    
main()