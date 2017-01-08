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
        self.__p_value = float(-1)
        self.__f_value = float(-1)
        self.__locus = locus
        self.__markers = markers
        
    def get_locus(self):
        return self.__locus
        
    def get_markers(self):
        return self.__markers
        
    def get_f_value(self):
        return self.__f_value
        
    def get_p_value(self):
        return self.__p_value
        
    def set_f_value(self, f_value):
        self.__f_value = f_value
        
    def set_p_value(self, p_value):
        self.__p_value = p_value

# Lees het locusbestand met markerdata
def read_file_locus(path):
    regel_nummer = 0
    genes, markers = [], []
    with open(path, "r") as bestand:
        for regel in bestand:
            if regel_nummer > 6:
                if ";" in regel:
                    locus = regel.split(";")[0].replace(" ","")
                else:
                    regel_markers = regel.replace(" ", "").replace("\n","").replace("\r","")
                    for marker in regel_markers:
                        markers.append(marker)
                if len(markers) == 162:
                    genes.append(Gene(locus, markers))
                    markers = []
            regel_nummer += 1
    bestand.close()
    return genes

# Lees het traitbestand met traitnames
def read_file_trait(path):
    regel_nummer = 0
    traits = []
    with open(path, "r") as bestand:
        for regel in bestand:
            if regel_nummer > 7 and "\t" in regel:
                traits.append(regel.split("\t")[1].replace("\n","").replace("\r","").replace(" ",""))
            regel_nummer += 1
    bestand.close()
    return traits

# Maak de a en de b groepen
def create_groups(traits, markers):
    a, b = [], []
    for i in range(0, len(markers)):
        if markers[i] == "a" and traits[i] != "*":
            a.append(float(traits[i]))
        if markers[i] == "b" and traits[i] != "*":
            b.append(float(traits[i]))
    return a, b

# Sla het eindbestand bestand op
def save_file(path, genes, sep="\t"):
    genes.sort(key=lambda gene: gene.get_p_value())
    with open(path, "w") as bestand:
        bestand.write("Locus{0}P-value{0}F-value".format(sep))
        for gene in genes:
            bestand.write("\n{1}{0}{2}{0}{3}".format(sep, gene.get_locus(), gene.get_p_value(), gene.get_f_value()))
    bestand.close()

# Main functie
def main():
    # Lees de bestanden
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
    # Sla de genen in een bestand op
    save_file("results_auto_anova.txt", genes)
    print("Done!")
    
main()