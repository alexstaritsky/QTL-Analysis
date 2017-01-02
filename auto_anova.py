"""
Automatische uitvoering van meerdere ANOVA's
Martijn Landman en Alex Staritsky
Datum: 28-12-2016
"""

from scipy import stats

# Lees het locusbestand
def read_file_locus(path):
    raise NotImplementedError

# Lees het treatbestand
def read_file_treat(path):
    raise NotImplementedError

# Maak de a en de b groep
def create_groups(locus_data, treat_data):
    raise NotImplementedError

# Doe een single factor ANOVA 
def single_factor_anova(group_a, group_b):
    f_value, p_value = stats.f_oneway(group_a, group_b)
    return f_value, p_value

# Sla het bestand op
def save_file(path, data, sep="\t"):
    raise NotImplementedError

# Haalt de test datasets voor de ANOVA op
# In Excel: p-value = 0.135202136
#           f-value = 2.254444988
# Waarden komen overeen!, dus we kunnen dit gebruiken!
def get_test_files(path_a="test_set_a.txt", path_b="test_set_b.txt"):
    group_a = []
    group_b = []
    try:
        with open(path_a, "r") as bestand_a:
            for regel in bestand_a:
                group_a.append(float(regel.replace("\n","").replace("\r","").replace("\t","").replace(" ","")))
        with open(path_b, "r") as bestand_b:
            for regel in bestand_b:
                group_b.append(float(regel.replace("\n","").replace("\r","").replace("\t","").replace(" ","")))
    except IOError:
        print("Er was een probleem bij het lezen van het bestand")
    except FileNotFoundError:
        print("Bestand niet gevonden!")
    except ValueError:
        print("Het bestand heeft waarschijnlijk een foute indeling")
    return group_a, group_b
    
# Main functie
def main():
    group_a, group_b = get_test_files()
    print(sum(group_a), len(group_a))
    print(sum(group_b), len(group_b))
    f_value, p_value = single_factor_anova(group_a, group_b)
    print("F-value = {0}\nP-value = {1}".format(f_value, p_value))
    
main()