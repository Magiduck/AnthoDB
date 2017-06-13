from flask import Flask, render_template, request
from Bio import Entrez
from Bio import Medline
from nltk.corpus import stopwords
from nltk.tokenize import word_tokenize
import re

app = Flask(__name__)


# Method voor het ophalen van het totaal aantal PubMed artikelen dat overeenkomt met de query.
# Input = mail, query(zoekterm). output = maxval, een int die het aantal hits aangeeft.
def get_number_of_hits(mail, query):
    Entrez.email = mail
    handle = Entrez.egquery(term=query)
    record = Entrez.read(handle)
    for row in record["eGQueryResult"]:
        if row["DbName"] == "pubmed":
            maxval = (row["Count"])

            return maxval


# Method voor het ophalen van abstracts uit PubMed die horen bij de query. input=mail, query, hits(maxval).
# output=records, de medline records die bestaan uit o.a title, abstract & PubMed ID.
def get_records(mail, query, hits):
    Entrez.email = mail
    handle = Entrez.esearch(db="pubmed", term=query, retmax=hits)
    record = Entrez.read(handle)
    idlist = record["IdList"]
    handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
    records = Medline.parse(handle)

    return records


# Method voor het filteren van pubmed abstracts. input(records) output(info_dict,no_abstract).
# info_dict: De key van deze dictionary is een accessiecode met als value een gefilterde abstract.
# no_abstract: lijst van accessiecodes die geen abstract bevatten.
def process_records(records):
    # stop_words is een lijst met meestgebruikte engelse woorden die uit het abstract gefilterd zullen worden.
    stop_words = set(stopwords.words("english"))
    # most_used_words is ook een lijst met veel gebruikte engelse woorden die uit de abstracts gefilter zullen worden.
    most_used_words = ["a", "about", "all", "also", "and", "as", "at", "be", "because", "but", "by", "can", "come",
                       "could", "day", "do", "even", "find", "first", "for", "from", "get", "give", "go", "have", "he",
                       "her", "here", "him", "his", "how", "I", "if", "in", "into", "it", "its", "just", "know", "like",
                       "look", "make", "man", "many", "me", "more", "my", "new", "no", "not", "now", "of", "on", "one",
                       "only", "or", "other", "our", "out", "people", "say", "see", "she", "so", "some", "take", "tell",
                       "than", "that", "the", "their", "them", "then", "there", "these", "they", "thing", "think",
                       "this", "those", "time", "to", "two", "up", "use", "very", "want", "way", "we", "well", "what",
                       "when", "which", "who", "will", "with", "would", "year", "you", "your", ",", ".", ":", ";",
                       "METHODS", "MATERIAL", "SCOPE", "CHEMICAL", "COMPOUNDS", "CONCLUSION", "CID", "AND", "RESULTS",
                       "BACKGROUND", "CONCLUSION", "CONCLUSIONS", "MAIN", ")", "(", "[", "]", "DNA", "RNA", "mRNA",
                       "RNAs", "DISCUSSION", "PubChem", "MATERIALS", "CONTEXT", "MESSAGE", "REVIEW", "KEY", "CO2",
                       "NaCl", "mM", "RNA-Seq", "PCR", "SNP", "QTLs"]
    info_dict = {}
    count = 0
    no_abstract = []

    # Deze for-loop zorgt dat iedere abstract wordt opgebroken in een lijst met woorden. Uit deze lijst worden de
    # woorden uit stop_words en most_used_words gefilterd. Het resultaat van deze for-loop is een gevulde no_abstract
    # en een gevulde info_dict.
    for record in records:
        unfiltered_ab = []
        filtered_ab = []
        count += 1
        try:
            unfiltered_ab = word_tokenize(record["AB"])
        except KeyError:
            no_abstract.append("PMID:" + str(record["PMID"]))

        for word in unfiltered_ab:
            if word not in stop_words and word not in most_used_words:
                filtered_ab.append(word)
        if filtered_ab:
            info_dict[record["PMID"]] = filtered_ab

    return [info_dict, no_abstract]


# Methode voor het vinden van gen namen middels een regular expression input=info_dict output=gene_dict, no_genes
# gene_dict heeft als key een PMID en als value een lijst met daarin een lijst met genen en een lijst met de gefilterde
# abstract. no_genes is een lijst met PubMed ID waarbij geen gen te vinden was in de abstract.
def search_genes(info_dict):
    regex = re.compile(r"[A-z]+\d*[A-Z]")
    no_genes = []
    gene_dict = {}

    for pmid, abstract in info_dict.items():
        genelist = []

        for word in abstract:
            if regex.match(word):
                if word not in genelist and word.startswith("Pub") is not True:
                    genelist.append(word)
        if genelist:
            gene_dict[pmid] = [genelist, abstract]
        else:
            no_genes.append("PMID:" + pmid)

    return [gene_dict, no_genes]


# Methode voor het vinden van stress omstandigheden middels een regular expression input=gene_dict output=result_dict,
# no_stress. result_dict heeft als key een PMID en als value een lijst met daarin een lijst met genen en een lijst met
# stress omstandigheden. no_stress is een lijst met PMIDs met abstracts waarin geen stress omstandigheden
# te vinden zijn.
def search_stress(gene_dict):
    no_stress = []
    result_dict = {}
    for pmid in gene_dict:
        stresslist = []
        for word in gene_dict[pmid][1]:
            if "stress" in word:
                stresslist.append(gene_dict[pmid][1][gene_dict[pmid][1].index(word) - 1] + " " + word)
        if stresslist:
            result_dict[pmid] = [gene_dict[pmid][0], stresslist]
        else:
            no_stress.append("PMID:" + pmid)

    return [result_dict, no_stress]


# Method die zorgt voor de datastructuur die gebruikt wordt voor de json_file_writer. Hiervoor wordt ook gekeken door
# deze method hoe vaak genen voor komen in combinatie met bepaalde stress omstandigheden. De input is result_dict.
# De output is completed_sunburst_dict, deze heeft als key een stress omstandigheid en als value een lijst met
# dictionaries die als key een gen hebben en als value een lijst met accessiecodes waarin de gen stress omstandigheid
# voorkomt.
def search_co_oc(result_dict):
    new_result_dict = {}
    dictionary = {}
    gene_stress_stringlist = []
    appearence_counter = []
    filtered_dict = {}
    final_dict = {}

    # Deze for-loop zorgt er voor dat een stress factor niet 2 keer in de streslijst voor kan komen.
    for PubMedID in result_dict:
        value = result_dict[PubMedID]
        genes = value[0]
        stresses = value[1]
        unique_stress_list = []
        for stress in stresses:
            if stress not in unique_stress_list:
                unique_stress_list.append(stress)
        # new _result_dict: dictionary key=pubmedID values= lijst van genen en lijst van unique stress omstandigheden
        new_result_dict[PubMedID] = [genes, unique_stress_list]

    # Deze for-loopt maakt strings in het formaat gen|-|stress deze string zal dienen als key voor een dictionary
    for PubMedID in new_result_dict:
        # values: een value uit new_result_dict bijv [[gen,gen1],[stres,stres1]]
        values = new_result_dict[PubMedID]
        genes = values[0]
        stresses = values[1]

        # deze geneste for-loop loopt over de gen lijst en stress lijst en maakt
        # een string door gen en stress aan elkaar te koppelen middels een |-| teken
        for gene in genes:
            for stress in stresses:
                gene_stress_stringlist_key = gene + "|-|" + stress
                # gene _stress_stringlist: een lijst me de verbanden die te vinden zijn tussen gen en stress
                gene_stress_stringlist.append(gene_stress_stringlist_key)

    # Deze for-loop maakt een dictionary in het formaat key = gen|-|stess value = hoevaak deze combinatie gevonden is
    for key in gene_stress_stringlist:
        if key not in appearence_counter:  # Als een key nog niet toegevoegd is aan dictionary
            appearence_counter.append(key)
            dictionary[key] = 1
        else:  # Als de key al toegevoegd is aan dictionary word de
            count = dictionary.get(key)  # value van deze opgehaald n word er 1 bij de value opgeteld
            count += 1
            dictionary[key] = count

    # In deze for-loop wordt de uiteindelijke count van een combinatie gen|-|stress opgehaald uit de dictionary
    for key in dictionary:
        count = dictionary[key]

        # In deze for-loop wordt er een dictionary gemaakt in het formaat key = gen|-|stress
        # value = [hoevaak de combinatie voorkomt, in welke abstracts de combinatie voorkomt]
        for PubMedID in result_dict:
            value = result_dict[PubMedID]
            genes = value[0]
            stresses = value[1]
            for gene in genes:
                for stress in stresses:
                    check = gene + "|-|" + stress
                    # Als de gen stress combinaite gelijk is aan de combinaite key in dictionary word het pubmed id(key)
                    # van result_dict toegevoegd aan een lijst, dit word zo gedaan dat er alleen unieke codes in de
                    # lijst aanwezig zijn.
                    if key == check:
                        if key not in final_dict:
                            pubmed_id_list = [PubMedID]
                            final_dict[key] = count, pubmed_id_list
                        if key in final_dict:
                            acc_codes = final_dict.get(key)
                            acc_code = acc_codes[1]
                            if PubMedID not in acc_code:
                                acc_code.append(PubMedID)
                                # Dictionary: key=gen|-|stress
                                # value= hoevaak de combinatie voorkomt, in welke abstracts de combinatie voorkomt.
                                final_dict[key] = count, acc_code
                            if PubMedID in acc_code:
                                final_dict[key] = count, acc_code
        # for-loop om eventueel relaties weg te filteren als deze maar in 1 keer in abstacts voorkomen om zo
        # significantere relaties over te houden. Om het filter te uit te zetten gebruik >=,
        # om het filter aan tezetten gebruik >, gewenste aantal voorkomen van een gen stress relatie.
        for relation in final_dict:
            appearances = final_dict[relation]
            if appearances[0] >= 1:
                filtered_dict[relation] = appearances

    stress_list = []
    sunburst_dict = {}

    # Deze for loopt maakt een dic die geschikt is om een json file van te maken voor de sunburst.
    # key=stess value=lijst van genen die een veband hebben met deze stress factor
    for gene_stress in filtered_dict:
        genes = gene_stress.split('|-|', 1)[0]
        stresses = gene_stress.split('|-|', 1)[1]
        stresses = stresses.lower()
        gene_list = []
        # Als stress nog niet in de stress_list voorkomt word deze een nieuwe key voor sunburst_dict
        if stresses not in stress_list:
            # Deze loop voegt dictionarys samen afhankelijk van de stress
            for gene_stress_2 in filtered_dict:
                genes_2 = gene_stress_2.split('|-|', 1)[0]
                stress_2 = gene_stress_2.split('|-|', 1)[1].lower()
                stress_list.append(stresses)
                # Als de stresomstandigheid van 2 verbanden uit gen_stress hetzelfde zijn
                if stresses == stress_2:
                    # Als gen 1 en 2 nog niet aanwezig zijn in genlist worden deze hieraan toegevoegd
                    if genes not in gene_list:
                        gene_list.append(genes)
                    if genes_2 not in gene_list:
                        gene_list.append(genes_2)
            # sunburst_dict: key= stressomstandigheid
            # value=lijst van unique genen die een veband hebben met de stres omstandigheid
            sunburst_dict[stresses] = gene_list

    completed_sunburst_dict = {}

    # deze for-loop maakt completed_sunburst_dict die geschikt is om een json file van te maken voor de sunburst.
    # Deze dictionary ziet er als volgt uit {stres1{gen{[acc,acc,acc]}gen{[acc,acc,acc]}},stres2{gen{[acc,acc,acc]}}}.
    # Completed_sunburst_dict wordt gemaakt door de accessiecodes van filtered_dict te koppelen aan een
    # bijhorend gen in sunburst_dict. De koppeling vindt plaats wanneer het gen en de stress uit filtered_dict
    # overeenkomen met het gen(value) en de stress(key) uit sunburst_dict.
    for stress_dict in sunburst_dict:
        gene_list = []
        for gene_stress in filtered_dict:  #
            genes = gene_stress.split('|-|', 1)[0]
            stresses = gene_stress.split('|-|', 1)[1]
            if stresses == stress_dict:
                acc_list = filtered_dict[gene_stress][1]
                for word in sunburst_dict[stress_dict]:
                    gen_code_dict = {}
                    if word == genes:
                        gen_code_dict[word] = acc_list
                        gene_list.append(gen_code_dict)
        completed_sunburst_dict[stress_dict] = gene_list

    return completed_sunburst_dict


# Deze method gebruikt search= zoekterm van de gebruiker en completed_sunburst_dict om een json file te schrijven.
# Dit wordt gedaan door over completed_sunburst te itereren en via het juiste format informatie weg te schrijven in
# flare.json. De zoekterm van de gebruiker is altijd de root van de sunburst. Als de json file volledig is geschreven
# wordt "done" geretourneerd.
def json_file_writer(search, completed_sunburst_dict):
    file_handler = open('templates/flare.json', "w")
    file_handler.write("{\n")
    file_handler.write('"name": "' + str(search) + '", \n')
    file_handler.write('"children": [ \n')
    count = 0

    for stress in completed_sunburst_dict:
        count += 1
        file_handler.write("\t{\n")
        file_handler.write('\t"name": "' + stress + '",\n')
        file_handler.write('\t"children": [\n')
        index = 0

        for gene_dict in completed_sunburst_dict[stress]:
            for gene, acc_code in gene_dict.items():
                index += 1
                if index == len(completed_sunburst_dict[stress]):
                    file_handler.write('\t\t{"name": "' + gene + '", "size": ' + str(len(acc_code)) + '}\n')
                else:
                    file_handler.write('\t\t{"name": "' + gene + '", "size": ' + str(len(acc_code)) + '},\n')
        file_handler.write('\t\t]\n')
        if count == len(completed_sunburst_dict):
            file_handler.write('\t}\n')
        else:
            file_handler.write('\t},\n')

    file_handler.write('\t]\n')
    file_handler.write('}')
    file_handler.close()

    return "done"


# Method om alle andere methods van de logica in de juiste volgorde aan te roepen. Deze main method neemt als input
# de zoekterm van de gebruiker en retourneerd done als alle berekeningen zijn uitgevoerd en er een json file is
# geschreven waarmee de resultaten gevisualiseerd kunnen worden in een sunburst.
def main(search):
    mail = "davidmaarn@gmail.com"
    query = search + " AND Anthocyanin"

    hits = get_number_of_hits(mail, query)
    records = get_records(mail, query, hits)
    info_dict = process_records(records)
    gene_dict = search_genes(info_dict[0])
    result_dict = search_stress(gene_dict[0])
    co_occurrence = search_co_oc(result_dict[0])
    done = json_file_writer(search, co_occurrence)

    return done


# App.route voor de home pagina van de website
@app.route('/', methods=["GET"])
def mainpage():
    return render_template('index.html')


# App.route voor de help pagina van de website
@app.route('/help', methods=["GET"])
def help_page():
    return render_template('help.html')


# App.route voor de contact pagina van de website
@app.route('/contact', methods=["GET"])
def contact_page():
    return render_template('contact.html')


# App.route voor de visualsatie pagina van de website. Deze route wordt gebruikt zo gauw op de visualise knop is
# gedrukt in textminedone.html
@app.route('/visualisatie', methods=["GET"])
def visualisation_page():
    return render_template('visualisatie.html')


# App.route voor het textmining gedeelte van de website. Deze route haalt de input van de gebruiker op en geeft hem
# mee aan main(). De berekeningen worden uitgevoerd en en json_file wordt geschreven. Zo gauw done == "done" wordt
# textminedone.html geretourneerd zo kan op de knop gedrukt worden voor het visualiseren van de resultaten in een
# sunburst.
@app.route('/textmine', methods=['GET', 'POST'])
def textmine():
    if request.method == "POST":
        query = request.form['query']
        done = main(query)
        if done == "done":
            return render_template("textminedone.html", query=query)
    return render_template("textmine.html")


if __name__ == '__main__':
    app.run()
